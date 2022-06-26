#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "timer.h"
#include <mpi.h>
#include <unistd.h>

#define MASTER 0

#define SOFTENING 1e-9f

typedef struct
{
    float x, y, z, vx, vy, vz;
} Body;

void determisticInitBodies(float *buf, int n)
{
    Body *p = (Body *)buf;
    for (int i = 0; i < n; i++)
    {
        p[i].x = 0.0 + i + 1;
        p[i].y = 0.0 + i + 1;
        p[i].z = 0.0 + i + 1;
        p[i].vx = 0.0f;
        p[i].vy = 0.0f;
        p[i].vz = 0.0f;
    }
}

void printBodies(float *buf, int n)
{
    Body *p = (Body *)buf;
    for (int i = 0; i < n; i++)
    {
        printf("%f ", p[i].x);
    }
    printf("\n");
}

void print_results(int rank, Body *p_recv, int nBodies)
{
    printf("proc %d: ", rank);
    printBodies((float *)p_recv, nBodies);
    printf("\n");
}

void prepere_send_buffer(Body *p, Body *p_recv, int own_portion, int start_own_portion)
{

    // ! init bodyforce for first call
    for (int i = 0; i < own_portion; i++)
    {
        p_recv[start_own_portion + i].x = p[start_own_portion + i].x;
        p_recv[start_own_portion + i].y = p[start_own_portion + i].y;
        p_recv[start_own_portion + i].z = p[start_own_portion + i].z;

        p_recv[start_own_portion + i].vx = p[start_own_portion + i].vx;
        p_recv[start_own_portion + i].vy = p[start_own_portion + i].vy;
        p_recv[start_own_portion + i].vz = p[start_own_portion + i].vz;
    }
}

void bodyForceSplit(Body *p, float dt, int own_portion, int start_own_portion, int start, int end, float Fx[], float Fy[], float Fz[])
{
    for (int i = 0; i < own_portion; i++)
    {
        for (int j = start; j < end; j++) // start = inizio nuova porzione da confrontare end = fine porzione da confrontare
        {
            float dx = p[j].x - p[start_own_portion + i].x;
            float dy = p[j].y - p[start_own_portion + i].y;
            float dz = p[j].z - p[start_own_portion + i].z;

            float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
            float invDist = 1.0f / sqrt(distSqr);
            float invDist3 = invDist * invDist * invDist;

            Fx[i] += dx * invDist3;
            Fy[i] += dy * invDist3;
            Fz[i] += dz * invDist3;
        }
    }
}

void integratePositionSplit(Body *p, float dt, int own_portion, int start_own_portion, float Fx[], float Fy[], float Fz[])
{
    for (int i = 0; i < own_portion; i++)
    {
        p[start_own_portion + i].vx += dt * Fx[i];
        p[start_own_portion + i].vy += dt * Fy[i];
        p[start_own_portion + i].vz += dt * Fz[i];

        p[start_own_portion + i].x += p[start_own_portion + i].vx * dt;
        p[start_own_portion + i].y += p[start_own_portion + i].vy * dt;
        p[start_own_portion + i].z += p[start_own_portion + i].vz * dt;
    }
}

void wait_some_work(MPI_Request bcast_recv[], int request_rank_indices[], int requests_ranks[], int worker_to_wait, int rank, Body *p, float Fx[], float Fy[], float Fz[], int proc_portion_size[], int proc_portion_start[])
{
    int ready_req = 0;
    int count = worker_to_wait;
    int start, end;
    int req_compl;
    const float dt = 0.01f;
    int own_portion = proc_portion_size[rank];
    int start_own_portion = proc_portion_start[rank];

    MPI_Waitsome(worker_to_wait, bcast_recv, &ready_req, request_rank_indices, MPI_STATUS_IGNORE);
    while (ready_req != MPI_UNDEFINED)
    {
        count -= ready_req;
        for (int i = 0; i < ready_req; i++)
        {
            req_compl = requests_ranks[request_rank_indices[i]];
            start = proc_portion_start[req_compl];
            end = start + proc_portion_size[req_compl];

            bodyForceSplit(p, dt, own_portion, start_own_portion, start, end, Fx, Fy, Fz);
        }
        MPI_Waitsome(worker_to_wait, bcast_recv, &ready_req, request_rank_indices, MPI_STATUS_IGNORE);
    }
    if (count > 0)
    {
        for (int i = 0; i < count; i++)
        {
            req_compl = requests_ranks[request_rank_indices[i]];
            start = proc_portion_start[req_compl];
            end = start + proc_portion_size[req_compl];

            bodyForceSplit(p, dt, own_portion, start_own_portion, start, end, Fx, Fy, Fz);
        }
    }
}

void determisticInitBodiesSplit(float *buf, int own_portion, int start_onw_portion)
{
    Body *p = (Body *)buf;
    for (int i = 0; i < own_portion; i++)
    {
        p[start_onw_portion + i].x = 0.0 + start_onw_portion + i + 1;
        p[start_onw_portion + i].y = 0.0 + start_onw_portion + i + 1;
        p[start_onw_portion + i].z = 0.0 + start_onw_portion + i + 1;
        p[start_onw_portion + i].vx = 0.0f;
        p[start_onw_portion + i].vy = 0.0f;
        p[start_onw_portion + i].vz = 0.0f;
    }
}

void calculate_portions(int proc_portion_size[], int proc_portion_start[], int n_workers, int portion, int offset)
{
    int i;
    int j;
    if (offset > 0)
    {
        for (i = 0; i < offset; i++, j++)
        {
            proc_portion_size[i] = portion + 1;
            proc_portion_start[i] = (portion * i) + i;
        }
        for (i = offset; i < n_workers; i++)
        {
            proc_portion_size[i] = portion;
            proc_portion_start[i] = portion * i + offset;
        }
    }
    else
    {
        for (int i = 0; i < n_workers; i++)
        {
            proc_portion_size[i] = portion;
            proc_portion_start[i] = portion * i;
        }
    }
}

void calculate_gatherv_portions(int recvcounts[], int displs[], int n_proc, int proc_portion_size[], int proc_portion_start[])
{
    displs[0] = 0;
    recvcounts[0] = 0;
    for (int i = 1, j = 0; i < n_proc; i++, j++)
    {
        recvcounts[i] = proc_portion_size[j] * 6;
        displs[i] = proc_portion_start[j] * 6;
    }
}

void initRequestsRanks(int requests_ranks[], int worker_rank, int n_workers)
{
    int i, j;

    for (i = 0, j = 0; i < worker_rank; i++, j++)
        requests_ranks[j] = i;

    for (i = worker_rank + 1; i < n_workers; i++, j++)
        requests_ranks[j] = i;
}

int main(int argc, char *argv[])
{
    //! MPI INIT
    MPI_Init(&argc, &argv);
    int comm_word_rank, n_proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_word_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    //! INIT INPUT PROGRAM
    int nBodies = 3000;
    if (argc > 1)
        nBodies = atoi(argv[1]);
    //! INIT TIMERS
    double totalTime = 0.0; // inizializzazione timer
    double end_time, start_time;
    double tElapsed;
    double avgTime;
    //! INIT PROGRAM VARIABLES AND CONSTANTS
    const float dt = 0.01f; // time step
    const int nIters = 10;  // simulation iterations
    int bytes = nBodies * sizeof(Body);
    float *buf = (float *)malloc(bytes);
    Body *p = (Body *)buf;
    //! DEFINE BUFFER FOR RECV RESULTS
    float *buf_recv = (float *)malloc(bytes);
    Body *p_recv = (Body *)buf_recv;
    //! CALCULATE WORKLOAD DISTRIBUCTION
    int portion = nBodies / (n_proc - 1);
    int offset = nBodies % (n_proc - 1);
    int proc_portion_size[n_proc - 1];
    int proc_portion_start[n_proc - 1];
    int recvcounts[n_proc];
    int displs[n_proc];
    //! DEFINE VARIABLES FOR WORKLOAD INTERVALS
    int own_portion = 0, start_own_portion = 0;
    int start = 0, end = 0;
    //! CREATE COMMUNICATOR for workers
    int n_workers;
    int worker_rank;
    MPI_Group orig_group, workers_group;
    MPI_Comm WORKERS_COMM;
    MPI_Comm_group(MPI_COMM_WORLD, &orig_group);
    int excl[1] = {0};
    MPI_Group_excl(orig_group, 1, excl, &workers_group);
    MPI_Comm_create(MPI_COMM_WORLD, workers_group, &WORKERS_COMM);
    calculate_portions(proc_portion_size, proc_portion_start, n_proc - 1, portion, offset);
    if (comm_word_rank != 0)
    {
        MPI_Comm_create_group(MPI_COMM_WORLD, workers_group, 1, &WORKERS_COMM);
        MPI_Comm_rank(WORKERS_COMM, &worker_rank);
        MPI_Comm_size(WORKERS_COMM, &n_workers);

        own_portion = proc_portion_size[worker_rank];
        start_own_portion = proc_portion_start[worker_rank];
        start = start_own_portion;
        end = start + own_portion;
    }
    float Fx[own_portion];
    float Fy[own_portion];
    float Fz[own_portion];

    calculate_gatherv_portions(recvcounts, displs, n_proc, proc_portion_size, proc_portion_start);

    //! START PROCESS WORK
    determisticInitBodiesSplit(buf, own_portion, start_own_portion);
    // determisticInitBodies(buf, nBodies);

    //! DEFINE PERSISTENT SEND REQUEST, RECV REQUEST, RESULT REQUEST

    int reqs_ranks[n_workers - 1];
    int reqs_ranks_indices[n_workers - 1];
    if (comm_word_rank != 0)
        initRequestsRanks(reqs_ranks, worker_rank, n_workers);

    //! DEFINE REQUESTS FOR PORTIONS EXCHANGES
    MPI_Request bcast_recv[n_workers];
    MPI_Request bsend_prec = MPI_REQUEST_NULL;
    MPI_Request bsend_next = MPI_REQUEST_NULL;
    //! DEFINE POINTER FOR PREC AND NEXT NONBLOCKING COMMUNICATION REQUESTS
    MPI_Request *p_bsend_next = &bsend_next;
    MPI_Request *p_bsend_prec = &bsend_prec;
    MPI_Request *p_tmp_bsend = &bsend_prec;
    //! DEFINE TMP VARIABLES
    Body *tmp_p_recv = p_recv;
    Body *tmp;
    int i, j = 0;
    MPI_Request igat = MPI_REQUEST_NULL;

    int recvcounts_all[n_workers];
    int displs_all[n_workers];

    if (comm_word_rank != 0)
    {
        for (int i = 1, j = 0; j < n_workers; i++, j++)
        {
            recvcounts_all[j] = proc_portion_size[j] * 6;
            displs_all[j] = proc_portion_start[j] * 6;
        }
    }

    if (comm_word_rank == MASTER)
        start_time = MPI_Wtime();

    for (int iter = 1; iter <= nIters; iter++)
    {
        if (comm_word_rank == MASTER)
        {
            StartTimer();
            MPI_Wait(&igat, MPI_STATUS_IGNORE);
            tmp_p_recv = p_recv;
            MPI_Igatherv(&p_recv, own_portion * 6, MPI_FLOAT, p_recv, recvcounts, displs, MPI_FLOAT, MASTER, MPI_COMM_WORLD, &igat);
        }
        else
        {
            for (i = 0, j = 0; i < worker_rank; i++, j++)
                MPI_Iallgatherv(p, proc_portion_size[i], MPI_FLOAT, p, recvcounts_all, displs_all, MPI_FLOAT, WORKERS_COMM, &bcast_recv[j]);

            MPI_Iallgatherv(p, proc_portion_size[i], MPI_FLOAT, p, recvcounts_all, displs_all, MPI_FLOAT, WORKERS_COMM, p_bsend_next);

            for (i = worker_rank + 1; i < n_workers; i++, j++)
                MPI_Iallgatherv(p, proc_portion_size[i], MPI_FLOAT, p, recvcounts_all, displs_all, MPI_FLOAT, WORKERS_COMM, &bcast_recv[j]);
            // INIT BODYFORCE FOR ITERATION i
            for (int i = 0; i < own_portion; i++) { Fx[i] = 0.0f; Fy[i] = 0.0f; Fz[i] = 0.0f; }
            // BODY FORCE COMPUTATION : CALCULATION OWN PORTION
            bodyForceSplit(p, dt, own_portion, start_own_portion, start, end, Fx, Fy, Fz);
            // BODY FORCE COMPUTATION : WAIT OTHER PORTIONS FROM OTHER PROCESSES FOR CONTINUE OWN COMPUTATION
            wait_some_work(bcast_recv, reqs_ranks_indices, reqs_ranks, n_workers - 1, worker_rank, p, Fx, Fy, Fz, proc_portion_size, proc_portion_start);
            // WAITS FOR THE PREVIOUS REQUEST TO WRITE IN THE SEND BUFFER
            MPI_Wait(p_bsend_prec, MPI_STATUS_IGNORE); // deve aspettare che invia a tutti i processi 
            MPI_Wait(&igat, MPI_STATUS_IGNORE); // deve aspettare che invia il risultato al padre 
            // quando ha finito di inviare la porzione allora può aggiornare il buffer
            prepere_send_buffer(p, p_recv, own_portion, start_own_portion);
            integratePositionSplit(p_recv, dt, own_portion, start_own_portion, Fx, Fy, Fz);
            
            MPI_Igatherv(&p_recv[start_own_portion], own_portion * 6, MPI_FLOAT, p_recv, recvcounts, displs, MPI_FLOAT, MASTER, MPI_COMM_WORLD, &igat);
        }
        //! SWAP BUFFER POINTERS FOR NEXT ITERATION
        tmp = p;
        p = p_recv;
        p_recv = tmp;
        //! SWAP REQEUST POINTER FOR NEXT ITERATION
        p_tmp_bsend = p_bsend_prec;
        p_bsend_prec = p_bsend_next;
        p_bsend_next = p_tmp_bsend;
        
        if (comm_word_rank == MASTER)
        {    
            const double tElapsed = GetTimer() / 1000.0;
            if(iter>1)
            {
                print_results(comm_word_rank, p_recv, nBodies);
                totalTime += tElapsed;
            }
            #ifndef SHMOO
            printf("Iteration %d: %.3f seconds\n", iter, tElapsed);
            #endif
        }
    }
    //! END PROCESES WORK AND PRINT TOTAL TIME FROM MASTER
    if (comm_word_rank == MASTER)
    {
        MPI_Wait(&igat, MPI_STATUS_IGNORE);
        print_results(comm_word_rank, tmp_p_recv, nBodies);
        end_time = MPI_Wtime();
        double avgTime = totalTime / (double)(nIters - 1);
        #ifdef SHMOO
        printf("%d, %0.3f\n", nBodies, 1e-9 * nBodies * nBodies / avgTime);
        #else
        printf("%d Bodies: average %0.3f Billion Interactions / second\n", nBodies, 1e-9 * nBodies * nBodies / avgTime);
        #endif
        printf("time : %f \n", (end_time - start_time));
    }

    // RESOURCES RELEASE
    free(buf);
    free(buf_recv);

    MPI_Finalize();
    return 0;
}
