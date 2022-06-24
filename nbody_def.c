
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "timer.h"
#include <mpi.h>

#define SOFTENING 1e-9f
#define MASTER 0

typedef struct { float x, y, z, vx, vy, vz; } Body;

FILE *output_file;
char file_name[50];

// PARALLEL PROGRAM FUNCTIONS

void determisticInitBodiesSplit(float *buf, int own_portion, int start_onw_portion)
{
    Body *p = (Body *)buf;
    for (int i = 0; i < own_portion; i++)
    {
        p[start_onw_portion + i].x = start_onw_portion + i + 1;
        p[start_onw_portion + i].y = start_onw_portion + i + 1;
        p[start_onw_portion + i].z = start_onw_portion + i + 1;
        p[start_onw_portion + i].vx = 0.0f;
        p[start_onw_portion + i].vy = 0.0f;
        p[start_onw_portion + i].vz = 0.0f;
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

void waitSomeWork(MPI_Request bcast_recv[], int request_rank_indices[], int requests_ranks[], int n_recv_req, int rank, Body *p, float Fx[], float Fy[], float Fz[], int proc_portion_size[], int proc_portion_start[])
{
    int ready_req = 0;
    int count = n_recv_req;
    int start, end;
    int req_compl;
    const float dt = 0.01f;
    int own_portion = proc_portion_size[rank];
    int start_own_portion = proc_portion_start[rank];

    MPI_Waitsome(n_recv_req, bcast_recv, &ready_req, request_rank_indices, MPI_STATUS_IGNORE);
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
        MPI_Waitsome(n_recv_req, bcast_recv, &ready_req, request_rank_indices, MPI_STATUS_IGNORE);
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

void prepereSendBuffer(Body *p, Body *p_recv, int own_portion, int start_own_portion)
{

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

void calculatePortions(int proc_portion_size[], int proc_portion_start[], int n_workers, int portion, int offset)
{
    proc_portion_size[0] = portion;
    proc_portion_start[0] = 0;

    int i, j;

    if (offset > 0)
    {
        int i, j;
        for (i = 1, j = 0; j < offset; i++, j++)
        {
            proc_portion_size[i] = portion + 1;
            proc_portion_start[i] = portion * i + j;
        }
        for (; i < n_workers; i++, j++)
        {
            proc_portion_size[i] = portion;
            proc_portion_start[i] = portion * i + offset;
        }
    }
    else
    {
        for (int i = 0, j = 0; i < n_workers; i++, j++)
        {
            proc_portion_size[i] = portion;
            proc_portion_start[i] = portion * j;
        }
    }
}

void calculateAllgathervPortions(int recvcounts[], int displs[], int n_workers, int proc_portion_size[], int proc_portion_start[])
{
    for (int i = 0; i < n_workers; i++)
    {
        recvcounts[i] = proc_portion_size[i] * 6;
        displs[i] = proc_portion_start[i] * 6;
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

void swapBuffers(Body *p, Body *p_send)
{
    Body* tmp_buf_swap = p;
    *p = *p_send;
    p_send = tmp_buf_swap;
}

void swapRequests(MPI_Request *allgath_pointer_prec_req, MPI_Request *allgath_pointer_next_req)
{
    MPI_Request *tmp_allgath_send_req_swap = allgath_pointer_prec_req;
    *allgath_pointer_prec_req = *allgath_pointer_next_req;
    allgath_pointer_next_req = tmp_allgath_send_req_swap;       
}

// SEQUENTIAL PROGRAM FUNCTION

void determisticInitBodies(float *buf, int n)
{
    Body *p = (Body *)buf;
    for (int i = 0; i < n; i++)
    {
        p[i].x = i + 1;
        p[i].y = i + 1;
        p[i].z = i + 1;
        p[i].vx = 0.0f;
        p[i].vy = 0.0f;
        p[i].vz = 0.0f;
    }
}

void bodyForce(Body *p, float dt, int n)
{
    for (int i = 0; i < n; i++)
    {
        float Fx = 0.0f;
        float Fy = 0.0f;
        float Fz = 0.0f;

        for (int j = 0; j < n; j++)
        {
            float dx = p[j].x - p[i].x;
            float dy = p[j].y - p[i].y;
            float dz = p[j].z - p[i].z;
            float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
            float invDist = 1.0f / sqrtf(distSqr);
            float invDist3 = invDist * invDist * invDist;

            Fx += dx * invDist3;
            Fy += dy * invDist3;
            Fz += dz * invDist3;
        }

        p[i].vx += dt * Fx;
        p[i].vy += dt * Fy;
        p[i].vz += dt * Fz;
    }
}

void integratePosition(Body *p, float dt, int nBodies)
{
    for (int i = 0; i < nBodies; i++) // integrate position
    {
        p[i].x += p[i].vx * dt;
        p[i].y += p[i].vy * dt;
        p[i].z += p[i].vz * dt;
    }
}

// COMMON FUNCTION FOR TESTING RESULTS
void printResults(Body *p, int nBodies)
{
    for (int i = 0; i < nBodies; i++)
        fprintf(output_file, "%f ", p[i].x);

    fprintf(output_file, "\n");
}

int main(int argc, char *argv[])
{
    // INIT TIMERS
    double totalTime = 0.0;
    double end_time, start_time;
    double tElapsed;
    double avgTime;
    // MPI INIT
    MPI_Init(&argc, &argv);
    int worker_rank, n_workers , print_res = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &worker_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_workers);
    // INIT INPUT PROGRAM
    int nBodies = 30000;
    if (argc > 1)
        nBodies = atoi(argv[1]);
    // INIT PROGRAM VARIABLES AND CONSTANTS
    const float dt = 0.01f; // time step
    const int nIters = 10;  // simulation iterations

    MPI_Barrier(MPI_COMM_WORLD);
    if (worker_rank == MASTER)
        start_time = MPI_Wtime();
    
    int bytes = nBodies * sizeof(Body);
    float *buf = (float *)malloc(bytes);
    Body *p = (Body *)buf;

    if (n_workers == 1) // SEQUENTIAL PROGRAM
    {
        if (argc > 2 && strcmp(argv[2],"-t") == 0 )
        {
            output_file = fopen("sequential_program", "w+");
            print_res = 1;
        }
        determisticInitBodies(buf, nBodies);
        for (int iter = 1; iter <= nIters; iter++)
        {
            StartTimer();
            bodyForce(p, dt, nBodies);
            integratePosition(p, dt, nBodies);
            tElapsed = GetTimer() / 1000.0;
            if (iter > 1)
                totalTime += tElapsed;
            printf("Iteration %d: %.3f seconds\n", iter, tElapsed);
            if (print_res == 1)
                printResults(p, nBodies);
        }
    }
    else // PARALLEL PROGRAM
    { 
        if (argc > 2 && strcmp(argv[2],"-t") == 0 )
        {
            sprintf(file_name, "parallel_%d", n_workers);
            output_file = fopen(file_name, "w+");
            print_res = 1;
        }
        // CALCULATE WORKLOAD DISTRIBUCTION
        int portion = nBodies / (n_workers);
        int rest = nBodies % (n_workers);
        int proc_portion_size[n_workers];
        int proc_portion_start[n_workers];
        calculatePortions(proc_portion_size, proc_portion_start, n_workers, portion, rest);
        int recvcounts[n_workers];
        int displs[n_workers];
        calculateAllgathervPortions(recvcounts, displs, n_workers, proc_portion_size, proc_portion_start);
        // DEFINE VARIABLES FOR WORKLOAD INTERVALS
        int own_portion = proc_portion_size[worker_rank];
        int start_own_portion = proc_portion_start[worker_rank];
        int start = start_own_portion;
        int end = start + own_portion;
        float Fx[own_portion];
        float Fy[own_portion];
        float Fz[own_portion];
        // START PROCESS WORK
        determisticInitBodiesSplit(buf, own_portion, start_own_portion);
        // DEFINE REQUESTS RANKS FOR MPI_WAIT SOME
        int reqs_ranks[n_workers - 1];
        int reqs_ranks_indices[n_workers - 1];
        initRequestsRanks(reqs_ranks, worker_rank, n_workers);
        // DEFINE REQUESTS FOR PORTIONS EXCHANGES
        MPI_Request allgath_recv[n_workers];
        MPI_Request allgath_send_prec = MPI_REQUEST_NULL;
        MPI_Request allgath_send_next = MPI_REQUEST_NULL;
        // DEFINE POINTER FOR PREC AND NEXT NONBLOCKING COMMUNICATION REQUESTS
        MPI_Request *allgath_pointer_next_req = &allgath_send_next;
        MPI_Request *allgath_pointer_prec_req = &allgath_send_prec;
        MPI_Request *tmp_allgath_send_req_swap = &allgath_send_prec;
        // DEFINE BUFFER FOR SEND AND NEW DATA WHILE PREC SEND REQUEST HAS NOT FINISHED
        float *buf_send = (float *)malloc(bytes);
        Body *p_send = (Body *)buf_send;
        // DEFINE TMP VARIABLES
        int i, j;
        Body *tmp_buf_swap;
        int proc_send_req_compltd;
        MPI_Request pending_allgath_req = MPI_REQUEST_NULL;
        MPI_Request proc_send_reqs[2];
        proc_send_reqs[0] = *allgath_pointer_prec_req;
        proc_send_reqs[1] = *allgath_pointer_next_req;
        int send_req_curr;
        int send_req_prec;
        // START i-th PROCES WORK
        for (int iter = 1; iter <= 10; iter++)
        {
            if (worker_rank == MASTER)
                StartTimer();
            // PORTION EXCHANGE BETWEEN PROCESSES: RECEIVE PORTION FROM OTHER PROCESSES
            for (i = 0, j = 0; i < worker_rank; i++, j++)
                MPI_Iallgatherv(p, proc_portion_size[i], MPI_FLOAT, p, recvcounts, displs, MPI_FLOAT, MPI_COMM_WORLD, &allgath_recv[j]);
            // PORTION EXCHANGE BETWEEN PROCESSES: SEND OWN PORTION TO OTHER PROCESSES
            MPI_Iallgatherv(p, proc_portion_size[i], MPI_FLOAT, p, recvcounts, displs, MPI_FLOAT, MPI_COMM_WORLD, allgath_pointer_next_req);
            // PORTION EXCHANGE BETWEEN PROCESSES: RECEIVE PORTION FROM OTHER PROCESSES
            for (i = worker_rank + 1; i < n_workers; i++, j++)
                MPI_Iallgatherv(p, proc_portion_size[i], MPI_FLOAT, p, recvcounts, displs, MPI_FLOAT, MPI_COMM_WORLD, &allgath_recv[j]);
            // INIT BODYFORCE FOR ITERATION i
            for (int i = 0; i < own_portion; i++) { Fx[i] = 0.0f; Fy[i] = 0.0f; Fz[i] = 0.0f; }
            // BODY FORCE COMPUTATION : CALCULATION OWN PORTION
            bodyForceSplit(p, dt, own_portion, start_own_portion, start, end, Fx, Fy, Fz);
            // BODY FORCE COMPUTATION : WAIT OTHER PORTIONS FROM OTHER PROCESSES FOR CONTINUE OWN COMPUTATION
            waitSomeWork(allgath_recv, reqs_ranks_indices, reqs_ranks, n_workers - 1, worker_rank, p, Fx, Fy, Fz, proc_portion_size, proc_portion_start);
            // MASTER PRINT INCOMING RESULTS
            if (worker_rank == MASTER && iter > 1 && print_res == 1)
                printResults(p, nBodies);
            // WAITS FOR THE PREVIOUS REQUEST TO WRITE IN THE SEND BUFFER
            MPI_Test(allgath_pointer_next_req, &send_req_curr, MPI_STATUS_IGNORE);
            MPI_Test(allgath_pointer_prec_req, &send_req_prec, MPI_STATUS_IGNORE);
            if ( (send_req_curr == 0 && send_req_prec == MPI_UNDEFINED) || (send_req_prec == 1 && send_req_curr == 0) )
            {
                prepereSendBuffer(p, p_send, own_portion, start_own_portion);
                integratePositionSplit(p_send, dt, own_portion, start_own_portion, Fx, Fy, Fz);
                // SWAP BUFFER POINTERS FOR NEXT ITERATION
                swapBuffers(p, p_send);
                // SWAP REQEUST POINTER FOR NEXT ITERATION
                swapRequests(allgath_pointer_prec_req, allgath_pointer_next_req);       
            }
            else if (send_req_curr == 0 && send_req_prec == 0)
            {
                MPI_Waitany(2, proc_send_reqs, &proc_send_req_compltd, MPI_STATUS_IGNORE);

                if(proc_send_req_compltd == 0)// la richiesta di inivio corrente NON è stata completata => usare l'altro buffer per l'invio per non aspettare
                {
                    prepereSendBuffer(p, p_send, own_portion, start_own_portion);
                    integratePositionSplit(p_send, dt, own_portion, start_own_portion, Fx, Fy, Fz);
                    // SWAP BUFFER POINTERS FOR NEXT ITERATION
                    swapBuffers(p, p_send);
                    // SWAP REQEUST POINTER FOR NEXT ITERATION
                    swapRequests(allgath_pointer_prec_req, allgath_pointer_next_req);                 
                }
                else // la richiesta di inivio corrente è stata completata => non c'è bisogno di usare l'altro buffer per l'invio
                    integratePositionSplit(p, dt, own_portion, start_own_portion, Fx, Fy, Fz);
            }
            else if( ( send_req_prec == 0 && send_req_curr == 1) || (send_req_curr == 1 && send_req_curr == 1) ) 
            {
                integratePositionSplit(p, dt, own_portion, start_own_portion, Fx, Fy, Fz);
            }
            // END i-th ITERATION
            if (worker_rank == MASTER)
            {
                tElapsed = GetTimer() / 1000.0;
                if (iter > 1)
                    totalTime += tElapsed;
                printf("Iteration %d: %.3f seconds\n", iter, tElapsed);
            }
        }
        // WAIT LAST PENDING REQUEST RESULTS COMPLETITION AND PROCESES SEND ONLY LAST RESULT TO THE MASTER
        MPI_Igatherv(&p[start_own_portion], own_portion * 6, MPI_FLOAT, p_send, recvcounts, displs, MPI_FLOAT, MASTER, MPI_COMM_WORLD, &pending_allgath_req);
        MPI_Wait(&pending_allgath_req, MPI_STATUS_IGNORE);

        if (worker_rank == MASTER && print_res == 1)
            printResults(p_send, nBodies); // stampiamo qui fuori perchè quando usciamo dal ciclo for la igather sta ancora lavorando e quindi qui dobbiamo stampare i risultati

        free(buf_send);
    }
    // END PROCESES WORK AND PRINT TOTAL TIME FROM MASTER
    if (worker_rank == MASTER)
    {
        end_time = MPI_Wtime();
        avgTime = totalTime / (double)(nIters - 1);
        printf("%d Bodies: average %0.3f Billion Interactions / second\n", nBodies, 1e-9 * nBodies * nBodies / avgTime);
        printf("time : %f \n", (end_time - start_time));
    }
    // RESOURCES RELEASE
    free(buf);

    if(worker_rank == MASTER)
    {
        FILE *programs_times;
        if(n_workers < 2)
            programs_times = fopen("programs_times.txt", "w");
        else
            programs_times = fopen("programs_times.txt", "a");

        fprintf(programs_times, "%f\n", (end_time - start_time));
        fclose(programs_times);
    }

    MPI_Finalize();
    return 0;
}