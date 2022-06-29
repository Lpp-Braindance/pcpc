#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "timer.h"
#include <mpi.h>

#define SOFTENING 1e-9f
#define MASTER 0

FILE *output_file;
char file_name[50];

// PARALLEL PROGRAM FUNCTIONS

typedef struct
{
    float x, y, z;
} BodyPosition;
typedef struct
{
    float vx, vy, vz;
} BodyVelocity;

void determisticInitBodiesSplit(BodyPosition *body_pos, BodyVelocity *body_vel, int own_portion, int start_own_portion)
{
    for (int i = 0; i < own_portion; i++)
    {
        body_pos[start_own_portion + i].x = start_own_portion + i + 1;
        body_pos[start_own_portion + i].y = start_own_portion + i + 1;
        body_pos[start_own_portion + i].z = start_own_portion + i + 1;

        body_vel[start_own_portion + i].vx = start_own_portion + i + 1;
        body_vel[start_own_portion + i].vy = start_own_portion + i + 1;
        body_vel[start_own_portion + i].vz = start_own_portion + i + 1;
    }
}

void bodyForceSplit(BodyPosition *p, float dt, int own_portion, int start_own_portion, int start, int end, float Fx[], float Fy[], float Fz[])
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

void waitSomeWork(MPI_Request bcast_recv[], int requests_ranks[], int n_recv_req, int rank, BodyPosition *p, float Fx[], float Fy[], float Fz[], int proc_portion_size[], int proc_portion_start[])
{
    int ready_req = 0;
    int count = n_recv_req;
    int start, end;
    int req_compl;
    const float dt = 0.01f;
    int own_portion = proc_portion_size[rank];
    int start_own_portion = proc_portion_start[rank];

    int index;
    while (count > 0)
    {
        MPI_Waitany(n_recv_req, bcast_recv, &index, MPI_STATUS_IGNORE);
        req_compl = requests_ranks[index];
        start = proc_portion_start[req_compl];
        end = start + proc_portion_size[req_compl];
        bodyForceSplit(p, dt, own_portion, start_own_portion, start, end, Fx, Fy, Fz);
        count--;
    }
}

void integratePositionSplit(BodyPosition *p, float dt, int own_portion, int start_own_portion, BodyVelocity *p_vel)
{
    for (int i = 0; i < own_portion; i++)
    {
        p[start_own_portion + i].x += p_vel[start_own_portion + i].vx * dt;
        p[start_own_portion + i].y += p_vel[start_own_portion + i].vy * dt;
        p[start_own_portion + i].z += p_vel[start_own_portion + i].vz * dt;
    }
}

void integrateVelocitySplit(BodyVelocity *p, float dt, int own_portion, int start_own_portion, float Fx[], float Fy[], float Fz[])
{
    for (int i = 0; i < own_portion; i++)
    {
        p[start_own_portion + i].vx += dt * Fx[i];
        p[start_own_portion + i].vy += dt * Fy[i];
        p[start_own_portion + i].vz += dt * Fz[i];
    }
}

void calculatePortions(int proc_portion_size[], int proc_portion_start[], int n_workers, int portion, int rest)
{
    proc_portion_size[0] = portion;
    proc_portion_start[0] = 0;

    int i, j;

    if (rest > 0)
    {
        int i, j;
        for (i = 1, j = 0; j < rest; i++, j++)
        {
            proc_portion_size[i] = portion + 1;
            proc_portion_start[i] = portion * i + j;
        }
        for (; i < n_workers; i++, j++)
        {
            proc_portion_size[i] = portion;
            proc_portion_start[i] = portion * i + rest;
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
        recvcounts[i] = proc_portion_size[i] * 3;
        displs[i] = proc_portion_start[i] * 3;
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

// SEQUENTIAL PROGRAM FUNCTION

typedef struct
{
    float x, y, z, vx, vy, vz;
} Body;

void determisticInitBodies(float *buf, int n)
{
    Body *p = (Body *)buf;
    for (int i = 0; i < n; i++)
    {
        p[i].x = i + 1;
        p[i].y = i + 1;
        p[i].z = i + 1;
        p[i].vx = i + 1;
        p[i].vy = i + 1;
        p[i].vz = i + 1;
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

// COMMON FUNCTION

void printIterTime(int iter, double tElapsed, double totalTime)
{
    tElapsed = GetTimer() / 1000.0;
    if (iter > 1)
        totalTime += tElapsed;

    printf("Iteration %d: %.3f seconds\n", iter, tElapsed);
}

// TESTING FUNCTION

void printResults(Body *p, int nBodies) // FOR SEQUENTIAL PROGRAM
{
    for (int i = 0; i < nBodies; i++)
        fprintf(output_file, "%0.3f %0.3f %0.3f %0.3f %0.3f %0.3f\n", p[i].x, p[i].y, p[i].z, p[i].vx, p[i].vy, p[i].vz);

    fprintf(output_file, "\n");
}

void printPositionsAndVelocitiesResults(BodyPosition *body_pos, BodyVelocity *body_vel, int nBodies) // FRO PARALLEL PROGRAM
{
    for (int i = 0; i < nBodies; i++)
        fprintf(output_file, "%0.3f %0.3f %0.3f %0.3f %0.3f %0.3f\n", body_pos[i].x, body_pos[i].y, body_pos[i].z, body_vel[i].vx, body_vel[i].vy, body_vel[i].vz);

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
    int worker_rank, n_workers, print_res = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &worker_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_workers);
    // FOR TESTING CORRECTNESS
    if (argc > 2 && strcmp(argv[2], "-t") == 0)
    {
        sprintf(file_name, "parallel_%d", n_workers);
        output_file = fopen(file_name, "w+");
        print_res = 1;
    }
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

    if (n_workers == 1) // SEQUENTIAL PROGRAM
    {
        int bytes = nBodies * sizeof(Body);
        float *buf = (float *)malloc(bytes);
        Body *p = (Body *)buf;

        determisticInitBodies(buf, nBodies);
        for (int iter = 1; iter <= nIters; iter++)
        {
            StartTimer();
            bodyForce(p, dt, nBodies);
            integratePosition(p, dt, nBodies);
            printIterTime(iter, tElapsed, totalTime);
            if (print_res == 1)
                printResults(p, nBodies);
        }
    }
    else // PARALLEL PROGRAM
    {
        // CALCULATE WORKLOAD DISTRIBUCTION
        int portion = nBodies / n_workers;
        int rest = nBodies % n_workers;
        int proc_portion_size[n_workers];
        int proc_portion_start[n_workers];
        calculatePortions(proc_portion_size, proc_portion_start, n_workers, portion, rest);

        int bytes = nBodies * sizeof(BodyPosition);
        float *bodies_positions = (float *)malloc(bytes);
        BodyPosition *body_pos = (BodyPosition *)bodies_positions;

        bytes = nBodies * sizeof(BodyPosition);
        float *bodies_positions2 = (float *)malloc(bytes);
        BodyPosition *p_send = (BodyPosition *)bodies_positions2;

        bytes = nBodies * sizeof(BodyVelocity);
        float *bodies_velocites = (float *)malloc(bytes);
        BodyVelocity *body_vel = (BodyVelocity *)bodies_velocites;

        int displs[n_workers];
        int recvcounts[n_workers];
        calculateAllgathervPortions(recvcounts, displs, n_workers, proc_portion_size, proc_portion_start);
        // DEFINE VARIABLES FOR WORKLOAD INTERVALS
        int own_portion = proc_portion_size[worker_rank];
        int start_own_portion = proc_portion_start[worker_rank];
        int start = start_own_portion;
        int end = start + own_portion;
        float Fx[own_portion];
        float Fy[own_portion];
        float Fz[own_portion];
        // INIT OWN BODY PORTION
        determisticInitBodiesSplit(body_pos, body_vel, own_portion, start_own_portion);
        // DEFINE REQUESTS RANKS FOR MPI_WAIT SOME
        int reqs_ranks[n_workers - 1];
        initRequestsRanks(reqs_ranks, worker_rank, n_workers);
        // DEFINE REQUESTS FOR PORTIONS EXCHANGES
        MPI_Request bcast_recv[n_workers];
        MPI_Request bcast_send_next = MPI_REQUEST_NULL;
        MPI_Request bcast_send_prec = MPI_REQUEST_NULL;

        // DEFINE POINTER FOR PREC AND NEXT NONBLOCKING COMMUNICATION REQUESTS
        MPI_Request *bcast_pointer_next_req = &bcast_send_next;
        MPI_Request *bcast_pointer_prec_req = &bcast_send_prec;
        MPI_Request *tmp_bcast_send_req_swap = &bcast_send_prec;

        // DEFINE TMP VARIABLES
        int i, j = 0, bcast_send_done;
        BodyPosition *tmp_buf_swap;
        MPI_Request gather_vels_req = MPI_REQUEST_NULL;
        // START i-th PROCES WORK
        for (int iter = 1; iter <= nIters; iter++)
        {
            if (worker_rank == MASTER)
                StartTimer();

            for (i = 0, j = 0; i < worker_rank; i++, j++)
                MPI_Iallgatherv(body_pos, proc_portion_size[i]*3, MPI_FLOAT, body_pos, recvcounts, displs, MPI_FLOAT, MPI_COMM_WORLD, &bcast_recv[j]);

            MPI_Iallgatherv(body_pos, proc_portion_size[i]*3, MPI_FLOAT, body_pos, recvcounts, displs, MPI_FLOAT, MPI_COMM_WORLD, bcast_pointer_next_req);
            
            for (i = worker_rank + 1; i < n_workers; i++, j++)
                MPI_Iallgatherv(body_pos, proc_portion_size[i]*3, MPI_FLOAT, body_pos, recvcounts, displs, MPI_FLOAT, MPI_COMM_WORLD, &bcast_recv[j]);

            for (int i = 0; i < own_portion; i++) { Fx[i] = 0.0f; Fy[i] = 0.0f; Fz[i] = 0.0f; }

            bodyForceSplit(body_pos, dt, own_portion, start_own_portion, start, end, Fx, Fy, Fz);
            waitSomeWork(bcast_recv, reqs_ranks, n_workers - 1, worker_rank, body_pos, Fx, Fy, Fz, proc_portion_size, proc_portion_start);
            // MASTER PRINT INCOMING RESULTS
            
            MPI_Wait(&gather_vels_req, MPI_STATUS_IGNORE);
            if (worker_rank == MASTER && iter > 1 && print_res == 1)
                printPositionsAndVelocitiesResults(body_pos, body_vel, nBodies);

            integrateVelocitySplit(body_vel, dt, own_portion, start_own_portion, Fx, Fy, Fz);
            MPI_Igatherv(&body_vel[start_own_portion], own_portion * 3, MPI_FLOAT, body_vel, recvcounts, displs, MPI_FLOAT, MASTER, MPI_COMM_WORLD, &gather_vels_req);

            MPI_Wait(bcast_pointer_next_req, MPI_STATUS_IGNORE);
            integratePositionSplit(body_pos, dt, own_portion, start_own_portion, body_vel);


            // MPI_Test(bcast_pointer_next_req, &bcast_send_done, MPI_STATUS_IGNORE);
            // if (bcast_send_done == 0)
            // {
            //     MPI_Wait(bcast_pointer_prec_req, MPI_STATUS_IGNORE);
                
            //     for (int i = 0; i < own_portion; i++)
            //     {
            //         p_send[start_own_portion + i].x = body_pos[start_own_portion + i].x;
            //         p_send[start_own_portion + i].y = body_pos[start_own_portion + i].y;
            //         p_send[start_own_portion + i].z = body_pos[start_own_portion + i].z;
            //     }
            //     // CALCULATION OF RESULT
            //     integratePositionSplit(p_send, dt, own_portion, start_own_portion,  body_vel);
            //     // SWAP BUFFER POINTERS FOR NEXT ITERATION
            //     tmp_buf_swap = body_pos;
            //     body_pos = p_send;
            //     p_send = tmp_buf_swap;
            //     // SWAP REQEUST POINTER FOR NEXT ITERATION
            //     tmp_bcast_send_req_swap = bcast_pointer_prec_req;
            //     bcast_pointer_prec_req = bcast_pointer_next_req;
            //     bcast_pointer_next_req = tmp_bcast_send_req_swap;
            // }
            // else
            //     integratePositionSplit(body_pos, dt, own_portion, start_own_portion, body_vel);

            // END i-th ITERATION
            if (worker_rank == MASTER)
                printIterTime(iter, tElapsed, totalTime);
        }
        // WAIT LAST PENDING REQUEST RESULTS COMPLETITION AND PROCESES SEND ONLY LAST RESULT TO THE MASTER
        MPI_Gatherv(&body_pos[start_own_portion], own_portion * 3, MPI_FLOAT, body_pos, recvcounts, displs, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
        MPI_Gatherv(&body_vel[start_own_portion], own_portion * 3, MPI_FLOAT, body_vel, recvcounts, displs, MPI_FLOAT, MASTER, MPI_COMM_WORLD);

        if (worker_rank == MASTER && print_res == 1)
            printPositionsAndVelocitiesResults(body_pos, body_vel, nBodies);
    }
    // END PROCESES WORK AND PRINT TOTAL TIME FROM MASTER
    if (worker_rank == MASTER)
    {
        end_time = MPI_Wtime();
        avgTime = totalTime / (double)(nIters - 1);
        printf("%d Bodies: average %0.3f Billion Interactions / second\n", nBodies, 1e-9 * nBodies * nBodies / avgTime);
        printf("time : %f \n\n", (end_time - start_time));
    }
    // RESOURCES RELEASE

    MPI_Finalize();
    return 0;
}