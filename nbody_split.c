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

typedef struct { float x, y, z; } BodyPosition;
typedef struct { float vx, vy, vz; } BodyVelocity;
typedef struct 
{ 
    int rank, own_portion, start_own_portion;
    int *procs_portions_sizes, *procs_portions_starts;
    float *Fx, *Fy, *Fz;

} ProcVariables;

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

void bodyForceSplit(BodyPosition *body_pos, float dt, ProcVariables proc, int start, int end)
{
    for (int i = 0; i < proc.own_portion; i++)
    {
        for (int j = start; j < end; j++) // start = inizio nuova porzione da confrontare end = fine porzione da confrontare
        {
            float dx = body_pos[j].x - body_pos[proc.start_own_portion + i].x;
            float dy = body_pos[j].y - body_pos[proc.start_own_portion + i].y;
            float dz = body_pos[j].z - body_pos[proc.start_own_portion + i].z;

            float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
            float invDist = 1.0f / sqrt(distSqr);
            float invDist3 = invDist * invDist * invDist;

            proc.Fx[i] += dx * invDist3;
            proc.Fy[i] += dy * invDist3;
            proc.Fz[i] += dz * invDist3;
        }
    }
}

void waitSomeWork(MPI_Request allgath_reqs[], int reqs_ranks[], int n_reqs, BodyPosition *body_pos, ProcVariables proc)
{
    const float dt = 0.01f;
    int count = n_reqs;
    int start, end;
    int index, req_compl;
    while (count > 0)
    {
        MPI_Waitany(n_reqs, allgath_reqs, &index, MPI_STATUS_IGNORE);
        req_compl = reqs_ranks[index];
        start = proc.procs_portions_starts[req_compl];
        end = start + proc.procs_portions_sizes[req_compl];
        bodyForceSplit(body_pos, dt, proc, start, end);
        count--;
    }
}

void integratePositionSplit(BodyPosition *body_pos, BodyVelocity *body_vel, float dt, int own_portion, int start_own_portion)
{
    for (int i = 0; i < own_portion; i++)
    {
        body_pos[start_own_portion + i].x += body_vel[start_own_portion + i].vx * dt;
        body_pos[start_own_portion + i].y += body_vel[start_own_portion + i].vy * dt;
        body_pos[start_own_portion + i].z += body_vel[start_own_portion + i].vz * dt;
    }
}

void integrateVelocitySplit(BodyVelocity *body_vel, float dt, int own_portion, int start_own_portion, float Fx[], float Fy[], float Fz[])
{
    for (int i = 0; i < own_portion; i++)
    {
        body_vel[start_own_portion + i].vx += dt * Fx[i];
        body_vel[start_own_portion + i].vy += dt * Fy[i];
        body_vel[start_own_portion + i].vz += dt * Fz[i];
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

void calculateAllgathervPortions(int recvcounts[], int displs[], int n_workers, int procs_portions_sizes[], int procs_portions_starts[])
{
    for (int i = 0; i < n_workers; i++)
    {
        recvcounts[i] = procs_portions_sizes[i] * 3;
        displs[i] = procs_portions_starts[i] * 3;
    }
}

void initRequestsRanks(int reqs_ranks[], int worker_rank, int n_workers)
{
    int i, j;

    for (i = 0, j = 0; i < worker_rank; i++, j++)
        reqs_ranks[j] = i;

    for (i = worker_rank + 1; i < n_workers; i++, j++)
        reqs_ranks[j] = i;
}

// SEQUENTIAL PROGRAM FUNCTION

typedef struct { float x, y, z, vx, vy, vz; } Body;

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
        // RESOURCES RELEASE
        free(buf);
    }
    else // PARALLEL PROGRAM
    {
        // CALCULATE WORKLOAD DISTRIBUCTION
        int portion = nBodies / n_workers;
        int rest = nBodies % n_workers;
        int procs_portions_sizes[n_workers];
        int procs_portions_starts[n_workers];
        calculatePortions(procs_portions_sizes, procs_portions_starts, n_workers, portion, rest);
        // DEFINE AND INIT BODIES BUFFERS
        int bytes = nBodies * sizeof(BodyPosition);
        float *bodies_positions = (float *)malloc(bytes);
        BodyPosition *body_pos = (BodyPosition *)bodies_positions;

        bytes = nBodies * sizeof(BodyVelocity);
        float *bodies_velocities = (float *)malloc(bytes);
        BodyVelocity *body_vel = (BodyVelocity *)bodies_velocities;

        int displs[n_workers];
        int recvcounts[n_workers];
        calculateAllgathervPortions(recvcounts, displs, n_workers, procs_portions_sizes, procs_portions_starts);
        // DEFINE PROCES VARIABLES
        ProcVariables proc;
        proc.rank = worker_rank;
        proc.own_portion = procs_portions_sizes[worker_rank];
        proc.start_own_portion = procs_portions_starts[worker_rank];
        float Fx[proc.own_portion];
        float Fy[proc.own_portion];
        float Fz[proc.own_portion];
        proc.Fx = Fx;  proc.Fy = Fy;  proc.Fz = Fz;
        proc.procs_portions_sizes = procs_portions_sizes;
        proc.procs_portions_starts = procs_portions_starts;
        int start = proc.start_own_portion;
        int end = start + proc.own_portion;
        // INIT OWN BODY PORTION
        determisticInitBodiesSplit(body_pos, body_vel, proc.own_portion, proc.start_own_portion);
        // DEFINE REQUESTS FOR PORTIONS EXCHANGES
        MPI_Request allgath_reqs[n_workers];
        MPI_Request allgath_send_req = MPI_REQUEST_NULL;
        MPI_Request gath_vels_reqs = MPI_REQUEST_NULL;
        int reqs_ranks[n_workers - 1];
        initRequestsRanks(reqs_ranks, worker_rank, n_workers);
        int i, j;
        // START i-th PROCES WORK
        for (int iter = 1; iter <= nIters; iter++)
        {
            if (worker_rank == MASTER)
                StartTimer();
            // PORTION EXCHANGE BETWEEN PROCESSES: RECEIVE PORTION FROM OTHER PROCESSES
            for (i = 0, j = 0; i < worker_rank; i++, j++)
                MPI_Iallgatherv(body_pos, proc.procs_portions_sizes[i] * 3, MPI_FLOAT, body_pos, recvcounts, displs, MPI_FLOAT, MPI_COMM_WORLD, &allgath_reqs[j]);
            // PORTION EXCHANGE BETWEEN PROCESSES: SEND OWN PORTION TO OTHER PROCESSES
            MPI_Iallgatherv(body_pos, proc.procs_portions_sizes[i] * 3, MPI_FLOAT, body_pos, recvcounts, displs, MPI_FLOAT, MPI_COMM_WORLD, &allgath_send_req);
            // PORTION EXCHANGE BETWEEN PROCESSES: RECEIVE PORTION FROM OTHER PROCESSES
            for (i = worker_rank + 1; i < n_workers; i++, j++)
                MPI_Iallgatherv(body_pos, proc.procs_portions_sizes[i] * 3, MPI_FLOAT, body_pos, recvcounts, displs, MPI_FLOAT, MPI_COMM_WORLD, &allgath_reqs[j]);
            // INIT INTERMEDIATE VALUES FOR BODYFORCE CALCULATION
            for (int i = 0; i < proc.own_portion; i++) {  proc.Fx[i] = 0.0f;  proc.Fy[i] = 0.0f;  proc.Fz[i] = 0.0f;  }
            bodyForceSplit(body_pos, dt, proc, start, end);
            waitSomeWork(allgath_reqs, reqs_ranks, n_workers - 1, body_pos, proc);
            // WAITS THE COMPLITION OF PREVIOUS REQUEST 
            MPI_Wait(&gath_vels_reqs, MPI_STATUS_IGNORE);
            // MASTER PRINT RESULTS RECEIVED
            if (worker_rank == MASTER && iter > 1 && print_res == 1)
                printPositionsAndVelocitiesResults(body_pos, body_vel, nBodies);
            // UPDATE BODIES VELOCITY AND SEND THEM
            integrateVelocitySplit(body_vel, dt, proc.own_portion, proc.start_own_portion, Fx, Fy, Fz);
            MPI_Igatherv(&body_vel[proc.start_own_portion], proc.own_portion * 3, MPI_FLOAT, body_vel, recvcounts, displs, MPI_FLOAT, MASTER, MPI_COMM_WORLD, &gath_vels_reqs);
            // WAIT SEND REQ COMPLITION BEFORE UPDATE BODIES POSITIONS 
            MPI_Wait(&allgath_send_req, MPI_STATUS_IGNORE);
            integratePositionSplit(body_pos, body_vel, dt, proc.own_portion, proc.start_own_portion);
            // END i-th ITERATION
            if (worker_rank == MASTER)
                printIterTime(iter, tElapsed, totalTime);
        }
        // SEND LAST RESULTS ONLY TO THE MASTER
        MPI_Gatherv(&body_pos[proc.start_own_portion], proc.own_portion * 3, MPI_FLOAT, body_pos, recvcounts, displs, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
        MPI_Gatherv(&body_vel[proc.start_own_portion], proc.own_portion * 3, MPI_FLOAT, body_vel, recvcounts, displs, MPI_FLOAT, MASTER, MPI_COMM_WORLD);

        if (worker_rank == MASTER && print_res == 1)
            printPositionsAndVelocitiesResults(body_pos, body_vel, nBodies);
        
        // RESOURCES RELEASE
        free(bodies_positions);
        free(bodies_velocities);
    }
    // END PROCESES WORK AND PRINT TOTAL TIME FROM MASTER
    if (worker_rank == MASTER)
    {
        end_time = MPI_Wtime();
        avgTime = totalTime / (double)(nIters - 1);
        printf("%d Bodies: average %0.3f Billion Interactions / second\n", nBodies, 1e-9 * nBodies * nBodies / avgTime);
        printf("time : %f \n\n", (end_time - start_time));
    }
    
    MPI_Finalize();
    return 0;
}