#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <mpi.h>
#include "timer.h"
#include "lib_sequential.h" // MY SEQUENTIAL PROGRAM FUNCTION
#include "lib_parallel.h"  // MY PARALLEL PROGRAM FUNCTIONS

#define MASTER 0
FILE *output_file;
char file_name[50];
const float dt = 0.01f; // time step
const int nIters = 10;  // simulation iterations

int main(int argc, char *argv[])
{
    // INIT TIMERS
    double end_time, start_time;
    // MPI INIT
    MPI_Init(&argc, &argv);
    int worker_rank, n_workers, print_res = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &n_workers);
    MPI_Comm_rank(MPI_COMM_WORLD, &worker_rank);
    // FOR TESTING CORRECTNESS
    if (argc > 2 && strcmp(argv[2], "-t") == 0)
    {
        // DEFINE FILE NAME
        sprintf(file_name, "parallel_%d", n_workers);
        // CREATE FILE FOR PRINT RESULTS
        output_file = fopen(file_name, "w+");
        print_res = 1;
    }
    // INIT INPUT PROGRAM
    int nBodies = 30000;
    if (argc > 1)
        nBodies = atoi(argv[1]);

    MPI_Barrier(MPI_COMM_WORLD);
    if (worker_rank == MASTER)
        start_time = MPI_Wtime();

    if (n_workers == 1) // SEQUENTIAL PROGRAM
    {
        Body p[nBodies];
        determisticInitBodies_sequential(p, nBodies);
        // RESET ACCUMULATORS
        float Fx[nBodies], Fy[nBodies], Fz[nBodies];
        // START WORK
        for (int iter = 1; iter <= nIters - 1; iter++)
        {
            StartTimer();
            // RESET ACCUMULATORS
            for (int i = 0; i < nBodies; i++) { Fx[i] = 0.0f; Fy[i] = 0.0f; Fz[i] = 0.0f; }
            // WORK ON OWN PORTION AND INCOMING PORTION 
            bodyForceSplit_sequential(p, dt, nBodies, 0, 0, nBodies, Fx, Fy, Fz);
            // UPDATE BODIES VALUES
            integratePosition_sequential(p, dt, nBodies, 0, Fx, Fy, Fz);
            // printIterTime(iter, tElapsed, totalTime);
            if (print_res == 1) printResults_sequential(p, nBodies, output_file);
                // printResults_sequential(p, nBodies, fh);
        }
    }
    else // PARALLEL PROGRAM
    {
        // DEFINE BODY TYPES
        MPI_Datatype BODY_POS, BODY_VEL;
        MPI_Datatype strcut_types[1];
        int blockcounts[1];
        MPI_Aint offsets[1], lb, extent;
        offsets[0] = 0;
        strcut_types[0] = MPI_FLOAT;
        blockcounts[0] = 3;

        MPI_Type_create_struct(1, blockcounts, offsets, strcut_types, &BODY_POS);
        MPI_Type_create_struct(1, blockcounts, offsets, strcut_types, &BODY_VEL);
        MPI_Type_commit(&BODY_POS);
        MPI_Type_commit(&BODY_VEL);
        // WORKLOADS DISTRIBUTION CALCULATION
        BodyPosition body_pos[nBodies];
        BodyVelocity body_vel[nBodies];

        int portion = nBodies / n_workers;
        int rest = nBodies % n_workers;
        int portions_sizes[n_workers], portions_starts[n_workers];
        calculatePortions(portions_sizes, portions_starts, n_workers, portion, rest);

        ProcVariables proc;
        proc.rank = worker_rank;
        proc.own_portion = portions_sizes[worker_rank];
        proc.start_own_portion = portions_starts[worker_rank];
        proc.end_own_portion = proc.start_own_portion + proc.own_portion;
        float Fx[proc.own_portion], Fy[proc.own_portion], Fz[proc.own_portion];
        proc.Fx = Fx;  proc.Fy = Fy;  proc.Fz = Fz;
        proc.procs_portions_sizes = portions_sizes;
        proc.procs_portions_starts = portions_starts;
        // INIT OWN BODY PORTION
        determisticInitBodiesSplit(body_pos, body_vel, proc);
        // DEFINE REQUESTS FOR PORTIONS EXCHANGES
        MPI_Request bcast_reqs[n_workers];
        MPI_Request bcast_reqs2[n_workers];
        int reqs_ranks[n_workers];
        int req_done_indexes[n_workers];
        // START i-th PROCES WORK
        for (int iter = 1; iter <= nIters; iter++)
        {
            // PORTION EXCHANGE
            for (int i = 0; i < n_workers; i++)
                MPI_Ibcast(&body_pos[portions_starts[i]], portions_sizes[i], BODY_POS, i, MPI_COMM_WORLD, &bcast_reqs[i]);
            // RESET ACCUMULATORS
            for (int i = 0; i < proc.own_portion; i++) {  Fx[i] = 0.0f;  Fy[i] = 0.0f;  Fz[i] = 0.0f;  }
            // WORK ON OWN PORTION AND INCOMING PORTIONS
            bodyForceSplit(body_pos, dt, proc, proc.start_own_portion, proc.end_own_portion);
            workOnIncomingPortions(n_workers, bcast_reqs, req_done_indexes, proc, body_pos, dt);
            // SLAVES SEND THEIR OWN VELOCITIES AND MASTER RECEIVES THEM
            if (iter > 1)
                MPI_Gatherv(&body_vel[proc.start_own_portion], proc.own_portion, BODY_VEL, &body_vel[proc.start_own_portion], portions_sizes, portions_starts, BODY_VEL, MASTER, MPI_COMM_WORLD);
            // MASTER PRINTS RESOULTS
            if (worker_rank == MASTER && iter > 1 && print_res == 1) 
                printResults(body_pos, body_vel, nBodies, output_file);
                // printResults(body_pos, body_vel, nBodies, fh);
            // UPDATE BODIES VALUES
            integrateVelocitySplit(body_vel, dt, proc);
            integratePositionSplit(body_pos, body_vel, dt, proc);
            // GO TO NEXT ITERATION
        }
        // FREE DATA TYPE
        MPI_Type_free(&BODY_POS);
        MPI_Type_free(&BODY_VEL);
    }
    // END PROCESES WORK AND PRINT TOTAL TIME FROM MASTER
    if (worker_rank == MASTER)
    {
        end_time = MPI_Wtime();
        fprintf(stderr, "time : %f \n\n", (end_time - start_time));
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
