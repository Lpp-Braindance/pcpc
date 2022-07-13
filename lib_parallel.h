#define SOFTENING 1e-9f

typedef struct { float x, y, z; } BodyPosition;
typedef struct { float vx, vy, vz; } BodyVelocity;

typedef struct
{
    int rank, own_portion, start_own_portion, end_own_portion;
    int *procs_portions_sizes, *procs_portions_starts;
    float *Fx, *Fy, *Fz;

} ProcVariables;

void calculatePortions(int procs_portions_sizes[], int procs_portions_starts[], int n_workers, int portion, int rest)
{
    for (int i = 0; i < n_workers; i++)
    {
        procs_portions_sizes[i] = portion;
        if (rest > 0)
        {
            procs_portions_sizes[i]++;
            rest--;
        }
    }

    procs_portions_starts[0] = 0;
    for (int i = 1; i < n_workers; i++)
        procs_portions_starts[i] = procs_portions_starts[i - 1] + procs_portions_sizes[i - 1];
}

void determisticInitBodiesSplit(BodyPosition *body_pos, BodyVelocity *body_vel, ProcVariables proc)
{
    for (int i = proc.start_own_portion; i < proc.end_own_portion; i++)
    {
        body_pos[i].x = i + 1;
        body_pos[i].y = i + 1;
        body_pos[i].z = i + 1;

        body_vel[i].vx = i + 1;
        body_vel[i].vy = i + 1;
        body_vel[i].vz = i + 1;
    }
}

void bodyForceSplit(BodyPosition *body_pos, float dt, ProcVariables proc, int start_new_portion, int end_new_portion)
{
    for (int proc_portion = proc.start_own_portion, own_accumulator = 0; proc_portion < proc.end_own_portion; proc_portion++, own_accumulator++)
    {
        for (int new_portion = start_new_portion; new_portion < end_new_portion; new_portion++)
        {
            float dx = body_pos[new_portion].x - body_pos[proc_portion].x;
            float dy = body_pos[new_portion].y - body_pos[proc_portion].y;
            float dz = body_pos[new_portion].z - body_pos[proc_portion].z;

            float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
            float invDist = 1.0f / sqrt(distSqr);
            float invDist3 = invDist * invDist * invDist;

            proc.Fx[own_accumulator] += dx * invDist3;
            proc.Fy[own_accumulator] += dy * invDist3;
            proc.Fz[own_accumulator] += dz * invDist3;
        }
    }
}

void workOnIncomingPortions(int n_req, MPI_Request bcast_reqs[], int req_done_indexes[], ProcVariables proc, BodyPosition body_pos[], const int dt)
{
    int num_req_compltd, start_new_portion, end_new_portion;
    while (1)
    {
        MPI_Waitsome(n_req, bcast_reqs, &num_req_compltd, req_done_indexes, MPI_STATUSES_IGNORE);
        if (num_req_compltd == MPI_UNDEFINED)
            break;

        for (int i = 0; i < num_req_compltd; i++)
        {
            start_new_portion = proc.procs_portions_starts[req_done_indexes[i]];
            end_new_portion = start_new_portion + proc.procs_portions_sizes[req_done_indexes[i]];
            bodyForceSplit(body_pos, dt, proc, start_new_portion, end_new_portion);
        }
    }
}

// void integratePositionSplit(BodyPosition *body_pos, BodyVelocity *body_vel, float dt, int own_portion, int start_own_portion)
// {
//     for (int i = 0; i < own_portion; i++)
//     {
//         body_pos[start_own_portion + i].x += body_vel[start_own_portion + i].vx * dt;
//         body_pos[start_own_portion + i].y += body_vel[start_own_portion + i].vy * dt;
//         body_pos[start_own_portion + i].z += body_vel[start_own_portion + i].vz * dt;
//     }
// }

// void integrateVelocitySplit(BodyVelocity *body_vel, float dt, int own_portion, int start_own_portion, float Fx[], float Fy[], float Fz[])
// {
//     for (int i = 0; i < own_portion; i++)
//     {
//         body_vel[start_own_portion + i].vx += dt * Fx[i];
//         body_vel[start_own_portion + i].vy += dt * Fy[i];
//         body_vel[start_own_portion + i].vz += dt * Fz[i];
//     }
// }



void integratePositionSplit(BodyPosition *body_pos, BodyVelocity *body_vel, float dt, ProcVariables proc)
{
    for (int i = proc.start_own_portion; i < proc.end_own_portion; i++)
    {
        body_pos[i].x += body_vel[i].vx * dt;
        body_pos[i].y += body_vel[i].vy * dt;
        body_pos[i].z += body_vel[i].vz * dt;
    }
}

void integrateVelocitySplit(BodyVelocity *body_vel, float dt, ProcVariables proc)
{
    for (int i = proc.start_own_portion, own_accumulator = 0; i < proc.end_own_portion; i++, own_accumulator++)
    {
        body_vel[i].vx += dt * proc.Fx[own_accumulator];
        body_vel[i].vy += dt * proc.Fy[own_accumulator];
        body_vel[i].vz += dt * proc.Fz[own_accumulator];
    }
}


void printResults(BodyPosition body_pos[], BodyVelocity body_vel[], int nBodies , FILE *output_file)
{
    for (int i = 0; i < nBodies; i++)
    {
        char line[100];
        fprintf(output_file, "%f %f %f %f %f %f\n", body_pos[i].x, body_pos[i].y, body_pos[i].z, body_vel[i].vx, body_vel[i].vy, body_vel[i].vz);
    }
}

