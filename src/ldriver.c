#include "stepper.h"
#include "shallow2d.h"

#define _OPENMP 1

#ifdef _OPENMP
#include <omp.h>
#elif defined SYSTIME
#include <sys/time.h>
#endif

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

//ldoc on
/**
 * # Driver code
 *
 * The driver code is where we put together the time stepper and
 * the physics routines to actually solve the equations and make
 * pretty pictures of the solutions.
 *
 * ## Diagnostics
 *
 * The numerical method is supposed to preserve (up to rounding
 * errors) the total volume of water in the domain and the total
 * momentum.  Ideally, we should also not see negative water heights,
 * since that will cause the system of equations to blow up.  For
 * debugging convenience, we'll plan to periodically print diagnostic
 * information about these conserved quantities (and about the range
 * of water heights).
 */

void solution_check(central2d_t* sim)
{
    int nx = sim->nx, ny = sim->ny;
    float* u = sim->u;
    float h_sum = 0, hu_sum = 0, hv_sum = 0;
    float hmin = u[central2d_offset(sim,0,0,0)];
    float hmax = hmin;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            float h = u[central2d_offset(sim,0,i,j)];
            h_sum += h;
            hu_sum += u[central2d_offset(sim,1,i,j)];
            hv_sum += u[central2d_offset(sim,2,i,j)];
            hmax = fmaxf(h, hmax);
            hmin = fminf(h, hmin);
        }
    float cell_area = sim->dx * sim->dy;
    h_sum *= cell_area;
    hu_sum *= cell_area;
    hv_sum *= cell_area;
    printf("-\n  Volume: %g\n  Momentum: (%g, %g)\n  Range: [%g, %g]\n",
           h_sum, hu_sum, hv_sum, hmin, hmax);
    assert(hmin > 0);
}

void solution_check_mpi(central2d_mpi_t* sim)
{
    int nx = sim->nx, ny = sim->ny;
    float* u = sim->u;
    float h_sum = 0, hu_sum = 0, hv_sum = 0;
    float hmin = u[central2d_mpi_offset(sim,0,0,0)];
    float hmax = hmin;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            float h = u[central2d_mpi_offset(sim,0,i,j)];
            h_sum += h;
            hu_sum += u[central2d_mpi_offset(sim,1,i,j)];
            hv_sum += u[central2d_mpi_offset(sim,2,i,j)];
            hmax = fmaxf(h, hmax);
            hmin = fminf(h, hmin);
        }
    float cell_area = sim->dx * sim->dy;
    h_sum *= cell_area;
    hu_sum *= cell_area;
    hv_sum *= cell_area;
    printf("-\n  Volume: %g\n  Momentum: (%g, %g)\n  Range: [%g, %g]\n",
           h_sum, hu_sum, hv_sum, hmin, hmax);
    assert(hmin > 0);
}

/**
 * ## I/O
 *
 * After finishing a run (or every several steps), we might want to
 * write out a data file for further processing by some other program
 * -- in this case, a Python visualizer.  The visualizer takes the
 * number of pixels in x and y in the first two entries, then raw
 * single-precision raster pictures.
 */

FILE* viz_open(const char* fname, central2d_t* sim, int vskip)
{
    FILE* fp = fopen(fname, "w");
    if (fp) {
        float xy[2] = {sim->nx/vskip, sim->ny/vskip};
        fwrite(xy, sizeof(float), 2, fp);
    }
    return fp;
}

void viz_close(FILE* fp)
{
    fclose(fp);
}

void viz_frame(FILE* fp, central2d_t* sim, int vskip)
{
    if (!fp)
        return;
    for (int iy = 0; iy < sim->ny; iy += vskip)
        for (int ix = 0; ix < sim->nx; ix += vskip)
            fwrite(sim->u + central2d_offset(sim,0,ix,iy),
                   sizeof(float), 1, fp);
}

/**
 * ## Lua driver routines
 *
 * A better way to manage simulation parameters is by a scripting
 * language.  Python is a popular choice, but I prefer Lua for many
 * things (not least because it is an easy build).  It's also quite
 * cheap to call a Lua function for every point in a mesh
 * (less so for Python, though it probably won't make much difference).
 *
 * ### Lua callback functions
 *
 * We specify the initial conditions by providing the simulator
 * with a callback function to be called at each cell center.
 * The callback function is assumed to be the `init` field of
 * a table at index 1.
 */

void lua_init_sim(lua_State* L, central2d_t* sim)
{
    lua_getfield(L, 1, "init");
    if (lua_type(L, -1) != LUA_TFUNCTION)
        luaL_error(L, "Expected init to be a string");

    int nx = sim->nx, ny = sim->ny, nfield = sim->nfield;
    float dx = sim->dx, dy = sim->dy;
    float* u = sim->u;

    for (int ix = 0; ix < nx; ++ix) {
        float x = (ix + 0.5) * dx;
        for (int iy = 0; iy < ny; ++iy) {
            float y = (iy + 0.5) * dy;
            lua_pushvalue(L, -1);
            lua_pushnumber(L, x);
            lua_pushnumber(L, y);
            lua_call(L, 2, nfield);
            for (int k = 0; k < nfield; ++k)
                u[central2d_offset(sim,k,ix,iy)] = lua_tonumber(L, k-nfield);
            lua_pop(L, nfield);
        }
    }

    lua_pop(L,1);
}


/**
 * ### Running the simulation
 *
 * The `run_sim` function looks a lot like the main routine of the
 * "ordinary" command line driver.  We specify the initial conditions
 * by providing the simulator with a callback function to be called at
 * each cell center.  Note that we have two different options for
 * timing the steps -- we can use the OpenMP timing routines
 * (preferable if OpenMP is available) or the POSIX `gettimeofday`
 * if the `SYSTIME` macro is defined.  If there's no OpenMP and
 * `SYSTIME` is undefined, we fall back to just printing the number
 * of steps without timing information.
 */

int run_sim(lua_State* L)
{
    // Get rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    printf("MPI worldrank: %d\n", world_rank);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // root processes print number of processes
    if (world_rank==0){
        // Get number of processes      
        printf("MPI worldsize: %d\n", world_size);
    };
    
    int n = lua_gettop(L);
    if (n != 1 || !lua_istable(L, 1))
        luaL_error(L, "Argument must be a table");

    lua_getfield(L, 1, "w");
    lua_getfield(L, 1, "h");
    lua_getfield(L, 1, "cfl");
    lua_getfield(L, 1, "ftime");
    lua_getfield(L, 1, "nx");
    lua_getfield(L, 1, "ny");
    lua_getfield(L, 1, "vskip");
    lua_getfield(L, 1, "frames");
    lua_getfield(L, 1, "batch");
    lua_getfield(L, 1, "block_n");
    lua_getfield(L, 1, "out");

    double w     = luaL_optnumber(L, 2, 2.0);
    double h     = luaL_optnumber(L, 3, w);
    double cfl   = luaL_optnumber(L, 4, 0.45);
    double ftime = luaL_optnumber(L, 5, 0.01);
    int nx_global = luaL_optinteger(L, 6, 200);
    int ny_global = luaL_optinteger(L, 7, nx_global);
    int vskip    = luaL_optinteger(L, 8, 1);
    int frames   = luaL_optinteger(L, 9, 50);
    int batch   = luaL_optinteger(L, 10, 1);
    int block_n   = luaL_optinteger(L, 11, 1);
    const char* fname = luaL_optstring(L, 12, "sim.out");
    lua_pop(L, 9);

    // Number of ghost cell
    int ng = 4 + 2 * (batch-1);
    central2d_t* sim_global = central2d_init(w, h, nx_global, ny_global, 3, shallow2d_flux, 
					     shallow2d_speed, cfl, ng);

    // Partition the x axis
    int npartx;
    int* offsets_x = alloc_partition(sim_global->nx, sim_global->ng, block_n, &npartx);

    // Partition the y axis
    int nparty;
    int* offsets_y = alloc_partition(sim_global->ny, sim_global->ng, block_n, &nparty);

    if (world_rank == 0){
        printf("batch size = %d, block size = %d \n", batch, block_n);
        printf("offsets_x: \n");
        for (int i = 0; i <= npartx; ++i)
            printf("%d, ", offsets_x[i]);
        printf("\n");
        
        printf("offsets_y: \n");
        for (int i = 0; i <= nparty; ++i)
            printf("%d, ", offsets_y[i]);
        printf("\n");
    }

    // set up sim_u_all to collect all u blocks from all processes
    // also set up relevant index info
    // all blocks of u matrices [0] [1]
    //                          [2] [3]
    int blockuindex[npartx*nparty];
    int blockunum[npartx*nparty];
    int elindex = 0;
    float* sim_u_all;
    int myel = 0;
    int nx, ny;
    blockuindex[0] = elindex;
    for (int j = 0; j < nparty; ++j){
        for (int i = 0; i < npartx; ++i){
            nx = offsets_x[i+1] - offsets_x[i];
            ny = offsets_y[j+1] - offsets_y[j];
            
            int nx_all = nx + 2*ng;
            int ny_all = ny + 2*ng;
            int nc = nx_all * ny_all;
            int N  = sim_global->nfield * nc;
            int blockindex = i + npartx * j;
            int curel = 4*N + 6*nx_all;
            if (world_rank == i + npartx * j){
                myel = curel;
            }
            elindex += curel;
            if (blockindex < npartx*nparty - 1){
                blockuindex[blockindex+1] = elindex;
            }
	        blockunum[blockindex] = curel;
        }
    }
    
    if (world_rank==0){
        // root process creates a matrix to store all the u elements
        sim_u_all = malloc(elindex * sizeof(float));
    }
    else   
        sim_u_all = NULL;

    // set up done_block that collects done information from all processes
    bool done_block[npartx*nparty];
    memset(done_block, false, npartx*nparty*sizeof(bool));

    // Create sim_local for each process
    central2d_mpi_t* sim_local = central2d_mpi_init(w/nx_global, h/ny_global, 3, shallow2d_flux, 
                                shallow2d_speed, cfl, ng, offsets_x, offsets_y, npartx, nparty);

    // Initial and periodic boundary condition of sim_global
    lua_init_sim(L, sim_global);
    central2d_periodic(sim_global->u, sim_global->nx, sim_global->ny, sim_global->ng, 3);

    // copy sim_global to sim_local for each process
    sub_mpi_copyin(sim_local, sim_global, offsets_x[sim_local->blocki], offsets_x[sim_local->blocki+1],
	           offsets_y[sim_local->blockj], offsets_y[sim_local->blockj+1]);

    FILE* viz;
    double tcompute;
    if (world_rank == 0){
        printf("%g %g %d %d %g %d %g\n", w, h, nx_global, ny_global, cfl, frames, ftime);
        viz = viz_open(fname, sim_global, vskip);
        solution_check(sim_global);
        viz_frame(viz, sim_global, vskip);
        tcompute = 0;
    }
    
    for (int i = 0; i < frames; ++i) {
#ifdef _OPENMP
        double t0 = omp_get_wtime();

        // Initial nstep and t arrays of all subdomains
        int nstep[npartx*nparty];
        float t[npartx*nparty];
        bool done = false;
        memset(nstep, 0, npartx*nparty*sizeof(int));
        memset(t, 0.0, npartx*nparty*sizeof(float));

        bool done_all = false;

        // Run
        while (!done_all){
            done_all = true;
            // run one batch time for each process
            central2d_mpi_batch_run(sim_local, ftime, batch, \
                    &nstep[sim_local->blocki + npartx * sim_local->blockj], \
                    &t[sim_local->blocki + npartx * sim_local->blockj], &done);

            //send and receive buffers left and right
            sub_start_sendrecv(sim_local, 0);
            sub_end_sendrecv(sim_local, 0);
            //send and receive buffers up and down
            sub_start_sendrecv(sim_local, 1);
            sub_end_sendrecv(sim_local, 1);

            // root process gather each done into done_block
            MPI_Gather(&done, 1, MPI_C_BOOL, done_block, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

            if (world_rank == 0){
                //check if is all done
                for (int index = 0; index < npartx*nparty; ++index)
                    if (done_block[index] == false) done_all = false;
            }

            // broadcast done_all to all processes
            MPI_Bcast(&done_all, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

        } // end while

        // Root process gather sim_local->u of all the processes into sim_u_all
        MPI_Gatherv(sim_local->u, myel, MPI_FLOAT, sim_u_all, blockunum, 
                    blockuindex, MPI_FLOAT, 0, MPI_COMM_WORLD);
        
        // Root process copy sim_u_all into sim_global
        if (world_rank == 0){
            for(int blockindexj = 0; blockindexj < nparty; ++blockindexj){
                for (int blockindexi = 0; blockindexi < npartx; ++blockindexi){
                    sub_field_copyout(sim_u_all + blockuindex[blockindexi + npartx * blockindexj], \
                                      sim_global->u, nx_global, ny_global, \
                                      ng, offsets_x[blockindexi], offsets_x[blockindexi+1], \
                                      offsets_y[blockindexj], offsets_y[blockindexj+1]);
                }
            }
        }

        double t1 = omp_get_wtime();
        double elapsed = t1-t0;
#elif defined SYSTIME
        struct timeval t0, t1;
        gettimeofday(&t0, NULL);
        int nstep = central2d_run(sim_global, ftime, batch);
        gettimeofday(&t1, NULL);
        double elapsed = (t1.tv_sec-t0.tv_sec) + (t1.tv_usec-t0.tv_usec)*1e-6;
#else   
        int nstep = central2d_run(sim_global, ftime, batch);
        double elapsed = 0;
#endif
        if (world_rank == 0){
            solution_check(sim_global);
            tcompute += elapsed;
            printf("  Time: %e (%e for %d steps)\n", elapsed, elapsed/nstep[0], nstep[0]);
            viz_frame(viz, sim_global, vskip);
	    }
    } // end frames
    if (world_rank==0)
    {
        printf("Total compute time: %e\n", tcompute);
        viz_close(viz);
    }
    central2d_free(sim_global);
    central2d_mpi_free(sim_local);

    return 0;
}


/**
 * ### Main
 *
 * The main routine has the usage pattern
 *
 *     lshallow tests.lua args
 *
 * where `tests.lua` has a call to the `simulate` function to run
 * the simulation.  The arguments after the Lua file name are passed
 * into the Lua script via a global array called `args`.
 */

int main(int argc, char** argv)
{   
    MPI_Init(&argc, &argv);
    if (argc < 2) {
        fprintf(stderr, "Usage: %s fname args\n", argv[0]);
        return -1;
    }

    lua_State* L = luaL_newstate();
    luaL_openlibs(L);
    lua_register(L, "simulate", run_sim);

    lua_newtable(L);
    for (int i = 2; i < argc; ++i) {
        lua_pushstring(L, argv[i]);
        lua_rawseti(L, 1, i-1);
    }
    lua_setglobal(L, "args");

    if (luaL_dofile(L, argv[1]))
        printf("%s\n", lua_tostring(L,-1));

    // Finalize MPI environment.
    MPI_Finalize();

    lua_close(L);
    return 0;
}
