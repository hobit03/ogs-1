
#include <iostream>
#include <stdio.h>
#include "mpi.h"
#include "lis.h"

#include "MathLib/LinAlg/LinearEquations/LisInterface.h"
#include "DiscreteLib/ogs5/par_ddc_group.h"

using namespace MathLib;

int nprocs = 0;
int my_rank = 0;

void initialize(int argc, char *argv[])
{
    LisSolver::initialize(argc, argv);

#ifdef USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
#else
    nprocs  = 1;
    my_rank = 0;
#endif
    if( my_rank==0 )
    {
        printf("\n");
        printf("number of processes = %d\n",nprocs);
    }
}

void finalize()
{
    LisSolver::finalize();
}

int testOGS5(int argc, char *argv[])
{
    OGS5::CPARDomainGroup *master;

    return 0;
}

int testLIS(int argc, char *argv[])
{
    // input parameters
    if( argc < 3 ) {
        if( my_rank==0 ) printf("Usage: test5 n gamma [options]\n");
        return -1;
    }

    int gn  = atoi(argv[1]);
    double gamma  = atof(argv[2]);

    // solve
    LisSolver lis;
#if 0
    SparseTableCRS<int>* crs = lis.createCRS(0, gn);
    int is = 0;
    int ie = 0;
    lis.getRange(is,ie);

    int k = 0;
    crs->row_ptr[0] = 0;
    for (int ii=is;ii<ie;ii++) {
        int jj=0;
        if( ii>1 )    { jj = ii - 2; crs->col_idx[k] = jj; crs->data[k++] = gamma;}
        if( ii<gn-1 ) { jj = ii + 1; crs->col_idx[k] = jj; crs->data[k++] = 1.0;}
        crs->col_idx[k] = ii; crs->data[k++] = 2.0;
        crs->row_ptr[ii-is+1] = k;
    }
#else
    lis.createDynamic(0, gn);
    int is = 0;
    int ie = 0;
    lis.getRange(is,ie);
    std::cout << my_rank << ": is=" << is << ", ie=" << ie << std::endl;
    for (int i=is;i<ie;i++) {
        if( i>1   )  lis.addA(i,i-2,gamma);
        if( i<gn-1 ) lis.addA(i,i+1,1.0);
        lis.addA(i,i,2.0);
    }
#endif

    lis.assembleMatrix();
    size_t iu = lis.createVector();
    lis.setVectorAll(iu, 1.0);
    lis.matvecToRHS(iu);
    lis.destroyVector(iu);
    lis.solve();

    return 0;
}

int main(int argc, char *argv[])
{
    initialize(argc, argv);
    testLIS(argc, argv);
    finalize();
#if 0
    MPI_Init(&argc, &argv);

    int size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank==0) {
        std::cout << "MPI program started with " << size << " processes " << std::endl;
        std::cout << std::endl;
        std::cout << "Who are you?" << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    std::cout << "Hello, I am process id " << rank  << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank==0) {
        std::cout << std::endl << "Talk each other!" << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    int say[] = {1, 2, 3, 4, 5, 6, 7};
    int hear[10] = {};
    if (rank%2==0) {
        int talk_to = (rank+1)%size;
        int hear_from = (rank+1)%size;
        int n_say = rank;
        int n_hear = rank+1;
        MPI_Send(say, n_say, MPI_INT, talk_to, 0, MPI_COMM_WORLD);
        std::cout << rank <<  " says hello to " << talk_to << std::endl;
        MPI_Status status;
        MPI_Recv(hear, n_hear, MPI_INT, talk_to, 0, MPI_COMM_WORLD, &status);
        std::cout << rank << " hear hello from " << talk_to << std::endl;
    } else {
        int talk_to = (rank-1)%size;
        int hear_from = (rank-1)%size;
        int n_say = rank;
        int n_hear = rank-1;
        MPI_Status status;
        MPI_Recv(hear, n_say, MPI_INT, talk_to, 0, MPI_COMM_WORLD, &status);
        std::cout << rank << " hear from " << talk_to << std::endl;
        MPI_Send(say, n_hear, MPI_INT, talk_to, 0, MPI_COMM_WORLD);
        std::cout << rank <<  " says hello to " << talk_to << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0) {
        std::cout << std::endl << "I distribute this to you!" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    int dist_v = 123;
    MPI_Bcast(&dist_v, 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::cout << rank <<  " received a number " << dist_v << std::endl;

    MPI_Finalize();
#endif
    return 0;
}

