
#include "LisInterface.h"

#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif


namespace MathLib
{

void CRSLisSolver::initialize()
{
    int argc=0;
    char **argv;
    lis_initialize(&argc, &argv);
}

void CRSLisSolver::finalize()
{
    lis_finalize();
}

CRSLisSolver::~CRSLisSolver()
{
    lis_vector_destroy(bb);
    lis_vector_destroy(xx);
}

void CRSLisSolver::setOption(const Base::Options &option)
{
    throw "LisSolver::setOption() is not implemented.";
}

void CRSLisSolver::solveEqs(CRSMatrix<double, signed> *A, double *b, double *x)
{
    long dimension = static_cast<long>(A->getNRows());

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "*** LIS solver computation" << std::endl;

    // Creating a matrix.
    LIS_SOLVER solver;
    LIS_MATRIX AA;
#ifndef USE_MPI
    int ierr = lis_matrix_create(0, &AA);
    ierr = lis_matrix_set_type(AA, LIS_MATRIX_CRS);
    ierr = lis_matrix_set_size(AA, 0, dimension);
#else
    lis_matrix_create(MPI_COMM_WORLD, &AA);
    lis_matrix_set_size(AA,dimension,0);
//    lis_matrix_get_size(AA, &_local_dim, &_global_dim);
#endif

    // Matrix solver and Precondition can be handled better way.
    const size_t MAX_ZEILE = 512;
    char solver_options[MAX_ZEILE], tol_option[MAX_ZEILE];

    int nthreads = omp_get_num_threads();
    //omp_set_num_threads (nthreads);

    sprintf(solver_options, "-i %d -p %d %s", _option.ls_method, _option.ls_precond, _option.ls_extra_arg.c_str()); 
    sprintf(tol_option, "-tol %e -maxiter %d -omp_num_threads %d -initx_zeros 0", _option.ls_error_tolerance, _option.ls_max_iterations, nthreads);

    int ierr;
    ierr = lis_matrix_set_crs(A->getNNZ(), (int*)A->getRowPtrArray(), (int*)A->getColIdxArray(), (double*)A->getEntryArray(), AA);
    ierr = lis_matrix_assemble(AA);

    // Assemble the vector, b, x
    ierr = lis_vector_duplicate(AA, &bb);
    ierr = lis_vector_duplicate(AA, &xx);
    #pragma omp parallel for
    for (long i=0; i < dimension; ++i)
    {
        ierr = lis_vector_set_value(LIS_INS_VALUE, i, x[i], xx);
        ierr = lis_vector_set_value(LIS_INS_VALUE, i, b[i], bb);
    }

    // Create solver
    ierr = lis_solver_create(&solver);

    ierr = lis_solver_set_option(solver_options, solver);
    ierr = lis_solver_set_option(tol_option, solver);
    ierr = lis_solver_set_option("-print mem", solver);
    
    ierr = lis_solve(AA, bb, xx, solver);
    int iter = 0;
    double resid = 0.0;
    ierr = lis_solver_get_iters(solver, &iter);
    ierr = lis_solver_get_residualnorm(solver, &resid);
    printf("\t iteration: %d/%d\n", iter, _option.ls_max_iterations);
    printf("\t residuals: %e\n", resid);
    //	lis_vector_print(xx);
    //	lis_vector_print(bb);

    // Update the solution (answer) into the x vector
    #pragma omp parallel for
    for(long i=0; i<dimension; ++i)
    {
        lis_vector_get_value(xx,i,&(x[i]));
    }

    // Clear memory
    //lis_matrix_destroy(AA); // CRSMatrix will free this memory
    lis_solver_destroy(solver);
    std::cout << "------------------------------------------------------------------" << std::endl;
}

void CRSLisSolver::gatherX(std::vector<double> &x)
{
    int local_n, global_n;
    lis_vector_get_size(xx, &local_n, &global_n);
    x.resize(global_n);
    lis_vector_gather(xx, &x[0]);
};




}
