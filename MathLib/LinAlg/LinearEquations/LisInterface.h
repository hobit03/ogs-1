
#pragma once

#include <string>

#include "lis.h"

#include "MathLib/LinAlg/Sparse/SparseTableCRS.h"
#include "ILinearEquations.h"
#include "SparseLinearEquationBase.h"

/*
LIS matrix solver options
CG -i {cg|1}
BiCG -i {bicg|2}
CGS -i {cgs|3}
BiCGSTAB -i {bicgstab|4}
BiCGSTAB(l) -i {bicgstabl|5} -ell [2] Value for l
GPBiCG -i {gpbicg|6}
TFQMR -i {tfqmr|7}
Orthomin(m) -i {orthomin|8} -restart [40] Value for Restart m
GMRES(m) -i {gmres|9} -restart [40] Value for Restart m
Jacobi -i {jacobi|10}
Gauss-Seidel -i {gs|11}
SOR -i {sor|12} -omega [1.9] Value for Relaxation Coefficient  (0 <  < 2)
BiCGSafe -i {bicgsafe|13}
CR -i {cr|14}
BiCR -i {bicr|15}
CRS -i {crs|16}
BiCRSTAB -i {bicrstab|17}
GPBiCR -i {gpbicr|18}
BiCRSafe -i {bicrsafe|19}
FGMRES(m) -i {fgmres|20} -restart [40] Value for Restart m
IDR(s) -i {idrs|21} -restart [40] Value for Restart s

Preconditioner Option Auxiliary Option
None -p {none|0}
Jacobi -p {jacobi|1}
ILU(k) -p {ilu|2} -ilu_fill [0] Fill level k
SSOR -p {ssor|3} -ssor_w [1.0] Relaxation Coefficient  (0 <  < 2)
Hybrid -p {hybrid|4} -hybrid_i [sor] Iterative method
-hybrid_maxiter [25] Maximum number of iterations
-hybrid_tol [1.0e-3] Convergence criteria
-hybrid_w [1.5] Relaxation Coefficient  for
the SOR method (0 <  < 2)
-hybrid_ell [2] Value for l of the BiCGSTAB(l) method
-hybrid_restart [40] Restart values for GMRES and Orthomin
I+S -p {is|5} -is_alpha [1.0] Parameter ?for preconditioner
of a I + ?(m) type
-is_m [3] Parameter m for preconditioner
of a I + ?(m) type
SAINV -p {sainv|6} -sainv_drop [0.05] Drop criteria
SA-AMG -p {saamg|7} -saamg_unsym [false] Selection of asymmetric version
Crout ILU -p {iluc|8} -iluc_drop [0.05] Drop criteria
-iluc_rate [5.0] Ratio of Maximum fill-in
ILUT -p {ilut|9} -ilut_drop [0.05] Drop criteria
-ilut_rate [5.0] Ratio of Maximum fill-in
additive Schwarz -adds true -adds_iter [1] Number of iterations
*/

namespace MathLib
{

struct LIS_option
{
    int ls_method;
    int ls_precond;
    std::string ls_extra_arg;
    long ls_max_iterations;
    double ls_error_tolerance;
    enum SolverType
    {
        CG = 1,
        BiCG = 2,
        CGS = 3,
        BiCGSTAB = 4,
        BiCGSTABl = 5,
        GPBiCG = 6,
        TFQMR = 7,
        Orthomin = 8,
        GMRES = 9
    };
    enum PreconType
    {
        NONE = 0,
        Jacobi = 1,
        ILU = 2
    };
    LIS_option()
    {
        ls_method = CG;
        ls_precond = NONE;
        ls_max_iterations = 500;
        ls_error_tolerance = 1.e-6;
    }
};

class CRSLisSolver : public CRSLinearEquationsBase<signed>
{
public:
    void initialize();
    void finalize();
    virtual ~CRSLisSolver();

    void setOption(const Base::Options &option);
    void setOption(const LIS_option &option)
    {
        _option = option;
    }
    LIS_option &getOption() 
    {
        return _option;
    }

    void gatherX(std::vector<double> &x);

protected:
    void solveEqs(CRSMatrix<double, signed> *A, double *rhs, double *x);

private:
    LIS_option _option;
    LIS_VECTOR bb,xx;
};

#if 1
class LisSolver // : public ILinearEquations
{
public:
    LisSolver() 
    {
        _global_dim = 0;
        _local_dim = 0;
        _dynamic = false;
    }
    virtual ~LisSolver();

    static void initialize(int argc, char* argv[]);
    static void finalize();

    void setOption(const LIS_option &option)
    {
        _option = option;
    }
    LIS_option &getOption() 
    {
        return _option;
    }

    void createDynamic(size_t local_n, size_t global_n);
    SparseTableCRS<int>* createCRS(size_t local_n, size_t global_n);
    void setOption(const Base::Options &option);
    void reset();

    size_t getDimension() const { return _global_dim; }
    size_t getLocalDimension() const { return _local_dim; }


    double getA(size_t rowId, size_t colId);
    void setA(size_t rowId, size_t colId, double v);
    void addA(size_t rowId, size_t colId, double v);
    void addA(std::vector<size_t> &vec_row_pos, std::vector<size_t> &vec_col_pos, MathLib::Matrix<double> &sub_matrix, double fkt=1.0);
    void addA(std::vector<size_t> &vec_pos, MathLib::Matrix<double> &sub_matrix, double fkt=1.0);

    double getRHS(size_t rowId);
    double* getRHS();
    void setRHS(size_t rowId, double v);
    void addRHS(std::vector<size_t> &vec_pos, double *sub_vector, double fkt=1.0);
    void addRHS(size_t rowId, double v);

    double* getX();

    void setKnownX(size_t row_id, double x);
    void setKnownX(const std::vector<size_t> &vec_id, const std::vector<double> &vec_x);
    void solve();

    void getRange(int &is, int &ie) const {is = _is; ie = _ie;};
    SparseTableCRS<int>* getCRS() {return &_crs;} ;
    void assembleMatrix();

    size_t createVector()
    {
        LIS_VECTOR u;
        int err = lis_vector_duplicate(_A,&u); CHKERR(err);
        _vec_u.push_back(u);
        return _vec_u.size()-1;
    }
    void destroyVector(size_t i)
    {
        LIS_VECTOR &u = _vec_u[i];
        lis_vector_destroy(u);
        _vec_u.erase(_vec_u.begin()+i);
    }
    void setVectorAll(size_t i, double v)
    {
        LIS_VECTOR &u = _vec_u[i];
        lis_vector_set_all(1.0, u);
    }
    void matvecToRHS(size_t i)
    {
        LIS_VECTOR &u = _vec_u[i];
        lis_matvec(_A,u,_b);
    }
private:
    LIS_option _option;
    LIS_MATRIX _A;
    LIS_VECTOR _b;
    LIS_VECTOR _x;
    int _local_dim;
    int _global_dim;
    std::vector<double> _tmp_b;
    std::vector<double> _tmp_x;
    int _is;
    int _ie;
    SparseTableCRS<int> _crs;
    std::vector<LIS_VECTOR> _vec_u;
    bool _dynamic;
};
#endif

}
