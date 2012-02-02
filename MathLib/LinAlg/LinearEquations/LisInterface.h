
#pragma once

#include <string>
#include "MathLib/LinAlg/Sparse/SparseTableCRS.h"
#include "ILinearEquations.h"
#include "LinearEquations.h"

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
	
typedef struct {
    int ls_method;
    int ls_precond;
    std::string ls_extra_arg;
    long ls_max_iterations;
    double ls_error_tolerance;
} LIS_option;

class LIS_Solver : public CRSLinearEquationsBase<signed>
{
public:
    void initialize();
    void finalize();

    void setOption(const Base::Options &option);
    void setOption(const LIS_option &option)
    {
        _option = option;
    }
    LIS_option &getOption() 
    {
        return _option;
    }

    void solve();
private:
    LIS_option _option;
};

}
