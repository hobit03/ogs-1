
#pragma once

#include <vector>
#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/Sparse/CRSMatrix.h"

namespace MathLib
{

/**
 * \brief 
 */
template<typename IDX_TYPE>
class CRSLinearEquationsBase : public ILinearEquations
{
public:
    void create(size_t length, RowMajorSparsity *sparsity)
    {
        SparseTableCRS<IDX_TYPE> *crs = convertRowMajorSparsityToCRS<IDX_TYPE>(*sparsity);
        assert (length == crs->dimension);
        _A = new CRSMatrix<double, IDX_TYPE>(crs->dimension, crs->row_ptr, crs->col_idx, crs->data);
        _b.resize(length);
        _x.resize(length);
    }

    void reset()
    {
        (*_A) = .0;
        _b.assign(_b.size(), .0);
        _x.assign(_x.size(), .0);
    }

    CRSMatrix<double, IDX_TYPE>* getA()
    {
        return _A;
    }

    double getA(size_t rowId, size_t colId)
    {
        return _A->getValue(rowId, colId);
    }

    void setA(size_t rowId, size_t colId, double v)
    {
        _A->setValue(rowId, colId, v);
    }

    void addA(size_t rowId, size_t colId, double v)
    {
        _A->addValue(rowId, colId, v);
    }

    void addA(std::vector<size_t> &vec_row_pos, std::vector<size_t> &vec_col_pos, MathLib::Matrix<double> &sub_matrix, double fkt=1.0)
    {
        for (size_t i=0; i<vec_row_pos.size(); i++) {
            const size_t rowId = vec_row_pos[i];
            for (size_t j=0; j<vec_col_pos.size(); j++) {
                const size_t colId = vec_row_pos[i];
                _A->addValue(rowId, colId, fkt*sub_matrix(i,j));
            }
        }
    }

    void addA(std::vector<size_t> &vec_pos, MathLib::Matrix<double> &sub_matrix, double fkt=1.0)
    {
        addA(vec_pos, vec_pos, sub_matrix, fkt);
    }

    double getRHS(size_t rowId)
    {
        return _b[rowId];
    }

    double* getRHS()
    {
        return &_b[0];
    }

    void setRHS(size_t rowId, double v)
    {
        _b[rowId] = v;
    }

    void addRHS(std::vector<size_t> &vec_row_pos, double *sub_vector, double fkt=1.0)
    {
        for (size_t i=0; i<vec_row_pos.size(); i++) {
            const size_t rowId = vec_row_pos[i];
            _b[rowId] += sub_vector[i];
        }
    }

    double* getX()
    {
        return &_x[0];
    }

    void setKnownX(size_t row_id, double x)
    {

    }

    void setKnownX(const std::vector<size_t> &vec_id, const std::vector<double> &vec_x)
    {
        for (size_t i=0; i<vec_id.size(); ++i)
            setKnownX(vec_id[i], vec_x[i]);
    }


private:
    CRSMatrix<double, IDX_TYPE> *_A;
    std::vector<double> _b;
    std::vector<double> _x;
};

/**
 * \brief 
 */
class SparseLinearEquations : public CRSLinearEquationsBase<unsigned>
{
public:
    enum SolverType
    {
        SolverCG,
        SolverBiCGStab,
        SolverGMRes
    };

    enum PreconditionerType
    {
        NONE,
        PreconDiag,
        PreconDiagScale
    };

    struct SpLinearOptions
    {
        SolverType solver_type;
        PreconditionerType precon_type;
        double error_tolerance;
        size_t max_iteration_step;

        SpLinearOptions()
        {
            solver_type = SolverCG;
            precon_type = NONE;
            error_tolerance = 1.e-10;
            max_iteration_step = 500;
        };
    };

    void initialize() {};
    void finalize() {};

    void setOption(const Base::Options &option);

    void setOption(const SpLinearOptions &option)
    {
        _option = option;
    }

    SpLinearOptions &getOption()
    {
        return _option;
    }


    void solve();

private:
    SpLinearOptions _option;

};

}