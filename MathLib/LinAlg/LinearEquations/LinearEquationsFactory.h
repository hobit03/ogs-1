
#pragma once

#include "ILinearEquations.h"
#include "LinearEquations.h"
#include "LisInterface.h"

namespace MathLib
{

struct LinearEquationsType
{
    enum type {
        DenseEquations,
        SparseEquations,
        LIS,
        PARDISO,
        PETSC
    };
};

class LinearEquationsFactory
{
public:
    static ILinearEquations* create(LinearEquationsType::type eq_type)
    {
        switch (eq_type) {
        case LinearEquationsType::DenseEquations:
            return new DenseLinearEquations();
        case LinearEquationsType::SparseEquations:
            return new SparseLinearEquations();
        case LinearEquationsType::LIS:
            return new LIS_Solver();
        }
        assert(false);
        return 0;
    }
};
}
