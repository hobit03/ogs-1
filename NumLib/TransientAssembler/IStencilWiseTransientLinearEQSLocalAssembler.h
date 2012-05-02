
#pragma once

#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"


namespace MeshLib
{
class IStencil;
};

namespace NumLib
{

class TimeStep;

class IStencilWiseTransientLinearEQSLocalAssembler
{
public:
	typedef MathLib::DenseLinearEquations LocalEquationType;
	typedef LocalEquationType::MatrixType LocalMatrixType;
	typedef LocalEquationType::VectorType LocalVectorType;

    virtual ~IStencilWiseTransientLinearEQSLocalAssembler() {};

    /// assemble a local linear equation for the given element
    /// @param time			time step
    /// @param e			element
    /// @param local_u_n1	guess of current time step value
    /// @param local_u_n	previous time step value
    /// @param eqs			local algebraic equation
    virtual void assembly(const TimeStep &time,  MeshLib::IStencil &nod, const LocalVectorType &local_u_n1, const LocalVectorType &local_u_n, LocalEquationType &eqs) = 0;
};

} //end
