
#pragma once

#include <vector>

#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"

#include "MeshLib/Core/IMesh.h"

#include "NumLib/Discrete/DoF.h"
#include "NumLib/Discrete/DiscreteSystemAssembler.h"

namespace NumLib
{

/** 
 * \brief Interface of discrete linear equations
 */
class IDiscreteLinearEquation
{
public:
    // construct discrete linear equations
    virtual void construct(IDiscreteSystemAssembler& assemler) = 0;
    // solve 
    virtual void solve() = 0;
    // get solution
    virtual void getX() = 0;
};

/**
 * \brief Mesh based discrete linear equations
 */
class MeshBasedDiscreteLinearEquation : public IDiscreteLinearEquation
{
private:
    MathLib::ILinearEquations* _eqs;
    MeshLib::IMesh* _msh;

public:
    MeshBasedDiscreteLinearEquation(MeshLib::IMesh &msh, MathLib::ILinearEquations &eqs)
    {
        _msh = &msh;
        _eqs = &eqs;
    }

    void construct(IDiscreteSystemAssembler& assemler)
    {
        DofMapManager *dofManager;
        std::vector<size_t> list_dofId;

        assemler.assembly(*_msh, *dofManager, list_dofId, *_eqs);
        // apply fixed boundary conditions
    }

    void solve()
    {
        _eqs->solve();
    }

    void getX()
    {

    }
};


}
