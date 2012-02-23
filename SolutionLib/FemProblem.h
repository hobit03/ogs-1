
#pragma once

#include "IProblem.h"

#include <vector>

#include "Base/CodingTools.h"
#include "FemLib/Function/FemFunction.h"
#include "FemLib/BC/FemDirichletBC.h"
#include "FemLib/BC/FemNeumannBC.h"


namespace SolutionLib
{

/**
 * \brief IVBV problems using FEM
 */
template<class T_ASSEMBLY>
class FemIVBVProblem : public AbstractMeshBasedDiscreteIVBVProblem
{
public:
    FemIVBVProblem(MeshLib::IMesh &msh, T_ASSEMBLY &user_assembly) : AbstractMeshBasedDiscreteIVBVProblem(msh), _user_assembly(user_assembly)
    {
        Base::zeroObject(_map_var, _map_ic);
    }

    virtual ~FemIVBVProblem()
    {
        Base::releaseObject(_map_var, _map_ic);
        Base::releaseObjectsInStdVector(_map_bc1);
        Base::releaseObjectsInStdVector(_map_bc2);
    }

    /// get the number of variables
    size_t getNumberOfVariables() const { return 1; }

    size_t createField(FemLib::PolynomialOrder::type order)
    {
        if (_map_var==0) {
            _map_var = new FemLib::FemNodalFunctionScalar(*getMesh(), order);
        }
        return 0;
    }

    FemLib::FemNodalFunctionScalar* getField(size_t) const
    {
        return _map_var;
    }

    void setIC(int, MathLib::IFunction<double, GeoLib::Point>& ic)
    {
        _map_ic = ic.clone();
    }

    MathLib::IFunction<double, GeoLib::Point>* getIC(int) const 
    {
        return _map_ic;
    };

    void addDirichletBC(int, GeoLib::GeoObject &geo, bool is_transient, MathLib::IFunction<double, GeoLib::Point>& bc1)
    {
        addDirichletBC(*new FemLib::FemDirichletBC<double>(_map_var, &geo, is_transient, &bc1, &FemLib::DiagonalizeMethod()));
    }


    size_t getNumberOfDirichletBC(int n=0) const {return _map_bc1.size();};

    MathLib::IFunction<double, GeoLib::Point>* getDirichletBC(int, int bc_id) const 
    {
        return _map_bc1[bc_id];
    };

    FemLib::FemDirichletBC<double>* getFemDirichletBC(int bc_id) const 
    {
        return _map_bc1[bc_id];
    };

    void addNeumannBC(int, GeoLib::GeoObject &geo, bool is_transient, MathLib::IFunction<double, GeoLib::Point>& bc2)
    {
        addNeumannBC(*new FemLib::FemNeumannBC<double, double>(_map_var, &geo, is_transient, &bc2));
    }

    size_t getNumberOfNeumannBC(int n=0) const {return _map_bc2.size();};

    MathLib::IFunction<double, GeoLib::Point>* getNeumannBC(int, int bc_id) const 
    {
        return _map_bc2[bc_id];
    };

    FemLib::FemNeumannBC<double, double>* getFemNeumannBC(int bc_id) const 
    {
        return _map_bc2[bc_id];
    };

    T_ASSEMBLY& getElementAssemlby()
    {
        return _user_assembly;
    }

private:
    T_ASSEMBLY _user_assembly;
    FemLib::FemNodalFunctionScalar* _map_var;
    MathLib::IFunction<double, GeoLib::Point>* _map_ic;
    std::vector<FemLib::FemDirichletBC<double>*> _map_bc1;
    std::vector<FemLib::FemNeumannBC<double, double>*> _map_bc2;


    void addDirichletBC(FemLib::FemDirichletBC<double>& bc1)
    {
        _map_bc1.push_back(&bc1);
    }

    void addNeumannBC(FemLib::FemNeumannBC<double, double>& bc2)
    {
        _map_bc2.push_back(&bc2);
    }

    DISALLOW_COPY_AND_ASSIGN(FemIVBVProblem);
};


/**
 * \brief IVBV problems using FEM
 */
template<class T_ASSEMBLY>
class FemIVBVProblemWithMultipleVariables : public AbstractMeshBasedDiscreteIVBVProblem
{
public:
    FemIVBVProblemWithMultipleVariables(MeshLib::IMesh &msh, T_ASSEMBLY &user_assembly) : AbstractMeshBasedDiscreteIVBVProblem(msh), _user_assembly(user_assembly)
    {
    }

    virtual ~FemIVBVProblemWithMultipleVariables()
    {
    }

    /// get the number of variables
    size_t getNumberOfVariables() const
    {
        return _map_var.size();
    }

    size_t createField(FemLib::PolynomialOrder::type order)
    {
        resizeAll(_map_var.size()+1);
        _map_var[_map_var.size()-1] = new FemLib::FemNodalFunctionScalar(getMesh(), order);
        return _map_var.size()-1;
    }

    FemLib::FemNodalFunctionScalar* getField(size_t i) const
    {
        return _map_var[i];
    }

    void setIC(int var_type, MathLib::IFunction<double, GeoLib::Point>& ic)
    {
        _map_ic[var_type] = &ic;
    }

    MathLib::IFunction<double, GeoLib::Point>* getIC(int var_type) const 
    {
        return _map_ic[var_type];
    };

    void addDirichletBC(int var_type, GeoLib::GeoObject &geo, MathLib::IFunction<double, GeoLib::Point>& bc1)
    {
        addDirichletBC(var_type, *new FemLib::FemDirichletBC<double>(_map_var[var_type], &geo, &bc1, &FemLib::DiagonalizeMethod()));
    }


    size_t getNumberOfDirichletBC(int var_type) const {return _map_bc1[var_type].size();};

    MathLib::IFunction<double, GeoLib::Point>* getDirichletBC(int var_type, int bc_id) const 
    {
        return _map_bc1[var_type][bc_id];
    };

    FemLib::FemDirichletBC<double>* getFemDirichletBC(int var_type, int bc_id) const 
    {
        return _map_bc1[var_type][bc_id];
    };

    void addNeumannBC(int var_type, GeoLib::GeoObject &geo, MathLib::IFunction<double, GeoLib::Point>& bc2)
    {
        addNeumannBC(var_type, *new FemLib::FemNeumannBC<double, double>(_map_var[var_type], &geo, &bc2));
    }

    size_t getNumberOfNeumannBC(int var_type) const {return _map_bc2[var_type].size();};

    MathLib::IFunction<double, GeoLib::Point>* getNeumannBC(int var_type, int bc_id) const 
    {
        return _map_bc2[var_type][bc_id];
    };

    FemLib::FemNeumannBC<double, double>* getFemNeumannBC(int var_type, int bc_id) const 
    {
        return _map_bc2[var_type][bc_id];
    };

    T_ASSEMBLY& getElementAssemlby()
    {
        return _user_assembly;
    }

private:
    T_ASSEMBLY _user_assembly;
    std::vector<FemLib::FemNodalFunctionScalar*> _map_var;
    std::vector<MathLib::IFunction<double, GeoLib::Point>*> _map_ic;
    std::vector<std::vector<FemLib::FemDirichletBC<double>*>> _map_bc1;
    std::vector<std::vector<FemLib::FemNeumannBC<double, double>*>> _map_bc2;


    void addDirichletBC(int var_type, FemLib::FemDirichletBC<double>& bc1)
    {
        _map_bc1[var_type].push_back(&bc1);
    }

    void addNeumannBC(int var_type, FemLib::FemNeumannBC<double, double>& bc2)
    {
        _map_bc2[var_type].push_back(&bc2);
    }

    void resizeAll(size_t n)
    {
        _map_var.resize(n);
        _map_ic.resize(n);
        _map_bc1.resize(n);
        _map_bc2.resize(n);
    }
    DISALLOW_COPY_AND_ASSIGN(FemIVBVProblemWithMultipleVariables);
};



} //end
