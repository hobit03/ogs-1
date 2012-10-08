/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NonLinearReactiveTransportJacobianLocalAssembler.h
 *
 * Created on 2012-09-24 by Haibing Shao
 */

/**
  * This file is same as the MassTransportTimeODELocalAssembler.h
  * The difference is, the compound molecular diffusion coefficient is disabled, 
  */

#ifndef NON_LINEAR_REACTIVE_TRANSPORT_JACOBIAN_LOCAL_ASSEMBLER_H
#define NON_LINEAR_REACTIVE_TRANSPORT_JACOBIAN_LOCAL_ASSEMBLER_H

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Compound.h"
#include "ChemLib/chemReductionKin.h"
#include "Concentrations.h"
#include "Ogs6FemData.h"

template <class T_NODAL_FUNCTION_SCALAR, class T_FUNCTION_DATA>
class NonLinearReactiveTransportJacobianLocalAssembler: public NumLib::IElementWiseTransientJacobianLocalAssembler
{
public:
	NonLinearReactiveTransportJacobianLocalAssembler(FemLib::LagrangianFeObjectContainer* feObjects, ogsChem::chemReductionKin* ReductionScheme, T_FUNCTION_DATA* concentrations)
        : _feObjects(*feObjects), _vel(NULL), _reductionKin(ReductionScheme), _concentrations(concentrations), _xi_mob_rates(NULL), _xi_immob_rates(NULL)
    {
    };

    virtual ~NonLinearReactiveTransportJacobianLocalAssembler() {};

    void setVelocity(const NumLib::ITXFunction *vel)
    {
        _vel = const_cast<NumLib::ITXFunction*>(vel);
    }

	void set_xi_mob_rates(   std::vector<T_NODAL_FUNCTION_SCALAR*> * xi_mob_rates )
	{
	    _xi_mob_rates = xi_mob_rates; 
	}
    
	void set_xi_immob_rates( std::vector<T_NODAL_FUNCTION_SCALAR*> * xi_immob_rates )
	{
	    _xi_immob_rates = xi_immob_rates; 
	}

	void set_xi_mob_drate_dxi( std::vector<T_NODAL_FUNCTION_SCALAR*> * drates_dxi )
	{
	    _drates_dxi = drates_dxi; 
	}

	T_FUNCTION_DATA* get_function_data(void)
	{
	    return _concentrations; 
	}

	void assembly(const NumLib::TimeStep &time, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &, const MathLib::LocalVector & u1, const MathLib::LocalVector & u0, MathLib::LocalMatrix & localJ)
    {
		size_t i, j, k, m, n, node_idx; 
  		size_t n_nodes = e.getNumberOfNodes(); 
		size_t n_xi_mob = _xi_mob_rates->size(); 
        size_t mat_id  = e.getGroupID(); 
 
        // clear the local Jacobian matrix
        localJ = MathLib::LocalMatrix::Zero(n_nodes*n_xi_mob, n_nodes*n_xi_mob);

		FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
        MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
        double cmp_mol_diffusion = .0;
        // _cmp->molecular_diffusion->eval(0, cmp_mol_diffusion);
        double dt = time.getTimeStepSize();  // time step size
        double theta = 1.0;

        MathLib::LocalMatrix matM(n_nodes, n_nodes);
        MathLib::LocalMatrix matDiff(n_nodes, n_nodes);
        MathLib::LocalMatrix matAdv(n_nodes, n_nodes);
        MathLib::LocalMatrix mat_dR(n_nodes, n_nodes); 
		MathLib::LocalVector vec_drates_dxi_value = MathLib::LocalVector::Zero(n_nodes);

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();

		size_t n_sp    = q->getNumberOfSamplingPoints();  // number of sampling points
		double gp_x[3], real_x[3];
        MathLib::LocalMatrix poro(1,1);
        MathLib::LocalMatrix d_poro(1,1);
        MathLib::LocalMatrix d_rate(1,1); 
        MathLib::LocalMatrix localJ_tmp = MathLib::LocalMatrix::Zero(n_nodes, n_nodes);
        NumLib::ITXFunction::DataType v;
		
        for (k=0; k<n_xi_mob; k++) {
			for (j=0; j<n_sp; j++) {
				q->getSamplingPoint(j, gp_x);
				fe->computeBasisFunctions(gp_x);
				fe->getRealCoordinates(real_x);

				pm->porosity->eval(real_x, poro);
				d_poro(0,0) = cmp_mol_diffusion * poro(0,0);
				_vel->eval(real_x, v);

				fe->integrateWxN(j, poro, matM);
				fe->integrateDWxDN(j, d_poro, matDiff);
				fe->integrateWxDN(j, v, matAdv); 
            
				// now dealing with the rate change terms
				// each location has n_xi_mob * n_xi_mob dR/dxi entries
				for (m=0; m<n_xi_mob; m++)
				{
					for (n=0; n<n_xi_mob; n++)
					{
						mat_dR = MathLib::LocalMatrix::Zero(n_nodes, n_nodes); 
						// loop over all the adjacent nodes to get the right drate/dxi value
               			for (k=0; k<n_nodes; k++)
             			{
							node_idx = e.getNodeID( k ); 
             				vec_drates_dxi_value( k ) = _drates_dxi->at(m*n_xi_mob+n)->getValue( node_idx ); 
						}  // end of for k<n_nodes
						// get the mean value
						d_rate(0,0) = vec_drates_dxi_value.mean(); 
						fe->integrateWxN(j, d_rate, mat_dR);
						// plugging mat_dR to the corresponding Jacobian matrix position
						localJ.block(n_nodes*m,n_nodes*n,n_nodes,n_nodes) -= mat_dR;
					}  // end of for n
				}  // end of for m
			}  // end of for j

			matM *= 1.0 / dt;
			matDiff *= theta;
			matAdv *= theta;
			localJ_tmp = matM;
			localJ_tmp += matDiff;
			localJ_tmp += matAdv;
            
            localJ.block(n_nodes*k,n_nodes*k,n_nodes,n_nodes) += localJ_tmp;

        }  // end of for k


        //std::cout << "M="; localM.write(std::cout); std::cout << std::endl;
        //std::cout << "L="; matDiff.write(std::cout); std::cout << std::endl;
        //std::cout << "A="; matAdv.write(std::cout); std::cout << std::endl;
    }  // end of function assembly

private:
    FemLib::LagrangianFeObjectContainer _feObjects;
    NumLib::ITXFunction* _vel;
	ogsChem::chemReductionKin* _reductionKin; 

	T_FUNCTION_DATA* _concentrations; 

	std::vector<T_NODAL_FUNCTION_SCALAR*> * _xi_mob_rates;
	std::vector<T_NODAL_FUNCTION_SCALAR*> * _xi_mob_rates_old;
    std::vector<T_NODAL_FUNCTION_SCALAR*> * _xi_immob_rates;
	std::vector<T_NODAL_FUNCTION_SCALAR*> * _drates_dxi;
};

#endif  // end of ifndef
