/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LocalProblem.h
 *
 * Created on 2013-05-29 by Reza Zolfaghari & Haibing Shao
 */
#ifndef LOCAL_PROBLEM_H
#define LOCAL_PROBLEM_H


#include   <functional>
#include   "ogs6THMC/ChemLib/chemReductionGIA.h"
#include "Local_ODE_Xi_immob_GIA.h"
#include "StepperBulischStoer.h"
//#include "chemconst.h"
//#include "chemcomp.h"
//#include "chemReactionEq.h"
//#include "BaseLib/OrderedMap.h"

class LocalProblem
{
public:
	/**
      * constructor of the class
      */
	LocalProblem(ogsChem::chemReductionGIA* ReductionGIA, MathLib::StepperBulischStoer<Local_ODE_Xi_immob_GIA>* sbs, Local_ODE_Xi_immob_GIA* local_ode_xi_immob_GIA);
	
	/**
      * destructor of the class
      */
	virtual ~LocalProblem(void)
    {
    	BaseLib::releaseObject(_ReductionGIA);
    	BaseLib::releaseObject(_local_ode_xi_immob_GIA);
    	BaseLib::releaseObject(_sbs);
    };

    /**
      * whether the reduction scheme has been initialized
      */
	bool IsInitialized(void) {return isInitialized;};

    /**
      * get the number of unknowns
      */
	size_t get_n_unknowns_local_problem(void) {return _n_unknowns; };

    /**
      * calculate the residual of the reaction system using 
      * the concentration vector of all components and
      * the total concentration constrain of the basis species
      */
    void calc_residual(double dt,
    				ogsChem::LocalVector & vec_unknowns,
    				ogsChem::LocalVector & vec_residual,
    				ogsChem::LocalVector & vec_Xi_Kin_bar);
    
    /**
      * calculate the Jacobi of the reaction system analytically using 
      * the concentration vector of all components, 
      * the total concentration constrain of the basis species, 
      * and an already calculated residual vector
      */
    void calc_Jacobian( double dt,
    					ogsChem::LocalVector & vec_unknowns,
    //				 ogsChem::LocalVector & vec_AI,
    				 ogsChem::LocalVector & vec_residual,
    				 ogsChem::LocalVector & vec_Xi_Kin_bar);

    /**
      * solve the equilibrium reaction system 
      * using the Newton iterations, input values are
      * the vector of concentration values
      * the convergence tolerance, 
      * the maximum number of iterations. 
      * and the status of convergence
      */
    void solve_LocalProblem_Newton_LineSearch(ogsChem::LocalVector & vec_unknowns,
    										  ogsChem::LocalVector & vec_tot_mass_constrain,
    										  std::size_t & node_idx,
    										  double  deltaT,
    										  const double iter_tol,
    										  const double rel_tol,
    										  const double max_iter);

    /**
      * solve the system J*dx = -b
      *
      * inputs are: 
      *    - idx_node is the index of node
      *    - J is Jacobi marix
      *    - b is residual vector
      * 
      * if rank of J is equal to number of b and dx, 
      * using direct solver provided by eigen
      * if rank of J is smaller than number of b, 
      * using minimization method
      * delta_x is the returned result. 
      */
    void Solv_Minimization(size_t      & idx_node,
                  ogsChem::LocalMatrix & J,
                  ogsChem::LocalVector & b,
                  ogsChem::LocalVector & delta_x);

    //TODO declare the residual functions


    void cal_ln_conc_vec(size_t    		        idx_size,
    					 ogsChem::LocalVector & conc_Mob,
    				     ogsChem::LocalVector & ln_conc_Mob);

    void cal_exp_conc_vec(size_t    		    idx_size,
    					  ogsChem::LocalVector & ln_conc,
    					  ogsChem::LocalVector & Conc);

    void ODE_solver(double dt,
    			    ogsChem::LocalVector & vec_unknowns,
    			    ogsChem::LocalVector & vec_Xi_Kin_bar);

private:
	/**
      * private flag indicating initialization
      */
	bool isInitialized; 

    /**
      * pointer to the reduction scheme.
      */
	ogsChem::chemReductionGIA* _ReductionGIA;

	MathLib::StepperBulischStoer<Local_ODE_Xi_immob_GIA>* _sbs;

	/**
	  * a list of all kinetic reactions
	  */
	std::vector<ogsChem::chemReactionKin*> & _list_kin_reactions;

    ogsChem::LocalMatrix _mat_Ald,_mat_c_mob_2_xi_mob,_mat_c_immob_2_xi_immob, _mat_c_mob_2_eta_mob, _mat_c_immob_2_eta_immob;

    /**
      * stoichiometric matrix S
      */	
	ogsChem::LocalMatrix _matStoi, _mat_S1min, _mat_Ssorp, _mat_A2kin, _mat_S1mob;

    /**
      * stoichiometric matrix from the input reactions
      */	
	ogsChem::LocalMatrix _matStoi_input;

    /**
      * vector of natural log K values
      */	
    ogsChem::LocalVector _logk_mob, _logk_min, _logk_sorp;


    /**
     * vector of total mass constrains ie. xi global and eta values
     */
    ogsChem::LocalVector _vec_eta , _vec_XiSorpTilde, _vec_XiMinTilde, _vec_etabar, _vec_rateKin, _vec_Xikin, _vec_XiBarKin, _vec_XiBarKin_old;

	/**
      * _I is the number of components and _J is the number of reactions
      */
	std::size_t const  _n_Comp;
	std::size_t _n_xi_global, _n_xi_local;

    //double  deltaT, theta_waterContent;

    /**
      * number of mobile, sorption and mineral components
      */
	size_t _I_mob, _I_sorp, _I_min, _I_NMin_bar;

	/**
	 * number of species
	 */
    size_t _n_xi_Mob,_n_xi_Min, _n_xi_Min_bar, _n_xi_Sorp_bar, _n_xi_Sorp_bar_li, _n_xi_Sorp_bar_ld, _J_tot_kin;

    /**
      * number of total unknowns
      */
	size_t _n_unknowns;

	/**
	 * number of linearly independent kinetic reactions
	 */
	size_t _n_xi_Kin_bar, _n_xi_Kin;


	//TODO specify each
	size_t _n_eta, _n_eta_bar, _n_xi_Sorp_tilde, _n_xi_Min_tilde, _n_xi_Sorp;


    /**
      * since jacobian matrix is unique we define it in the scope of the class
      * its size is equal to the _n_unknowns by _n_unknowns
      */
    ogsChem::LocalMatrix _mat_Jacobian;

    /**
      * the nested ode local problem
      */
    Local_ODE_Xi_immob_GIA*                       _local_ode_xi_immob_GIA;


    void reaction_rates(ogsChem::LocalVector & conc,
    					ogsChem::LocalVector & vec_rateKin);

    /**
      * increment the unknown with a damping factor
      * to prevent the concentrations falling into negative values
      */    
    void increment_unknown(ogsChem::LocalVector & x_old,
                           ogsChem::LocalVector & delta_x,
                           ogsChem::LocalVector & x_new);
    
    /**
      * update the concentration of minerals
      */
    void update_minerals_conc_AI(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalVector & vec_AI);
    
    /**
      * calcuate one particular mineral concentration 
      * by the amount of basis concentration and total mass constrain
      */
    double cal_cbarmin_by_constrain(size_t idx_min,
    								ogsChem::LocalVector & conc_Mob,
    								ogsChem::LocalVector & conc_NonMin_bar,
    								ogsChem::LocalVector & conc_Min_ba);

    /**
      * calculate the residual of mobile species
      */




    // Eq. 3.55
    void residual_conc_Mob			(ogsChem::LocalVector & conc_Mob,
    								 ogsChem::LocalVector & vec_residual);
    // Eq. 3.56
    void residual_Eta				(ogsChem::LocalVector & conc_Mob,
    								 ogsChem::LocalVector & vec_residual);
    // Eq. 3.57 - 58
    void residual_xi_Sorp_tilde	(ogsChem::LocalVector & conc_Mob,
    								 ogsChem::LocalVector & conc_NonMin_bar,
    								 ogsChem::LocalVector & conc_Min_bar,
    								 ogsChem::LocalVector & vec_residual);
    // Eq. 3.59
    void residual_xi_Min_tilde		(ogsChem::LocalVector & conc_Mob,
    							     ogsChem::LocalVector & conc_NonMin_bar,
    							     ogsChem::LocalVector & conc_Min_bar,
    							     ogsChem::LocalVector & vec_residual);
    // Eq. 3.60
    void residual_xi_Kin			(ogsChem::LocalVector & conc_Mob,
			 	 	 	 	 	 	 ogsChem::LocalVector & vec_residual);
    // Eq. 3.61
    void residual_conc_Sorp		(ogsChem::LocalVector & conc_Mob,
    								ogsChem::LocalVector & conc_NonMin_bar,
			 	 	 	 	 	 	 ogsChem::LocalVector & vec_residual);
    // Eq. 3.62
    void residual_conc_Min			(ogsChem::LocalVector & conc_Mob,
								     ogsChem::LocalVector & vec_AI,
								     ogsChem::LocalVector & vec_residual);
    // Eq. 3.63
    void residual_Eta_bar			(ogsChem::LocalVector & conc_NonMin_bar,
		     	 	 	 	 	 	 ogsChem::LocalVector & conc_Min_bar,
		     	 	 	 	 	 	 ogsChem::LocalVector & vec_residual);
    // Eq. 3.64
    void residual_xi_KinBar_Eq		(ogsChem::LocalVector & conc_NonMin_bar,
    								 ogsChem::LocalVector & conc_Min_bar,
    								 ogsChem::LocalVector & Xi_Kin_bar,
    								 ogsChem::LocalVector & vec_residual);
    // Eq. 3.65
    void residual_xi_KinBar_Kin	(   double deltaT,
    								ogsChem::LocalVector & conc_Mob,
    								ogsChem::LocalVector & conc_NonMin_bar,
    								 ogsChem::LocalVector & conc_Min_bar,
    								 ogsChem::LocalVector & Xi_Kin_bar,
    								 ogsChem::LocalVector & vec_residual);


};

#endif
