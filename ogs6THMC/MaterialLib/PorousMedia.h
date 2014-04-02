/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PorousMedia.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "IMedium.h"
#include "NumLib/Function/FunctionLinear.h"
#include "logog.hpp"
#include "MathLib/MathTools.h"


namespace NumLib
{
class ITXFunction;
class FunctionLinear1D;
}

namespace MaterialLib
{

struct PorousMedia : public IMedium
{
    NumLib::ITXFunction* hydraulic_conductivity;
    NumLib::ITXFunction* permeability;
    NumLib::ITXFunction* porosity;
    NumLib::ITXFunction* storage;
    NumLib::ITXFunction* geo_area;
    NumLib::ITXFunction* dispersivity_long;
    NumLib::ITXFunction* dispersivity_trans; 
	double res_saturation; 
    double max_saturation;      
    //double exp_saturation;
	double Pb;
	double saturation_exponent;
    NumLib::FunctionLinear1D* capp_sat_curve; 
    std::vector<size_t> perm_saturation_model;    
    std::vector<NumLib::FunctionLinear1D*> perm_saturation_curve; 
    size_t capp_sat_model;          
	double minimum_relative_permeability;  

    PorousMedia()
    {
        BaseLib::zeroObject(
                hydraulic_conductivity,
                permeability,
                porosity,
                storage,
                geo_area
               // res_saturation, 
                //max_saturation, 
				//exp_saturation
                );
        BaseLib::zeroObject(
                capp_sat_curve,
                dispersivity_long,
                dispersivity_trans
                );
        capp_sat_model = 0;
        minimum_relative_permeability = .0;
    }

    virtual ~PorousMedia()
    {
        BaseLib::releaseObject(
                hydraulic_conductivity,
                permeability,
                porosity,
                storage,
                geo_area
               // res_saturation,
				//max_saturation,
                //exp_saturation
                );
        perm_saturation_model.clear();
        BaseLib::releaseObjectsInStdVector(perm_saturation_curve);
        BaseLib::releaseObject(capp_sat_curve);

    }

    virtual MediumType::type getMediumType() const {return MediumType::PorousMedium;};


	virtual double getKrelbySw(const double Sw/*wetting saturation*/, size_t idx_phase)
    {
        double kr = 0.0, Sr, Sm, ex, Sw2, se, m, epsilon;
		epsilon = std::numeric_limits<double>::epsilon();
        //bool phase_shift = false;
        size_t model; 
        model = this->perm_saturation_model[idx_phase];
		minimum_relative_permeability =  1.0e-9; // Could be set as default at a beter place

        switch(model)
        {
	     default: // Error occurs! 
		    ERR("ERROR in PermeabilitySaturationFunction(): Unrecognized relative permeability method.\n");
            exit(0);
	        break;
	    case 0: // CURVE
			this->perm_saturation_curve[idx_phase]->eval( Sw, kr ); 
            if( kr < minimum_relative_permeability )
	        kr = minimum_relative_permeability;
            break;
		case 4: // Van Genuchten
			{
			Sr = res_saturation;
			Sm = max_saturation;
			Sw2 = MathLib::MRange(Sr+epsilon, Sw, Sm);
			se = (Sw2-Sr)/(Sm-Sr);
			m   = saturation_exponent;
			kr = sqrt(se) * pow(1.0-pow(1.0-pow(se,1.0/m),m),2);
			 if( kr < minimum_relative_permeability )
	        kr = minimum_relative_permeability;
			}
			break;
		case 6: // Brooks-Corey
			{
			//double Sr, Sm, ex, Sw2, epsilon;
			Sr = res_saturation;
			Sm = max_saturation;
			ex = (2+3*saturation_exponent)/saturation_exponent;
			//epsilon = std::numeric_limits<double>::epsilon();
			Sw2 = MathLib::MRange(Sr+epsilon, Sw, Sm); // because Sw is defined and used as const double in case 0, Sw2 as double is needed, WH
			kr = pow((Sw2-Sr)/(Sm-Sr),ex);
			 if( kr < minimum_relative_permeability )
	        kr = minimum_relative_permeability;
			}
			break;
        }
        return kr; 
    }

    /**
    * return the water saturation value by capillary pressure
    */
    virtual double getSwbyPc(double Pc, double Density)
    {	
        double Sw = 0.0, lamda, Sm, Sr, epsilon, m, se, Pb1;
		epsilon = std::numeric_limits<double>::epsilon();
        switch ( this->capp_sat_model )
        {
        case 0:  // curve value
            // get the value from curve. 
            capp_sat_curve->eval( Pc, Sw );
            break;
		case 4: // Van Genuchten
			lamda = saturation_exponent;
			Sr = res_saturation;
			Sm = max_saturation;
			Pb1 = Density*9.81/Pb;
			m   = saturation_exponent;
			if(Pc < 0.0) Pc = 0.0;
			se = pow(Pc/Pb1, 1.0/(1.0-m)) + 1.0;
			se = pow(se,-m);
			Sw = se*(Sm-Sr) + Sr;			
			Sw = MathLib::MRange (Sr+epsilon, Sw, Sm);
			break;
		case 6: // Brooks-corey
			//double lamda, Sm, Sr, epsilon;
			lamda = saturation_exponent;
			Sr = res_saturation;
			Sm = max_saturation;
			if (Pc < Pb) Pc = Pb;
			Sw = pow(Pb/Pc,lamda)*(Sm-Sr)+Sr;
			//epsilon = std::numeric_limits<double>::epsilon();
			Sw = MathLib::MRange (Sr+epsilon, Sw, Sm);
			break;
        default: 
            ERR("No valid capilary pressure vs water saturation model! "); 
            exit(0); 
            break;
        } // end of switch capp_sat_model
        return Sw; 
    }

	virtual double PorousMedia::getdSwdPc (double Pc, double Sw, double Density)
    {
        double dSwdPc = 0.0, lamda, Sm, Sr, v1, v2, m, epsilon, Pb1;
		epsilon = std::numeric_limits<double>::epsilon();

        switch ( this->capp_sat_model )
        {
        case 0:  // curve value
            // get the value from curve. 
            capp_sat_curve->eval_slope( Pc, dSwdPc, Sw);
            break;
		case 4: // Van Genuchten
			Sr = res_saturation;
			Sm = max_saturation;
			Pb1 = Density*9.81/Pb;
			m = saturation_exponent;
			if(Pc < 0.0) Pc = epsilon;
			v1 = pow((Pc/Pb1),(1.0/(1.0-m)));
			v2 = pow((1.0+v1),(-1.0-m));
			dSwdPc = (m*v1*v2*(Sm-Sr)) / ((m-1.0)*Pc);
			break;
		case 6: // Brooks-corey
			lamda = saturation_exponent;
			Sr = res_saturation;
			Sm = max_saturation;
			v1 = pow((Pc/Pb),-lamda);
			if (Pc < Pb) Pc = Pb;
			dSwdPc = (lamda*v1*(Sr - Sm)) / Pc; 
			break;
        default: 
            ERR("Error in getSwbyPc: No valid capilary pressure vs water saturation model! "); 
            exit(0); 
            break;
        } // end of switch capp_sat_model
        return dSwdPc; 
    }


};

} //end
 
