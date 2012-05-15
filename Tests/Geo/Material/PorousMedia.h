
#pragma once

#include "Base/CodingTools.h"

#include "Solid.h"

namespace NumLib
{
class ITXFunction;
}

namespace Geo
{

struct PorousMedia
{
	NumLib::ITXFunction* hydraulic_conductivity;
	NumLib::ITXFunction* porosity;
	NumLib::ITXFunction* storage;
	Solid* solidphase;

	PorousMedia()
	{
		Base::zeroObject(hydraulic_conductivity, porosity, storage, solidphase);
	}
	~PorousMedia()
	{
		Base::releaseObject(hydraulic_conductivity, porosity, storage);
	}
};

} //end
 
