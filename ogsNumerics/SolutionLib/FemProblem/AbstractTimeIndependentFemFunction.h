
#pragma once

#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"

namespace SolutionLib
{

template <
	size_t N_IN_PARAMETER,
	size_t N_OUT_PARAMETER
	>
class AbstractTimeIndependentFemFunction
: public NumLib::TemplateTransientMonolithicSystem
{
public:
	AbstractTimeIndependentFemFunction()
    {
        TemplateTransientMonolithicSystem::resizeInputParameter(N_IN_PARAMETER);
        TemplateTransientMonolithicSystem::resizeOutputParameter(N_OUT_PARAMETER);
    }

    double suggestNext(const NumLib::TimeStep &/*time_current*/) { return .0; }

    bool isAwake(const NumLib::TimeStep &/*time*/) { return true;  }

    virtual void accept(const NumLib::TimeStep &/*time*/) {};

};

}
