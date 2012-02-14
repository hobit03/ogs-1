
#pragma once

#include "TimeStep.h"

namespace NumLib
{

/**
 * \brief Interface class of transient coupling problems
 */
class ITransientProblem
{
public:
    virtual double suggestNext(const TimeStep &time_current) = 0;
    virtual bool isAwake(const TimeStep &time) = 0;
    virtual int solveTimeStep(const TimeStep &time) = 0;

    virtual int solve(const TimeStep &time) = 0;

    const TimeStep& getCurrentTime() const {return _current_time;};
    void setCurrentTime(const TimeStep &t) {_current_time = t;};
private:
    //int solve()
    //{
    //    return solveTimeStep(_current_time);
    //}
    TimeStep _current_time;
};




}
