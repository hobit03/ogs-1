
#pragma once

#include "BaseLib/Options.h"
#include "TimeStep.h"

namespace NumLib
{

/**
 * \brief Interface of any transient systems to be solved
 *
 * This class defines interface of any transient systems which work under time stepping concept.
 * Purpose of this class is to communicate with time stepping controllers regarding,
 * - to inform in which time step this system is active
 * - to update state of the system, e.g. solve
 */
class ITransientSystem
{
public:
    /// 
    virtual ~ITransientSystem() {};

    /// initialize
    virtual bool initialize(const BaseLib::Options &/*op*/) {return true;};
    /// finalize
    virtual void finalize() {};

    /// suggest the next time step
    virtual double suggestNext(const TimeStep &time_current) = 0;
    /// return if this system is awake with the given time
    virtual bool isAwake(const TimeStep &time) = 0;
    /// solve this system with the given time step
    virtual int solveTimeStep(const TimeStep &time) = 0;
    /// accept the given time step
    virtual void accept(const TimeStep &time) = 0;

    /// return current time
    const TimeStep& getCurrentTime() const {return *_current_time;};
    /// set current time
    void setCurrentTime(const TimeStep &t) {_current_time = const_cast<TimeStep*>(&t);};

private:
    TimeStep* _current_time;
};




}
