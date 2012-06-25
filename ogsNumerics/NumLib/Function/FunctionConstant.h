
#pragma once

#include "IFunction.h"

namespace NumLib
{

template<typename Tpos, typename Tval>
class FunctionConstant : public TemplateFunction<Tpos,Tval>
{
public:
    FunctionConstant(const Tval &v)
    {
        _v = v;
    };

    virtual ~FunctionConstant() {};

    virtual void eval(const Tpos&, Tval &v)
    {
        v = _v;
    };

    virtual void eval(const Tpos&, Tval &v) const
    {
        v = _v;
    };

    virtual void eval(Tval &v) const
    {
        v = _v;
    };

    virtual TemplateFunction<Tpos,Tval>* clone() const
    {
    	return new FunctionConstant(_v);
    }
private:
    Tval _v;
};

} //end
