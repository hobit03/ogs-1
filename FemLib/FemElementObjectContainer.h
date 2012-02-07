
#pragma once

#include <map>
#include "Base/MemoryTools.h"
#include "FemLib/Core/IFemElement.h"
#include "FemLib/Core/Element/FemElementFactory.h"

namespace FemLib
{

class FeObjectCachePerFeType
{
public:
    virtual ~FeObjectCachePerFeType()
    {
        Base::destroyStdMapWithPointers(_mapFeObj);
    }

    IFiniteElement* getFeObject(FiniteElementType::type fe_type, MeshLib::IMesh* msh)
    {
        IFiniteElement *fe = 0;
        if (_mapFeObj.count(fe_type)==0) {
            fe = FemElementFactory::create(fe_type, msh);
            _mapFeObj[fe_type] = fe;
        } else {
            fe = _mapFeObj[fe_type];
        }
        return fe;
    }

private:
    std::map<FiniteElementType::type, IFiniteElement*> _mapFeObj;
};


class IFeObjectContainer
{
public:
    virtual IFiniteElement* getFeObject(const MeshLib::IElement &e, MeshLib::IMesh *msh) = 0;
};

class FeObjectContainerPerElement : public IFeObjectContainer
{
public:
    FeObjectContainerPerElement(size_t ele_size)
    {
        _vec_fem.resize(ele_size);
    }
    virtual ~FeObjectContainerPerElement()
    {
        Base::destroyStdVectorWithPointers(_vec_fem);
    }

    void addFiniteElement(size_t i, FiniteElementType::type fe_type, MeshLib::IMesh* msh)
    {
        _vec_fem[i] = FemElementFactory::create(fe_type, msh);
    }

    virtual IFiniteElement* getFeObject(const MeshLib::IElement &e, MeshLib::IMesh *msh) 
    {
        return _vec_fem[e.getID()];
    }

private:
    std::vector<IFiniteElement*> _vec_fem;
};



class LagrangianFeObjectContainer : public FeObjectCachePerFeType, public IFeObjectContainer
{
public:
    LagrangianFeObjectContainer() 
    {
        _order = 1;
    }

    void setPolynomialOrder(size_t order) 
    {
        _order = order;
    }

    virtual IFiniteElement* getFeObject(const MeshLib::IElement &e, MeshLib::IMesh *msh)
    {
        FiniteElementType::type fe_type = getFeType(e.getElementType(), _order);
        return FeObjectCachePerFeType::getFeObject(fe_type, msh);
    }
private:
    size_t _order;

    FiniteElementType::type getFeType(MeshLib::ElementType::type ele_type, size_t order)
    {
        switch (ele_type)
        {
            case MeshLib::ElementType::LINE:
                return (order==1) ? FiniteElementType::LINE2 : FiniteElementType::LINE3;
            case MeshLib::ElementType::QUAD:
                return (order==1) ? FiniteElementType::QUAD4 : FiniteElementType::QUAD9;
        }
        return FiniteElementType::INVALID;
    };
};



}
