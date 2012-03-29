
#pragma once

#include <valarray>
#include <cmath>

namespace FemLib
{
/**
 * \brief Template class for FEM integration point-based functions
 */
template<typename Tvalue>
class TemplateFEMIntegrationPointFunction : public MathLib::TemplateFunction<GeoLib::Point,Tvalue>
{
public:
    typedef DiscreteLib::DiscreteVector<std::valarray<Tvalue> > DiscreteVectorType;

    TemplateFEMIntegrationPointFunction(DiscreteLib::DiscreteSystem &dis, MeshLib::IMesh &msh)
    {
        initialize(dis, msh, msh.getNumberOfElements());
    };

    TemplateFEMIntegrationPointFunction(const TemplateFEMIntegrationPointFunction &src)
    {
        initialize(*src._discrete_system, *src._msh, src._values->size());
        (*this->_values) = (*src._values);
    };

    MathLib::TemplateFunction<GeoLib::Point, Tvalue>* clone() const
    {
        TemplateFEMIntegrationPointFunction<Tvalue> *obj = new TemplateFEMIntegrationPointFunction<Tvalue>(*this);
        return obj;
    };

    const MeshLib::IMesh* getMesh() const
    {
        return _msh;
    }

    void eval(const GeoLib::Point &pt, Tvalue &v)
    {
        throw "The method or operation is not implemented.";
    };

    void setIntegrationPointValue( size_t i_e, size_t ip, Tvalue &q )
    {
        assert(ip<(*_values)[i_e].size());
        (*_values)[i_e][ip] = q;
    }

    void setNumberOfIntegationPoints(size_t i_e, size_t n)
    {
        (*_values)[i_e].resize(n);
    }

    const std::valarray<Tvalue>& getIntegrationPointValues(size_t i_e) const
    {
        return (*_values)[i_e];
    }

    double norm_diff(const TemplateFEMIntegrationPointFunction<Tvalue> &ref) const
    {
    	const size_t n = _values->size();
    	if (n!=ref._values->size()) {
    		std::cout << "***Warning in TemplateFEMIntegrationPointFunction::norm_diff(): size of two vectors is not same." << std::endl;
    		return .0;
    	}

    	double mnorm = .0;
    	for (size_t i=0; i<n; ++i) {
    		const std::valarray<Tvalue> &val1 = (*_values)[i];
    		const std::valarray<Tvalue> &val2 = (*ref._values)[i];
    		const size_t n_gp = val1.size();
        	if (n_gp!=val2.size()) {
        		std::cout << "***Warning in TemplateFEMIntegrationPointFunction::norm_diff(): size of two vectors is not same." << std::endl;
        		return .0;
        	} else if (n_gp==0) {
                std::cout << "***Warning in TemplateFEMIntegrationPointFunction::norm_diff(): size of two vectors is zero." << std::endl;
                return .0;
            }
        	std::valarray<Tvalue> val_diff = val1 - val2;

        	val_diff = std::abs(val_diff);
			mnorm = std::max(mnorm, val_diff.max());
    	}

    	return mnorm;
    }

    const DiscreteVectorType* getNodalValues() const
    {
        return _values;
    }

    void printouf() const
    {
    	std::cout << "integration_pt_values = ";
    	for (size_t i=_values->getRangeBegin(); i<_values->getRangeEnd(); ++i) {
    		const std::valarray<Tvalue> &val1 = (*_values)[i];
    		std::cout << "(";
        	for (size_t j=0; j<val1.size(); ++j) std::cout << val1[j] << " ";
    		std::cout << ") ";
    	}
    	std::cout << std::endl;
    }

private:
    MeshLib::IMesh* _msh;
    //std::vector<std::vector<Tvalue> >* _values;
    DiscreteVectorType* _values;
    DiscreteLib::DiscreteSystem* _discrete_system;

    void initialize(DiscreteLib::DiscreteSystem &dis, MeshLib::IMesh &msh, size_t n)
    {
        _msh = &msh;
        _discrete_system = &dis;
//        _values = new std::vector<std::vector<Tvalue> >(n);
        _values = _discrete_system->createVector<std::valarray<Tvalue> >(n);
    }


};

typedef TemplateFEMIntegrationPointFunction<double> FEMIntegrationPointFunctionScalar2d;
typedef TemplateFEMIntegrationPointFunction<MathLib::Vector2D> FEMIntegrationPointFunctionVector2d;

} //end
