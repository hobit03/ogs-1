
#pragma once

#include <string>
#include <vector>

#include "BaseLib/CodingTools.h"
#include "BaseLib/OrderedMap.h"
#include "MeshLib/Core/IMesh.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"

#include "IOutputTiming.h"

/**
 * 
 */
struct OutputVariableInfo
{
    enum DataType
    {
        Char,
        Int,
        Real
    };
    enum ObjectType
    {
        Node,
        Element
    };

    std::string name;
    ObjectType object_type;
    DataType data_type;
    size_t nr_of_components;
    NumLib::ITXFunction* value;

    OutputVariableInfo() 
        : name(""), object_type(Node), data_type(Char), nr_of_components(1), value(NULL)
    {};

    OutputVariableInfo(const std::string &s, ObjectType o, DataType d, size_t n, NumLib::ITXFunction* v) 
    : name(s), object_type(o), data_type(d), nr_of_components(n), value(v)
    {
    };
};

/**
 * \brief Interface of Output class
 */
class IOutput
{
public:
    ///
    IOutput() : _output_timing(0), _geo_obj(0) {};

    ///
    virtual ~IOutput()
    {
        BaseLib::releaseObject(_output_timing);
    };

    ///
    void setOutputPath(const std::string &dir, const std::string &base_name)
    {
        _output_dir_path = dir;
        _base_name = base_name;
    }

    ///
    std::string getOutputDir() const {return _output_dir_path;};
    
    ///
    std::string getOutputBaseName() const {return _base_name;};

    ///
    void addVariable(const std::string &var_name)
    {
        _list_var_name.push_back(var_name);
    };

    ///
    void addVariable(const std::vector<std::string> &vec_var_name)
    {
        _list_var_name.insert(_list_var_name.end(), vec_var_name.begin(), vec_var_name.end());
    };
    
    ///
    std::vector<std::string>& getListOfVariables() {return _list_var_name;};
    
    ///
    bool hasVariable(const std::string &var_name) const
    {
        return (std::find(_list_var_name.begin(), _list_var_name.end(), var_name) != _list_var_name.end());
    }

    ///
    void setGeometry(const GeoLib::GeoObject* geo_obj) {_geo_obj = geo_obj;};
    ///
    const GeoLib::GeoObject* getGeometry() const {return _geo_obj;};

    void setMesh(MeshLib::IMesh* msh) {_msh = msh;};
    MeshLib::IMesh* getMesh() const {return _msh;};

    ///
    void setOutputTiming(IOutputTiming* timing) {_output_timing = timing;};
    ///
    bool isActive(const NumLib::TimeStep &current_time) const {return _output_timing->isActive(current_time);};

    ///
    virtual void write( const NumLib::TimeStep &current_time, 
                        BaseLib::OrderedMap<std::string, OutputVariableInfo> &data) = 0;

private:
    std::string _output_dir_path;
    std::string _base_name;
    IOutputTiming* _output_timing;
    std::vector<std::string> _list_var_name;
    const GeoLib::GeoObject* _geo_obj;
    MeshLib::IMesh* _msh;
};

