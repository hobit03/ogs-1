
#include "SimulationInfo.h"

#include <cstdio>
#include <iostream>

#include "logog/include/logog.hpp"

#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

#include "Configure.h"
#include "BuildInfo.h"

namespace ogs6
{

void SimulationInfo::outputHeader ( void )
{
//    char buf[128];
//
//    const int len = 47;

    INFO("");
    INFO("          ###################################################");
    INFO("          ##                                               ##");
    INFO("          ##              OpenGeoSys-Project 6             ##");
#ifdef USE_LIS
    INFO("          ## %s ##", BaseLib::bothPadding("powered by LIS",45).c_str());
#endif
    INFO("          ##                                               ##");
    INFO("          ##   Contributors                                ##");
    INFO("          ##   * Helmholtz Centre for Environmental        ##");
    INFO("          ##     Research - UFZ                            ##");
    INFO("          ##   * TU Dresden                                ##");
    INFO("          ##   * University of Kiel                        ##");
    INFO("          ##   * University of Edinburgh                   ##");
    INFO("          ##   * University of Tuebingen (ZAG)             ##");
    INFO("          ##   * Federal Institute for Geosciences         ##");
    INFO("          ##     and Natural Resources (BGR)               ##");
    INFO("          ##   * Helmholtz Centre Potsdam GFZ              ##");
    INFO("          ##     German Research Centre for Geosciences    ##");
    INFO("          ##                                               ##");
    INFO("          ##   Program version                             ##");
    INFO("          ##   * Version: %s ##", BaseLib::rightPadding(OGS_VERSION, 32).c_str());
    INFO("          ##   * Date   : %s ##", BaseLib::rightPadding(OGS_DATE, 32).c_str());
    INFO("          ##   * Rev.   :                                  ##");
    INFO("          ##     %s ##", BaseLib::rightPadding(GIT_COMMIT_INFO, 41).c_str());
    INFO("          ##                                               ##");
    INFO("          ###################################################");
    INFO("");
    INFO("");
}

SimulationInfo::SimulationInfo(const std::string &project_path)
{
    this->setProjectPath(project_path);
}


void SimulationInfo::setProjectPath(const std::string& path)
{
    _project_path = path;
    _project_dir = BaseLib::getFileDirecotryPath(path);
    _project_name = BaseLib::getFileBaseName(path);
};

} //end
