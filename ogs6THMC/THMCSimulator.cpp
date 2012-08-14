/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file THMCSimulator.cpp
 *
 * Created on 2012-07-17 by Norihiro Watanabe
 */

#include "THMCSimulator.h"

#include <iostream>

// external library
#include "logog.hpp"
#include "tclap/CmdLine.h"
#ifdef USE_LIS
#include "lis.h"
#endif

// internal library
#include "BaseLib/CodingTools.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Options.h"
#include "BaseLib/OptionsXMLReader.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/RunTime.h"
#include "NumLib/TransientCoupling/TransientCouplingStructureBuilder.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "FemIO/ogs5/Ogs5FemIO.h"

// this module
#include "SimulationInfo.h"
#include "FormatterCustom.h"
#include "Ogs6FemData.h"
#include "Ogs5ToOgs6.h"
#include "MyConvergenceCheckerFactory.h"
#include "TimeSteppingControllerWithOutput.h"

namespace ogs6
{

////////////////////////////////////////////////////////////////////////////////
// Variables
////////////////////////////////////////////////////////////////////////////////
static logog::Cout *logogCout;
static logog::LogFile *logog_file;
static FormatterCustom *custom_format;
static bool isOgsInitCalled = false;
static bool isOgsExitCalled = false;
////////////////////////////////////////////////////////////////////////////////

void ogsInit(int argc, char* argv[])
{
    if (isOgsInitCalled) return;
    isOgsInitCalled = true;

    LOGOG_INITIALIZE();
    custom_format = new FormatterCustom();
    logogCout = new logog::Cout();
    logogCout->SetFormatter(*custom_format);
    logog_file = NULL;

#ifdef USE_LIS
    lis_initialize((LIS_INT*)&argc, &argv);
#endif
}

void ogsExit()
{
    if (isOgsExitCalled) return;
    isOgsExitCalled = true;

#ifdef USE_LIS
    lis_finalize();
#endif

    INFO("exit ogs6.");
    BaseLib::releaseObject(custom_format, logogCout, logog_file);
    LOGOG_SHUTDOWN();
}

THMCSimulator::THMCSimulator(int argc, char* argv[])
: _sim_info(NULL), _cpl_system(NULL)
{
    try {
        // Command line parser
        TCLAP::CmdLine cmd("ogs6", ' ', "0.1");
        TCLAP::ValueArg<std::string> input_arg("i", "input", "input file", false, "", "string");
        cmd.add( input_arg );
        TCLAP::ValueArg<std::string> output_dir_arg("o", "output", "output directory", false, "", "string");
        cmd.add( output_dir_arg );
        TCLAP::ValueArg<std::string> logfile_arg("l", "log", "log file", false, "", "string");
        cmd.add( logfile_arg );
        TCLAP::ValueArg<unsigned> verbosity_arg("v", "verbose", "level of verbosity [0 very low information, 1 much information]", false, 0, "number");
        cmd.add( verbosity_arg );
        TCLAP::ValueArg<unsigned> pcs_arg("m", "modules", "list available modules [0 off, 1 on]", false, 1, "number");
        cmd.add( pcs_arg );
        cmd.parse( argc, argv ); // process can exit in this function

        // initialize
        ogsInit(argc, argv);

        // get parsed data
        // log file
        if (! logfile_arg.getValue().empty()) {
            if (!logog_file) delete logog_file;
            std::string log_file = logfile_arg.getValue();
            BaseLib::truncateFile(log_file); // do this not to append log into an existing file
            logog_file = new logog::LogFile(log_file.c_str());
            logog_file->SetFormatter( *custom_format );
        }


        SimulationInfo::outputHeader();
        // list modules
        unsigned flag_list_modules (pcs_arg.getValue());
        if (flag_list_modules!=0) {
            ProcessBuilder::getInstance()->output();
        }

        INFO("->Parsing input arguments");
        if (!input_arg.getValue().empty()) {
            INFO("* Project path    : %s", input_arg.getValue().c_str());
        }
        if (! logfile_arg.getValue().empty()) {
            INFO("* Log file path   : %s", logfile_arg.getValue().c_str());
        }
        
        // data output directory
        std::string output_dir_path = "";
        if (! output_dir_arg.getValue().empty()) {
            output_dir_path = output_dir_arg.getValue();
            INFO("* Output directory: %s", output_dir_path.c_str());
        }

        if (! input_arg.getValue().empty()) {
            const std::string proj_path = input_arg.getValue();
            if (checkInputFiles(proj_path)) {
                _sim_info = new SimulationInfo(proj_path, output_dir_path);
            } else {
                ERR("***Error: Cannot find a project - %s", proj_path.c_str());
            }
        }

    } catch (TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }

}

THMCSimulator::~THMCSimulator()
{
    BaseLib::releaseObject(_sim_info, _cpl_system);
    ogsExit();
}

bool THMCSimulator::checkInputFiles(const std::string& proj_path)
{
    // meanwhile OGS5 files are default
    if(!BaseLib::IsFileExisting(proj_path+".pcs"))
    {
        ERR("Cannot find a PCS file - %s.pcs", proj_path.c_str());
        return false;
    }
    if(!BaseLib::IsFileExisting(proj_path+ ".prop"))
    {
        ERR("Cannot find a property file - %s.prop", proj_path.c_str());
        return false;
    }

    return true;
}

int THMCSimulator::execute()
{
    if (!_sim_info) return 0;

    BaseLib::Options op;
    Ogs6FemData* ogs6fem = Ogs6FemData::getInstance();
    const std::string proj_path = _sim_info->getProjectPath();
    ogs6fem->project_name = _sim_info->getProjectName();
    ogs6fem->output_dir = _sim_info->getOutputDirPath();

    //-------------------------------------------------------------------------
    // Read files
    //-------------------------------------------------------------------------
    INFO("->Reading input files...");
    // ogs5fem
    ogs5::Ogs5FemData ogs5femdata;
    ogs5::Ogs5FemIO::read(proj_path, ogs5femdata);
    Ogs5ToOgs6::convert(ogs5femdata, *ogs6fem, op);

    // coupling
    BaseLib::addXMLtoOptions(proj_path+".prop", op);

    // ddc

    //-------------------------------------------------------------------------
    // Setup simulation
    //-------------------------------------------------------------------------
    //INFO("Setting up simulation...");

    // construct mesh
    INFO("->Constructing meshes... %d mesh loaded", ogs6fem->list_dis_sys.size());
    for (size_t i=0; i<ogs6fem->list_dis_sys.size(); i++) {
        MeshLib::IMesh* msh = ogs6fem->list_dis_sys[i]->getMesh();
        msh->constructGeometricProperty();
        INFO("->mesh id %d: dim=%d, nodes=%d, elements=%d", i, msh->getDimension(), msh->getNumberOfNodes(), msh->getNumberOfElements());
    }

    // construct coupling system
    INFO("->Generating coupling system...");
    MyConvergenceCheckerFactory checkFac;
    NumLib::TransientCoulplingStrucutreBuilder cpl_builder;
    if (_cpl_system!=NULL) delete _cpl_system;
    _cpl_system = cpl_builder.build(&op, *GeoProcessBuilder::getInstance(), checkFac);

    if (!_cpl_system->check()) {
        ERR("***Error while checking coupled system");
        return 0;
    }

    // list up monolithic processes
    INFO("->Initializing all processes...");
    std::vector<NumLib::AbstractTransientMonolithicSystem*> &list_mono_system = cpl_builder.getListOfMonolithicSystem();
    std::vector<std::string> &list_mono_system_name = cpl_builder.getListOfMonolithicSystemName();
    for (size_t i=0; i<list_mono_system.size(); i++) {
        std::string &pcs_name = list_mono_system_name[i];
        ProcessLib::Process* pcs = list_mono_system[i];
        ogs6fem->list_pcs.insert(pcs_name, pcs);
        const BaseLib::Options* opPCS = op.getSubGroup("ProcessData")->getSubGroup(pcs_name);
        bool isPcsReady = pcs->initialize(opPCS!=NULL ? *opPCS : op);
        if (!isPcsReady) {
            ERR("***Error while setting up processes");
            return 0;
        }
    }

    INFO("->Setting time stepping...");
    TimeSteppingControllerWithOutput timestepping(&ogs6fem->outController);
    timestepping.setTransientSystem(*_cpl_system);

    //TODO the following calculation should be done in TimeSteppingController
    double t_start = std::numeric_limits<double>::max();
    double t_end = -1 * std::numeric_limits<double>::max();
    for (size_t i=0; i<ogs6fem->list_tim.size(); i++) {
        t_start = std::min(t_start, ogs6fem->list_tim[i]->getBeginning());
        t_end = std::max(t_end, ogs6fem->list_tim[i]->getEnd());
    }

    INFO("->Outputting the initial values...");
    ogs6fem->outController.outputData(NumLib::TimeStep(t_start));

    //-------------------------------------------------------------------------
    // Run simulation
    //-------------------------------------------------------------------------
    INFO("->Simulation is ready! start=%f, end=%f", t_start, t_end);
    BaseLib::RunTime runTime;
    runTime.start();

    timestepping.setBeginning(t_start); //TODO really need this? start, end is already given in timestep function
    size_t n_timesteps = timestepping.solve(t_end);

    runTime.stop();

    INFO("");
    INFO("->Simulation is finished.\n");
    INFO("#############################################################");
    INFO("*** Summary of this simulation");
    INFO("total time step : %d", n_timesteps);
    INFO("elapsed time   : %g sec", runTime.elapsed());
    INFO("#############################################################");
    INFO("");

    return 0;
}

} //end ogs6
