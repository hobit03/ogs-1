/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Ogs5FemData.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <vector>

#include "GeoLib/GEOObjects.h"
#include "MeshLib/Core/IMesh.h"
#include "rf_pcs.h"
#include "rf_mfp_new.h"
#include "rfmat_cp.h"
#include "rf_bc_new.h"
#include "rf_st_new.h"
#include "rf_ic_new.h"
#include "rf_out_new.h"
#include "rf_tim_new.h"
#include "rf_msp_new.h"
#include "rf_mmp_new.h"
#include "rf_num_new.h"
#include "Output.h"

namespace ogs5
{

class Ogs5FemData
{

public:
    Ogs5FemData() {};
    ~Ogs5FemData();

public:
    std::vector<CRFProcess*> pcs_vector;
    std::vector<CFluidProperties*> mfp_vector;
    std::vector<CompProperties*> cp_vector;
    std::vector<CBoundaryCondition*> bc_vector;
    std::vector<CSourceTerm*> st_vector;
    std::vector<CInitialCondition*> ic_vector;
    std::vector<COutput*> out_vector;
    std::vector<CTimeDiscretization*> time_vector;
    std::vector<CSolidProperties*> msp_vector;
    std::vector<CMediumProperties*> mmp_vector;
    std::vector<CNumerics*>num_vector;
    GeoLib::GEOObjects* geo_obj;
    std::string geo_unique_name;
    std::vector<MeshLib::IMesh*> list_mesh;
};

}
