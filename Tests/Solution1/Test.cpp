
#include <cmath>
#include <gtest/gtest.h>

#include "Base/CodingTools.h"

#include "MathLib/Vector.h"
#include "MathLib/LinAlg/LinearEquations/LisInterface.h"
#include "GeoLib/Shape/Line.h"
#include "MeshLib/Tools/MeshGenerator.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TransientCoupling/AsyncPartitionedSystem.h"
#include "NumLib/Nonlinear/TemplateDiscreteNonlinearSolver.h"
#include "NumLib/Coupling/Algorithm/TransientPartitionedAlgorithm.h"
#include "FemLib/Function/FemNorm.h"
#include "Tests/Geo/Model/Head.h"
#include "Tests/Geo/Model/Velocity.h"
#include "Tests/Geo/Model/Concentration.h"
#include "Tests/Geo/Model/Displacement.h"

#include "Tests/ExactSolution/OgataBank.h"

#include "TestUtil.h"


class FemFunctionConvergenceCheck : public IConvergenceCheck
{
	FemLib::NormOfFemNodalFunction<double> _norm;
public:
	explicit FemFunctionConvergenceCheck(DiscreteLib::DiscreteSystem *dis) : _norm(dis)
	{

	}

	bool isConverged(UnnamedParameterSet& vars_prev, UnnamedParameterSet& vars_current, double eps, double &v_diff)
	{
	    for (size_t i=0; i<vars_prev.size(); i++) {
#if 1
	    	if (vars_prev.getName(i).compare("h")==0 || vars_prev.getName(i).compare("c")==0) {
		        const FemNodalFunctionScalar* f_fem_prev = vars_prev.get<FemNodalFunctionScalar>(i);
		        const FemNodalFunctionScalar* f_fem_cur = vars_current.get<FemNodalFunctionScalar>(i);
	    		//v_diff = f_fem_cur->norm_diff(*f_fem_prev);
	    		v_diff = _norm(*f_fem_prev, *f_fem_cur);
	    	} else if (vars_prev.getName(i).compare("v")==0) {
	    		const FEMIntegrationPointFunctionVector* f_fem_prev = vars_prev.get<FEMIntegrationPointFunctionVector>(i);
	    		const FEMIntegrationPointFunctionVector* f_fem_cur = vars_current.get<FEMIntegrationPointFunctionVector>(i);
	    		//v_diff = f_fem_cur->norm_diff(*f_fem_prev);
	    		v_diff = _norm(*f_fem_prev, *f_fem_cur);
	    	}
#endif
	        if (v_diff>eps) {
	            return false;
	        }
	    }
	    return true;
	}
};

Geo::GWFemProblem* createGWProblem(DiscreteSystem &dis, Geo::PorousMedia &pm)
{
    LagrangianFeObjectContainer* _feObjects = new LagrangianFeObjectContainer(*dis.getMesh());
    //equations
    Geo::GWFemEquation::LinearAssemblerType* linear_assembler = new Geo::GWFemEquation::LinearAssemblerType(*_feObjects, pm);
    Geo::GWFemEquation::ResidualAssemblerType* r_assembler = new Geo::GWFemEquation::ResidualAssemblerType(*_feObjects, pm);
    Geo::GWFemEquation::JacobianAssemblerType* j_eqs = new Geo::GWFemEquation::JacobianAssemblerType(*_feObjects, pm);
    Geo::GWFemEquation* eqs = new  Geo::GWFemEquation(linear_assembler, r_assembler, j_eqs);
    //IVBV problem
    Geo::GWFemProblem* _problem = new Geo::GWFemProblem(&dis);
    _problem->setEquation(eqs);
    return _problem;
}


Geo::GWFemProblem* defineGWProblem(DiscreteSystem &dis, Rectangle &_rec, Geo::PorousMedia &pm)
{
    Geo::GWFemProblem* _problem = createGWProblem(dis, pm);
    // var
    FemVariable* head = _problem->addVariable("head");
    // IC
    FemNodalFunctionScalar* h0 = new FemNodalFunctionScalar(dis, PolynomialOrder::Linear, 0);
    head->setIC(h0);
    //BC
    Polyline* poly_left = _rec.getLeft();
    Polyline* poly_right = _rec.getRight();
    head->addDirichletBC(new FemDirichletBC(dis.getMesh(), poly_right, new NumLib::TXFunctionConstant(.0)));
    head->addNeumannBC(new FemNeumannBC(dis.getMesh(), h0->getFeObjectContainer(), poly_left, new NumLib::TXFunctionConstant(-1e-5)));

    return _problem;
}

Geo::GWFemProblem* defineGWProblem1D(DiscreteSystem &dis, GeoLib::Line &line, Geo::PorousMedia &pm)
{
    Geo::GWFemProblem* _problem = createGWProblem(dis, pm);
    // var
    FemVariable* head = _problem->addVariable("head");
    // IC
    FemNodalFunctionScalar* h0 = new FemNodalFunctionScalar(dis, PolynomialOrder::Linear, 0);
    head->setIC(h0);
    //BC
    head->addDirichletBC(new FemDirichletBC(dis.getMesh(), line.getPoint2(), new NumLib::TXFunctionConstant(.0)));
    head->addNeumannBC(new FemNeumannBC(dis.getMesh(), h0->getFeObjectContainer(), line.getPoint1(), new NumLib::TXFunctionConstant(-1e-5)));

    return _problem;
}

Geo::MassFemProblem* defineMassTransportProblem(DiscreteSystem &dis, Rectangle &_rec, Geo::PorousMedia &pm, Geo::Compound &comp)
{
    LagrangianFeObjectContainer* _feObjects = new LagrangianFeObjectContainer(*dis.getMesh());
    //equations
    Geo::MassFemEquation::LinearAssemblerType* linear_assembler = new Geo::MassFemEquation::LinearAssemblerType(*_feObjects, pm, comp);
    Geo::MassFemEquation::ResidualAssemblerType* r_assembler = new Geo::MassFemEquation::ResidualAssemblerType(*_feObjects, pm, comp);
    Geo::MassFemEquation::JacobianAssemblerType* j_eqs = new Geo::MassFemEquation::JacobianAssemblerType(*_feObjects, pm, comp);
    Geo::MassFemEquation* eqs = new  Geo::MassFemEquation(linear_assembler, r_assembler, j_eqs);
    //IVBV problem
    Geo::MassFemProblem* _problem = new Geo::MassFemProblem(&dis);
    _problem->setEquation(eqs);
    // var
    FemVariable* c = _problem->addVariable("c");
    // IC
    FemNodalFunctionScalar* c0 = new FemNodalFunctionScalar(dis, PolynomialOrder::Linear, 0);
    c->setIC(c0);
    //BC
    Polyline* poly_left = _rec.getLeft();
    c->addDirichletBC(new FemDirichletBC(dis.getMesh(), poly_left, new NumLib::TXFunctionConstant(1.0)));

    return _problem;
}


static void getGWExpectedHead(std::vector<double> &expected)
{
    expected.resize(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) expected[i] = 2.e+6;
        if (i%3==1) expected[i] = 1.e+6;
        if (i%3==2) expected[i] = 0.e+6;
    }
}

typedef Geo::FunctionHead<CRSLisSolver> MyFunctionHead;
typedef Geo::FunctionVelocity MyFunctionVelocity;
typedef Geo::FunctionConcentration<CRSLisSolver> MyFunctionConcentration;
//typedef Geo::FunctionDisplacement<CRSLisSolver> MyFunctionDisplacement;


TEST(Solution, CouplingFem2D)
{
	try {
	    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
	    Rectangle* _rec = new Rectangle(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
	    Geo::PorousMedia pm;
	    pm.hydraulic_conductivity = new NumLib::TXFunctionConstant(1.e-11);
	    pm.porosity = new NumLib::TXFunctionConstant(0.2);
	    DiscreteSystem dis(*msh);
	    Geo::GWFemProblem* pGW = defineGWProblem(dis, *_rec, pm);
        TimeStepFunctionConstant tim(.0, 100.0, 10.0);
        pGW->setTimeSteppingFunction(tim);
	    // options
	    Base::Options options;
	    Base::Options* op_lis = options.addSubGroup("Lis");
	    op_lis->addOption("solver_type", "CG");
	    op_lis->addOption("precon_type", "NONE");
	    op_lis->addOptionAsNum("error_tolerance", 1e-10);
	    op_lis->addOptionAsNum("max_iteration_step", 500);
	    Base::Options* op_nl = options.addSubGroup("Nonlinear");
	    op_nl->addOption("solver_type", "Picard");
	    op_nl->addOptionAsNum("error_tolerance", 1e-6);
	    op_nl->addOptionAsNum("max_iteration_step", 500);

		MyFunctionHead f_head;
		f_head.define(&dis, pGW, options);
        f_head.setOutputParameterName(0, "h");
		MyFunctionVelocity f_vel;
		f_vel.define(dis, pm);
        f_vel.setInputParameterName(0, "h");
        f_vel.setOutputParameterName(0, "v");

        FemFunctionConvergenceCheck checker(&dis);
	    SerialStaggeredMethod method(checker, 1e-5, 100);
	    AsyncPartitionedSystem apart1;
        apart1.setAlgorithm(method);
        apart1.resizeOutputParameter(2);
        apart1.setOutputParameterName(0, "h");
        apart1.setOutputParameterName(1, "v");
        apart1.addProblem(f_head);
        apart1.addProblem(f_vel);
        apart1.connectParameters();

	    TimeSteppingController timestepping;
	    timestepping.addTransientSystem(apart1);

	    //const double epsilon = 1.e-3;
	    timestepping.setBeginning(.0);
	    timestepping.solve(1.0);

	    const FemNodalFunctionScalar* r_f_head = apart1.getOutput<FemNodalFunctionScalar>(apart1.getOutputParameterID("h"));
	    const FEMIntegrationPointFunctionVector* r_f_v = apart1.getOutput<FEMIntegrationPointFunctionVector>(apart1.getOutputParameterID("v"));
	    const IDiscreteVector<double>* vec_h = r_f_head->getNodalValues();
	    //const FEMIntegrationPointFunctionVector2d::DiscreteVectorType* vec_v = r_f_v->getNodalValues();

	    r_f_head->printout();
	    r_f_v->printout();

	    std::vector<double> expected;
	    getGWExpectedHead(expected);

	    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], &(*vec_h)[0], vec_h->size());

	} catch (const char* e) {
		std::cout << "***Exception caught! " << e << std::endl;
	}

}

TEST(FEM, line)
{
    try {
        const double len = 2.0;
        const size_t div = 2;
        const double h = len / div;
        MeshLib::IMesh *msh = MeshGenerator::generateLineMesh(len, div, .0, .0, .0);
        GeoLib::Line* line = new GeoLib::Line(Point(0.0, 0.0, 0.0),  Point(2.0, 0.0, 0.0));
        Geo::PorousMedia pm;
	    pm.hydraulic_conductivity = new NumLib::TXFunctionConstant(1.e-11);
	    pm.porosity = new NumLib::TXFunctionConstant(0.2);
        DiscreteSystem dis(*msh);
        Geo::GWFemProblem* pGW = defineGWProblem1D(dis, *line, pm);
        TimeStepFunctionConstant tim(.0, 10.0, 10.0);
        pGW->setTimeSteppingFunction(tim);
        // options
        Base::Options options;
        Base::Options* op_lis = options.addSubGroup("Lis");
        op_lis->addOption("solver_type", "CG");
        op_lis->addOption("precon_type", "NONE");
        op_lis->addOptionAsNum("error_tolerance", 1e-10);
        op_lis->addOptionAsNum("max_iteration_step", 500);

        MyFunctionHead f_head;
        f_head.define(&dis, pGW, options);
        f_head.setOutputParameterName(0, "h");
        MyFunctionVelocity f_vel;
        f_vel.define(dis, pm);
        f_vel.setInputParameterName(0, "h");
        f_vel.setOutputParameterName(0, "v");

        FemFunctionConvergenceCheck checker(&dis);
        SerialStaggeredMethod method(checker, 1e-5, 100);
        AsyncPartitionedSystem apart1;
        apart1.setAlgorithm(method);
        apart1.resizeOutputParameter(2);
        apart1.setOutputParameterName(0, "h");
        apart1.setOutputParameterName(1, "v");
        apart1.addProblem(f_head);
        apart1.addProblem(f_vel);
        apart1.connectParameters();

        TimeSteppingController timestepping;
        timestepping.addTransientSystem(apart1);

        //const double epsilon = 1.e-3;
        timestepping.setBeginning(.0);
        timestepping.solve(1.0);

        const FemNodalFunctionScalar* r_f_head = apart1.getOutput<FemNodalFunctionScalar>(apart1.getOutputParameterID("h"));
        const FEMIntegrationPointFunctionVector* r_f_v = apart1.getOutput<FEMIntegrationPointFunctionVector>(apart1.getOutputParameterID("v"));
        const IDiscreteVector<double>* vec_h = r_f_head->getNodalValues();
        //const FEMIntegrationPointFunctionVector2d::DiscreteVectorType* vec_v = r_f_v->getNodalValues();

        r_f_head->printout();
        r_f_v->printout();

        std::vector<double> expected;
        expected.resize(div+1);
        const double p_left = 2.e+6;
        const double p_right = .0;
        const double x_len = 2.0;
        for (size_t i=0; i<expected.size(); i++) {
            double x = i*h;
            expected[i] = (p_right-p_left) / x_len * x + p_left;
        }

        ASSERT_DOUBLE_ARRAY_EQ(&expected[0], &(*vec_h)[0], vec_h->size());

    } catch (const char* e) {
        std::cout << "***Exception caught! " << e << std::endl;
    }

}

TEST(Solution, CouplingFem2)
{
	try {
        //space
	    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 20, .0, .0, .0);
	    Rectangle* _rec = new Rectangle(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
        //time
        //TimeStepFunctionConstant tim(.0, 1e+3, 1e+3);
        TimeStepFunctionConstant tim(.0, 1e+4, 1e+3);
        //material
	    Geo::PorousMedia pm;
	    pm.hydraulic_conductivity = new NumLib::TXFunctionConstant(1.e-11);
	    pm.porosity = new NumLib::TXFunctionConstant(1.0);
	    Geo::Compound tracer;
	    tracer.molecular_diffusion = new NumLib::TXFunctionConstant(1.e-6);
        //problems
	    DiscreteSystem dis(*msh);
	    Geo::GWFemProblem* pGW = defineGWProblem(dis, *_rec, pm);
	    Geo::MassFemProblem* pMass = defineMassTransportProblem(dis, *_rec, pm, tracer);
        pGW->setTimeSteppingFunction(tim);
        pMass->setTimeSteppingFunction(tim);
	    //options
	    Base::Options optionsGW;
	    Base::Options* op_lis = optionsGW.addSubGroup("Lis");
	    op_lis->addOption("solver_type", "CG");
	    op_lis->addOption("precon_type", "NONE");
	    op_lis->addOptionAsNum("error_tolerance", 1e-10);
	    op_lis->addOptionAsNum("max_iteration_step", 1000);
	    Base::Options* op_nl = optionsGW.addSubGroup("Nonlinear");
	    op_nl->addOption("solver_type", "Picard");
	    op_nl->addOptionAsNum("error_tolerance", 1e-6);
	    op_nl->addOptionAsNum("max_iteration_step", 500);
        Base::Options optionsMT;
        op_lis = optionsMT.addSubGroup("Lis");
        op_lis->addOption("solver_type", "BiCG");
        op_lis->addOption("precon_type", "NONE");
        op_lis->addOptionAsNum("error_tolerance", 1e-10);
        op_lis->addOptionAsNum("max_iteration_step", 1000);
	    Base::Options* op_nl2 = optionsMT.addSubGroup("Nonlinear");
	    op_nl2->addOption("solver_type", "Picard");
	    op_nl2->addOptionAsNum("error_tolerance", 1e-6);
	    op_nl2->addOptionAsNum("max_iteration_step", 500);
        // unknowns
		MyFunctionHead f_head;
		f_head.define(&dis, pGW, optionsGW);
        f_head.setOutputParameterName(0, "h");
		MyFunctionVelocity f_vel;
		f_vel.define(dis, pm);
        f_vel.setInputParameterName(0, "h");
        f_vel.setOutputParameterName(0, "v");
		MyFunctionConcentration f_c;
		f_c.define(&dis, pMass, optionsMT);
        f_c.setInputParameterName(0, "v");
        f_c.setOutputParameterName(0, "c");


        FemFunctionConvergenceCheck checker(&dis);
        SerialStaggeredMethod method(checker, 1e-5, 100);
	    AsyncPartitionedSystem apart1;
        apart1.setAlgorithm(method);
        apart1.resizeOutputParameter(3);
        apart1.setOutputParameterName(0, "h");
        apart1.setOutputParameterName(1, "v");
        apart1.setOutputParameterName(2, "c");
        apart1.addProblem(f_head);
        apart1.addProblem(f_vel);
        apart1.addProblem(f_c);
        apart1.connectParameters();

	    TimeSteppingController timestepping;
	    timestepping.addTransientSystem(apart1);

	    //const double epsilon = 1.e-3;
	    timestepping.setBeginning(.0);
	    timestepping.solve(tim.getEnd());

	    const FemNodalFunctionScalar* r_f_head = apart1.getOutput<FemNodalFunctionScalar>(apart1.getOutputParameterID("h"));
	    const FEMIntegrationPointFunctionVector* r_f_v = apart1.getOutput<FEMIntegrationPointFunctionVector>(apart1.getOutputParameterID("v"));
	    const FemNodalFunctionScalar* r_f_c = apart1.getOutput<FemNodalFunctionScalar>(apart1.getOutputParameterID("c"));
	    const IDiscreteVector<double>* vec_h = r_f_head->getNodalValues();
	    //const FEMIntegrationPointFunctionVector2d::DiscreteVectorType* vec_v = r_f_v->getNodalValues();
	    const IDiscreteVector<double>* vec_c = r_f_c->getNodalValues();

	    //r_f_head->printout();
	    //r_f_v->printout();
//#undef OUTPUT_C
#define OUTPUT_C
#ifdef OUTPUT_C
        r_f_c->printout();
#endif
	    std::vector<double> expectedHead(21);
    	const double p_left = 2.e+6;
    	const double p_right = .0;
    	const double x_len = 2.0;
#ifdef OUTPUT_C
        std::cout << std::endl << "expected p=";
#endif
        for (size_t i=0; i<expectedHead.size(); i++) {
        	double x = i*0.1;
        	expectedHead[i] = (p_right-p_left) / x_len * x + p_left;
#ifdef OUTPUT_C
            std::cout << expectedHead[i] << " ";
#endif
        }
//	    getGWExpectedHead(expectedHead);
        ASSERT_DOUBLE_ARRAY_EQ(&expectedHead[0], &(*vec_h)[0], expectedHead.size(), 1e-5);


        std::vector<double> expectedC;
        expectedC.resize(21);

#ifdef OUTPUT_C
        std::cout << std::endl << "expected C=";
#endif
        for (size_t i=0; i<expectedC.size(); i++) {
            expectedC[i] = analyticalOgataBank(i*0.1, tim.getEnd(), 1e-5/1.0, 1e-6);
#ifdef OUTPUT_C
            std::cout << expectedC[i] << " ";
#endif
        }
#ifdef OUTPUT_C
        std::cout << std::endl;
#endif
	    ASSERT_DOUBLE_ARRAY_EQ(&expectedC[0], &(*vec_c)[0], expectedC.size(), 5e-2);


	} catch (const char* e) {
		std::cout << "***Exception caught! " << e << std::endl;
	}

}

#if 0
TEST(Fem, LinearElastic2D)
{
	try {
	    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
	    Rectangle* _rec = new Rectangle(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
	    Geo::PorousMedia pm;
	    pm.hydraulic_conductivity = new NumLib::TXFunctionConstant(1.e-11);
	    pm.porosity = new NumLib::TXFunctionConstant(0.2);
	    DiscreteSystem dis(*msh);
	    Geo::FemLinearElasticProblem* pGW = 0; //defineGWProblem(dis, *_rec, pm);
        TimeStepFunctionConstant tim(.0, 100.0, 10.0);
        pGW->setTimeSteppingFunction(tim);
	    // options
	    Base::Options options;
	    Base::Options* op_lis = options.addSubGroup("Lis");
	    op_lis->addOption("solver_type", "CG");
	    op_lis->addOption("precon_type", "NONE");
	    op_lis->addOptionAsNum("error_tolerance", 1e-10);
	    op_lis->addOptionAsNum("max_iteration_step", 500);
	    Base::Options* op_nl = options.addSubGroup("Nonlinear");
	    op_nl->addOption("solver_type", "Picard");
	    op_nl->addOptionAsNum("error_tolerance", 1e-6);
	    op_nl->addOptionAsNum("max_iteration_step", 500);

		MyFunctionDisplacement f_u;
		f_u.define(&dis, pGW, options);
		f_u.setOutputParameterName(0, "u");
//		MyFunctionVelocity f_vel;
//		f_vel.define(dis, pm);
//        f_vel.setInputParameterName(0, "h");
//        f_vel.setOutputParameterName(0, "v");
//
//        FemFunctionConvergenceCheck checker;
//	    SerialStaggeredMethod method(checker, 1e-5, 100);
//	    AsyncPartitionedSystem apart1;
//        apart1.setAlgorithm(method);
//        apart1.resizeOutputParameter(2);
//        apart1.setOutputParameterName(0, "h");
//        apart1.setOutputParameterName(1, "v");
//        apart1.addProblem(f_head);
//        apart1.addProblem(f_vel);
//        apart1.connectParameters();
//
//	    TimeSteppingController timestepping;
//	    timestepping.addTransientSystem(apart1);
//
//	    //const double epsilon = 1.e-3;
//	    timestepping.setBeginning(.0);
//	    timestepping.solve(1.0);
//
//	    const FemNodalFunctionScalar* r_f_head = apart1.getOutput<FemNodalFunctionScalar>(apart1.getOutputParameterID("h"));
//	    const FEMIntegrationPointFunctionVector2d* r_f_v = apart1.getOutput<FEMIntegrationPointFunctionVector2d>(apart1.getOutputParameterID("v"));
//	    const DiscreteVector<double>* vec_h = r_f_head->getNodalValues();
//	    //const FEMIntegrationPointFunctionVector2d::DiscreteVectorType* vec_v = r_f_v->getNodalValues();
//
//	    r_f_head->printout();
//	    r_f_v->printout();
//
//	    std::vector<double> expected;
//	    getGWExpectedHead(expected);
//
//	    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], &(*vec_h)[0], vec_h->size());

	} catch (const char* e) {
		std::cout << "***Exception caught! " << e << std::endl;
	}

}
#endif
