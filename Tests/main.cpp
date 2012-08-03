/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file main.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include <gtest/gtest.h>

//int add (int x, int y) {return x+y;};
//
//TEST(Dummy, Test1)
//{
//    ASSERT_EQ(2, add(1, 1));
//}


int main(int argc, char *argv[])
{
#if 1
    argc = 2;
    //argv[1] = "--gtest_filter=Math.Nonlinear*";
    //argv[1] = "--gtest_filter=Num.Discrete*:FEM.*";
    //argv[1] = "--gtest_filter=Discrete.NDDC*";
    //argv[1] = "--gtest_filter=Solution.Coupling*";
    //argv[1] = "--gtest_filter=Math.Matrix*";
    //argv[1] = "--gtest_filter=Coupling.*";
    //argv[1] = "--gtest_filter=Math.SystemOfEqs*";
    argv[1] = "--gtest_filter=Fem.LinearElastic2D";
    //argv[1] = "--gtest_filter=Solution.CouplingF*";
    //argv[1] = "--gtest_filter=*";
#endif
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

