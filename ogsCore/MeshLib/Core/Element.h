/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Element.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma  once

#include "ElementTopology.h"
#include "TemplateElement.h"

namespace MeshLib
{

// elements
typedef TemplateElement<LineTopology> Line;
typedef TemplateElement<TriangleTopology> Triangle;
typedef TemplateElement<Quad9Topology> Quadrirateral;
typedef TemplateElement<TetraTopology> Tetrahedron;
typedef TemplateElement<PyramidTopology> Pyramid;
typedef TemplateElement<PrismTopology> Prism;
typedef TemplateElement<HexTopology> Hexahedron;

} // end namespace

