
#include "Tools.h"

#include <algorithm>
#include <exception>

#include "MeshLib/Core/ElementFactory.h"
#include "MeshLib/Tools/MeshNodesAlongPolyline.h"
#include "MeshLib/Topology/Topology.h"

namespace MeshLib
{

/// seeks all nodes located on a given polyline and returns a list of the found node pointers.
void findNodesOnPolyline(IMesh const* msh, GeoLib::Polyline const* poly, std::vector<size_t> *vec_nodes)
{
    MeshNodesAlongPolyline obj(poly, msh);
    std::vector<size_t> vec_node_id = obj.getNodeIDs();
    vec_nodes->assign(vec_node_id.begin(), vec_node_id.end());
};

///
void findNodesOnGeometry(IMesh const* msh, GeoLib::GeoObject const* obj, std::vector<size_t> *vec_nodes)
{
    switch (obj->getGeoType()) {
        case GeoLib::GeoObjType::POLYLINE:
            findNodesOnPolyline(msh, static_cast<GeoLib::Polyline const*>(obj), vec_nodes);
            break;
        default:
            throw std::exception("This geo type is not supported in MeshLib::findNodesOnGeometry()");
            break;
    }
};

/// get a list of elements connected to one of give nodes
void findConnectedElements(IMesh const* msh, const std::vector<size_t> &nodes, std::vector<size_t> &connected_elements)
{
    for (size_t i=0; i<msh->getNumberOfElements(); i++) {
        IElement* e = msh->getElemenet(i);
        for (size_t j=0; j<e->getNumberOfNodes(); j++) {
            if (std::find(nodes.begin(), nodes.end(), e->getNodeID(j))!=nodes.end()) {
                connected_elements.push_back(e->getID());
                break;
            }
        }
    }
};

/// find edge elements on a given polyline
void findEdgeElementsOnPolyline(IMesh * msh, GeoLib::Polyline const* poly, std::vector<IElement*> *vec_edges_on_poly)
{
    // get a list of nodes on the polyline
    std::vector<size_t> nodes_on_poly;
    findNodesOnPolyline(msh, poly, &nodes_on_poly);
    // get a list of elements having the nodes
    std::vector<size_t> elements_near_poly;
    findConnectedElements(msh, nodes_on_poly, elements_near_poly);
    // get a list of edges made of the nodes
    std::vector<IElement*> selected_edges;
    createEdgeElements(msh, elements_near_poly, selected_edges);
    const size_t nr_edges = selected_edges.size();
    for (size_t i=0; i<nr_edges; i++) {
        IElement *edge_e = selected_edges[i];
        // check
        size_t cnt_match = 0;
        for (size_t j=0; j<edge_e->getNumberOfNodes(); j++) {
            //const INode* edge_nod = msh->getNode(edge_e->getNodeID(j));
            if (std::find(nodes_on_poly.begin(), nodes_on_poly.end(), edge_e->getNodeID(j)) != nodes_on_poly.end())
                cnt_match++;
            else
                break;
        }
        // update the list
        if (cnt_match==edge_e->getNumberOfNodes())
            vec_edges_on_poly->push_back(edge_e);
    }
};

/// find boundary elements located on given geometries
void findBoundaryElementsOnGeometry(IMesh * msh, GeoLib::GeoObject const* obj, std::vector<IElement*> *vec_eles)
{
    switch (obj->getGeoType()) {
        case GeoLib::GeoObjType::POLYLINE:
            findEdgeElementsOnPolyline(msh, static_cast<GeoLib::Polyline const*>(obj), vec_eles);
            break;
        default:
            throw std::exception("This geo type is not supported in MeshLib::findNodesOnGeometry()");
            break;
    }
};

IElement* createEdgeElement(IMesh& msh, IElement &e, size_t j)
{
    std::vector<size_t> vec_edge_nodes;
    IElement* edge = ElemenetFactory::createNewElement(e.getEdgeElementType(j)); 
    e.getNodeIDsOfEdgeElement(j, vec_edge_nodes);
    for (size_t k=0; k<vec_edge_nodes.size(); k++) {
        edge->setNodeID(k, vec_edge_nodes[k]);
    }
    msh.addEdgeElement(edge);
    return edge;
}

IElement* findEdgeElement(const std::vector<IElement*> &edges, std::vector<size_t> &vec_edge_nodes)
{
    std::sort(vec_edge_nodes.begin(), vec_edge_nodes.end());

    IElement *e_edge = 0;
    for (std::vector<IElement*>::const_iterator itr=edges.begin(); itr!=edges.end(); itr++) {
        if ((*itr)->hasNodeIds(vec_edge_nodes)) {
            e_edge = *itr;
            break;
        }
    }
    return e_edge;
}

void findEdgeElements(IMesh& msh, IElement &e, std::vector<IElement*> &edges)
{
    std::vector<size_t> vec_edge_nodes;
    const size_t n_edges = e.getNumberOfEdges();
    for (size_t j=0; j<n_edges; j++) {
        //check if already exists
        e.getNodeIDsOfEdgeElement(j, vec_edge_nodes);
        IElement* edge = findEdgeElement(edges, vec_edge_nodes);
        if (edge!=0) {
            if (e.getEdgeElement(j)==0)
                e.setEdgeElement(j, edge);
            continue;
        }

        edge = e.getEdgeElement(j);
        if (edge==0) {
            //if new, create 
            if (edge==0) {
                edge = createEdgeElement(msh, e, j);
            }
            e.setEdgeElement(j, edge);
        }
        edges.push_back(edge);

    }
}

/// get a list of edge elements of given elements
void createEdgeElements(IMesh * msh, const std::vector<size_t> &selected_ele, std::vector<IElement*> &edges)
{
    for (size_t i=0; i<selected_ele.size(); i++) {
        IElement *e = msh->getElemenet(selected_ele[i]);
        findEdgeElements(*msh, *e, edges);
    }
};

void createEdgeElements(IMesh * msh)
{
    std::vector<IElement*> edges;
    for (size_t i=0; i<msh->getNumberOfElements(); i++) {
        IElement *e = msh->getElemenet(i);
        findEdgeElements(*msh, *e, edges);
    }
};


/**
 * generate higher order mesh
 */
void generateHigherOrderUnstrucuredMesh(UnstructuredMesh &msh, size_t order)
{
    assert(order<3);

    TopologySequentialNodes2Elements nod2ele(msh);

    // make all edges higher order
    for (size_t i=0; i<msh.getNumberOfEdges(); i++) {
        IElement* e = msh.getEdgeElement(i);
        double const* const pnt0(msh.getNodeCoordinatesRef(e->getNodeID(0))->getData());
        double const* const pnt1(msh.getNodeCoordinatesRef(e->getNodeID(1))->getData());
        e->setMaximumOrder(order);
        if (order==2) {
            // add midpoint
            size_t nod_id = msh.addNode(GeoLib::Point(.5 * (pnt0[0] + pnt1[0]), 0.5 * (pnt0[1] + pnt1[1]), 0.5 * (pnt0[2] + pnt1[2])), 2);
            e->setNodeID(2, nod_id);
        } else {
            //
        }
    }

    // set the new node ids to all elements
    for (size_t i=0; i<msh.getNumberOfElements(); i++) {
        IElement* e = msh.getElemenet(i);
        size_t e_nnodes1 = e->getNumberOfNodes(1);
        e->setMaximumOrder(order);
        // for each edge
        for (size_t j=0; j<e->getNumberOfEdges(); j++) {
            IElement *edge = e->getEdgeElement(j);
            size_t new_nod = edge->getNodeID(2);
            size_t local_id = e_nnodes1 + j;
            e->setNodeID(local_id, new_nod);
        }

        // Quad 9
        if (e->getShapeType() == ElementShape::QUAD || e->getNumberOfNodes(2)==9)
        {
            double x0, y0, z0;
            x0 = y0 = z0 = 0.0;
            std::vector<size_t> list_e_nodes;
            e->getNodeIDList(1, list_e_nodes);
            for (size_t i = 0; i < list_e_nodes.size(); i++) // Nodes
            {
                const GeoLib::Point* pt = msh.getNodeCoordinatesRef(list_e_nodes[i]);
                double const* const pnt_i(pt->getData());
                x0 += pnt_i[0];
                y0 += pnt_i[1];
                z0 += pnt_i[2];
            }
            x0 /= (double) e_nnodes1;
            y0 /= (double) e_nnodes1;
            z0 /= (double) e_nnodes1;
            size_t nodid = msh.addNode(GeoLib::Point(x0,y0,z0),2);
            e->setNodeID(8, nodid);
        }
    }
}

} // end namespace