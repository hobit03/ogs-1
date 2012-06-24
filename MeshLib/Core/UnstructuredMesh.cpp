
#include "UnstructuredMesh.h"

#include "BaseLib/CodingTools.h"


namespace MeshLib
{

UnstructuredMesh::~UnstructuredMesh()
{
    BaseLib::releaseObjectsInStdVector(_list_nodes);
    BaseLib::releaseObjectsInStdVector(_list_elements);
    BaseLib::releaseObjectsInStdVector(_list_edge_elements);
}

size_t UnstructuredMesh::addNode( GeoLib::Point &x, size_t order ) 
{
    size_t new_node_id = _list_nodes.size();
    this->_list_nodes.push_back(new Node(new_node_id, x));
    this->_map_order_nnodes[order]++;
    return  new_node_id;
}

size_t UnstructuredMesh::setNodeCoordinates( size_t node_id, GeoLib::Point &x ) 
{
    size_t new_node_id = node_id;
    if (node_id<_list_nodes.size()) {
        INode *node = this->_list_nodes.at(node_id);
        node->setX(x);
    } else {
        new_node_id = addNode(x);
    }
    return new_node_id;
};

void UnstructuredMesh::getListOfNodeCoordinates(const std::vector<size_t> &vec_node_id, std::vector<GeoLib::Point> &vec_pt) const 
{
    for (size_t i=0; i<vec_node_id.size(); i++)
        vec_pt.push_back(*this->getNodeCoordinatesRef(vec_node_id[i]));
}

size_t UnstructuredMesh::getNumberOfElements(ElementShape::type ele_type) const
{
    std::map<ElementShape::type, size_t>::const_iterator itr = _n_eles_type.find(ele_type);
    if (itr==_n_eles_type.end()) return 0;
    else return itr->second;
}

size_t UnstructuredMesh::addElement( IElement *e) 
{
    assert(e!=0);
    e->setID(_list_elements.size());
    _list_elements.push_back(e);
    _n_eles_type[e->getShapeType()]++;

    return e->getID();
};

void UnstructuredMesh::addEdgeElement(IElement *e) 
{
    assert(e!=0);
    e->setID(_list_edge_elements.size());
    _list_edge_elements.push_back(e);
}

size_t UnstructuredMesh::getNumberOfEdges() const {
    return _list_edge_elements.size();
}

IElement* UnstructuredMesh::getEdgeElement(size_t edge_id) {
    if (edge_id<_list_edge_elements.size())
        return _list_edge_elements[edge_id];
    return 0;
}

} // end namespace
