
#pragma  once

#include <vector>

#include "GeoLib/Core/Point.h"

#include "MeshLib/Core/INode.h"
//#include "CoordinateSystem.h"

namespace MeshLib
{

/// List of available element types
struct ElementShape
{
    enum type {
        LINE = 1,
        QUAD = 2,
        HEXAHEDRON = 3,
        TRIANGLE = 4,
        TETRAHEDRON = 5,
        PRISM = 6,
        PYRAMID = 7,
        INVALID = -1
    };
};

struct ElementType
{
    enum type {
        LINE2 = 1,
        QUAD4 = 2,
        HEXAHEDRON8 = 3,
        TRIANGLE3 = 4,
        TETRAHEDRON4 = 5,
        PRISM6 = 6,
        PYRAMID5 = 7,
        LINE3 = 8,
        QUAD8 = 9,
        QUAD9 = 10,
        HEXAHEDRON20 = 11,
        TRIANGLE6 = 12,
        TETRAHEDRON10 = 13,
        INVALID = -1
    };
};

class IElementCoordinatesMapping;

/**
 * \brief Interface of element classes
 *
 */
class IElement
{
public:
    virtual ~IElement() {};
    ///
    virtual void initialize() = 0;
    virtual IElement* clone() const = 0;

    /// return if this element has the given list of nodes 
    virtual bool hasNodeIds(const std::vector<size_t> &node_ids) const = 0;

    /// return this element id
    virtual size_t getID() const = 0;
    /// set this element id
    virtual void setID(size_t id) = 0;

    /// return the group id of this element
    virtual size_t getGroupID() const = 0;
    /// set the group if of this element
    virtual void setGroupID(size_t id) = 0;

    /// return the shape type
    virtual ElementShape::type getShapeType() const = 0;
    /// return intrinsic dimensions of this element
    virtual size_t getDimension() const = 0;

    /// return the number of nodes under current order
    virtual size_t getNumberOfNodes() const = 0;
    /// set node id
    virtual void setNodeID(const size_t &local_node_id, const size_t &node_id) = 0;
    /// return node id
    virtual size_t getNodeID(size_t local_node_id) const = 0;
    /// get a list of node ids
    virtual void getNodeIDList( std::vector<size_t> &e_node_id_list ) const = 0;

    /// get the number of element faces
    virtual size_t getNumberOfFaces() const = 0;

    /// get the number of element edges
    virtual size_t getNumberOfEdges() const = 0;
    /// get an edge element
    virtual IElement* getEdgeElement(size_t edge_id) const = 0;
    /// set an edge element
    virtual void setEdgeElement(size_t edge_id, IElement* e) = 0;
    /// get node ids of the edge element
    virtual void getNodeIDsOfEdgeElement(size_t edge_id, std::vector<size_t> &vec_node_ids) const = 0;
    /// get edge element type
    virtual ElementShape::type getEdgeElementType(size_t edge_id) const = 0;

    /// set mapped coordinates
    virtual void setMappedCoordinates(IElementCoordinatesMapping* mapping) = 0;
    /// get mapped coordinates
    virtual IElementCoordinatesMapping* getMappedCoordinates() = 0;

    /// set max. order
    virtual void setMaximumOrder(size_t) = 0;
    /// get max. order
    virtual size_t getMaximumOrder() const = 0;
    /// set current order
    virtual void setCurrentOrder(size_t) = 0;
    /// get current order
    virtual size_t getCurrentOrder() const = 0;
    /// return the number of nodes under the given order
    virtual size_t getNumberOfNodes(size_t order) const = 0;
    /// get a list of node id
    virtual void getNodeIDList( size_t order, std::vector<size_t> &e_node_id_list ) const = 0;
    /// get a list of node id of the specified edge
    virtual void getNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_node_ids) const = 0;
    /// get a list of node size with orders
    virtual void getListOfNumberOfNodesForAllOrders(std::vector<size_t> &vec) const = 0;
};

}
