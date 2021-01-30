#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <memory>

enum GeometryTypes {
    Convolution =   0,
    Graph =         1,
    Contours =      2,
    MarchingSquare =   3,
    SurfaceMeshQMeshPatch =   4,
    SurfaceMeshOpenMesh =   5,
};

class Geometry
{
public:
    Geometry();
    ~Geometry() {}

    virtual GeometryTypes getGeometryType() = 0;
    virtual std::shared_ptr<Geometry> getEntities2Slice(bool isIncremental, double height) = 0;
};

#endif // GEOMETRY_H
