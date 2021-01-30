// QSurfaceMesh.h: interface for the QMeshPatch class.
//
//////////////////////////////////////////////////////////////////////

#ifndef QSurfaceMesh_h
#define QSurfaceMesh_h

#include <vector>
#include <algorithm>
#include <memory>

#include "Geometry.h"
#include "QMesh/QMeshPatch.h"

class QSurfaceMesh : public QMeshPatch, public Geometry {
public:
    QSurfaceMesh() {QSurfaceMesh::QMeshPatch();}
    ~QSurfaceMesh() {}

    GeometryTypes getGeometryType() {return SurfaceMeshQMeshPatch;}
    std::shared_ptr<Geometry> getEntities2Slice(bool isIncremental, double height);
};

#endif
