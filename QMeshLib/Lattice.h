// Lattice.h: interface for the QMeshPatch class.
//
//////////////////////////////////////////////////////////////////////

#ifndef Lattice_h
#define Lattice_h

#include <vector>
#include <algorithm>
#include <memory>

#include "Geometry.h"
#include "../GLKLib/GLKObList.h"




class Lattice : public Geometry, public GLKObject {
public:
    Lattice() {}
    ~Lattice() {}


};

#endif
