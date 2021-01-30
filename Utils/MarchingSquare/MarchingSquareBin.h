#pragma once

#include "../GLKLib/GLKObList.h"
#include <unordered_map>
#include <string>

class QMeshNode;
class QMeshFace;
class QMeshEdge;
class QMeshPatch;

class MarchingSquareBin
{
public:
    MarchingSquareBin(std::unordered_map<std::string, double> &m_image_, size_t xres_, size_t zres_,
                      double xzgridwidth_, double height_, double* imageorigin_);
    ~MarchingSquareBin();

    void doContouringWithOptim();
    void doContouring();
    void writeContours(const std::string &path);
    QMeshPatch *getContours(bool isDeleteData);

private:
    void FindAllSticks();
    void basicContouring();
    void BiLaplacianSmoothOp(double smootherrorcoeff, double movementcoeffict);
    double PerformVSAOnContour(int iter, double paradistterror, int simpratio);
    void convert2XYPlane();

private:
    double imageorigin[2];
    size_t xres;
    size_t zres;
    double xzgridwidth;
    double height;
    std::unordered_map<std::string, double>* m_image; // string stores "rowidx" + wp + "colidx"; if a grid point can be found, it is an inner point
    GLKArray *stickindexx;			//each grid edge has an index which is encoded as
    GLKArray *stickindexz;
    QMeshPatch *m_contour;
};

