/*************************************************************************
*
* ImplicitSlicer is developed by Qiang Zou for research use. All rights about
* the program (esp. surface reconstruction) are reserved by Qiang Zou. This C++
* source codes are available only to a primary user for academic purposes. No
* secondary use, such as copy, distribution, diversion, business purpose, etc.,
* is allowed. In no event shall the author be liable to any party for direct,
* indirect, special, incidental, or consequential damage arising out of the use
* of this program. ImplicitSlicer is self-contained.
* Qiang Zou, Aug. 2020, qzou.code@gmail.com
*/


#ifndef LatticeModeler_H
#define LatticeModeler_H

#include <list>
#include <string>
#include <memory>
#include <deque>
#include <set>


#include "QMeshLib/Geometry.h"
#include "QMeshLib/QSurfaceMesh.h"

#include "Eigen/Core"


struct LatticeEdge;
struct LatticeNode{
    double p[3]; // coordinates
    std::vector<LatticeEdge *> adjEdges; // adjacent edges
    size_t idx/* = std::numeric_limits<size_t>::max()*/;
    void addAdjEdge(std::shared_ptr<LatticeEdge> e) {adjEdges.push_back(e.get());}
    auto& getAdjEdges() {return adjEdges;}
};

struct LatticeEdge{
    LatticeNode* start;
    LatticeNode* end;
    double weight; // for storing radius
    size_t idx;
    LatticeNode* getOtherPoint(std::shared_ptr<LatticeNode> n) {
        return n.get() == start ? end : start;
    }
};

struct LatticeGraph{
    std::vector<std::shared_ptr<LatticeNode>> nodeArray;
    std::vector<std::shared_ptr<LatticeEdge>> edgeArray;
    void addNode(std::shared_ptr<LatticeNode> n) {n->idx = nodeArray.size(); nodeArray.push_back(n);}
    void addEdge(std::shared_ptr<LatticeEdge> e) {e->idx = edgeArray.size(); edgeArray.push_back(e);}
    bool addEdge(size_t idxStart, size_t idxEnd, double weight = 1);
    auto& getNodes() {return nodeArray;}
    auto& getNode(size_t i) {return nodeArray[i];}
    auto& getEdges() {return edgeArray;}
    auto& getEdge(size_t i) {return edgeArray[i];}
    void clear() {edgeArray.clear(); nodeArray.clear();}
    bool isEmpty() {return nodeArray.empty() || edgeArray.empty();}
};

class LatticeModeler
{
public:
    LatticeModeler();
    ~LatticeModeler();

    // lattice geometry and file
    void clearLatticeModel();
    bool writeLatticeModel(bool isHybrid, string path);
    bool readLatticeModel(std::string path); // hybrid format
    void getBBox(Eigen::Vector3d& left, Eigen::Vector3d& right);
    template <typename T> void getBBox(Eigen::Vector3d& left, Eigen::Vector3d& right, T &edges, bool hasIdx = false);

    // 3d printing
    bool sliceModelImplicitStream(std::string path);

    void save2Pic(std::vector<std::vector<double>> &digitImage, std::string path, int dpi);

	// g-code utilities
	bool writeCLIFile(std::string path, std::vector<std::shared_ptr<QMeshPatch>>& layers);
	bool writeCLIFileBin(string path, std::vector<std::shared_ptr<QMeshPatch> > &layers);

private:
    std::string surfaceMeshFilePath;
    std::string volumeMeshFilePath;
    std::unique_ptr<LatticeGraph> ptrLatticeModel;
    double tolerence;
    double avgLength;
    double m_radiusConvol;
    double m_thresholdConvol;
	std::vector<std::shared_ptr<QMeshPatch>> layers;
};


#endif // SLICER_H
