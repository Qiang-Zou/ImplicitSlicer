#include "LatticeModeler.h"
#include "Utils/MarchingSquare/MarchingSquareBin.h"

#include "QMeshLib/QMesh/QMeshNode.h"
#include "QMeshLib/QMesh/QMeshEdge.h"
#include "QMeshLib/QMesh/QMeshFace.h"
#include "QMeshLib/QMesh/QMeshPatch.h"

//#include "Utils/tetgen/tetgen.h"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <experimental/unordered_map>

#include <QImage>

bool LatticeGraph::addEdge(size_t idxStart, size_t idxEnd, double weight)
{
    if (idxStart >= nodeArray.size() || idxEnd >= nodeArray.size()) return false;
    // add edge
    std::shared_ptr<LatticeEdge> e(new LatticeEdge);
    e->start = nodeArray[idxStart].get();
    e->end = nodeArray[idxEnd].get();
    e->weight = weight;
    e->idx = edgeArray.size();
    edgeArray.push_back(e);
    return true;
}

LatticeModeler::LatticeModeler() :
    surfaceMeshFilePath("../../Results/in_surface.off"),
    volumeMeshFilePath("../../Results/out_volume"),
    tolerence(1e-6),
    ptrLatticeModel(new LatticeGraph)
{

}

LatticeModeler::~LatticeModeler()
{
    // delete files
    // if(remove("../results/in_surface.off") != 0)
    //    std::cerr << "LatticeModeler > Error deleting files" << endl;
}

void LatticeModeler::clearLatticeModel()
{
    ptrLatticeModel.reset();
    std::cout << "LatticeModeler > successfully deleting the graph model" << endl;
}

bool LatticeModeler::writeLatticeModel(bool isHybrid, string path)
{
    if (ptrLatticeModel->isEmpty()) {
        std::cerr << "LatticeMode∫ler > error writing latticeModel." << endl;
        return false;
    }

    if (!isHybrid) { // save points and edges in two sections
        std::ofstream outfile(path);
        if (!outfile.is_open()) {
            std::cerr << "LatticeModeler > error writing latticeModel." << endl;
            return false;
        }

        outfile<< ptrLatticeModel->nodeArray.size() << " " << ptrLatticeModel->edgeArray.size();
        // output points
        for (size_t i(0); i < ptrLatticeModel->nodeArray.size(); ++i) {
            auto& xyz = ptrLatticeModel->nodeArray[i]->p;
            outfile << endl << xyz[0] << " " << xyz[1] << " " << xyz[2];
        }
        // output edges
        for (size_t i(0); i < ptrLatticeModel->edgeArray.size(); ++i) {
            auto e = ptrLatticeModel->edgeArray[i];
            outfile << endl << e->start->idx << " " << e->end->idx << " " << e->weight;
        }
        outfile.close();
        std::cout << "LatticeModeler > successfully write latticeModel." << endl;
    }
    else { // save points and edges together
        std::ofstream outfile(path);
        if (!outfile.is_open()) {
            std::cerr << "LatticeModeler > error writing latticeModel." << endl;
            return false;
        }

        outfile<< ptrLatticeModel->nodeArray.size() << " " << ptrLatticeModel->edgeArray.size();
        // 2-n lines: p1 p2 adjacent-edges
        for (size_t i(0); i < ptrLatticeModel->getEdges().size(); ++i) {
            auto e = ptrLatticeModel->edgeArray[i];
            auto& p1 = e->start->p;
            auto& p2 = e->end->p;
            outfile << endl << p1[0] << " " << p1[1] << " " << p1[2] << " "
                    << p2[0] << " " << p2[1] << " " << p2[2] << " " << e->weight;
        }
        outfile.close();
        std::cout << "LatticeModeler > successfully write latticeModel." << endl;
    }

    return true;
}

bool LatticeModeler::readLatticeModel(string path)
{
    ptrLatticeModel->clear();

    // open file
    std::ifstream ifs(path.c_str(), ifstream::in);
    std::string line;
    if (!(ifs.good() && !ifs.eof())) {
        std::cerr << "LatticeModeler > error openning file" << endl;
        return false;
    }

    int numOfEdges;
    double radiusEdge;

    std::getline(ifs, line); // first line
    stringstream str(line);
    str >> numOfEdges >> radiusEdge >> m_radiusConvol >> m_thresholdConvol;
    //cout << numOfEdges << endl;

    ptrLatticeModel->nodeArray.reserve(numOfEdges * 2);
    ptrLatticeModel->edgeArray.reserve(numOfEdges);
    for (size_t i(0); i < numOfEdges; ++i) {
        // add nodes
        std::shared_ptr<LatticeNode> n1(new LatticeNode);
        std::shared_ptr<LatticeNode> n2(new LatticeNode);
        double radius;

        std::getline(ifs, line);
        str = stringstream(line);
        str >> n1->p[0] >> n1->p[1] >> n1->p[2];
        str >> n2->p[0] >> n2->p[1] >> n2->p[2];
        str >> radius;
        double scale = 1.4;
        n1->p[0] *= scale; n1->p[1] *= scale; n1->p[2] *= scale;// scale
        n2->p[0] *= scale; n2->p[1] *= scale; n2->p[2] *= scale;
        ptrLatticeModel->addNode(n1);
        ptrLatticeModel->addNode(n2);

        // add edges
        ptrLatticeModel->addEdge(n1->idx, n2->idx, radius);
    }
    ifs.close();

    return true;
}

void LatticeModeler::getBBox(Eigen::Vector3d &left, Eigen::Vector3d &right)
{
    if (ptrLatticeModel->isEmpty()) {
        left = right = Eigen::Vector3d(0, 0, 0);
        return;
    }

    left[0] = left[1] = left[2] = numeric_limits<double>::max();
    right[0] = right[1] = right[2] = numeric_limits<double>::min();
    for (auto iter(ptrLatticeModel->nodeArray.begin()); iter != ptrLatticeModel->nodeArray.end(); ++iter) {
        auto& p = (*iter)->p;
        left[0] = std::min(left[0], p[0]);
        left[1] = std::min(left[1], p[1]);
        left[2] = std::min(left[2], p[2]);
        right[0] = std::max(right[0], p[0]);
        right[1] = std::max(right[1], p[1]);
        right[2] = std::max(right[2], p[2]);
    }
}

template <typename T>
void LatticeModeler::getBBox(Eigen::Vector3d &left, Eigen::Vector3d &right, T &edges, bool hasIdx)
{
    if (edges.empty()) {
        left = right = Eigen::Vector3d(0, 0, 0);
        return;
    }

    left[0] = left[1] = left[2] = numeric_limits<double>::max();
    right[0] = right[1] = right[2] = numeric_limits<double>::min();
    size_t prefix = hasIdx ? 1 : 0;
    for (auto& e : edges) {
        for (size_t i(0); i < 2; ++i) {
            Eigen::Vector3d p(e[i * 3 + 0 + prefix], e[i * 3 + 1 + prefix], e[i * 3 + 2 + prefix]);
            left[0] = std::min(left[0], p[0] - e[6 + prefix]); // p[6] stores raidus
            left[1] = std::min(left[1], p[1] - e[6 + prefix]);
            left[2] = std::min(left[2], p[2] - e[6 + prefix]);
            right[0] = std::max(right[0], p[0] + e[6 + prefix]);
            right[1] = std::max(right[1], p[1] + e[6 + prefix]);
            right[2] = std::max(right[2], p[2] + e[6 + prefix]);
        }
    }
}

bool LatticeModeler::sliceModelImplicitStream(std::string path)
{
    typedef Eigen::Matrix<double, 7, 1> Vector7d; // store edge

    //  read model
    if (readLatticeModel(path)) std::cout << "LatticeModeler > successfully read the model." << endl;
    if (ptrLatticeModel->isEmpty()) { std::cerr << "LatticeModeler > can't find a model" << endl; return false;}

    // basic setup
    double convolRadius(m_radiusConvol), convolThreshold(m_thresholdConvol); // todo: add user input
    Eigen::Vector3d left, right;
    getBBox(left, right);
    std::cout << "bounding box: " << left.transpose() << ", " << right.transpose()<< std::endl;
    // padding bounding box to attain robustness
    left += Eigen::Vector3d(-convolRadius, -convolRadius, -convolRadius);
    right += Eigen::Vector3d(convolRadius, convolRadius, convolRadius);

    // sort edges w.r.t. z-coordinate, and write to file
    long time = clock();
    std::sort(ptrLatticeModel->edgeArray.begin(), ptrLatticeModel->edgeArray.end(), [](auto& lhs, auto& rhs){
        double z1 = std::min(lhs->start->p[2], lhs->end->p[2]);
        double z2 = std::min(rhs->start->p[2], rhs->end->p[2]);
        return z1 < z2;
    });
    std::string pathLatticeGraph(volumeMeshFilePath + "_sortedGraph.txt");
    writeLatticeModel(true, pathLatticeGraph);
    clearLatticeModel(); // release memory
    std::cout << "LatticeModeler > Edge Sorting and writing Time: " << clock() - time << endl;


    // get layer thick
    double layerThick;
    std::cout << "Choose the thickness of layers: ";
    std::cin >> layerThick;

    // get resolution
    size_t resX, resY;
    std::cout << "Choose the resolution: ";
    std::cin >> resX >> resY;

    time = clock();

    // a functor for slicing at each iteration
    auto intersectEdgesPlane = [&](std::vector<Vector7d>& edges, double layerLevel, size_t layerNum){
        if (edges.empty()) return;
        // edges: 0-2: first point coordinates; 3-5: second point coordinates; 6: radius
        // setup the digital image size, resolution, and convolutional radius
        Eigen::Vector3d leftLocal, rightLocal;
        getBBox(leftLocal, rightLocal, edges);
        leftLocal += Eigen::Vector3d(-convolRadius, -convolRadius, 0);
        rightLocal += Eigen::Vector3d(convolRadius, convolRadius, 0);
        leftLocal[2] = rightLocal[2] = layerLevel;

        double cellSize(std::min((rightLocal[0] - leftLocal[0]) / min(resX, resY), (rightLocal[1] - leftLocal[1]) / min(resX, resY))); // todo: input from the user, or determined by the bar radius
        size_t xRes(std::ceil((rightLocal[0] - leftLocal[0]) / cellSize));
        size_t yRes(std::ceil((rightLocal[1] - leftLocal[1]) / cellSize));
        if (xRes < yRes) yRes = max(yRes, max(resX, resY));
        else xRes = max(xRes, max(resX, resY));
        std::cout << "LatticeModeler > digital image resoltuion at "<< layerNum << "-th layer :" << xRes << "*" << yRes << endl;
        std::unordered_map<std::string, double> digitImage; // string: "rowIdx+space+colIdx"

        // do convolution
        for (auto& e : edges) {
            double xMin, yMin, xMax, yMax;
            Eigen::Vector3d p0(e.segment(0,3)), p1(e.segment(3,3));
            // clamp the edge in case it's too long
            if (abs(p0[2] - p1[2]) > 1e-6) {
                Eigen::Vector3d p3(p0), p4(p1);
                double t = (layerLevel + convolRadius - p0[2]) / (p0[2] - p1[2]);
                if (t >= 0 && t <= 1) p4 = p0 + t * (p1 - p0);
                t = (layerLevel - convolRadius - p0[2]) / (p0[2] - p1[2]);
                if (t >= 0 && t <= 1) p3 = p0 + t * (p1 - p0);
                p0 = p3;
                p1 = p4;
            }

            double radius(e[6]);
            xMin = std::min(p0[0], p1[0]) - convolRadius;
            xMax = std::max(p0[0], p1[0]) + convolRadius;
            yMin = std::min(p0[1], p1[1]) - convolRadius;
            yMax = std::max(p0[1], p1[1]) + convolRadius;
            int gridIdxMinX(std::floor((xMin - leftLocal[0]) / cellSize));
            int gridIdxMinY(std::floor((yMin - leftLocal[1]) / cellSize));
            int gridIdxMaxX(std::ceil((xMax - leftLocal[0]) / cellSize));
            int gridIdxMaxY(std::ceil((yMax - leftLocal[1]) / cellSize));

            for (auto rowIdx(gridIdxMinX); rowIdx <= gridIdxMaxX; ++rowIdx) {
                for (auto colIdx(gridIdxMinY); colIdx <= gridIdxMaxY; ++colIdx) {
                    // compute convolution value for this point with regards to edge e
                    Eigen::Vector3d p(leftLocal[0] + cellSize * rowIdx, leftLocal[1] + cellSize * colIdx, layerLevel);
                    double l((p1 - p0).norm()), a((p0 - p).dot(p1 - p0)), b((p0 - p).norm());
                    // compute line-sphere intersection
                    l *= l; a *= 2; b = b * b - convolRadius * convolRadius;
                    double delta = a * a - 4 * l * b;
                    if (delta >= 1e-6) {
                        double intersect1 = (- a - std::sqrt(delta)) / (2 * l);
                        double intersect2 = (- a + std::sqrt(delta)) / (2 * l);
                        double start = std::max(0., std::min(intersect1, intersect2));
                        double end = std::min(1., std::max(intersect1, intersect2));
                        if (start <= 1 && end >= 0) {
                            l = std::sqrt(l); a /= -2; b += convolRadius * convolRadius; b = std::sqrt(b);
                            double f = radius / (15 * std::pow(convolRadius, 4)) *
                                    (3 * std::pow(l, 4) * (std::pow(end, 5) - std::pow(start, 5)) -
                                     15 * a * l * l * (std::pow(end, 4) - std::pow(start, 4)) +
                                     20 * a * a * (std::pow(end, 3) - std::pow(start, 3)));
                            auto key = std::to_string(rowIdx) + " " + std::to_string(colIdx);
                            auto iter = digitImage.find(key);
                            if (iter == digitImage.end()) digitImage.insert(std::make_pair(key, f)); // do convolution
                            else iter->second += f;
                        }
                    }
                    // end of convolution
                }
            }
        }

        // create a bindary image: exclude non-solid parts
        for (auto iter = digitImage.begin(), last = digitImage.end(); iter != last; ) {
            if (iter->second <= convolThreshold) iter = digitImage.erase(iter);
            else ++iter;
        }
        if(digitImage.empty()) return;

        // save the bindary image to file
        size_t width = xRes, height = yRes;
        QImage img(width, height, QImage::Format_RGB32);
        for (size_t i(0); i < width; ++i) {
            for (size_t j(0); j < height; ++j) {
                auto iter = digitImage.find(std::to_string(i) + " " + std::to_string(j));
                if (iter != digitImage.end()) img.setPixel(i, j, qRgb(0, 0, 0));
                else img.setPixel(i, j, qRgb(255, 255, 255));
            }
        }
        std::string path = volumeMeshFilePath + "." + std::to_string(layerNum) + ".slice.jpg";
        if (!img.save(QString::fromStdString(path), "JPEG", 50)) std::cerr << "LatticeModeler > error writing picture at layer: " << layerNum << endl;
        else std::cout << "LatticeModeler > successfully writing picture at layer: " << layerNum << endl;

        // do contouring with optimization
        MarchingSquareBin msb(digitImage, xRes, yRes, cellSize, layerLevel, leftLocal.data()); // binImage will be moved, not copy, into MarchingSquareBin
        msb.doContouring();
        auto result = msb.getContours(false);
        //.......................
        // do things with "result"

        //........................
        std::cout << "LatticeModeler > successfully generating contours at layer: " << layerNum << endl;
    };

    // steaming from file
    std::vector<Vector7d> edgeList; // to be used as a heap
    auto cmp = [](auto& left, auto& right) {return std::max(left[2], left[5]) > std::max(right[2], right[5]);}; // for max heap

    // setup slicing parameters
    double zMin(left[2]), zMax(right[2]);
    size_t layerNum(1);
    double layerLevel(zMin + layerThick);
    zMax -= layerThick / 10;

    std::ifstream ifs(pathLatticeGraph.c_str(), ifstream::in);
    std::string line; std::getline(ifs, line); // exclude the first line
    while (layerLevel < zMax) {
        // get new edges
        while (ifs.good() && !ifs.eof() && std::getline(ifs, line)) {
            Vector7d edge;
            stringstream str(line);
            str >> edge[0] >> edge[1] >> edge[2]; // point 1
            str >> edge[3] >> edge[4] >> edge[5]; // point 2
            str >> edge[6]; // radius

            // push into heap if the edge's zMin is below the current layerLevel
            double lowerBound = std::min(edge[2] - convolRadius, edge[5] - convolRadius);
            double upperBound =  std::max(edge[2] + convolRadius, edge[5] + convolRadius);
            if (upperBound >= layerLevel - tolerence) {
                edgeList.push_back(edge);
                std::push_heap(edgeList.begin(), edgeList.end(), cmp);
            }
            if (lowerBound >= layerLevel + tolerence) break;
        }

        // remove if the edge's zMax is below the current layerLevel
        while (!edgeList.empty()) {
            double upperBound = std::max(edgeList.front()[2] + convolRadius, edgeList.front()[5] + convolRadius);
            if (upperBound < layerLevel + tolerence)
            {std::pop_heap(edgeList.begin(),edgeList.end(), cmp); edgeList.pop_back();}
            else break;
        }
        if (edgeList.empty()) { layerLevel += layerThick; continue; }

        // do slicing
        time = clock();
        std::cout << "LatticeModeler > Slicing edge num: " << edgeList.size() << endl;
        intersectEdgesPlane(edgeList, layerLevel, layerNum);
        std::cout << "LatticeModeler > Slicing Time (micro second) for "<< layerNum << "-th layer:" << double(clock() - time) / CLOCKS_PER_SEC << endl;

        ++layerNum;
        layerLevel += layerThick;
    }
    ifs.close();
    std::cout << "LatticeModeler > successfully slicing the model" << std::endl;

    return true;
}

void LatticeModeler::save2Pic(std::vector<std::vector<double> > &digitImage, string path, int dpi)
{
    int width = digitImage.size();
    int height = digitImage[0].size();
    std::shared_ptr<QImage> img(new QImage(width, height, QImage::Format_RGB32));
    for (size_t i(0); i < width; ++i) {
        for (size_t j(0); j < height; ++j) {
            int v = 255 * digitImage[i][j];
            img->setPixel(i, j, qRgb(v, v, v));
        }
    }

    // save
    if (!img->save(QString::fromStdString(path), "JPG", 100)) {
        std::cerr << "LatticeModeler > error writing picture" << endl;
        return;
    }
}















