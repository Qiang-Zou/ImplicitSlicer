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
#include <vector>
#include<stack>
#include <experimental/unordered_map>

#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/remove.h>
#include <thrust/sequence.h>
#include <thrust/fill.h>
#include <thrust/count.h>
#include <thrust/functional.h>
#include <thrust/sequence.h>
#include <thrust/iterator/reverse_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/unique.h>
#include <cstdlib>

#include <QImage>
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cutil.h"
#include "LDNIcudaOperation.h"
#include "LDNIcudaSolid.h"
#include "call_cuda.h"
#define BLOCKSIZE	64
//#define VERSION2
#define VERSION3
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
    surfaceMeshFilePath("Results/in_surface.off"),
    volumeMeshFilePath("Results/out_volume"),
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
        ptrLatticeModel->addEdge(n1->idx, n2->idx, 0.01);
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
#ifdef VERSION3
bool LatticeModeler::sliceModelImplicitStream(std::string path)
{
	typedef Eigen::Matrix<double, 7, 1> Vector7d; // store edge

	//  read model
	if (readLatticeModel(path)) std::cout << "LatticeModeler > successfully read the model." << endl;
	if (ptrLatticeModel->isEmpty()) { std::cerr << "LatticeModeler > can't find a model" << endl; return false; }

	// basic setup
	double convolRadius(m_radiusConvol), convolThreshold(m_thresholdConvol); // todo: add user input
	Eigen::Vector3d left, right;
	getBBox(left, right);
	std::cout << "bounding box: " << left.transpose() << ", " << right.transpose() << std::endl;
	// padding bounding box to attain robustness
	left += Eigen::Vector3d(-3 * convolRadius, -3 * convolRadius, -3 * convolRadius);
	right += Eigen::Vector3d(3 * convolRadius, 3 * convolRadius, 3 * convolRadius);

	// sort edges w.r.t. z-coordinate, and write to file
	long time = clock();
	std::sort(ptrLatticeModel->edgeArray.begin(), ptrLatticeModel->edgeArray.end(), [](auto& lhs, auto& rhs) {
		double z1 = std::max(lhs->start->p[2], lhs->end->p[2]);
		double z2 = std::max(rhs->start->p[2], rhs->end->p[2]);
		return z1 > z2;
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
	// Support structure
	bool *gridNodes;
	//vector<vector<int>> binaryNodes[40];
	int imageSize[3] = { 0,0,0 };


	double cellSize = min((right[0] - left[0]) / min(resX, resY), (right[1] - left[1]) / min(resX, resY));
	resX = (right[0] - left[0]) / cellSize;
	resY = (right[1] - left[1]) / cellSize;
	imageSize[0] = resX; //cout << (right[0] - left[0]) << endl;
	imageSize[2] = resY; //cout << (right[1] - left[1]) << endl;

	int nodeNum = imageSize[0] * imageSize[2];
	bool* upperLayer;
	bool* lowerLayer;
	CUDA_SAFE_CALL(cudaMalloc((void**)&(upperLayer), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)upperLayer, false, nodeNum * sizeof(bool)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&(lowerLayer), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)lowerLayer, false, nodeNum * sizeof(bool)));

	bool* upperSupt;
	bool* lowerSupt;
	CUDA_SAFE_CALL(cudaMalloc((void**)&(upperSupt), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)upperSupt, false, nodeNum * sizeof(bool)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&(lowerSupt), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)lowerSupt, false, nodeNum * sizeof(bool)));
	int*** binaryNodes;
	binaryNodes = new int**[100];
	for (int i = 0; i < 100; i++)
	{
		binaryNodes[i] = new int*[resX];
		for (int j = 0; j < resX; j++)
		{
			binaryNodes[i][j] = new int[resY];
		}
	}
	
	//User specified variables
	double anchorR = 2;
	double thres = 0.2;
	double nSampleWidth = 0.005;
	double cylinderRadius = 6.0;
	double patternThickness = 3;

	bool *suptNodes;
	bool *tempImage;
	bool *tmpImg;
	bool *targetImage;
	bool *assistImage;
	bool *temp3D;


	double anchorRadius = anchorR;
	double threshold = thres;

	int suptRadius = (int)floor(anchorRadius*1.414 / nSampleWidth);
	int suptGridResX = (imageSize[0] - 1) / suptRadius + 1;
	int suptGridResZ = (imageSize[2] - 1) / suptRadius + 1;

	int3 imgRes = make_int3(imageSize[0], imageSize[1], imageSize[2]);
	int3 suptimgRes = make_int3(imageSize[0], imageSize[1] - 1, imageSize[2]);


	long ti[10] = { 0 };

	/*CUDA_SAFE_CALL(cudaMalloc((void**)&(temp3D), imageSize[0] * (imageSize[1] - 1)*imageSize[2] * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)temp3D, 0, imageSize[0] * (imageSize[1] - 1)*imageSize[2] * sizeof(bool)));


	CUDA_SAFE_CALL(cudaMalloc((void**)&(suptNodes), imageSize[0] * (imageSize[1] - 1)*imageSize[2] * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)suptNodes, false, imageSize[0] * (imageSize[1] - 1)*imageSize[2] * sizeof(bool)));*/

	CUDA_SAFE_CALL(cudaMalloc((void**)&(tempImage), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)tempImage, false, nodeNum * sizeof(bool)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&(targetImage), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)targetImage, false, nodeNum * sizeof(bool)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&(assistImage), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)assistImage, false, nodeNum * sizeof(bool)));


	short2 **disTextures = (short2 **)malloc(2 * sizeof(short2 *));

	int disTexSize = max(imageSize[0], imageSize[2]);
	int factor = ceil((float)disTexSize / BLOCKSIZE);
	disTexSize = BLOCKSIZE * factor;

	int disMemSize = disTexSize * disTexSize * sizeof(short2);

	// Allocate 2 textures

	cudaMalloc((void **)&disTextures[0], disMemSize);
	cudaMalloc((void **)&disTextures[1], disMemSize);

	unsigned int *LinkIndex;
	CUDA_SAFE_CALL(cudaMalloc((void **)&LinkIndex, (nodeNum + 1) * sizeof(unsigned int)));
	CUDA_SAFE_CALL(cudaMemset((void*)LinkIndex, 0, (nodeNum + 1) * sizeof(unsigned int)));

	long t;
	int flag = 3;
	stack<std::shared_ptr<QMeshPatch>> stk;

	unsigned int LinkNum;

	short *linkLayerD;
	short2 *linkID;
	unsigned int *linkLayerC;
	unsigned int *temp2D;

	CUDA_SAFE_CALL(cudaMalloc((void **)&temp2D, nodeNum * sizeof(unsigned int)));
	CUDA_SAFE_CALL(cudaMemset((void*)temp2D, 0, nodeNum * sizeof(unsigned int)));
	
	// a functor for slicing at each iteration
	auto intersectEdgesPlane = [&](std::vector<Vector7d>& edges, double layerLevel, size_t layerNum) {
		if (edges.empty()) return;
		// edges: 0-2: first point coordinates; 3-5: second point coordinates; 6: radius
		// setup the digital image size, resolution, and convolutional radius
		/*Eigen::Vector3d leftLocal, rightLocal;
		getBBox(leftLocal, rightLocal, edges);
		leftLocal += Eigen::Vector3d(-convolRadius, -convolRadius, 0);
		rightLocal += Eigen::Vector3d(convolRadius, convolRadius, 0);
		leftLocal[2] = rightLocal[2] = layerLevel;

		double cellSize(std::min((rightLocal[0] - leftLocal[0]) / min(resX, resY), (rightLocal[1] - leftLocal[1]) / min(resX, resY))); // todo: input from the user, or determined by the bar radius
		size_t xRes(std::ceil((rightLocal[0] - leftLocal[0]) / cellSize));
		size_t yRes(std::ceil((rightLocal[1] - leftLocal[1]) / cellSize));
		if (xRes < yRes) yRes = max(yRes, max(resX, resY));
		else xRes = max(xRes, max(resX, resY));
		std::cout << "LatticeModeler > digital image resoltuion at "<< layerNum << "-th layer :" << xRes << "*" << yRes << endl;*/
		//std::unordered_map<std::string, double> digitImage; // string: "rowIdx+space+colIdx"
		SCALAR* digitImage;
		CUDA_SAFE_CALL(cudaMalloc((void **)&digitImage, nodeNum * sizeof(SCALAR)));
		CUDA_SAFE_CALL(cudaMemset((void*)digitImage, 0.0, nodeNum * sizeof(SCALAR)));

		long time = clock();
		// do convolution
		for (auto& e : edges) {
			double xMin, yMin, xMax, yMax;
			Eigen::Vector3d p0(e.segment(0, 3)), p1(e.segment(3, 3));
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
			/*int gridIdxMinX(std::floor((xMin - leftLocal[0]) / cellSize));
			int gridIdxMinY(std::floor((yMin - leftLocal[1]) / cellSize));
			int gridIdxMaxX(std::ceil((xMax - leftLocal[0]) / cellSize));
			int gridIdxMaxY(std::ceil((yMax - leftLocal[1]) / cellSize));*/

			int gridIdxMinX(std::floor((xMin - left[0]) / cellSize));
			int gridIdxMinY(std::floor((yMin - left[1]) / cellSize));
			int gridIdxMaxX(std::ceil((xMax - left[0]) / cellSize));
			int gridIdxMaxY(std::ceil((yMax - left[1]) / cellSize));

			int pNum = (gridIdxMaxX - gridIdxMinX)*(gridIdxMaxY- gridIdxMinY);
			int2 XY = make_int2((gridIdxMaxX - gridIdxMinX), (gridIdxMaxY - gridIdxMinY));
			call_DoConvolution(digitImage, gridIdxMinX, gridIdxMinY, cellSize, pNum, XY, imgRes, left[0], left[1], p0, p1, convolRadius, radius,layerLevel);
			/*for (auto rowIdx(gridIdxMinX); rowIdx <= gridIdxMaxX; ++rowIdx) {
				for (auto colIdx(gridIdxMinY); colIdx <= gridIdxMaxY; ++colIdx) {
					// compute convolution value for this point with regards to edge e
					Eigen::Vector3d p(left[0] + cellSize * rowIdx, left[1] + cellSize * colIdx, layerLevel);
					double l((p1 - p0).norm()), a((p0 - p).dot(p1 - p0)), b((p0 - p).norm());
					// compute line-sphere intersection
					l *= l; a *= 2; b = b * b - convolRadius * convolRadius;
					double delta = a * a - 4 * l * b;
					if (delta >= 1e-6) {
						double intersect1 = (-a - std::sqrt(delta)) / (2 * l);
						double intersect2 = (-a + std::sqrt(delta)) / (2 * l);
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
							else iter->second += f; //summarization for convolution
						}
					}
					// end of convolution
				}
			}*/
		}

		// create a bindary image: exclude non-solid parts
		/*for (auto iter = digitImage.begin(), last = digitImage.end(); iter != last; ) {
			if (iter->second <= convolThreshold) iter = digitImage.erase(iter);
			else ++iter;
		}
		if (digitImage.empty()) return;*/
		std::cout << "LatticeModeler > DoConvolution Time (micro second) for " << layerNum << "-th layer:" << double(clock() - time) / CLOCKS_PER_SEC << endl;
		// Support structure
		/*for (size_t i = 0; i < imageSize[0]; i++)
		{

			for (size_t j = 0; j < imageSize[2]; j++)
			{
				auto iter = digitImage.find(std::to_string(i) + " " + std::to_string(j));
				if (iter != digitImage.end())
				{
					binaryNodes[layerNum][i][j] = 1;
					//binaryNodes[layerNum].push_back(1);
				}
				else
				{
					binaryNodes[layerNum][i][j] = 0;
					//binaryNodes[layerNum].push_back(0);
				}
			}
		}*/

		// save the bindary image to file
		SCALAR* hostDigit = new SCALAR[nodeNum];
		CUDA_SAFE_CALL(cudaMemcpy((void*)hostDigit, (void*)digitImage, nodeNum * sizeof(SCALAR), cudaMemcpyDeviceToHost));

		bool* tmp = new bool[nodeNum];
		int cnt = 0;
		size_t width = resX, height = resY;
		QImage img(width, height, QImage::Format_RGB32);
		for (size_t j(0); j < height; ++j) {
			for (size_t i(0); i < width; ++i) {
				//auto iter = digitImage.find(std::to_string(i) + " " + std::to_string(j));
				if (hostDigit[i+j*imgRes.x]> convolThreshold)
				{
					img.setPixel(i, j, qRgb(0, 0, 0));
					//binaryNodes[layerNum][i][j] = 1;
					tmp[cnt++] = 1;
				}
				else
				{
					img.setPixel(i, j, qRgb(255, 255, 255));
					//binaryNodes[layerNum][i][j] = 0;
					tmp[cnt++] = 0;
				}
			}
		}
		std::string path = volumeMeshFilePath + "." + std::to_string(layerNum) + ".slice.jpg";
		if (!img.save(QString::fromStdString(path), "JPEG", 50)) std::cerr << "LatticeModeler > error writing picture at layer: " << layerNum << endl;
		else std::cout << "LatticeModeler > successfully writing picture at layer: " << layerNum << endl;

		// do contouring with optimization
		/*MarchingSquareBin msb(digitImage, xRes, yRes, cellSize, layerLevel, leftLocal.data()); // binImage will be moved, not copy, into MarchingSquareBin
		msb.doContouring();
		auto result = msb.getContours(false);*/
		//.......................
		// do things with "result"
		if (flag == 2)
		{
			long time = clock();
			//........................
			std::unordered_map<std::string, double> DImage;
			for (int i = 0; i < resX; i++)
			{
				for (int j = 0; j < resY; j++)
				{
					double f = hostDigit[i+j*imgRes.x];
					if (f > 0)
					{
						auto key = std::to_string(i) + " " + std::to_string(j);
						DImage.insert(std::make_pair(key, f));
					}
				}
			}
			MarchingSquareBin msb(DImage, resX, resY, cellSize, layerLevel, left.data());
			msb.doContouring();
			path = volumeMeshFilePath + "." + std::to_string(layerNum) + ".contours.off";
			//msb.writeContours(path);
			auto result = msb.getContours(true);
			stk.push(std::shared_ptr<QMeshPatch>(result));
			//layers.push_back(std::shared_ptr<QMeshPatch>(result));
			//auto result2 = msb.getContours(true);
			//std::string clipath = "Results/Contours/out_volume.contours1.cli";
			//writeCLIFile(clipath, layers);
			std::cout << "LatticeModeler > successfully generating contours at layer: " << layerNum << endl;
			std::cout << "LatticeModeler > MarchingSquare Time (micro second) for " << layerNum << "-th layer:" << double(clock() - time) / CLOCKS_PER_SEC << endl;
		}
		

		if (layerNum == 1)
		{
			CUDA_SAFE_CALL(cudaMemcpy((void*)upperLayer, (void*)tmp, nodeNum * sizeof(bool), cudaMemcpyHostToDevice));
			CUDA_SAFE_CALL(cudaMemcpy((void*)lowerLayer, (void*)tmp, nodeNum * sizeof(bool), cudaMemcpyHostToDevice));
		}
		else
		{
			CUDA_SAFE_CALL(cudaMemcpy((void*)upperLayer, (void*)lowerLayer, nodeNum * sizeof(bool), cudaMemcpyDeviceToDevice));
			CUDA_SAFE_CALL(cudaMemcpy((void*)lowerLayer, (void*)tmp, nodeNum * sizeof(bool), cudaMemcpyHostToDevice));
		}
		if (layerNum >= 2)
		{
			t = clock();
			//Region Subtraction
			//if(flag==1)
			//call_krSLAContouring_Initialization(tempImage, targetImage, gridNodes, nodeNum, imgRes, imageSize[1]-layerNum);
			//else
			call_RegionSubtractionSLA(tempImage, targetImage, upperLayer, lowerLayer, nodeNum, imgRes);
			CUDA_SAFE_CALL(cudaMemcpy(assistImage, targetImage, nodeNum * sizeof(bool), cudaMemcpyDeviceToDevice));
			ti[0] += clock() - t;

			t = clock();
			//Growth and Swallow for Si'
			LDNIcudaOperation::LDNISLAContouring_GrowthAndSwallow(threshold, targetImage, tempImage, layerNum, imageSize, nSampleWidth, disTextures[0], disTextures[1], disTexSize);

			ti[1] += clock() - t;
			//add new support cylinder if necessary
			CUDA_SAFE_CALL(cudaMemset((void*)tempImage, false, nodeNum * sizeof(bool)));
			CUDA_SAFE_CALL(cudaMemset((void*)assistImage, false, nodeNum * sizeof(bool)));
			call_krSLAContouring_Filter1(assistImage, tempImage, targetImage, suptGridResX*suptGridResZ, make_int2(imageSize[0], imageSize[2]), suptRadius, layerNum);

			//first step: support region growth in first class cylinders
			LDNIcudaOperation::LDNISLAContouring_GrowthAndSwallow(anchorRadius, targetImage, assistImage, layerNum, imageSize, nSampleWidth, disTextures[0], disTextures[1], disTexSize);

			//second step: prepare second class cylinders and perform support region growth in second class cylinders
			CUDA_SAFE_CALL(cudaMemset((void*)assistImage, false, nodeNum * sizeof(bool)));
			call_krSLAContouring_OrthoSearchRemainAnchorZ(assistImage, tempImage, targetImage, suptGridResX,
				make_int2(suptRadius, imageSize[2]), make_int2(imageSize[0], imageSize[2]), layerNum);

			call_krSLAContouring_OrthoSearchRemainAnchorX(assistImage, tempImage, targetImage, suptGridResZ,
				make_int2(suptRadius, imageSize[0]), make_int2(imageSize[0], imageSize[2]), layerNum);



			LDNIcudaOperation::LDNISLAContouring_GrowthAndSwallow(anchorRadius, targetImage, assistImage, layerNum, imageSize, nSampleWidth, disTextures[0], disTextures[1], disTexSize);
			long time = clock();

			//third step: prepare third class cylinders and support region growth in all third class cylinders
			CUDA_SAFE_CALL(cudaMemset((void*)assistImage, false, nodeNum * sizeof(bool)));
			LDNIcudaOperation::LDNISLAContouring_ThirdClassCylinder(anchorRadius, targetImage, assistImage, tempImage, make_int2(imageSize[0], imageSize[2]), nSampleWidth, layerNum, disTextures[0], disTextures[1], disTexSize);
			
			
			
			//generate support structure map for this layer and update the cylinder position information arrays
			/*if (flag == 1)
			{
				call_krSLAContouring_Filter5(gridNodes, tempImage, suptNodes, LinkIndex, nodeNum, imgRes, imageSize[1] - layerNum);
			}*/
			call_Filter5forSLA(upperLayer, lowerLayer, tempImage, upperSupt, lowerSupt, LinkIndex, nodeNum, imgRes, layerNum,flag);
			/*if (flag == 1)
			{
				call_krFDMContouring_CopyNodesrom2Dto3D(tempImage, temp3D, nodeNum, suptimgRes, imageSize[1] - layerNum);
			}*/
			if(flag==1)
			{
				//call_krFDMContouring_CopyNodesrom2Dto3D(tempImage, temp3D, nodeNum, suptimgRes, i);

				call_myFindAllLinks(LinkIndex, linkLayerC, linkLayerD, linkID, tempImage, temp2D, nodeNum, suptimgRes, imageSize[1]-layerNum);
				//call_krSLAContouring_FindAllLinks(LinkIndex, linkLayerC, linkLayerD, linkID, temp3D, temp2D, suptimgRes.x*suptimgRes.y*suptimgRes.z, suptimgRes);
			
				call_myRelateAllLinksBetweenLayers(LinkIndex, linkLayerC, lowerLayer, upperSupt, nodeNum, imgRes, layerNum, imageSize[1] - layerNum);
				//call_krSLAContouring_RelateAllLinksBetweenLayers(LinkIndex, linkLayerC, gridNodes, suptNodes, suptimgRes.x*suptimgRes.y*suptimgRes.z, imgRes);

			}
			bool *uppertmp;
			CUDA_SAFE_CALL(cudaMalloc((void**)&(uppertmp), nodeNum * sizeof(bool)));
			CUDA_SAFE_CALL(cudaMemset((void*)uppertmp, false, nodeNum * sizeof(bool)));
			CUDA_SAFE_CALL(cudaMemcpy(uppertmp, upperSupt, nodeNum * sizeof(bool), cudaMemcpyDeviceToDevice));

			CUDA_SAFE_CALL(cudaMemcpy((void*)upperSupt, (void*)lowerSupt, nodeNum * sizeof(bool), cudaMemcpyDeviceToDevice));

			if (flag == 0)
			{

				double realThreshold = (cylinderRadius*nSampleWidth - nSampleWidth) / nSampleWidth;
				int gridRadius = (int)floor(realThreshold);
				CUDA_SAFE_CALL(cudaMalloc((void**)&(tmpImg), nodeNum * sizeof(bool)));
				CUDA_SAFE_CALL(cudaMemset((void*)tmpImg, false, nodeNum * sizeof(bool)));
				CUDA_SAFE_CALL(cudaMemcpy(tmpImg, lowerSupt, nodeNum * sizeof(bool), cudaMemcpyDeviceToDevice));
				call_krFDMContouring_Dilation(lowerSupt, tmpImg, nodeNum, imgRes, realThreshold, gridRadius, layerNum);
				CUDA_SAFE_CALL(cudaMemcpy(lowerSupt, tmpImg, nodeNum * sizeof(bool), cudaMemcpyDeviceToDevice));

				CUDA_SAFE_CALL(cudaMemset((void*)tmpImg, false, nodeNum * sizeof(bool)));
				CUDA_SAFE_CALL(cudaMemcpy(tmpImg, uppertmp, nodeNum * sizeof(bool), cudaMemcpyDeviceToDevice));
				call_krFDMContouring_Dilation(uppertmp, tmpImg, nodeNum, imgRes, realThreshold, gridRadius, layerNum);
				CUDA_SAFE_CALL(cudaMemcpy(uppertmp, tmpImg, nodeNum * sizeof(bool), cudaMemcpyDeviceToDevice));
				if (layerNum >= 3)
				{
					call_myVerticalSpptPxlProp(lowerLayer, uppertmp, lowerSupt, nodeNum, imgRes);
					//call_krFDMContouring_VerticalSpptPxlProp(gridNodes, suptNodes, nodeNum, imgRes, i);
				}
				int linkThreshold = (int)(0.4 / nSampleWidth);
				int lengthofLayer = linkThreshold / 8;
				int furtherStepLength = lengthofLayer / 2;

				LDNIcudaOperation::LDNISLAContouring_GenerateConnectionforCylinders2(linkLayerC, linkLayerD, linkID, upperLayer,lowerLayer,lowerSupt, imageSize, linkThreshold,
					lengthofLayer, furtherStepLength, LinkNum, nSampleWidth,imageSize[1]-layerNum);

				realThreshold = (patternThickness*nSampleWidth - nSampleWidth) / nSampleWidth;
				gridRadius = (int)floor(realThreshold);

				CUDA_SAFE_CALL(cudaMemset((void*)tmpImg, false, nodeNum * sizeof(bool)));
				CUDA_SAFE_CALL(cudaMemcpy(tmpImg, lowerSupt, nodeNum * sizeof(bool), cudaMemcpyDeviceToDevice));
				call_krFDMContouring_Dilation(lowerSupt, tmpImg, nodeNum, imgRes, realThreshold, gridRadius, layerNum);
				CUDA_SAFE_CALL(cudaMemcpy(lowerSupt, tmpImg, nodeNum * sizeof(bool), cudaMemcpyDeviceToDevice));
			}

			
			//std::cout << "Generate support map time for " << i + 1 << "-th layer:" << double(clock() - time) / CLOCKS_PER_SEC << endl;
			if (flag == 0)
			{
				bool *supttmp = new bool[nodeNum];
				CUDA_SAFE_CALL(cudaMemcpy((void*)supttmp, (void*)lowerSupt, nodeNum * sizeof(bool), cudaMemcpyDeviceToHost));
				cnt = 0;
				QImage suptimg(width, height, QImage::Format_RGB32);
				for (size_t j(0); j < height; ++j)
				{
					for (size_t i(0); i < width; ++i)
					{
						if (supttmp[cnt++])suptimg.setPixel(i, j, qRgb(255, 0, 0));
						else suptimg.setPixel(i, j, qRgb(255, 255, 255));
					}
				}
				path = volumeMeshFilePath + "." + std::to_string(layerNum-1) + ".support.jpg";
				if (!suptimg.save(QString::fromStdString(path), "JPEG", 50)) std::cerr << "Support > error writing picture at layer: " << layerNum-1 << endl;
				else std::cout << "Support > successfully writing picture at layer: " << layerNum-1 << endl;
			}
			
		}
		

	};

	// steaming from file
	std::vector<Vector7d> edgeList; // to be used as a heap
	auto cmp = [](auto& left, auto& right) {return std::min(left[2], left[5]) < std::min(right[2], right[5]); }; // for min heap

	// setup slicing parameters
	double zMin(left[2]), zMax(right[2]);
	size_t layerNum(1);
	double layerLevel(zMax - layerThick);
	zMin += layerThick / 10;

	imageSize[1] = (int)((zMax - zMin) / layerThick);
	std::cout << "layer number: " << imageSize[1] << endl;
	imgRes = make_int3(imageSize[0], imageSize[1], imageSize[2]);
	suptimgRes = make_int3(imageSize[0], imageSize[1] - 1, imageSize[2]);

	/*CUDA_SAFE_CALL(cudaMalloc((void**)&(temp3D), imageSize[0] * (imageSize[1] - 1)*imageSize[2] * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)temp3D, 0, imageSize[0] * (imageSize[1] - 1)*imageSize[2] * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(gridNodes), imageSize[0] * imageSize[1] * imageSize[2] * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)gridNodes, 0.0, imageSize[0] * imageSize[1] * imageSize[2] * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(suptNodes), imageSize[0] * (imageSize[1] - 1)*imageSize[2] * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)suptNodes, false, imageSize[0] * (imageSize[1] - 1)*imageSize[2] * sizeof(bool)));*/
	while (flag--)
	{
		std::ifstream ifs(pathLatticeGraph.c_str(), ifstream::in);
		std::string line; std::getline(ifs, line); // exclude the first line
		while (layerLevel > zMin) {
			// get new edges
			while (ifs.good() && !ifs.eof() && std::getline(ifs, line)) {
				Vector7d edge;
				stringstream str(line);
				str >> edge[0] >> edge[1] >> edge[2]; // point 1
				str >> edge[3] >> edge[4] >> edge[5]; // point 2
				str >> edge[6]; // radius

				// push into heap if the edge's zMin is below the current layerLevel
				double lowerBound = std::min(edge[2] - convolRadius, edge[5] - convolRadius);
				double upperBound = std::max(edge[2] + convolRadius, edge[5] + convolRadius);
				if (lowerBound <= layerLevel + tolerence) {
					edgeList.push_back(edge);
					std::push_heap(edgeList.begin(), edgeList.end(), cmp);
				}
				if (upperBound <= layerLevel - tolerence) break;
			}

			// remove if the edge's zMax is below the current layerLevel
			while (!edgeList.empty()) {
				double lowerBound = std::min(edgeList.front()[2] - convolRadius, edgeList.front()[5] - convolRadius);
				if (lowerBound > layerLevel - tolerence)
				{
					std::pop_heap(edgeList.begin(), edgeList.end(), cmp); edgeList.pop_back();
				}
				else break;
			}
			if (edgeList.empty()) { layerLevel -= layerThick; continue; }

			// do slicing
			time = clock();
			std::cout << "LatticeModeler > Slicing edge num: " << edgeList.size() << endl;
			intersectEdgesPlane(edgeList, layerLevel, layerNum);
			std::cout << "LatticeModeler > Slicing and MarchingSquare Time (micro second) for " << layerNum << "-th layer:" << double(clock() - time) / CLOCKS_PER_SEC << endl;

			++layerNum;
			layerLevel -= layerThick;
		}
		ifs.close();
		//imageSize[1] = layerNum - 1;
		
		/*if (flag == 2)
		{
			

			int idx = 0;
			bool* tmp2 = new bool[imageSize[0] * imageSize[1] * imageSize[2]];
			for (int j = 0; j < imageSize[2]; j++)
			{
				for (int k = imageSize[1]; k >=1; k--)
				{

					for (int i = 0; i < imageSize[0]; i++)
					{
						//gridNodes[idx] = binaryNodes[i][j];
						if (binaryNodes[k][i][j] == 1)
							tmp2[idx] = true;
						else
							tmp2[idx] = false;
						idx++;
					}

				}
			}
			CUDA_SAFE_CALL(cudaMemcpy((void*)gridNodes, (void*)tmp2, imageSize[0] * imageSize[1] * imageSize[2] * sizeof(bool), cudaMemcpyHostToDevice));
		}*/
		if (flag == 2)
		{
			call_func::setdev_ptr(LinkIndex, nodeNum, LinkNum);
			CUDA_SAFE_CALL(cudaMalloc((void **)&linkLayerC, LinkNum * sizeof(unsigned int)));
			//CUDA_SAFE_CALL(cudaMemset((void*)linkLayerC, imgRes.y+1, LinkNum*sizeof(unsigned int) ) );
			CUDA_SAFE_CALL(cudaMalloc((void **)&linkLayerD, LinkNum * sizeof(short)));
			CUDA_SAFE_CALL(cudaMemset((void*)linkLayerD, 0, LinkNum * sizeof(short)));
			CUDA_SAFE_CALL(cudaMalloc((void **)&linkID, LinkNum * sizeof(short2)));
			CUDA_SAFE_CALL(cudaMemset((void*)linkID, 0, LinkNum * sizeof(short2)));
			call_func::setdev_ptr2(linkLayerC, LinkNum, imgRes);
			//call_krSLAContouring_FindAllLinks(LinkIndex, linkLayerC, linkLayerD, linkID, temp3D, temp2D, suptimgRes.x*suptimgRes.y*suptimgRes.z, suptimgRes);
			//call_krSLAContouring_RelateAllLinksBetweenLayers(LinkIndex, linkLayerC, gridNodes, suptNodes, suptimgRes.x*suptimgRes.y*suptimgRes.z, imgRes);
		}

		edgeList.clear();
		layerNum = 1;
		layerLevel = zMax - layerThick;
	}
	
	std::cout << "LatticeModeler > successfully slicing the model" << std::endl;
	// write to .cli files
	int stkSize = stk.size();
	for (int i = 0; i < stkSize; i++)
	{
		std::shared_ptr<QMeshPatch> result = stk.top();
		layers.push_back(result);
		stk.pop();
	}
	string clipath = "Results/CLIFileforPart/out_volume.contours.cli";
	//writeCLIFileBin(path, nullptr, 0, 3); // 3: finish writing the cli file
	//    readCLIFileBin(path);
	writeCLIFileBin(clipath, layers);
	clipath = "Results/CLIFileforPart/out_volume.contours1.cli";
	writeCLIFile(clipath, layers);

	// Support structure
	/*imageSize[1] = layerNum - 1;
	CUDA_SAFE_CALL(cudaMalloc((void**)&(gridNodes), imageSize[0] * imageSize[1] * imageSize[2] * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)gridNodes, 0.0, imageSize[0] * imageSize[1] * imageSize[2] * sizeof(bool)));

	int idx = 0;
	bool* tmp = new bool[imageSize[0] * imageSize[1] * imageSize[2]];
	for (int j = 0; j < imageSize[2]; j++)
	{
		for (int k = 1; k <= imageSize[1]; k++)
		{

			for (int i = 0; i < imageSize[0]; i++)
			{
				//gridNodes[idx] = binaryNodes[i][j];
				if (binaryNodes[k][i][j] == 1)
					tmp[idx] = true;
				else
					tmp[idx] = false;
				idx++;
			}

		}
	}
	CUDA_SAFE_CALL(cudaMemcpy((void*)gridNodes, (void*)tmp, imageSize[0] * imageSize[1] * imageSize[2] * sizeof(bool), cudaMemcpyHostToDevice));
	

	for (int i = imageSize[1] - 2; i > -1; i--)
	{
		
	}

	cudaFree(disTextures[0]);
	cudaFree(disTextures[1]);
	free(disTextures);
	cudaFree(assistImage);
	cudaFree(targetImage);
	cudaFree(tempImage);

	
	


	


	

	

	cudaFree(temp3D);
	cudaFree(temp2D);



	

	cudaFree(LinkIndex);
	std::cout << "Relate all links time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;



	CUDA_SAFE_CALL(cudaMalloc((void**)&(tempImage), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)tempImage, false, nodeNum * sizeof(bool)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&(targetImage), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)targetImage, false, nodeNum * sizeof(bool)));


	double realThreshold = (cylinderRadius*nSampleWidth - nSampleWidth) / nSampleWidth;
	int gridRadius = (int)floor(realThreshold);
	time = clock();
	for (int i = imageSize[1] - 2; i > -1; i--)
	{

		call_krFDMContouring_CopyNodesrom3Dto2D(targetImage, suptNodes, nodeNum, suptimgRes, i);
		CUDA_SAFE_CALL(cudaMemcpy(tempImage, targetImage, nodeNum * sizeof(bool), cudaMemcpyDeviceToDevice));
		call_krFDMContouring_Dilation(targetImage, tempImage, nodeNum, imgRes, realThreshold, gridRadius, i);


		call_krFDMContouring_CopyNodesrom2Dto3D(tempImage, suptNodes, nodeNum, suptimgRes, i);

	}
	std::cout << "Dilation time:" << double(clock() - time) / CLOCKS_PER_SEC << endl;
	cudaFree(targetImage);
	cudaFree(tempImage);
	time = clock();
	for (int i = imageSize[1] - 3; i > -1; i--)
		call_krFDMContouring_VerticalSpptPxlProp(gridNodes, suptNodes, nodeNum, imgRes, i);
	std::cout << "VerticalSpptPxlProp time:" << double(clock() - time) / CLOCKS_PER_SEC << endl;

	int linkThreshold = (int)(0.4 / nSampleWidth);
	int lengthofLayer = linkThreshold / 8;
	int furtherStepLength = lengthofLayer / 2;



	time = clock();
	LDNIcudaOperation::LDNISLAContouring_GenerateConnectionforCylinders(linkLayerC, linkLayerD, linkID, gridNodes, suptNodes, imageSize, linkThreshold,
		lengthofLayer, furtherStepLength, LinkNum, nSampleWidth);
	std::cout << "Generate Connection for Cylinders time:" << double(clock() - time) / CLOCKS_PER_SEC << endl;
	cudaFree(linkLayerC);
	cudaFree(linkLayerD);
	cudaFree(linkID);
	//cudaFree(gridNodes);


	realThreshold = (patternThickness*nSampleWidth - nSampleWidth) / nSampleWidth;
	gridRadius = (int)floor(realThreshold);

	CUDA_SAFE_CALL(cudaMalloc((void**)&(tempImage), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)tempImage, false, nodeNum * sizeof(bool)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&(targetImage), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)targetImage, false, nodeNum * sizeof(bool)));
	time = clock();
	for (int i = imageSize[1] - 2; i > -1; i--)
	{

		call_krFDMContouring_CopyNodesrom3Dto2D(targetImage, suptNodes, nodeNum, suptimgRes, i);
		CUDA_SAFE_CALL(cudaMemcpy(tempImage, targetImage, nodeNum * sizeof(bool), cudaMemcpyDeviceToDevice));
		call_krFDMContouring_Dilation(targetImage, tempImage, nodeNum, imgRes, realThreshold, gridRadius, i);
		call_krFDMContouring_CopyNodesrom2Dto3D(tempImage, suptNodes, nodeNum, suptimgRes, i);

	}
	std::cout << "Dilation time:" << double(clock() - time) / CLOCKS_PER_SEC << endl;
	cudaFree(targetImage);
	cudaFree(tempImage);


	bool* gridtmp = new bool[imageSize[0] * imageSize[1] * imageSize[2]];
	CUDA_SAFE_CALL(cudaMemcpy((void*)gridtmp, (void*)gridNodes, imageSize[0] * imageSize[1] * imageSize[2] * sizeof(bool), cudaMemcpyDeviceToHost));
	int cnt = 0;
	size_t width = imageSize[0], height = imageSize[2], slicenum = imageSize[1];
	QImage *img = new QImage[slicenum];
	for (int i = 0; i < slicenum; i++)
	{
		QImage tmpimg(width, height, QImage::Format_RGB32);
		img[i] = tmpimg;
	}
	for (int j = 0; j < imageSize[2]; j++)
	{

		for (int k = 0; k < imageSize[1]; k++)
		{

			for (int i = 0; i < imageSize[0]; i++)
			{
				if (gridtmp[cnt++])img[k].setPixel(i, j, qRgb(0, 0, 0));
				else img[k].setPixel(i, j, qRgb(255, 255, 255));

			}
		}

	}
	for (int i = 1; i <= imageSize[1]; i++)
	{
		std::string path = "Results/Support/out_volume." + std::to_string(i) + ".origin.jpg";
		if (!img[i - 1].save(QString::fromStdString(path), "JPEG", 50)) std::cerr << "Origin > error writing picture at layer: " << i << endl;
		else std::cout << "Origin > successfully writing picture at layer: " << i << endl;
	}


	std::unordered_map<std::string, double> *digitImage = new std::unordered_map<std::string, double>[slicenum - 1];//generate contours for support structure
	std::vector<std::shared_ptr<QMeshPatch>> suptlayers;
	layerLevel = zMin + layerThick;
	bool* supttmp = new bool[imageSize[0] * (imageSize[1] - 1)*imageSize[2]];
	CUDA_SAFE_CALL(cudaMemcpy((void*)supttmp, (void*)suptNodes, imageSize[0] * (imageSize[1] - 1) * imageSize[2] * sizeof(bool), cudaMemcpyDeviceToHost));
	cnt = 0;
	QImage *suptimg = new QImage[slicenum - 1];
	for (int i = 0; i < slicenum - 1; i++)
	{
		QImage tmpimg(width, height, QImage::Format_RGB32);
		suptimg[i] = tmpimg;
	}
	for (int j = 0; j < imageSize[2]; j++)
	{
		for (int k = 0; k < imageSize[1] - 1; k++)
		{

			for (int i = 0; i < imageSize[0]; i++)
			{
				if (supttmp[cnt++])
				{
					suptimg[k].setPixel(i, j, qRgb(255, 0, 0));
					auto key = std::to_string(i) + " " + std::to_string(j);
					digitImage[k].insert(std::make_pair(key, 1));
				}
				else suptimg[k].setPixel(i, j, qRgb(255, 255, 255));

			}
		}

	}
	for (int i = 1; i <= imageSize[1] - 1; i++)
	{
		std::string path = "Results/Support/out_volume." + std::to_string(i) + ".support.jpg";
		if (!suptimg[i - 1].save(QString::fromStdString(path), "JPEG", 50)) std::cerr << "Support > error writing picture at layer: " << i << endl;
		else std::cout << "Support > successfully writing picture at layer: " << i << endl;

		long time = clock();
		MarchingSquareBin msb(digitImage[i - 1], resX, resY, cellSize, layerLevel, left.data());
		msb.doContouring();
		path = "Results/Support/out_volume." + std::to_string(i) + ".contours.off";
		//msb.writeContours(path);
		auto result = msb.getContours(true);
		suptlayers.push_back(std::shared_ptr<QMeshPatch>(result));
		layerLevel += layerThick;
		std::cout << "Support >  MarchingSquare Time (micro second) for " << i << "-th layer:" << double(clock() - time) / CLOCKS_PER_SEC << endl;
	}
	// write to .cli files
	string suptclipath = "Results/Support/CLIFileforSupt/out_volume.contours.cli";
	writeCLIFileBin(suptclipath, suptlayers);
	suptclipath = "Results/Support/CLIFileforSupt/out_volume.contours1.cli";
	writeCLIFile(suptclipath, suptlayers);*/
	return true;
}
#endif

#ifdef VERSION2
bool LatticeModeler::sliceModelImplicitStream(std::string path)
{
	typedef Eigen::Matrix<double, 7, 1> Vector7d; // store edge

	//  read model
	if (readLatticeModel(path)) std::cout << "LatticeModeler > successfully read the model." << endl;
	if (ptrLatticeModel->isEmpty()) { std::cerr << "LatticeModeler > can't find a model" << endl; return false; }

	// basic setup
	double convolRadius(m_radiusConvol), convolThreshold(m_thresholdConvol); // todo: add user input
	Eigen::Vector3d left, right;
	getBBox(left, right);
	std::cout << "bounding box: " << left.transpose() << ", " << right.transpose() << std::endl;
	// padding bounding box to attain robustness
	left += Eigen::Vector3d(-3 * convolRadius, -3 * convolRadius, -3 * convolRadius);
	right += Eigen::Vector3d(3 * convolRadius, 3 * convolRadius, 3 * convolRadius);

	// sort edges w.r.t. z-coordinate, and write to file
	long time = clock();
	std::sort(ptrLatticeModel->edgeArray.begin(), ptrLatticeModel->edgeArray.end(), [](auto& lhs, auto& rhs) {
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
	// Support structure
	bool *gridNodes;
	//vector<vector<int>> binaryNodes[40];
	int imageSize[3] = { 0,0,0 };


	double cellSize = min((right[0] - left[0]) / min(resX, resY), (right[1] - left[1]) / min(resX, resY));
	resX = (right[0] - left[0]) / cellSize;
	resY = (right[1] - left[1]) / cellSize;
	imageSize[0] = resX; //cout << (right[0] - left[0]) << endl;
	imageSize[2] = resY; //cout << (right[1] - left[1]) << endl;
	int*** binaryNodes;
	binaryNodes = new int**[100];
	for (int i = 0; i < 100; i++)
	{
		binaryNodes[i] = new int*[resX];
		for (int j = 0; j < resX; j++)
		{
			binaryNodes[i][j] = new int[resY];
		}
	}
	// a functor for slicing at each iteration
	auto intersectEdgesPlane = [&](std::vector<Vector7d>& edges, double layerLevel, size_t layerNum) {
		if (edges.empty()) return;
		// edges: 0-2: first point coordinates; 3-5: second point coordinates; 6: radius
		// setup the digital image size, resolution, and convolutional radius
		/*Eigen::Vector3d leftLocal, rightLocal;
		getBBox(leftLocal, rightLocal, edges);
		leftLocal += Eigen::Vector3d(-convolRadius, -convolRadius, 0);
		rightLocal += Eigen::Vector3d(convolRadius, convolRadius, 0);
		leftLocal[2] = rightLocal[2] = layerLevel;

		double cellSize(std::min((rightLocal[0] - leftLocal[0]) / min(resX, resY), (rightLocal[1] - leftLocal[1]) / min(resX, resY))); // todo: input from the user, or determined by the bar radius
		size_t xRes(std::ceil((rightLocal[0] - leftLocal[0]) / cellSize));
		size_t yRes(std::ceil((rightLocal[1] - leftLocal[1]) / cellSize));
		if (xRes < yRes) yRes = max(yRes, max(resX, resY));
		else xRes = max(xRes, max(resX, resY));
		std::cout << "LatticeModeler > digital image resoltuion at "<< layerNum << "-th layer :" << xRes << "*" << yRes << endl;*/
		std::unordered_map<std::string, double> digitImage; // string: "rowIdx+space+colIdx"


		// do convolution
		for (auto& e : edges) {
			double xMin, yMin, xMax, yMax;
			Eigen::Vector3d p0(e.segment(0, 3)), p1(e.segment(3, 3));
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
			/*int gridIdxMinX(std::floor((xMin - leftLocal[0]) / cellSize));
			int gridIdxMinY(std::floor((yMin - leftLocal[1]) / cellSize));
			int gridIdxMaxX(std::ceil((xMax - leftLocal[0]) / cellSize));
			int gridIdxMaxY(std::ceil((yMax - leftLocal[1]) / cellSize));*/

			int gridIdxMinX(std::floor((xMin - left[0]) / cellSize));
			int gridIdxMinY(std::floor((yMin - left[1]) / cellSize));
			int gridIdxMaxX(std::ceil((xMax - left[0]) / cellSize));
			int gridIdxMaxY(std::ceil((yMax - left[1]) / cellSize));

			for (auto rowIdx(gridIdxMinX); rowIdx <= gridIdxMaxX; ++rowIdx) {
				for (auto colIdx(gridIdxMinY); colIdx <= gridIdxMaxY; ++colIdx) {
					// compute convolution value for this point with regards to edge e
					Eigen::Vector3d p(left[0] + cellSize * rowIdx, left[1] + cellSize * colIdx, layerLevel);
					double l((p1 - p0).norm()), a((p0 - p).dot(p1 - p0)), b((p0 - p).norm());
					// compute line-sphere intersection
					l *= l; a *= 2; b = b * b - convolRadius * convolRadius;
					double delta = a * a - 4 * l * b;
					if (delta >= 1e-6) {
						double intersect1 = (-a - std::sqrt(delta)) / (2 * l);
						double intersect2 = (-a + std::sqrt(delta)) / (2 * l);
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
							else iter->second += f; //summarization for convolution
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
		if (digitImage.empty()) return;

		// Support structure
		/*for (size_t i = 0; i < imageSize[0]; i++)
		{
			
			for (size_t j = 0; j < imageSize[2]; j++)
			{
				auto iter = digitImage.find(std::to_string(i) + " " + std::to_string(j));
				if (iter != digitImage.end())
				{
					binaryNodes[layerNum][i][j] = 1;
					//binaryNodes[layerNum].push_back(1);
				}
				else
				{
					binaryNodes[layerNum][i][j] = 0;
					//binaryNodes[layerNum].push_back(0);
				}
			}
		}*/

		// save the bindary image to file
		size_t width = resX, height = resY;
		QImage img(width, height, QImage::Format_RGB32);
		for (size_t i(0); i < width; ++i) {
			for (size_t j(0); j < height; ++j) {
				auto iter = digitImage.find(std::to_string(i) + " " + std::to_string(j));
				if (iter != digitImage.end())
				{
					img.setPixel(i, j, qRgb(0, 0, 0));
					binaryNodes[layerNum][i][j] = 1;
				}
				else
				{
					img.setPixel(i, j, qRgb(255, 255, 255));
					binaryNodes[layerNum][i][j] = 0;
				}
			}
		}
		std::string path = volumeMeshFilePath + "." + std::to_string(layerNum) + ".slice.jpg";
		if (!img.save(QString::fromStdString(path), "JPEG", 50)) std::cerr << "LatticeModeler > error writing picture at layer: " << layerNum << endl;
		else std::cout << "LatticeModeler > successfully writing picture at layer: " << layerNum << endl;

		// do contouring with optimization
		/*MarchingSquareBin msb(digitImage, xRes, yRes, cellSize, layerLevel, leftLocal.data()); // binImage will be moved, not copy, into MarchingSquareBin
		msb.doContouring();
		auto result = msb.getContours(false);*/
		//.......................
		// do things with "result"
		long time = clock();
		//........................
		MarchingSquareBin msb(digitImage, resX, resY, cellSize, layerLevel, left.data());
		msb.doContouring();
		path = volumeMeshFilePath + "." + std::to_string(layerNum) + ".contours.off";
		//msb.writeContours(path);
		auto result = msb.getContours(true);
		layers.push_back(std::shared_ptr<QMeshPatch>(result));
		//auto result2 = msb.getContours(true);
		//std::string clipath = "Results/Contours/out_volume.contours1.cli";
		//writeCLIFile(clipath, layers);
		std::cout << "LatticeModeler > successfully generating contours at layer: " << layerNum << endl;
		std::cout << "LatticeModeler > MarchingSquare Time (micro second) for " << layerNum << "-th layer:" << double(clock() - time) / CLOCKS_PER_SEC << endl;
	};

	// steaming from file
	std::vector<Vector7d> edgeList; // to be used as a heap
	auto cmp = [](auto& left, auto& right) {return std::max(left[2], left[5]) > std::max(right[2], right[5]); }; // for max heap

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
			double upperBound = std::max(edge[2] + convolRadius, edge[5] + convolRadius);
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
			{
				std::pop_heap(edgeList.begin(), edgeList.end(), cmp); edgeList.pop_back();
			}
			else break;
		}
		if (edgeList.empty()) { layerLevel += layerThick; continue; }

		// do slicing
		time = clock();
		std::cout << "LatticeModeler > Slicing edge num: " << edgeList.size() << endl;
		intersectEdgesPlane(edgeList, layerLevel, layerNum);
		std::cout << "LatticeModeler > Slicing and MarchingSquare Time (micro second) for " << layerNum << "-th layer:" << double(clock() - time) / CLOCKS_PER_SEC << endl;

		++layerNum;
		layerLevel += layerThick;
	}
	ifs.close();
	std::cout << "LatticeModeler > successfully slicing the model" << std::endl;
	// write to .cli files
	string clipath = "Results/CLIFileforPart/out_volume.contours.cli";
	//writeCLIFileBin(path, nullptr, 0, 3); // 3: finish writing the cli file
	//    readCLIFileBin(path);
	writeCLIFileBin(clipath, layers);
	clipath = "Results/CLIFileforPart/out_volume.contours1.cli";
	writeCLIFile(clipath, layers);

	// Support structure
	imageSize[1] = layerNum - 1;
	CUDA_SAFE_CALL(cudaMalloc((void**)&(gridNodes), imageSize[0] * imageSize[1] * imageSize[2] * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)gridNodes, 0.0, imageSize[0] * imageSize[1] * imageSize[2] * sizeof(bool)));

	int idx = 0;
	bool* tmp = new bool[imageSize[0] * imageSize[1] * imageSize[2]];
	for (int j = 0; j < imageSize[2]; j++)
	{
		for (int k = 1; k <= imageSize[1]; k++)
		{

			for (int i = 0; i < imageSize[0]; i++)
			{
				//gridNodes[idx] = binaryNodes[i][j];
				if (binaryNodes[k][i][j] == 1)
					tmp[idx] = true;
				else
					tmp[idx] = false;
				idx++;
			}

		}
	}
	CUDA_SAFE_CALL(cudaMemcpy((void*)gridNodes, (void*)tmp, imageSize[0] * imageSize[1] * imageSize[2] * sizeof(bool), cudaMemcpyHostToDevice));
	//User specified variables
	double anchorR = 2;
	double thres = 0.2;
	double nSampleWidth = 0.005;
	double cylinderRadius = 6.0;
	double patternThickness = 3;
	
	bool *suptNodes;
	bool *tempImage;
	bool *targetImage;
	bool *assistImage;
	bool *temp3D;


	double anchorRadius = anchorR;
	double threshold = thres;

	int suptRadius = (int)floor(anchorRadius*1.414 / nSampleWidth);
	int suptGridResX = (imageSize[0] - 1) / suptRadius + 1;
	int suptGridResZ = (imageSize[2] - 1) / suptRadius + 1;

	int3 imgRes = make_int3(imageSize[0], imageSize[1], imageSize[2]);
	int3 suptimgRes = make_int3(imageSize[0], imageSize[1] - 1, imageSize[2]);
	int nodeNum = imageSize[0] * imageSize[2];

	long ti[10] = { 0 };

	CUDA_SAFE_CALL(cudaMalloc((void**)&(temp3D), imageSize[0] * (imageSize[1] - 1)*imageSize[2] * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)temp3D, 0, imageSize[0] * (imageSize[1] - 1)*imageSize[2] * sizeof(bool)));


	CUDA_SAFE_CALL(cudaMalloc((void**)&(suptNodes), imageSize[0] * (imageSize[1] - 1)*imageSize[2] * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)suptNodes, false, imageSize[0] * (imageSize[1] - 1)*imageSize[2] * sizeof(bool)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&(tempImage), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)tempImage, false, nodeNum * sizeof(bool)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&(targetImage), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)targetImage, false, nodeNum * sizeof(bool)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&(assistImage), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)assistImage, false, nodeNum * sizeof(bool)));


	short2 **disTextures = (short2 **)malloc(2 * sizeof(short2 *));

	int disTexSize = max(imageSize[0], imageSize[2]);
	int factor = ceil((float)disTexSize / BLOCKSIZE);
	disTexSize = BLOCKSIZE * factor;

	int disMemSize = disTexSize * disTexSize * sizeof(short2);

	// Allocate 2 textures

	cudaMalloc((void **)&disTextures[0], disMemSize);
	cudaMalloc((void **)&disTextures[1], disMemSize);

	unsigned int *LinkIndex;
	CUDA_SAFE_CALL(cudaMalloc((void **)&LinkIndex, (nodeNum + 1) * sizeof(unsigned int)));
	CUDA_SAFE_CALL(cudaMemset((void*)LinkIndex, 0, (nodeNum + 1) * sizeof(unsigned int)));

	long t;
	
	for (int i = imageSize[1] - 2; i > -1; i--)
	{
		
		t = clock();
		call_krSLAContouring_Initialization(tempImage, targetImage, gridNodes, nodeNum, imgRes, i);
		CUDA_SAFE_CALL(cudaMemcpy(assistImage, targetImage, nodeNum * sizeof(bool), cudaMemcpyDeviceToDevice));
		ti[0] += clock() - t;

		t = clock();
		LDNIcudaOperation::LDNISLAContouring_GrowthAndSwallow(threshold, targetImage, tempImage, i, imageSize, nSampleWidth, disTextures[0], disTextures[1], disTexSize);
		

		ti[1] += clock() - t;
		//add new support cylinder if necessary
		CUDA_SAFE_CALL(cudaMemset((void*)tempImage, false, nodeNum * sizeof(bool)));
		CUDA_SAFE_CALL(cudaMemset((void*)assistImage, false, nodeNum * sizeof(bool)));
		call_krSLAContouring_Filter1(assistImage, tempImage, targetImage, suptGridResX*suptGridResZ, make_int2(imageSize[0], imageSize[2]), suptRadius, i);


		//first step: support region growth in first class cylinders
		LDNIcudaOperation::LDNISLAContouring_GrowthAndSwallow(anchorRadius, targetImage, assistImage, i, imageSize, nSampleWidth, disTextures[0], disTextures[1], disTexSize);


		//second step: prepare second class cylinders and perform support region growth in second class cylinders
		CUDA_SAFE_CALL(cudaMemset((void*)assistImage, false, nodeNum * sizeof(bool)));
		call_krSLAContouring_OrthoSearchRemainAnchorZ(assistImage, tempImage, targetImage, suptGridResX,
			make_int2(suptRadius, imageSize[2]), make_int2(imageSize[0], imageSize[2]), i);



		call_krSLAContouring_OrthoSearchRemainAnchorX(assistImage, tempImage, targetImage, suptGridResZ,
			make_int2(suptRadius, imageSize[0]), make_int2(imageSize[0], imageSize[2]), i);



		LDNIcudaOperation::LDNISLAContouring_GrowthAndSwallow(anchorRadius, targetImage, assistImage, i, imageSize, nSampleWidth, disTextures[0], disTextures[1], disTexSize);
		long time = clock();

		//third step: prepare third class cylinders and support region growth in all third class cylinders
		CUDA_SAFE_CALL(cudaMemset((void*)assistImage, false, nodeNum * sizeof(bool)));
		LDNIcudaOperation::LDNISLAContouring_ThirdClassCylinder(anchorRadius, targetImage, assistImage, tempImage, make_int2(imageSize[0], imageSize[2]), nSampleWidth, i, disTextures[0], disTextures[1], disTexSize);


		//generate support structure map for this layer and update the cylinder position information arrays

		call_krSLAContouring_Filter5(gridNodes, tempImage, suptNodes, LinkIndex, nodeNum, imgRes, i);




		call_krFDMContouring_CopyNodesrom2Dto3D(tempImage, temp3D, nodeNum, suptimgRes, i);

		//std::cout << "Generate support map time for " << i + 1 << "-th layer:" << double(clock() - time) / CLOCKS_PER_SEC << endl;
	}

	cudaFree(disTextures[0]);
	cudaFree(disTextures[1]);
	free(disTextures);
	cudaFree(assistImage);
	cudaFree(targetImage);
	cudaFree(tempImage);

	unsigned int LinkNum;
	call_func::setdev_ptr(LinkIndex, nodeNum, LinkNum);


	short *linkLayerD;
	short2 *linkID;
	unsigned int *linkLayerC;
	unsigned int *temp2D;

	CUDA_SAFE_CALL(cudaMalloc((void **)&temp2D, nodeNum * sizeof(unsigned int)));
	CUDA_SAFE_CALL(cudaMemset((void*)temp2D, 0, nodeNum * sizeof(unsigned int)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&linkLayerC, LinkNum * sizeof(unsigned int)));
	//CUDA_SAFE_CALL(cudaMemset((void*)linkLayerC, imgRes.y+1, LinkNum*sizeof(unsigned int) ) );
	CUDA_SAFE_CALL(cudaMalloc((void **)&linkLayerD, LinkNum * sizeof(short)));
	CUDA_SAFE_CALL(cudaMemset((void*)linkLayerD, 0, LinkNum * sizeof(short)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&linkID, LinkNum * sizeof(short2)));
	CUDA_SAFE_CALL(cudaMemset((void*)linkID, 0, LinkNum * sizeof(short2)));


	call_func::setdev_ptr2(linkLayerC, LinkNum, imgRes);

	time = clock();
	call_krSLAContouring_FindAllLinks(LinkIndex, linkLayerC, linkLayerD, linkID, temp3D, temp2D, suptimgRes.x*suptimgRes.y*suptimgRes.z, suptimgRes);

	cudaFree(temp3D);
	cudaFree(temp2D);



	call_krSLAContouring_RelateAllLinksBetweenLayers(LinkIndex, linkLayerC, gridNodes, suptNodes, suptimgRes.x*suptimgRes.y*suptimgRes.z, imgRes);

	cudaFree(LinkIndex);
	std::cout<< "Relate all links time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;



	CUDA_SAFE_CALL(cudaMalloc((void**)&(tempImage), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)tempImage, false, nodeNum * sizeof(bool)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&(targetImage), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)targetImage, false, nodeNum * sizeof(bool)));


	double realThreshold = (cylinderRadius*nSampleWidth - nSampleWidth) / nSampleWidth;
	int gridRadius = (int)floor(realThreshold);
	time = clock();
	for (int i = imageSize[1] - 2; i > -1; i--)
	{
		
		call_krFDMContouring_CopyNodesrom3Dto2D(targetImage, suptNodes, nodeNum, suptimgRes, i);
		CUDA_SAFE_CALL(cudaMemcpy(tempImage, targetImage, nodeNum * sizeof(bool), cudaMemcpyDeviceToDevice));
		call_krFDMContouring_Dilation(targetImage, tempImage, nodeNum, imgRes, realThreshold, gridRadius, i);


		call_krFDMContouring_CopyNodesrom2Dto3D(tempImage, suptNodes, nodeNum, suptimgRes, i);
		
	}
	std::cout << "Dilation time:" << double(clock() - time) / CLOCKS_PER_SEC << endl;
	cudaFree(targetImage);
	cudaFree(tempImage);
	time = clock();
	for (int i = imageSize[1] - 3; i > -1; i--)
		call_krFDMContouring_VerticalSpptPxlProp(gridNodes, suptNodes, nodeNum, imgRes, i);
	std::cout<< "VerticalSpptPxlProp time:"<< double(clock() - time) / CLOCKS_PER_SEC << endl;

	int linkThreshold = (int)(0.4 / nSampleWidth);
	int lengthofLayer = linkThreshold / 8;
	int furtherStepLength = lengthofLayer / 2;



	time = clock();
	LDNIcudaOperation::LDNISLAContouring_GenerateConnectionforCylinders(linkLayerC, linkLayerD, linkID, gridNodes, suptNodes, imageSize, linkThreshold,
		lengthofLayer, furtherStepLength, LinkNum, nSampleWidth);
	std::cout<< "Generate Connection for Cylinders time:" << double(clock() - time) / CLOCKS_PER_SEC << endl;
	cudaFree(linkLayerC);
	cudaFree(linkLayerD);
	cudaFree(linkID);
	//cudaFree(gridNodes);


	realThreshold = (patternThickness*nSampleWidth - nSampleWidth) / nSampleWidth;
	gridRadius = (int)floor(realThreshold);

	CUDA_SAFE_CALL(cudaMalloc((void**)&(tempImage), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)tempImage, false, nodeNum * sizeof(bool)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&(targetImage), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)targetImage, false, nodeNum * sizeof(bool)));
	time = clock();
	for (int i = imageSize[1] - 2; i > -1; i--)
	{
		
		call_krFDMContouring_CopyNodesrom3Dto2D(targetImage, suptNodes, nodeNum, suptimgRes, i);
		CUDA_SAFE_CALL(cudaMemcpy(tempImage, targetImage, nodeNum * sizeof(bool), cudaMemcpyDeviceToDevice));
		call_krFDMContouring_Dilation(targetImage, tempImage, nodeNum, imgRes, realThreshold, gridRadius, i);
		call_krFDMContouring_CopyNodesrom2Dto3D(tempImage, suptNodes, nodeNum, suptimgRes, i);
		
	}
	std::cout << "Dilation time:" << double(clock() - time) / CLOCKS_PER_SEC << endl;
	cudaFree(targetImage);
	cudaFree(tempImage);


	bool* gridtmp = new bool[imageSize[0] * imageSize[1] * imageSize[2]];
	CUDA_SAFE_CALL(cudaMemcpy((void*)gridtmp, (void*)gridNodes, imageSize[0] * imageSize[1] * imageSize[2] * sizeof(bool), cudaMemcpyDeviceToHost));
	int cnt = 0;
	size_t width = imageSize[0], height = imageSize[2], slicenum = imageSize[1];
	QImage *img = new QImage[slicenum];
	for (int i = 0; i < slicenum; i++)
	{
		QImage tmpimg(width, height, QImage::Format_RGB32);
		img[i] = tmpimg;
	}
	for (int j = 0; j < imageSize[2]; j++)
	{

		for (int k = 0; k < imageSize[1]; k++)
		{

			for (int i = 0; i < imageSize[0]; i++)
			{
				if (gridtmp[cnt++])img[k].setPixel(i, j, qRgb(0, 0, 0));
				else img[k].setPixel(i, j, qRgb(255, 255, 255));

			}
		}

	}
	for (int i = 1; i <= imageSize[1]; i++)
	{
		std::string path = "Results/Support/out_volume." + std::to_string(i) + ".origin.jpg";
		if (!img[i - 1].save(QString::fromStdString(path), "JPEG", 50)) std::cerr << "Origin > error writing picture at layer: " << i << endl;
		else std::cout << "Origin > successfully writing picture at layer: " << i << endl;
	}


	std::unordered_map<std::string, double> *digitImage = new std::unordered_map<std::string, double>[slicenum - 1];//generate contours for support structure
	std::vector<std::shared_ptr<QMeshPatch>> suptlayers;
	layerLevel = zMin + layerThick;
	bool* supttmp = new bool[imageSize[0] * (imageSize[1] - 1)*imageSize[2]];
	CUDA_SAFE_CALL(cudaMemcpy((void*)supttmp, (void*)suptNodes, imageSize[0] * (imageSize[1] - 1) * imageSize[2] * sizeof(bool), cudaMemcpyDeviceToHost));
	cnt = 0;
	QImage *suptimg = new QImage[slicenum - 1];
	for (int i = 0; i < slicenum - 1; i++)
	{
		QImage tmpimg(width, height, QImage::Format_RGB32);
		suptimg[i] = tmpimg;
	}
	for (int j = 0; j < imageSize[2]; j++)
	{
		for (int k = 0; k < imageSize[1] - 1; k++)
		{

			for (int i = 0; i < imageSize[0]; i++)
			{
				if (supttmp[cnt++])
				{
					suptimg[k].setPixel(i, j, qRgb(255, 0, 0));
					auto key = std::to_string(i) + " " + std::to_string(j);
					digitImage[k].insert(std::make_pair(key, 1));
				}
				else suptimg[k].setPixel(i, j, qRgb(255, 255, 255));

			}
		}

	}
	for (int i = 1; i <= imageSize[1] - 1; i++)
	{
		std::string path = "Results/Support/out_volume." + std::to_string(i) + ".support.jpg";
		if (!suptimg[i - 1].save(QString::fromStdString(path), "JPEG", 50)) std::cerr << "Support > error writing picture at layer: " << i << endl;
		else std::cout << "Support > successfully writing picture at layer: " << i << endl;

		long time = clock();
		MarchingSquareBin msb(digitImage[i - 1], resX, resY, cellSize, layerLevel, left.data());
		msb.doContouring();
		path = "Results/Support/out_volume." + std::to_string(i) + ".contours.off";
		//msb.writeContours(path);
		auto result = msb.getContours(true);
		suptlayers.push_back(std::shared_ptr<QMeshPatch>(result));
		layerLevel += layerThick;
		std::cout<< "Support >  MarchingSquare Time (micro second) for " << i << "-th layer:" << double(clock() - time) / CLOCKS_PER_SEC << endl;
	}
	// write to .cli files
	string suptclipath = "Results/Support/CLIFileforSupt/out_volume.contours.cli";
	writeCLIFileBin(suptclipath, suptlayers);
	suptclipath = "Results/Support/CLIFileforSupt/out_volume.contours1.cli";
	writeCLIFile(suptclipath, suptlayers);
	return true;
}
#endif
/*bool LatticeModeler::sliceModelImplicitStream(std::string path)
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
    left += Eigen::Vector3d(-3*convolRadius, -3 * convolRadius, -3 * convolRadius);
    right += Eigen::Vector3d(3 * convolRadius, 3 * convolRadius, 3 * convolRadius);

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
    // Support structure
    bool *gridNodes;
    vector<vector<int>> binaryNodes[40];
    int imageSize[3] = {0,0,0};
    

	double cellSize = min((right[0]-left[0])/min(resX,resY),(right[1]-left[1])/min(resX,resY));
	resX = (right[0] - left[0]) / cellSize;
	resY = (right[1] - left[1]) / cellSize;
	imageSize[0] = resX; //cout << (right[0] - left[0]) << endl;
	imageSize[2] = resY; //cout << (right[1] - left[1]) << endl;
    // a functor for slicing at each iteration
    auto intersectEdgesPlane = [&](std::vector<Vector7d>& edges, double layerLevel, size_t layerNum){
        if (edges.empty()) return;
        // edges: 0-2: first point coordinates; 3-5: second point coordinates; 6: radius
        // setup the digital image size, resolution, and convolutional radius
        /*Eigen::Vector3d leftLocal, rightLocal;
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
            /*int gridIdxMinX(std::floor((xMin - leftLocal[0]) / cellSize));
            int gridIdxMinY(std::floor((yMin - leftLocal[1]) / cellSize));
            int gridIdxMaxX(std::ceil((xMax - leftLocal[0]) / cellSize));
            int gridIdxMaxY(std::ceil((yMax - leftLocal[1]) / cellSize));

			int gridIdxMinX(std::floor((xMin - left[0]) / cellSize));
			int gridIdxMinY(std::floor((yMin - left[1]) / cellSize));
			int gridIdxMaxX(std::ceil((xMax - left[0]) / cellSize));
			int gridIdxMaxY(std::ceil((yMax - left[1]) / cellSize));

            for (auto rowIdx(gridIdxMinX); rowIdx <= gridIdxMaxX; ++rowIdx) {
                for (auto colIdx(gridIdxMinY); colIdx <= gridIdxMaxY; ++colIdx) {
                    // compute convolution value for this point with regards to edge e
                    Eigen::Vector3d p(left[0] + cellSize * rowIdx, left[1] + cellSize * colIdx, layerLevel);
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
                            else iter->second += f; //summarization for convolution
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

        // Support structure
        for(size_t i = 0; i < imageSize[0]; i++)
        {
			vector<int> v;
			binaryNodes[layerNum].push_back(v);
            for(size_t j = 0; j < imageSize[2]; j++)
            {
                auto iter = digitImage.find(std::to_string(i) + " " + std::to_string(j));
                if (iter != digitImage.end())
                {
					binaryNodes[layerNum][i].push_back(1);
                    //binaryNodes[layerNum].push_back(1);
                }
                else
                {
					binaryNodes[layerNum][i].push_back(0);
                    //binaryNodes[layerNum].push_back(0);
                }
            }
        }

        // save the bindary image to file
        size_t width = resX, height = resY;
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
        /*MarchingSquareBin msb(digitImage, xRes, yRes, cellSize, layerLevel, leftLocal.data()); // binImage will be moved, not copy, into MarchingSquareBin
        msb.doContouring();
        auto result = msb.getContours(false);
        //.......................
        // do things with "result"

        //........................
		MarchingSquareBin msb(digitImage, resX, resY, cellSize, layerLevel, left.data());
		msb.doContouring();
		path = volumeMeshFilePath + "." + std::to_string(layerNum) + ".contours.off";
		//msb.writeContours(path);
		auto result = msb.getContours(true);
		layers.push_back(std::shared_ptr<QMeshPatch>(result));
		//auto result2 = msb.getContours(true);
		//std::string clipath = "Results/Contours/out_volume.contours1.cli";
		//writeCLIFile(clipath, layers);
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
	// write to .cli files
	string clipath = "Results/CLIFileforPart/out_volume.contours.cli";
	//writeCLIFileBin(path, nullptr, 0, 3); // 3: finish writing the cli file
	//    readCLIFileBin(path);
	writeCLIFileBin(clipath, layers);
	clipath = "Results/CLIFileforPart/out_volume.contours1.cli";
	writeCLIFile(clipath, layers);

    // Support structure
    imageSize[1] = layerNum-1;
    CUDA_SAFE_CALL( cudaMalloc( (void**)&(gridNodes), imageSize[0]*imageSize[1]*imageSize[2]*sizeof(bool) ) );
    CUDA_SAFE_CALL( cudaMemset( (void*)gridNodes, 0.0, imageSize[0]*imageSize[1]*imageSize[2]*sizeof(bool) ) );

    int idx = 0;
	bool* tmp = new bool[imageSize[0] * imageSize[1] * imageSize[2]];
    for(int j=0;j<imageSize[2];j++)
    {		
		for (int k = 1; k <= imageSize[1]; k++)
        {
			
			for (int i = 0; i < imageSize[0]; i++)
			{
				//gridNodes[idx] = binaryNodes[i][j];
				if (binaryNodes[k][i][j] == 1)
					tmp[idx] = true;
				else
					tmp[idx] = false;
				idx++;
			}
            
        }
    }
	CUDA_SAFE_CALL(cudaMemcpy((void*)gridNodes, (void*)tmp, imageSize[0] * imageSize[1] * imageSize[2]*sizeof(bool), cudaMemcpyHostToDevice));
    bool *suptNodes;
    bool *solidNodes;
    bool *suptTemp;
    bool *grossImage;
    bool *test1;
    bool *test2;
    double t = 0.025;
    double nSampleWidth = 0.005;
    int3 imgRes = make_int3(imageSize[0], imageSize[1], imageSize[2]);
    int3 suptimgRes = make_int3(imageSize[0], imageSize[1]-1, imageSize[2]);
    int nodeNum = imageSize[0]*imageSize[2];
    CUDA_SAFE_CALL( cudaMalloc( (void**)&(suptNodes), imageSize[0]*(imageSize[1]-1)*imageSize[2]*sizeof(bool) ) );
    CUDA_SAFE_CALL( cudaMemset( (void*)suptNodes, false, imageSize[0]*(imageSize[1]-1)*imageSize[2]*sizeof(bool) ) );

    CUDA_SAFE_CALL( cudaMalloc( (void**)&(solidNodes), nodeNum*sizeof(bool) ) );
    CUDA_SAFE_CALL( cudaMemset( (void*)solidNodes, false, nodeNum*sizeof(bool) ) );

    CUDA_SAFE_CALL( cudaMalloc( (void**)&(suptTemp), nodeNum*sizeof(bool) ) );
    CUDA_SAFE_CALL( cudaMemset( (void*)suptTemp, false, nodeNum*sizeof(bool) ) );

    CUDA_SAFE_CALL( cudaMalloc( (void**)&(grossImage), nodeNum*sizeof(bool) ) );
    CUDA_SAFE_CALL( cudaMemset( (void*)grossImage, false, nodeNum*sizeof(bool) ) );

    CUDA_SAFE_CALL( cudaMalloc( (void**)&(test1), nodeNum*sizeof(bool) ) );
    CUDA_SAFE_CALL( cudaMemset( (void*)test1, false, nodeNum*sizeof(bool) ) );

    CUDA_SAFE_CALL( cudaMalloc( (void**)&(test2), nodeNum*sizeof(bool) ) );
    CUDA_SAFE_CALL( cudaMemset( (void*)test2, false, nodeNum*sizeof(bool) ) );

    short2 **disTextures = (short2 **) malloc(2 * sizeof(short2 *));
    int disTexSize = max(imageSize[0],imageSize[2]);
    int factor = ceil((float)disTexSize / 64);
    disTexSize = 64*factor;
    int disMemSize = disTexSize * disTexSize * sizeof(short2);
    // Allocate 2 textures
    cudaMalloc((void **) &disTextures[0], disMemSize);
    cudaMalloc((void **) &disTextures[1], disMemSize);

    long ti[10] = {0,0,0,0,0,0,0,0,0,0};
    long te;
    for(int i=imageSize[1]-2; i>-1; i--)
    {
        te = clock();
        CUDA_SAFE_CALL( cudaMemset( (void*)solidNodes, false, nodeNum*sizeof(bool) ) );
        call_krFDMContouring_SubtractSolidRegion(gridNodes, suptNodes, solidNodes, nodeNum, imgRes, i);
        cudaThreadSynchronize();
        ti[0] += clock()-te;

        te = clock();
        call_krFDMContouring_CopyNodesrom3Dto2D(suptTemp, suptNodes, nodeNum, suptimgRes, i);
        cudaThreadSynchronize();
        ti[1] += clock()-te;

        te = clock();
        LDNIcudaOperation::LDNIFDMContouring_GrowthAndSwallow(gridNodes, i, imageSize, t, nSampleWidth, suptTemp, solidNodes, disTextures[0], disTextures[1], disTexSize, test1, test2);
        cudaThreadSynchronize();
        ti[2] += clock()-te;

        te = clock();
        call_krFDMContouring_Filter4(solidNodes, suptNodes, suptTemp, nodeNum, suptimgRes, i);
        cudaThreadSynchronize();
        ti[3] += clock()-te;

        te = clock();
        call_krFDMContouring_integrateImageintoGrossImage(grossImage, suptNodes, nodeNum, suptimgRes, i);
        cudaThreadSynchronize();
        ti[4] += clock()-te;

        te = clock();
        call_krFDMContouring_CopyNodesrom2Dto3D(grossImage, suptNodes, nodeNum, suptimgRes, i);
        cudaThreadSynchronize();
        ti[5] += clock() -te;

		te = clock();
		LDNIcudaOperation::LDNIFDMContouring_Closing(imageSize, 2.0*t, nSampleWidth, suptNodes, i, test1, test2);
		cudaThreadSynchronize();
		ti[6] += clock() - te;

		te = clock();
		call_krFDMContouring_Filter5(suptNodes, solidNodes, gridNodes, nodeNum, i, imgRes);
		cudaThreadSynchronize();
		ti[7] += clock() - te;
    }

	bool* gridtmp = new bool[imageSize[0] * imageSize[1] *imageSize[2]];
	CUDA_SAFE_CALL(cudaMemcpy((void*)gridtmp, (void*)gridNodes, imageSize[0] * imageSize[1] * imageSize[2] * sizeof(bool), cudaMemcpyDeviceToHost));
	int cnt = 0;
	size_t width = imageSize[0], height = imageSize[2], slicenum = imageSize[1];
	QImage *img = new QImage[slicenum];
	for (int i = 0; i < slicenum; i++)
	{
		QImage tmpimg(width, height, QImage::Format_RGB32);
		img[i] = tmpimg;
	}
	for (int j = 0; j < imageSize[2]; j++)
	{
		
		for (int k = 0; k < imageSize[1]; k++)
		{
			
			for (int i = 0; i < imageSize[0]; i++)
			{
				if (gridtmp[cnt++])img[k].setPixel(i, j, qRgb(0, 0, 0));
				else img[k].setPixel(i, j, qRgb(255, 255, 255));

			}
		}
		
	}
	for (int i = 1; i <= imageSize[1]; i++)
	{
		std::string path = "Results/Support/out_volume." + std::to_string(i) + ".origin.jpg";
		if (!img[i-1].save(QString::fromStdString(path), "JPEG", 50)) std::cerr << "Origin > error writing picture at layer: " << i << endl;
		else std::cout << "Origin > successfully writing picture at layer: " << i << endl;
	}


	std::unordered_map<std::string, double> *digitImage = new std::unordered_map<std::string, double>[slicenum-1];//generate contours for support structure
	std::vector<std::shared_ptr<QMeshPatch>> suptlayers;
	layerLevel = zMin + layerThick;
	bool* supttmp = new bool[imageSize[0] * (imageSize[1] - 1)*imageSize[2]];
	CUDA_SAFE_CALL(cudaMemcpy((void*)supttmp, (void*)suptNodes, imageSize[0] * (imageSize[1]-1) * imageSize[2] * sizeof(bool), cudaMemcpyDeviceToHost));
    cnt = 0;
	QImage *suptimg = new QImage[slicenum-1];
	for (int i = 0; i < slicenum-1; i++)
	{
		QImage tmpimg(width, height, QImage::Format_RGB32);
		suptimg[i] = tmpimg;
	}
    for(int j=0; j<imageSize[2];j++)
    {
		for (int k = 0; k < imageSize[1] - 1; k++)
        {
            
			for (int i = 0; i < imageSize[0]; i++)
            {
				if (supttmp[cnt++])
				{
					suptimg[k].setPixel(i, j, qRgb(255, 0, 0));
					auto key = std::to_string(i) + " " + std::to_string(j);
					digitImage[k].insert(std::make_pair(key, 1));
				}
                else suptimg[k].setPixel(i, j, qRgb(255, 255, 255));

            }
        }
        
    }
	for (int i = 1; i <= imageSize[1] - 1; i++)
	{
		std::string path = "Results/Support/out_volume." + std::to_string(i) + ".support.jpg";
		if (!suptimg[i-1].save(QString::fromStdString(path), "JPEG", 50)) std::cerr << "Support > error writing picture at layer: " << i << endl;
		else std::cout << "Support > successfully writing picture at layer: " << i << endl;

		MarchingSquareBin msb(digitImage[i-1], resX, resY, cellSize, layerLevel, left.data());
		msb.doContouring();
		path = "Results/Support/out_volume." + std::to_string(i) + ".contours.off";
		//msb.writeContours(path);
		auto result = msb.getContours(true);
		suptlayers.push_back(std::shared_ptr<QMeshPatch>(result));
		layerLevel += layerThick;
	}
	// write to .cli files
	string suptclipath = "Results/Support/CLIFileforSupt/out_volume.contours.cli";
	writeCLIFileBin(suptclipath, suptlayers);
	suptclipath = "Results/Support/CLIFileforSupt/out_volume.contours1.cli";
	writeCLIFile(suptclipath, suptlayers);
    return true;
}*/

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

bool LatticeModeler::writeCLIFile(string path, std::vector<std::shared_ptr<QMeshPatch> > &layers)
{
	// bounding box
	Eigen::Vector3d left(1e8 * Eigen::Vector3d::Ones()), right(-1e8 * Eigen::Vector3d::Ones()), p;
	for (auto& l : layers) {
		for (auto iter = l->FaceBegin(); iter != l->FaceEnd(); iter = iter->Next()) {
			auto face = l->GetFaceAt(iter);
			//int tmp = face->GetEdgeNum(); cout << endl << tmp;
			for (size_t i = 0; i < face->GetEdgeNum(); i++) {
				face->GetNodeRecordPtr(i)->GetCoord3D(p[0], p[1], p[2]);
				left[0] = std::min(left[0], p[0]);
				left[1] = std::min(left[1], p[1]);
				left[2] = std::min(left[2], p[2]);
				right[0] = std::max(right[0], p[0]);
				right[1] = std::max(right[1], p[1]);
				right[2] = std::max(right[2], p[2]);
			}
		}
	}

	// shift the bottom to the xoy plane
	vector<double> layerHeights;
	Eigen::Vector3d offset(left[0] - 0.01, left[1] - 0.01, left[2] - 0.01);
	for (auto& l : layers) {
		for (auto iter = l->FaceBegin(); iter != l->FaceEnd(); iter = iter->Next()) {
			auto face = l->GetFaceAt(iter);
			for (size_t i = 0; i < face->GetEdgeNum(); i++) {
				face->GetNodeRecordPtr(i)->GetCoord3D(p[0], p[1], p[2]);
				p[0] -= offset[0]; p[1] -= offset[1]; p[2] -= offset[2];
				face->GetNodeRecordPtr(i)->SetCoord3D(p[0], p[1], p[2]);
			}
		}
		layerHeights.push_back(p[2]);
	}
	left -= offset; right -= offset;

	// write to the file
	// a functor for determining ccw
	auto areaCal = [](QMeshFace* ptrFace) {
		auto num = ptrFace->GetEdgeNum();
		double area(0);
		double p1[3], p2[3];
		for (size_t i(0); i < num; ++i) {
			ptrFace->GetNodeRecordPtr(i)->GetCoord3D(p1[0], p1[1], p1[2]);
			ptrFace->GetNodeRecordPtr((i + 1) % num)->GetCoord3D(p2[0], p2[1], p2[2]);
			area += p1[0] * p2[1] - p1[1] * p2[0];
		}
		return area / 2;
	};

	std::ofstream outfile;
	outfile.open(path, ios::out | ios::binary);
	if (!outfile.is_open()) {
		std::cerr << "LatticeModeler > error openning CLI file" << endl;
		return false;
	}
	// write header
	outfile << "$$HEADERSTART" << std::endl;
	outfile << "$$BINARY" << std::endl << "$$UNITS\/00000000.010000" << std::endl << "$$VERSION\/100" << std::endl;
	outfile << "$$LABEL\/1,part1" << std::endl << "$$DATE\/070920" << std::endl;
	outfile << "$$DIMENSION\/" << std::fixed << std::setprecision(6) << std::setw(8) << std::setfill('0') <<
		left[0] << "," << left[1] << "," << left[2] << "," << right[0] << "," << right[1] << "," << right[2] << endl;
	outfile << "$$LAYERS\/" << layers.size() << endl;
	outfile << "$$HEADEREND";

	// write geometry
	for (int i(0); i < layers.size(); ++i) {
		auto ptrContours = layers[i];
		double layerThick = layerHeights[i];
		//        if (i == 0) layerThick = layerHeights[i];
		//        else layerThick = layerHeights[i] - layerHeights[i - 1];
		outfile << "$$LAYER\/" << layerThick << std::endl;
		size_t dir, num;
		double area;
		double p[3];
		for (auto iter = ptrContours->FaceBegin(); iter != ptrContours->FaceEnd(); iter = iter->Next()) {
			auto face = ptrContours->GetFaceAt(iter);
			num = face->GetEdgeNum();
			if (num < 3) continue;
			area = areaCal(face);
			if (abs(area) < 1e-6) continue;
			dir = area > 0 ? 1 : 0;
			outfile << "$$POLYLINE\/1," << dir << "," << num + 1 << ",";
			for (size_t i = 0; i < num; i++) {
				face->GetNodeRecordPtr(i)->GetCoord3D(p[0], p[1], p[2]);
				outfile /*<< std::fixed*/ << p[0] << "," << p[1] << ",";
			}
			// the last/first point
			face->GetNodeRecordPtr(0)->GetCoord3D(p[0], p[1], p[2]);
			outfile /*<< std::fixed*/ << p[0] << "," << p[1] << std::endl;
		}
	}

	outfile.close();
	return true;
}

bool LatticeModeler::writeCLIFileBin(string path, std::vector<std::shared_ptr<QMeshPatch> > &layers)
{
	// bounding box
	Eigen::Vector3d left(1e8 * Eigen::Vector3d::Ones()), right(-1e8 * Eigen::Vector3d::Ones()), p;
	for (auto& l : layers) {
		for (auto iter = l->FaceBegin(); iter != l->FaceEnd(); iter = iter->Next()) {
			auto face = l->GetFaceAt(iter);
			for (size_t i = 0; i < face->GetEdgeNum(); i++) {
				face->GetNodeRecordPtr(i)->GetCoord3D(p[0], p[1], p[2]);
				left[0] = std::min(left[0], p[0]);
				left[1] = std::min(left[1], p[1]);
				left[2] = std::min(left[2], p[2]);
				right[0] = std::max(right[0], p[0]);
				right[1] = std::max(right[1], p[1]);
				right[2] = std::max(right[2], p[2]);
			}
		}
	}

	// shift the bottom to the xoy plane
	vector<double> layerHeights;
	Eigen::Vector3d offset(left[0] - 0.01, left[1] - 0.01, left[2] - 0.01);
	for (auto& l : layers) {
		for (auto iter = l->FaceBegin(); iter != l->FaceEnd(); iter = iter->Next()) {
			auto face = l->GetFaceAt(iter);
			for (size_t i = 0; i < face->GetEdgeNum(); i++) {
				face->GetNodeRecordPtr(i)->GetCoord3D(p[0], p[1], p[2]);
				p[0] -= offset[0]; p[1] -= offset[1]; p[2] -= offset[2];
				face->GetNodeRecordPtr(i)->SetCoord3D(p[0], p[1], p[2]);
			}
		}
		layerHeights.push_back(p[2]);
	}
	left -= offset; right -= offset;

	// write to the file
	// a functor for determining ccw
	auto areaCal = [](QMeshFace* ptrFace) {
		auto num = ptrFace->GetEdgeNum();
		double area(0);
		double p1[3], p2[3];
		for (size_t i(0); i < num; ++i) {
			ptrFace->GetNodeRecordPtr(i)->GetCoord3D(p1[0], p1[1], p1[2]);
			ptrFace->GetNodeRecordPtr((i + 1) % num)->GetCoord3D(p2[0], p2[1], p2[2]);
			area += p1[0] * p2[1] - p1[1] * p2[0];
		}
		return area / 2;
	};

	auto writeUInt = [](uint16_t i, ofstream& out) {
		union u { uint16_t i; unsigned char c[2]; } ic;
		ic.i = i;
		out.write((char*)ic.c, 2);
	};

	std::ofstream outfile;
	outfile.open(path, ios::out | ios::binary);
	if (!outfile.is_open()) {
		std::cerr << "LatticeModeler > error openning CLI file" << endl;
		return false;
	}
	// write header
	outfile << "$$HEADERSTART" << std::endl;
	outfile << "$$BINARY" << std::endl << "$$UNITS\/00000000.010000" << std::endl << "$$VERSION\/100" << std::endl;
	outfile << "$$LABEL\/1,part1" << std::endl << "$$DATE\/070920" << std::endl;
	outfile << "$$DIMENSION\/" << std::fixed << std::setprecision(6) << std::setfill('0') << std::setw(15) <<
		left[0] << "," << std::setfill('0') << std::setw(15) << left[1] << "," << std::setfill('0') << std::setw(15) << left[2] << "," << std::setfill('0') << std::setw(15) <<
		right[0] << "," << std::setfill('0') << std::setw(15) << right[1] << "," << std::setfill('0') << std::setw(15) << right[2] << endl;
	outfile << "$$LAYERS\/" << std::setfill('0') << std::setw(6) << layers.size() << endl;
	outfile << "$$HEADEREND";

	// write geometry
	for (int i(0); i < layers.size(); ++i) {
		uint16_t command = 128; // layer start
		writeUInt(command, outfile);
		//        double layerThick;
		//        if (i == 0) layerThick = layerHeights[i];
		//        else layerThick = layerHeights[i] - layerHeights[i - 1];
		double layerThick = layerHeights[i];
		writeUInt(layerThick * 100, outfile);
		uint16_t dir, num, temp(1);
		double area, p[3];

		auto ptrContours = layers[i];
		for (auto iter = ptrContours->FaceBegin(); iter != ptrContours->FaceEnd(); iter = iter->Next()) {
			auto face = ptrContours->GetFaceAt(iter);
			num = uint16_t(face->GetEdgeNum());
			if (num < 3) continue;
			area = areaCal(face);
			if (abs(area) < 1e-6) continue; // super small regions
			dir = area > 0 ? 1 : 0;
			command = 129;
			writeUInt(command, outfile);
			writeUInt(temp, outfile);
			writeUInt(dir, outfile);
			writeUInt(num + 1, outfile); // +1 due to p0 = pn

			uint16_t prevP[2] = { 10000, 10000 };
			for (size_t i = 0; i < num; i++) {
				face->GetNodeRecordPtr(i)->GetCoord3D(p[0], p[1], p[2]);
				// just in case
				if (p[0] < 0 || p[1] < 0) {
					cout << p[0] << "  " << p[1] << endl;
					std::cerr << "LatticeModeler > error writing CLI file" << endl;
					return false;
				}
				if (prevP[0] == uint16_t(p[0] * 100) && prevP[1] == uint16_t(p[1] * 100)) continue;
				else { prevP[0] = uint16_t(p[0] * 100); prevP[1] = uint16_t(p[1] * 100); }

				//cout << uint16_t(p[0] * 100) << "  " << uint16_t(p[1] * 100) << endl;
				writeUInt(p[0] * 100, outfile);
				writeUInt(p[1] * 100, outfile);
			}
			// the last/first point
			face->GetNodeRecordPtr(0)->GetCoord3D(p[0], p[1], p[2]);
			//cout << uint16_t(p[0] * 100) << "  " << uint16_t(p[1] * 100) << endl;
			writeUInt(p[0] * 100, outfile);
			writeUInt(p[1] * 100, outfile);
		}
	}

	outfile.close();
	return true;
}













