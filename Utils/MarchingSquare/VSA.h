// VSA.h: interface for the Variational Shape Approximation class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _HP_VSA
#define _HP_VSA

#include "../GLKLib/GLKHeap.h"

class QMeshNode;
class QMeshFace;
class QMeshEdge;
class QMeshPatch;

class GLKObList;

class VSAHeapNode : public GLKHeapNode
{
public:
	VSAHeapNode();
	virtual ~VSAHeapNode();

public:
	int whichproxyagainst;
};

class VSANode
{
public:
	VSANode();
	~VSANode();

public:
	int ProxyIndicator;
	double centerpnt[3];		//note: for generality, we just use 3d coord to represent 2d point with the coordinate in y direction being 0 
	double area;
	void *meshobj;
};

class VSA
{
public:
	VSA();
	~VSA();

	void InitializeVSA(QMeshPatch *parapatch, int pararegionnum, short paradim);

	void PerformVSA2D(int iternum, double desireddistterror, bool needcontourbdcheck);		//the parameter needcontourbdcheck is for whether you are dealing with closed or open contour
	float EvaluateDistortionError2D(double *proxy, double v1, double v2, double u1, double u2, double length);
	void SimplifyMeshBasedOnVSARegions2D();		//note: this function will not release memory for the old patch
	void SimplifyMeshBasedOnVSARegions2D_LineIntsct();		//note: this function will not release memory for the old patch

	//these functions are just for binary image VSA application in order to keep every black sampling point in the contour and every white sampling point out of the contour  
	void BinaryImageInOutCorrection(GLKArray *stickindexx, GLKArray *stickindexz, double *biorigin, double bigridwidth, int baseedgenumber);
	bool _IsTwoSegmentIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);
	int _RecursiveLocalVSA(QMeshPatch *parapatch, QMeshNode *startnode, QMeshNode *endnode, int &currregionnum,
		GLKArray *stickindexx, GLKArray *stickindexz, double *biorigin, double bigridwidth);

	//these function are for verification for distortion error, they are also in the partten of recursion like the above block of functions
	double DistortionErrorCorrection(double desireddistterror);
	int _RecursiveLocalVSAForDistterror(QMeshPatch *parapatch, QMeshNode *startnode, QMeshNode *endnode, int &currregionnum, double desireddistterror, double &largerdistterror);

private:
	double VSACoreIterations(int maxiter, GLKHeap *heap, double ***covmatrixarray, double *totalareaforregions, double *regioncenterx, double *regioncenterz, double *mindisttforregions,
							double *regiondistortionerrors, VSANode **tempvsanodelist, int &maxdistterrorregion, int &maxdisttvsanodeind, bool needcontourbdcheck);

public:
	int m_regionnum;
	int m_facenum;
	short m_dim;
	QMeshPatch *m_meshpatch;
	VSANode **m_vsanodes;
	double **m_proxies;

	int fakeregionnum;
};



#endif
