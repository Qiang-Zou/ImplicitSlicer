// QMeshNode.h: interface for the QMeshNode class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _CW_QMESHNODE
#define _CW_QMESHNODE

#include "../GLKLib/GLKObList.h"
#include "QMeshPatch.h"
#include "QMeshEdge.h"
#include "QMeshFace.h"

class QMeshPatch;
class QMeshFace;
class QMeshEdge;

class QMeshNode : public GLKObject  
{
public:
	QMeshNode();
	virtual ~QMeshNode();

public:
	int GetIndexNo();		//from 1 to n
	void SetIndexNo( const int _index = 1 );

	bool GetAttribFlag( const int whichBit );
	void SetAttribFlag( const int whichBit, const bool toBe = true );

	void GetCoord2D( double &x, double &y );
	void SetCoord2D( double x, double y );

	void GetCoord3D( double &x, double &y, double &z );
	void SetCoord3D( double x, double y, double z );

	void GetCoord3D_last( double &x, double &y, double &z );
	void SetCoord3D_last( double x, double y, double z );

	void SetMeanCurvatureNormalVector(double kHx, double kHy, double kHz);
	void GetMeanCurvatureNormalVector(double &kHx, double &kHy, double &kHz);
	
	void SetGaussianCurvature(double kG);
	double GetGaussianCurvature();
	
	void SetPMaxCurvature(double k1);
	double GetPMaxCurvature();

	void SetPMinCurvature(double k2);
	double GetPMinCurvature();

    void CalNormal(double normal_[]);
    void CalNormal() {CalNormal(normal);}
    void SetNormal(double nx, double ny, double nz );
    void GetNormal(double &nx, double &ny, double &nz);

	void SetBoundaryDis(double dist);
	double GetBoundaryDis();

	void SetDensityFuncValue(double density) {m_densityFuncValue=density;};
	double GetDensityFuncValue() {return m_densityFuncValue;};

	void SetMeshPatchPtr(QMeshPatch* _mesh);
	QMeshPatch* GetMeshPatchPtr();

	void AddFace(QMeshFace *_face);
	int GetFaceNumber();
	QMeshFace* GetFaceRecordPtr(int No);	//from 1 to n
    GLKObList& GetFaceList();

	void AddEdge(QMeshEdge *_edge);
	int GetEdgeNumber();
	QMeshEdge* GetEdgeRecordPtr(int No);	//from 1 to n
    GLKObList& GetEdgeList();

	void AddNode(QMeshNode *_node);
	int GetNodeNumber();
	QMeshNode* GetNodeRecordPtr(int No);	//from 1 to n
    GLKObList& GetNodeList();
	bool IsNodeInNodeList(QMeshNode *_node);

	void SetMinCurvatureVector(double vx, double vy, double vz);
	void GetMinCurvatureVector(double &vx, double &vy, double &vz);

	void SetMaxCurvatureVector(double vx, double vy, double vz);
	void GetMaxCurvatureVector(double &vx, double &vy, double &vz);

	GLKObList attachedList;

	double m_trackingPos[3];	int m_collideConstraintIndex;
	QMeshFace* m_trackingFace;

	void* attachedPointer;

private:
	int indexno;
	bool flags[8];
				// bit 0 -- True for boundary points
				// bit 1 -- True for points on coarse mesh
				// bit 2 -- True for points on interpolation curve 
				// bit 3 -- True for points on hand (temp use) (or points adjacent to boundary face)
				// bit 4 -- True for points on arm (temp use) (or branch vertices)
				// bit 5 -- True for sharp-feature vertex (or vertex cannot be moved)
				// bit 6 -- True for sharp-feature-region vertex
				// bit 7 -- True for points moved (temp use or newly created)

	double coord2D[2];
				// 2D space coordinate
	double coord3D[3];
				// 3D space coordinate
	double coord3D_last[3];  // just for reset sewing operation for one step
                // 3D space coordinate in the last sewing step
    double normal[3]; // Point normal on the mesh surface

	double m_meanCurvatureNormalVector[3], m_gaussianCurvature, m_pMaxCurvature, m_pMinCurvature;
	double m_minCurvatureVector[3], m_maxCurvatureVector[3];
	double m_boundaryDist, m_densityFuncValue;

	QMeshPatch *meshSurface;		// QMesh contains this node

	GLKObList faceList;	// a list of faces contained this node (QMeshFace)
	GLKObList edgeList;	// a list of edges contained this node (QMeshEdge)
	GLKObList nodeList;	// a list of nodes coincident with this node (QMeshNode)
};

#endif
