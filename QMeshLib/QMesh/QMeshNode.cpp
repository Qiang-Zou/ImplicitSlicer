// QMeshNode.cpp: implementation of the QMeshNode class.
//
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include "QMeshPatch.h"
#include "QMeshFace.h"
#include "QMeshEdge.h"
#include "QMeshNode.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

QMeshNode::QMeshNode()
{
	indexno=0;		m_trackingFace=NULL;
	m_meanCurvatureNormalVector[0]=0.0;	
	m_meanCurvatureNormalVector[1]=0.0;	
	m_meanCurvatureNormalVector[2]=0.0;	
	m_gaussianCurvature=0.0;	m_pMaxCurvature=0.0;	m_pMinCurvature=0.0;
	m_boundaryDist=0.0;
	faceList.RemoveAll();
	edgeList.RemoveAll();
	nodeList.RemoveAll();
	for(int i=0;i<8;i++) flags[i]=false;

	m_trackingPos[0]=0.0;
	m_trackingPos[1]=0.0;
	m_trackingPos[2]=0.0;
	attachedPointer=NULL;
}

QMeshNode::~QMeshNode()
{
	faceList.RemoveAll();
	edgeList.RemoveAll();
	nodeList.RemoveAll();
}

//////////////////////////////////////////////////////////////////////
// Implementation
//////////////////////////////////////////////////////////////////////

int QMeshNode::GetIndexNo() 
{
	return indexno;
}
	
void QMeshNode::SetIndexNo( const int _index )
{
	indexno=_index;
}

bool QMeshNode::GetAttribFlag( const int whichBit )
{
	return flags[whichBit];
}

void QMeshNode::SetAttribFlag( const int whichBit, const bool toBe )
{
	flags[whichBit]=toBe;
}

void QMeshNode::GetCoord2D( double &x, double &y )
{
	x=coord2D[0];	y=coord2D[1];
}

void QMeshNode::SetCoord2D( double x, double y )
{
	coord2D[0]=x;	coord2D[1]=y;
}

void QMeshNode::GetCoord3D( double &x, double &y, double &z )
{
	x=coord3D[0];	y=coord3D[1];	z=coord3D[2];
}

void QMeshNode::SetCoord3D( double x, double y, double z )
{
	coord3D[0]=x;	coord3D[1]=y;	coord3D[2]=z;
}

void QMeshNode::GetCoord3D_last( double &x, double &y, double &z )
{
	x=coord3D_last[0];	y=coord3D_last[1];	z=coord3D_last[2];
}

void QMeshNode::SetCoord3D_last( double x, double y, double z )
{
	coord3D_last[0]=x;	coord3D_last[1]=y;	coord3D_last[2]=z;
}

void QMeshNode::SetMeanCurvatureNormalVector(double kHx, double kHy, double kHz)
{
	m_meanCurvatureNormalVector[0]=kHx;
	m_meanCurvatureNormalVector[1]=kHy;
	m_meanCurvatureNormalVector[2]=kHz;
}

void QMeshNode::GetMeanCurvatureNormalVector(double &kHx, double &kHy, double &kHz)
{
	kHx=m_meanCurvatureNormalVector[0];
	kHy=m_meanCurvatureNormalVector[1];
	kHz=m_meanCurvatureNormalVector[2];
}

void QMeshNode::SetGaussianCurvature(double kG)
{
	m_gaussianCurvature=kG;
}

double QMeshNode::GetGaussianCurvature()
{
	return m_gaussianCurvature;
}
	
void QMeshNode::SetPMaxCurvature(double k1)
{
	m_pMaxCurvature=k1;
}

double QMeshNode::GetPMaxCurvature()
{
	return m_pMaxCurvature;
}

void QMeshNode::SetPMinCurvature(double k2)
{
	m_pMinCurvature=k2;
}

double QMeshNode::GetPMinCurvature()
{
	return m_pMinCurvature;
}

void QMeshNode::SetMinCurvatureVector(double vx, double vy, double vz)
{
	m_minCurvatureVector[0]=vx;	m_minCurvatureVector[1]=vy;	m_minCurvatureVector[2]=vz;
}

void QMeshNode::GetMinCurvatureVector(double &vx, double &vy, double &vz)
{
	vx=m_minCurvatureVector[0];	vy=m_minCurvatureVector[1];	vz=m_minCurvatureVector[2];
}

void QMeshNode::SetMaxCurvatureVector(double vx, double vy, double vz)
{
	m_maxCurvatureVector[0]=vx;	m_maxCurvatureVector[1]=vy;	m_maxCurvatureVector[2]=vz;
}

void QMeshNode::GetMaxCurvatureVector(double &vx, double &vy, double &vz)
{
	vx=m_maxCurvatureVector[0];	vy=m_maxCurvatureVector[1];	vz=m_maxCurvatureVector[2];
}

void QMeshNode::SetBoundaryDis(double dist)
{
	m_boundaryDist=dist;
}

double QMeshNode::GetBoundaryDis()
{
	return m_boundaryDist;
}

void QMeshNode::CalNormal(double normal_[])
{
	double nx,ny,nz,tt;
	nx=0.0;	ny=0.0;	nz=0.0;

	GLKPOSITION Pos;
	for(Pos=faceList.GetHeadPosition();Pos!=NULL;)
	{
		double a,b,c,d;

		QMeshFace *temp=(QMeshFace *)(faceList.GetNext(Pos));
		temp->GetPlaneEquation(a,b,c,d);
		nx+=a;	ny+=b;	nz+=c;
	}
	tt=nx*nx+ny*ny+nz*nz;
	tt=sqrt(tt);

    normal_[0]=(double)(nx/tt);	normal_[1]=(double)(ny/tt);	normal_[2]=(double)(nz/tt);
}

void QMeshNode::GetNormal(double &nx, double &ny, double &nz)
{
    nx=normal[0];	ny=normal[1];	nz=normal[2];
}

void QMeshNode::SetNormal(double nx, double ny, double nz)
{
    normal[0]=nx;	normal[1]=ny;	normal[2]=nz;
}

void QMeshNode::SetMeshPatchPtr(QMeshPatch* _mesh)
{
	meshSurface=_mesh;
}

QMeshPatch* QMeshNode::GetMeshPatchPtr()
{
	return meshSurface;
}

void QMeshNode::AddFace(QMeshFace *_face)
{
	faceList.AddTail(_face);
}

int QMeshNode::GetFaceNumber() 
{
	return faceList.GetCount();
}

QMeshFace* QMeshNode::GetFaceRecordPtr(int No) 	//from 1 to n
{
	if( (No < 1) || (No > faceList.GetCount()))    return  NULL;
    return (QMeshFace *)faceList.GetAt(faceList.FindIndex(No-1));
}

GLKObList& QMeshNode::GetFaceList()
{
	return faceList;
}

void QMeshNode::AddEdge(QMeshEdge *_edge)
{
	edgeList.AddTail(_edge);
}

int QMeshNode::GetEdgeNumber() 
{
	return edgeList.GetCount();
}

QMeshEdge* QMeshNode::GetEdgeRecordPtr(int No) 	//from 1 to n
{
	if( (No < 1) || (No > edgeList.GetCount()))    return  NULL;
    return (QMeshEdge *)edgeList.GetAt(edgeList.FindIndex(No-1));
}

GLKObList& QMeshNode::GetEdgeList()
{
	return edgeList;
}

void QMeshNode::AddNode(QMeshNode *_node)
{
	nodeList.AddTail(_node);
}

int QMeshNode::GetNodeNumber()
{
	return nodeList.GetCount();
}

bool QMeshNode::IsNodeInNodeList(QMeshNode *_node)
{
	GLKPOSITION Pos;

	for(Pos=nodeList.GetHeadPosition();Pos!=NULL;) {
		QMeshNode *tempnode=(QMeshNode *)(nodeList.GetNext(Pos));
		if (tempnode==_node) return true;
	}

	return false;
}

QMeshNode* QMeshNode::GetNodeRecordPtr(int No)	//from 1 to n
{
	if ((No < 1) || (No > nodeList.GetCount())) return  NULL;
    return (QMeshNode *)nodeList.GetAt(nodeList.FindIndex(No-1));
}

GLKObList& QMeshNode::GetNodeList()
{
	return nodeList;
}
