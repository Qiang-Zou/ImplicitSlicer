// QMeshEdge.h: interface for the QMeshEdge class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _CW_QMESHEDGE
#define _CW_QMESHEDGE

#include "../GLKLib/GLKObList.h"
#include "QMeshPatch.h"
#include "QMeshFace.h"
#include "QMeshNode.h"

class QMeshPatch;
class QMeshNode;
class QMeshFace;

class QMeshEdge : public GLKObject  
{
public:
	QMeshEdge();
	virtual ~QMeshEdge();

public:
	bool GetAttribFlag( const int whichBit );
	void SetAttribFlag( const int whichBit, const bool toBe = true );

	int GetIndexNo();		//from 1 to n
	void SetIndexNo( const int _index = 1 );

	bool IsBoundaryEdge();

	QMeshNode * GetStartPoint();
	void SetStartPoint( QMeshNode * _pStartPoint = NULL );

	QMeshNode * GetEndPoint();
	void SetEndPoint( QMeshNode * _pEndPoint = NULL );

	QMeshFace * GetLeftFace();
	void SetLeftFace( QMeshFace * _pLeftFace = NULL );

	QMeshFace * GetRightFace();
	void SetRightFace( QMeshFace * _pRightFace = NULL );

	void SetMeshPatchPtr(QMeshPatch* _mesh);
	QMeshPatch* GetMeshPatchPtr();

	void CalNormal(double normal[]);

	void SetSharpFactor(int factor);
	int GetSharpFactor();

	double CalLength();
	double GetLength() {return m_edgeLength;};
	double Cal2DLength();
	double Get2DLength() {return m_edge2DLength;};

    GLKObList& GetAttachedList() {return attachedList;};

	////******************************************************************
	////	Interpolate Points List
	//int attrPntNum;
	//float **attrPnt;		//	attrPnt[attrPntNum][3]
	//float *attrLength;		//	attrLength[attrPntNum-1]
	//float attrTotalLength;	

	void *attachedPointer;

    int regionind; // for storing info. for the variational Shape Approximation class

private:
	int indexno;
	bool flags[8];
				// bit 0 -- TRUE for boundary edge
				// bit 1 -- TRUE for sharp-feature edge	(or edges on key feature curves)
				// bit 2 -- TRUE for sharp-feature-region edge (or edges on accessory feature curves)
				// bit 4 -- TRUE for the edge needs to be splitted
				// bit 6 -- TRUE for edges need to be cut off
				// bit 7 -- temp use

                                    //*** Edge vector definition ***
                                    //                             *
		                            //         end point           *
	QMeshNode * pStartPoint;		//           /|\               *
	QMeshNode * pEndPoint;			//            |                *
                                    //  left face | right face     *
	QMeshFace * pLeftFace;			//            |                *
	QMeshFace * pRightFace;			//            |                *
		                            //       start point           *
		                            //                             *
		                            //******************************

	QMeshPatch *meshSurface;		// QMesh contain this edge

	GLKObList attachedList;			// a list of attached object

	int m_sharpFactor;	
	double m_edgeLength,m_edge2DLength;
};

#endif
