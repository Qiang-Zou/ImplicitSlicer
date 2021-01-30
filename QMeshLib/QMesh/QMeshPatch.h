// QMeshPatch.h: interface for the QMeshPatch class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _CW_QMESHPATCH
#define _CW_QMESHPATCH

#include <vector>
#include <algorithm>
using namespace std;

#include "../GLKLib/GLKObList.h"

#include "QMeshNode.h"
#include "QMeshEdge.h"
#include "QMeshFace.h"

class QMeshFace;
class QMeshEdge;
class QMeshNode;

class QMeshPatch : public GLKObject
{
public:
    QMeshPatch();
    virtual ~QMeshPatch();

public:
    void ClearAll();

    bool GetAttribFlag( const int whichBit );
    void SetAttribFlag( const int whichBit, const bool toBe = true );

    int GetIndexNo();		//from 1 to n
    void SetIndexNo( const int _index = 1 );

    int GetFaceNumber();
    QMeshFace* GetFaceRecordPtr(int No);	//from 1 to n
    GLKObList& GetFaceList();
    QMeshFace* GetFaceAt(GLKPOSITION Pos) {auto t = faceList.GetAt(Pos);
                                           auto tt = (QMeshFace*)t;
        return (QMeshFace*)(faceList.GetAt(Pos));}
    GLKPOSITION FaceBegin() {return faceList.GetHeadPosition();}
    GLKPOSITION FaceEnd() {return NULL;}

    int GetEdgeNumber();
    QMeshEdge* GetEdgeRecordPtr(int No);	//from 1 to n
    GLKObList& GetEdgeList();
    QMeshEdge* GetEdgeAt(GLKPOSITION Pos) {return (QMeshEdge*)(edgeList.GetAt(Pos));}
    GLKPOSITION EdgeBegin() {return edgeList.GetHeadPosition();}
    GLKPOSITION EdgeEnd() {return NULL;}

    int GetNodeNumber();
    QMeshNode* GetNodeRecordPtr(int No);	//from 1 to n
    GLKObList& GetNodeList();
    QMeshNode* GetVertAt(GLKPOSITION Pos) {return (QMeshNode*)(nodeList.GetAt(Pos));}
    GLKPOSITION VertBegin() {return nodeList.GetHeadPosition();}
    GLKPOSITION VertEnd() {return NULL;}

    void SetMaterial(bool bDir, int material);
    int GetMaterial(bool bDir);

    bool inputOBJFile(char* filename);
    bool inputMFile(char* filename);
    bool inputPLY2File(char* filename);
    bool inputOFFFile(char* filename);

    void outputOBJFile(char* filename, bool b2D=false);
    void outputTrglOBJFile(char* filename);
    void outputOFFFile(char* filename);

    void InverseOrientation();

    void constructionFromVerFaceTable(int nodeNum, float *nodeTable, int faceNum, unsigned int* faceTable);

    vector<double> getBBox(); // 0-2 elements in the returned vector store left min, and the rest right max
    double CompAverageEdgeLength();

private:
    int indexno;			// start from 1 to n

    bool flags[8];			// bit 0 -- TRUE for displaying the valence on nodes
    //			FALSE for NOT displaying the valence on nodes
    // bit 1 -- TRUE for displaying the tensile energy on edges
    //			FALSE for NOT displaying the tensile energy on edges
    // bit 2 -- TRUE for the 2D pattern has been determined
    //			FALSE for the 2D pattern has NOT been determined

    int m_materialPositiveDir,m_materialNegativeDir;

    GLKObList faceList;		// a list of mesh's faces (QMeshFace)
    GLKObList edgeList;		// a list of mesh's edges (QMeshEdge)
    GLKObList nodeList;		// a list of mesh's nodes (QMeshNode)
};


#endif
