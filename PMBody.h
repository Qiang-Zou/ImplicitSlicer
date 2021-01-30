#ifndef _CCL_POLYMESH_BODY
#define _CCL_POLYMESH_BODY

#include "GLKLib/GLK.h"
#include "GLKLib/GLKObList.h"

#define SHARPEDGE_NORMAL_THRESHOLD		0.85

class QMeshFace;
class QMeshEdge;
class QMeshNode;

class PMBody : public GLKEntity  
{
public:
	PMBody();
	virtual ~PMBody();

	void DeleteGLList();
	void BuildGLList(bool bVertexNormalShading, bool bDisplayBndDistContour=false, int activeRegion=0, int activePatchID=0);

	virtual void drawShade();
	virtual void drawMesh();
	virtual void drawProfile();
	virtual void drawPreMesh();
	virtual void drawHighLight();
	virtual float getRange() {return m_range;}

	void drawBox(float xx, float yy, float zz, float r);

	void ClearAll();
	void computeRange();
	
	GLKObList &GetMeshList() {return meshList;};

	void ExportPolytopeFile(char *filename);
	void ExportRegionOBJFile(char *filename, int regionID);

	int m_materialTypeNum;

	void FlipModel(short nDir);
	void CompBoundingBox(double boundingBox[]);
	void Transformation(double dx, double dy, double dz);

	double CompAverageEdgeLength();

public:
	GLKObList meshList;
	float m_range;
	int m_drawListID;

	void _buildDrawShadeList(bool bVertexNormalShading, bool bDisplayBndDistContour);
	void _buildDrawShadeListMultiMaterial(bool bVertexNormalShading);
	void _buildDrawMeshList();
	void _buildDrawProfileList();
	void _buildDrawPreMeshList();
	void _changeValueToColor(int nType, float & nRed, float & nGreen, float & nBlue);
	void _changeValueToColor(double maxValue, double minValue, double Value, 
								 float & nRed, float & nGreen, float & nBlue);

	void _buildStripObjColorNormalDisplayList();
	void _extendTrglStrip(QMeshFace *face, QMeshEdge *edge, bool orderFlag);

	void _buildDrawShadeListWithNormalProcessed();

	bool _compLocalNormalForSharpVertex(QMeshNode* node, double **normalx, double **normaly, double **normalz);

private:
	int activeMaterialRegion,activePatch;
};

#endif
