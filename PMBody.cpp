// PMBody.cpp: implementation of the PMBody class.
//
//////////////////////////////////////////////////////////////////////
#define _CRT_SECURE_NO_DEPRECATE

#include <math.h>
#include <memory.h>

#include "PMBody.h"

#include "QMeshLib/QMesh/QMeshPatch.h"
#include "QMeshLib/QMesh/QMeshFace.h"
#include "QMeshLib/QMesh/QMeshEdge.h"
#include "QMeshLib/QMesh/QMeshNode.h"


extern GLK _pGLK;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

PMBody::PMBody()
{
	meshList.RemoveAll();
	m_drawListID=-1;
	m_materialTypeNum=1;
	activeMaterialRegion=activePatch=0;
}

PMBody::~PMBody()
{
	ClearAll();
	if (m_drawListID!=-1) glDeleteLists(m_drawListID, 4);
}

//////////////////////////////////////////////////////////////////////
// Implementation
//////////////////////////////////////////////////////////////////////

void PMBody::FlipModel(short nDir)
{
	GLKPOSITION PosMesh;	
	GLKPOSITION Pos;
	double bndBox[6],cx,cy,cz,xx,yy,zz;

	CompBoundingBox(bndBox);	
	cx=(bndBox[0]+bndBox[1])*0.5;
	cy=(bndBox[2]+bndBox[3])*0.5;
	cz=(bndBox[4]+bndBox[5])*0.5;

	for(PosMesh=meshList.GetHeadPosition();PosMesh!=NULL;) {
		QMeshPatch *mesh=(QMeshPatch *)(meshList.GetNext(PosMesh));
		for(Pos=mesh->GetNodeList().GetHeadPosition();Pos!=NULL;) {
			QMeshNode *node=(QMeshNode *)(mesh->GetNodeList().GetNext(Pos));
			node->GetCoord3D(xx,yy,zz);

			switch(nDir) {
			case 0:{xx=cx+(cx-xx);   }break;
			case 1:{yy=cy+(cy-yy);   }break;
			case 2:{zz=cz+(cz-zz);   }break;
			}
			node->SetCoord3D(xx,yy,zz);
		}

		mesh->InverseOrientation();
	}
}

double PMBody::CompAverageEdgeLength()
{
	double div, avgLength;
	GLKPOSITION PosMesh;
	GLKPOSITION Pos;

	div = avgLength = 0.0;
	for (PosMesh = meshList.GetHeadPosition(); PosMesh != NULL;) {
		QMeshPatch *mesh = (QMeshPatch *)(meshList.GetNext(PosMesh));
		for (Pos = mesh->GetEdgeList().GetHeadPosition(); Pos != NULL;) {
			QMeshEdge *edge = (QMeshEdge *)(mesh->GetEdgeList().GetNext(Pos));
            div += 1.0;
            avgLength += edge->CalLength();
		}
	}
	avgLength = avgLength / div;
	printf("Average Length = %lf \n", avgLength);

	return avgLength;
}

void PMBody::CompBoundingBox(double boundingBox[])
{
	GLKPOSITION PosMesh;	
	GLKPOSITION Pos;
	double xx,yy,zz;

	boundingBox[0]=boundingBox[2]=boundingBox[4]=1.0e+32;	
	boundingBox[1]=boundingBox[3]=boundingBox[5]=-1.0e+32;	

	for(PosMesh=meshList.GetHeadPosition();PosMesh!=NULL;) {
		QMeshPatch *mesh=(QMeshPatch *)(meshList.GetNext(PosMesh));
		for(Pos=mesh->GetNodeList().GetHeadPosition();Pos!=NULL;) {
			QMeshNode *node=(QMeshNode *)(mesh->GetNodeList().GetNext(Pos));
			node->GetCoord3D(xx,yy,zz);

			if (xx<boundingBox[0]) boundingBox[0]=xx;
			if (xx>boundingBox[1]) boundingBox[1]=xx;
			if (yy<boundingBox[2]) boundingBox[2]=yy;
			if (yy>boundingBox[3]) boundingBox[3]=yy;
			if (zz<boundingBox[4]) boundingBox[4]=zz;
			if (zz>boundingBox[5]) boundingBox[5]=zz;
		}
	}
}

void PMBody::Transformation(double dx, double dy, double dz)
{
	GLKPOSITION PosMesh;	
	GLKPOSITION Pos;
	double xx,yy,zz;

	for(PosMesh=meshList.GetHeadPosition();PosMesh!=NULL;) {
		QMeshPatch *mesh=(QMeshPatch *)(meshList.GetNext(PosMesh));
		for(Pos=mesh->GetNodeList().GetHeadPosition();Pos!=NULL;) {
			QMeshNode *node=(QMeshNode *)(mesh->GetNodeList().GetNext(Pos));
			node->GetCoord3D(xx,yy,zz);
			node->SetCoord3D(xx+dx,yy+dy,zz+dz);
		}
	}
}

void PMBody::ExportPolytopeFile(char *filename)
{
	FILE *fp;	GLKPOSITION Pos;	GLKPOSITION PosMesh;	GLKPOSITION PosFace;
	QMeshNode *node;	QMeshFace *face,*otherFace;	
	QMeshEdge *edge;	QMeshPatch *mesh;
	GLKObList faceList,newFaceList;
	int i,j,eNum,*nodeNum,*faceNum; bool *bVisited;
	double xx,yy,zz;

	//-----------------------------------------------------------------------------
	//	Step 1: search the polytopes
	for(PosMesh=meshList.GetHeadPosition();PosMesh!=NULL;) {
		mesh=(QMeshPatch *)(meshList.GetNext(PosMesh));
		for(Pos=mesh->GetFaceList().GetHeadPosition();Pos!=NULL;) {
			face=(QMeshFace *)(mesh->GetFaceList().GetNext(Pos));
			face->m_nIdentifiedPatchIndex=-1;
		}
	}
	//-----------------------------------------------------------------------------
	//	Flooding based labeling of faces and find the total number of polytopes
	int nPolytopeNum=0;		QMeshFace *seedFace;
	do{
		seedFace=NULL;	faceList.RemoveAll();
		for(PosMesh=meshList.GetHeadPosition();PosMesh!=NULL;) {
			mesh=(QMeshPatch *)(meshList.GetNext(PosMesh));
			for(Pos=mesh->GetFaceList().GetHeadPosition();Pos!=NULL;) {
				face=(QMeshFace *)(mesh->GetFaceList().GetNext(Pos));
				if (face->m_nIdentifiedPatchIndex<0) {seedFace=face; break;}
			}
			if (seedFace!=NULL) break;
		}
		if (seedFace!=NULL) {
			nPolytopeNum++; 
			seedFace->m_nIdentifiedPatchIndex=nPolytopeNum; 
			faceList.AddTail(seedFace);
		}

		//-------------------------------------------------------------------------
		// flooding across the edge of faces
		bool bNewFaceFound;
		do{
			bNewFaceFound=false;	newFaceList.RemoveAll();
			for(Pos=faceList.GetHeadPosition();Pos!=NULL;) {
				face=(QMeshFace *)(faceList.GetNext(Pos));
				eNum=face->GetEdgeNum();
				for(i=0;i<eNum;i++) {
					edge=face->GetEdgeRecordPtr(i);
					if (edge->GetLeftFace()==face)
						otherFace=edge->GetRightFace();
					else
						otherFace=edge->GetLeftFace();
					if (otherFace!=NULL && otherFace->m_nIdentifiedPatchIndex<0) {
						otherFace->m_nIdentifiedPatchIndex=nPolytopeNum;
						newFaceList.AddTail(otherFace);
						bNewFaceFound=true;
					}
				}
			}
			faceList.RemoveAll();	faceList.AddTail(&newFaceList);
		}while(bNewFaceFound);

	}while(seedFace!=NULL);
	printf("Total %d polytopes are found.\n",nPolytopeNum);
	//-----------------------------------------------------------------------------
	if (nPolytopeNum>0) {
		nodeNum=new int[nPolytopeNum]; faceNum=new int[nPolytopeNum]; bVisited=new bool[nPolytopeNum];
	}

	//-----------------------------------------------------------------------------
	//	Step 2: count the number of nodes and number of triangles in each polytope
	memset(nodeNum,0,nPolytopeNum*sizeof(int));
	memset(faceNum,0,nPolytopeNum*sizeof(int));
	for(PosMesh=meshList.GetHeadPosition();PosMesh!=NULL;) {
		mesh=(QMeshPatch *)(meshList.GetNext(PosMesh));
		for(Pos=mesh->GetNodeList().GetHeadPosition();Pos!=NULL;) {
			node=(QMeshNode *)(mesh->GetNodeList().GetNext(Pos));
			memset(bVisited,false,nPolytopeNum*sizeof(bool));
			for(PosFace=node->GetFaceList().GetHeadPosition();PosFace!=NULL;) {
				face=(QMeshFace *)(node->GetFaceList().GetNext(PosFace));
				i=face->m_nIdentifiedPatchIndex-1;
				if (bVisited[i]) continue;
				nodeNum[i]++;	bVisited[i]=true;
			}
		}
	}
	//-----------------------------------------------------------------------------
	for(PosMesh=meshList.GetHeadPosition();PosMesh!=NULL;) {
		mesh=(QMeshPatch *)(meshList.GetNext(PosMesh));
		for(Pos=mesh->GetFaceList().GetHeadPosition();Pos!=NULL;) {
			face=(QMeshFace *)(mesh->GetFaceList().GetNext(Pos));
			eNum=face->GetEdgeNum();
			(faceNum[face->m_nIdentifiedPatchIndex-1])+=(eNum-2);
		}
	}
	//-----------------------------------------------------------------------------
//	for(i=0;i<nPolytopeNum;i++) {printf("node_Num: %d    face_Num: %d\n", nodeNum[i], faceNum[i]);}

	//-----------------------------------------------------------------------------
	//	Step 3: save the polytopes into the file one by one
	fp = fopen(filename, "w");
    if(!fp) {
		printf("===============================================");  
		printf("Can not open the data file - PLT File Export!");
		printf("===============================================");  return;
	}
	//-----------------------------------------------------------------------------
	fprintf(fp,"Total_polytope_num: %d\n",nPolytopeNum);
	for(i=0;i<nPolytopeNum;i++) {
		fprintf(fp,"---------------------------------------\n");
		fprintf(fp,"Polytope_ID: %d\n",i+1);
		fprintf(fp,"node_Num: %d    face_Num: %d\n", nodeNum[i], faceNum[i]);
		int nodeIndex=1;
		for(PosMesh=meshList.GetHeadPosition();PosMesh!=NULL;) {
			mesh=(QMeshPatch *)(meshList.GetNext(PosMesh));
			for(Pos=mesh->GetNodeList().GetHeadPosition();Pos!=NULL;) {
				node=(QMeshNode *)(mesh->GetNodeList().GetNext(Pos));
				for(PosFace=node->GetFaceList().GetHeadPosition();PosFace!=NULL;) {
					face=(QMeshFace *)(node->GetFaceList().GetNext(PosFace));
					if ((face->m_nIdentifiedPatchIndex-1)==i) {
						node->SetIndexNo(nodeIndex);	nodeIndex++;
						node->GetCoord3D(xx,yy,zz);
						fprintf(fp,"v %lf %lf %lf\n",xx,yy,zz);
						break;
					}
				}
			}
		}
		for(PosMesh=meshList.GetHeadPosition();PosMesh!=NULL;) {
			mesh=(QMeshPatch *)(meshList.GetNext(PosMesh));
			for(Pos=mesh->GetFaceList().GetHeadPosition();Pos!=NULL;) {
				face=(QMeshFace *)(mesh->GetFaceList().GetNext(Pos));
				if ((face->m_nIdentifiedPatchIndex-1)!=i) continue;
				eNum=face->GetEdgeNum();
				for(j=0;j<(eNum-2);j++) {
					fprintf(fp,"f %d %d %d\n",
						face->GetNodeRecordPtr(0)->GetIndexNo(),
						face->GetNodeRecordPtr(j+1)->GetIndexNo(),
						face->GetNodeRecordPtr(j+2)->GetIndexNo());
				}
			}
		}
	}
	//-----------------------------------------------------------------------------
	fclose(fp);
	printf("--------------------------------------------\n");
	printf("File output is completed.\n");

	//-----------------------------------------------------------------------------
	//	Step 4: release the memory
	if (nPolytopeNum>0) {delete []nodeNum;  delete []faceNum;  delete []bVisited;}
}

void PMBody::ExportRegionOBJFile(char *filename, int regionID)
{
	FILE *fp;	bool bDir;	int i,num,index;
	GLKPOSITION Pos;
	GLKPOSITION PosMesh;	
	QMeshNode *node;	QMeshFace *face;	QMeshPatch *mesh;
	double xx,yy,zz;

	fp = fopen(filename, "w");
    if(!fp)
	{ 
		printf("===============================================");
	    printf("Can not open the data file - OBJ File Export!");
		printf("===============================================");
	    return;
	}

	i=1;
	for(PosMesh=meshList.GetHeadPosition();PosMesh!=NULL;) {
		mesh=(QMeshPatch *)(meshList.GetNext(PosMesh));
		if ((mesh->GetMaterial(true)!=regionID) && (mesh->GetMaterial(false)!=regionID)) continue;
		for(Pos=mesh->GetNodeList().GetHeadPosition();Pos!=NULL;i++) {
			node=(QMeshNode *)(mesh->GetNodeList().GetNext(Pos));
			node->SetIndexNo(i);
			node->GetCoord3D(xx,yy,zz);
			fprintf(fp,"v %.12f %.12f %.12f\n",xx,yy,zz);
		}
//		break;
	}

	fprintf(fp,"\n");

	for(PosMesh=meshList.GetHeadPosition();PosMesh!=NULL;) {
		mesh=(QMeshPatch *)(meshList.GetNext(PosMesh));
		if ((mesh->GetMaterial(true)!=regionID) && (mesh->GetMaterial(false)!=regionID)) continue;
		if (mesh->GetMaterial(false)==regionID) bDir=false; else bDir=true;
		for(Pos=mesh->GetFaceList().GetHeadPosition();Pos!=NULL;) {
			face=(QMeshFace *)(mesh->GetFaceList().GetNext(Pos));
			num=face->GetEdgeNum();
			fprintf(fp,"f ");
//			bDir=true;
			if (bDir)
			for(i=0;i<num;i++) {
				index=face->GetNodeRecordPtr(i)->GetIndexNo();
				fprintf(fp,"%d ",index);
			}
			else
			for(i=num-1;i>=0;i--) {
				index=face->GetNodeRecordPtr(i)->GetIndexNo();
				fprintf(fp,"%d ",index);
			}
			fprintf(fp,"\n");
		}
//		break;
	}

	fclose(fp);
}

void PMBody::DeleteGLList()
{
	if (m_drawListID!=-1) {
		glDeleteLists(m_drawListID, 4);
		m_drawListID=-1;
	}
}

void PMBody::BuildGLList(bool bVertexNormalShading, bool bDisplayBndDistContour, int activeRegion, int activePatchID)
{
	if (activePatchID!=0) {activePatch=activePatchID;	activeRegion=0;}

	if (m_drawListID!=-1) glDeleteLists(m_drawListID, 4);
	m_drawListID = glGenLists(4);	
//	printf("PMBody GL_ID = %d\n",m_drawListID);

	//if (bVertexNormalShading) 
	//	_buildDrawShadeListWithNormalProcessed(); 
	//else
	activeMaterialRegion=activeRegion;
	if (m_materialTypeNum==1) {
		_buildDrawShadeList(bVertexNormalShading, bDisplayBndDistContour);
	}
	else {
		_buildDrawShadeListMultiMaterial(bVertexNormalShading);
	}

//	_buildStripObjColorNormalDisplayList();
	_buildDrawMeshList();
	_buildDrawPreMeshList();
	_buildDrawProfileList();
}

void PMBody::_buildDrawShadeList(bool bVertexNormalShading, bool bDisplayBndDistContour)
{
	GLKPOSITION Pos;
	GLKPOSITION PosFace;
	QMeshFace *face;
	QMeshNode *node;
	QMeshPatch *mesh;
	double xx,yy,zz,dd;		float rr,gg,bb;
	int k,i,num,meshIndex;

	double maxValue=0.0;
	if (bDisplayBndDistContour) {
		for(Pos=meshList.GetHeadPosition();Pos!=NULL;) {
			mesh=(QMeshPatch *)(meshList.GetNext(Pos));
			for(PosFace=mesh->GetNodeList().GetHeadPosition();PosFace!=NULL;) {
				QMeshNode *node=(QMeshNode *)(mesh->GetNodeList().GetNext(PosFace));
				if (node->GetBoundaryDis()>maxValue) maxValue=node->GetBoundaryDis();
			}
		}
//		printf("maxValue=%lf\n",maxValue);
	}

	glNewList(m_drawListID, GL_COMPILE);

	glEnable(GL_NORMALIZE);
	glEnable(GL_LIGHTING);

	meshIndex=0;
	for(Pos=meshList.GetHeadPosition();Pos!=NULL;meshIndex++) {
		mesh=(QMeshPatch *)(meshList.GetNext(Pos));

		glBegin(GL_TRIANGLES);
		for(PosFace=(mesh->GetFaceList()).GetHeadPosition();PosFace!=NULL;)
		{
			face=(QMeshFace *)((mesh->GetFaceList()).GetNext(PosFace));

			if (face->m_nIdentifiedPatchIndex<0) {
				if (meshIndex==0) {rr=gg=bb=0.8f;}
				else if (meshIndex==1) {rr=gg=0.8f; bb=0.5f;}
				else {rr=gg=bb=0.8f;}
				glColor3f(rr,gg,bb);
				if (face->GetAttribFlag(4)) glColor3f(1.0f,1.0f,0.0f);
				if (face->GetAttribFlag(5)) glColor3f(0.0f,0.0f,0.0f);
			}
			else {
				_changeValueToColor(face->m_nIdentifiedPatchIndex,rr,gg,bb);
				glColor3f(rr,gg,bb);
			}
			if (face->GetAttribFlag(7)) glColor3f(0,1,0);

			num=face->GetEdgeNum();
			for(k=0;k<num-2;k++) {
				for(i=0;i<3;i++) {
					if (i<1)
						node=face->GetNodeRecordPtr(i);
					else
						node=face->GetNodeRecordPtr(i+k);

					if (bVertexNormalShading) { // && face->m_nIdentifiedPatchIndex>0) {
						double normal[3];
						node->CalNormal(normal);
						glNormal3dv(normal);
					}
					else {
						face->GetPlaneEquation(xx,yy,zz,dd);
						glNormal3d(xx,yy,zz);
					}
					if (bDisplayBndDistContour) {
						_changeValueToColor(maxValue,0.0,node->GetBoundaryDis(),rr,gg,bb);
						glColor3f(rr,gg,bb);
					}
					node->GetCoord3D(xx,yy,zz);
					glVertex3d(xx,yy,zz);
				}
			}
		}
		glEnd();

	//_pGLK.GetScale(scale);
	//range=_pGLK.GetRange();
	//_pGLK.GetSize(sx,sy);	width=(sx>sy)?sx:sy;
	//gwidth=m_solid->GetGridWidth();
	//glPointSize((float)(gwidth*1.0*width*0.5/range*scale)*0.866f);
		float size=0.002,range;
		range=_pGLK.GetRange();		size=size*range;
		for(PosFace=mesh->GetNodeList().GetHeadPosition();PosFace!=NULL;) {
			QMeshNode *node=(QMeshNode *)(mesh->GetNodeList().GetNext(PosFace));
			if (node->GetAttribFlag(3)) {
				glColor3d(1,0,0);
				node->GetCoord3D(xx,yy,zz);
				drawBox(xx,yy,zz,size);
			}
			if (node->GetAttribFlag(4)) {
				glColor3d(0,0,1);
				node->GetCoord3D(xx,yy,zz);
				drawBox(xx,yy,zz,size);
			}
		}
	}

	glEndList();
}

void PMBody::_buildDrawShadeListMultiMaterial(bool bVertexNormalShading)
{
	GLKPOSITION Pos;
	GLKPOSITION PosFace;
	QMeshFace *face;
	QMeshNode *node; 
	QMeshPatch *mesh;
	double xx,yy,zz,dd;		float rr,gg,bb;
	int k,i,num,meshIndex;

	glNewList(m_drawListID, GL_COMPILE);

	glEnable(GL_NORMALIZE);
	glEnable(GL_LIGHTING);

	meshIndex=0;
	for(Pos=meshList.GetHeadPosition();Pos!=NULL;meshIndex++) {
		mesh=(QMeshPatch *)(meshList.GetNext(Pos));
//		if (meshIndex==0) continue;
		if ((activePatch!=0) && (meshIndex!=activePatch-1)) continue;
		if ((activeMaterialRegion!=0) && (mesh->GetMaterial(false)!=activeMaterialRegion)
			 && (mesh->GetMaterial(true)!=activeMaterialRegion)) continue;

		if (activeMaterialRegion==0) 
			_changeValueToColor(mesh->GetMaterial(false),rr,gg,bb);
		else
			_changeValueToColor(activeMaterialRegion,rr,gg,bb);

		glColor3f(rr,gg,bb);
		glBegin(GL_TRIANGLES);
		for(PosFace=(mesh->GetFaceList()).GetHeadPosition();PosFace!=NULL;)
		{
			face=(QMeshFace *)((mesh->GetFaceList()).GetNext(PosFace));

			num=face->GetEdgeNum();
			for(k=0;k<num-2;k++) {
				for(i=0;i<3;i++) {
					if (i<1)
						node=face->GetNodeRecordPtr(i);
					else
						node=face->GetNodeRecordPtr(i+k);

					if (bVertexNormalShading) {
						double normal[3];
						node->CalNormal(normal);
						glNormal3dv(normal);
					}
					else {
						face->GetPlaneEquation(xx,yy,zz,dd);
						glNormal3d(xx,yy,zz);
					}
					node->GetCoord3D(xx,yy,zz);
					glVertex3d(xx,yy,zz);
				}
			}
		}
		glEnd();
	}

	glEndList();
}

void PMBody::_buildStripObjColorNormalDisplayList()
{
	GLKPOSITION Pos;
	GLKPOSITION PosFace;
	QMeshFace *face;		QMeshPatch *mesh;
	double xx,yy,zz,dd;		int i,meshIndex;	float rr,gg,bb;

	glNewList(m_drawListID, GL_COMPILE);

	glEnable(GL_NORMALIZE);
	glEnable(GL_LIGHTING);
	glShadeModel(GL_FLAT);

	meshIndex=0;
	for(Pos=meshList.GetHeadPosition();Pos!=NULL;meshIndex++) {
		mesh=(QMeshPatch *)(meshList.GetNext(Pos));

		i=0;
		for(PosFace=(mesh->GetFaceList()).GetHeadPosition();PosFace!=NULL;i++) {
			face=(QMeshFace *)((mesh->GetFaceList()).GetNext(PosFace));
			face->SetAttribFlag(7,false);
			face->SetIndexNo(i);
			if (face->GetEdgeNum()!=3) return;
		}

		for(PosFace=(mesh->GetFaceList()).GetHeadPosition();PosFace!=NULL;) {
			face=(QMeshFace *)((mesh->GetFaceList()).GetNext(PosFace));
			if (face->GetAttribFlag(7)) continue;

			QMeshFace *otherFace,*dirFace=NULL;	short nEdgeIndex;
			for(i=0;i<3;i++) {
				QMeshEdge *edge=face->GetEdgeRecordPtr(i);
				if (edge->GetLeftFace()==face)
					otherFace=edge->GetRightFace();
				else
					otherFace=edge->GetLeftFace();
				if (otherFace==NULL) continue;
				if (otherFace->GetAttribFlag(7)) continue;
				dirFace=otherFace;	nEdgeIndex=i;	break;
			}

			if (dirFace==NULL) continue;

			glBegin(GL_TRIANGLE_STRIP);

			face->GetNodePos(nEdgeIndex,xx,yy,zz);
			glVertex3d(xx,yy,zz);
			face->GetNodePos((nEdgeIndex+1)%3,xx,yy,zz);
			glVertex3d(xx,yy,zz);

			face->GetPlaneEquation(xx,yy,zz,dd);
			glNormal3d(xx,yy,zz);
			_changeValueToColor(face->GetIndexNo(),rr,gg,bb);
			glColor3f(rr,gg,bb);
			face->SetAttribFlag(7,true);
			_extendTrglStrip(face,face->GetEdgeRecordPtr(nEdgeIndex),true);

			glEnd();
		}

		glBegin(GL_TRIANGLES);
		for(PosFace=(mesh->GetFaceList()).GetHeadPosition();PosFace!=NULL;) {
			face=(QMeshFace *)((mesh->GetFaceList()).GetNext(PosFace));
			if (face->GetAttribFlag(7)) continue;

			face->GetPlaneEquation(xx,yy,zz,dd);
			glNormal3d(xx,yy,zz);
			_changeValueToColor(face->GetIndexNo(),rr,gg,bb);
			glColor3f(rr,gg,bb);

			face->GetNodePos(0,xx,yy,zz);	glVertex3d(xx,yy,zz);
			face->GetNodePos(1,xx,yy,zz);	glVertex3d(xx,yy,zz);
			face->GetNodePos(2,xx,yy,zz);	glVertex3d(xx,yy,zz);
		}
		glEnd();
	}

	glEndList();
}

void PMBody::_extendTrglStrip(QMeshFace *face, QMeshEdge *edge, bool orderFlag)
{
	short i;		double xx,yy,zz,dd;
	QMeshNode *node;	QMeshFace *otherFace;		QMeshEdge *nextEdge;

	for(i=0;i<3;i++) {
		QMeshEdge *temp=face->GetEdgeRecordPtr(i);
		if (temp==edge) break;
	}
	node=face->GetNodeRecordPtr((i+2)%3);
	node->GetCoord3D(xx,yy,zz);
	glVertex3d(xx,yy,zz);

	if (orderFlag)
		nextEdge=face->GetEdgeRecordPtr((i+1)%3);
	else
		nextEdge=face->GetEdgeRecordPtr((i+2)%3);
	if (nextEdge->GetLeftFace()==face)
		otherFace=nextEdge->GetRightFace();
	else
		otherFace=nextEdge->GetLeftFace();
	if (otherFace==NULL) return;
	if (otherFace->GetAttribFlag(7)) return;

	otherFace->GetPlaneEquation(xx,yy,zz,dd);
	glNormal3d(xx,yy,zz);
	float rr,gg,bb;
	_changeValueToColor(otherFace->GetIndexNo(),rr,gg,bb);
	glColor3f(rr,gg,bb);
	otherFace->SetAttribFlag(7,true);
	_extendTrglStrip(otherFace,nextEdge,!orderFlag);
}

void PMBody::_changeValueToColor(int nType, float & nRed, float & nGreen, float & nBlue)
{
	float color[][3]={
		{255,255,255},
		{255,0,128},
		{0,255,255},
		{128,255,0},
		{128,128,64},
		{255,0,0},{0,255,0},{0,0,255},
		{128,128,192},
		{255,255,128},
		{255,128,0},
		{255,128,255},
		{255,214,202},
		{128,128,192},
		{255,165,0}, //orange
		{255,128,192},
//		{39, 64, 139},//RoyalBlue
		{128,128,64},
		{0,255,255},
		{238,130,238},//violet
		{220,220,220},//gainsboro
		{188, 143, 143}, // rosy brown 
		{46, 139, 87},//sea green
		{210, 105, 30 },//chocolate
		{100, 149, 237}//cornflower blue 
	};

//	printf("%d ",nType);
	nRed=color[nType%22][0]/255.0f;
	nGreen=color[nType%22][1]/255.0f;
	nBlue=color[nType%22][2]/255.0f;
}

void PMBody::_buildDrawShadeListWithNormalProcessed()
{
	GLKPOSITION Pos;
	GLKPOSITION PosFace;
	GLKPOSITION PosEdge;
	GLKPOSITION PosNode;
	QMeshFace *face;
	QMeshEdge *edge;
	QMeshNode *node;
	QMeshPatch *mesh;
	double **normalx,**normaly,**normalz;
	int id,i,k,num;	double xx,yy,zz,dd,nx,ny,nz;	float rr,gg,bb;
	double normal[3];
	if (meshList.GetCount()==0) return;

	glNewList(m_drawListID, GL_COMPILE);

	glEnable(GL_NORMALIZE);
	glEnable(GL_LIGHTING);	int meshIndex=0;

	for(Pos=meshList.GetHeadPosition();Pos!=NULL;meshIndex++) {
		mesh=(QMeshPatch *)(meshList.GetNext(Pos));

		//--------------------------------------------------------------------------------------
		//	Detecting sharp edges
		for(PosNode=(mesh->GetNodeList()).GetHeadPosition();PosNode!=NULL;) {
			node=(QMeshNode *)((mesh->GetNodeList()).GetNext(PosNode));
			node->SetAttribFlag(7,false);
		}
		for(PosEdge=(mesh->GetEdgeList()).GetHeadPosition();PosEdge!=NULL;) {
			edge=(QMeshEdge*)((mesh->GetEdgeList()).GetNext(PosEdge));
			edge->SetAttribFlag(1,false);
			if (!(edge->GetLeftFace() && edge->GetRightFace())) continue;
			edge->GetLeftFace()->GetPlaneEquation(xx,yy,zz,dd);
			edge->GetRightFace()->GetPlaneEquation(nx,ny,nz,dd);
			dd=nx*xx+ny*yy+nz*zz;
			if (dd<SHARPEDGE_NORMAL_THRESHOLD) {
				edge->SetAttribFlag(1,true);
				edge->GetStartPoint()->SetAttribFlag(7,true);
				edge->GetEndPoint()->SetAttribFlag(7,true);
			}
		}

		//--------------------------------------------------------------------------------------
		//	Computing vertex normals in each triangle or polygon
		//--------------------------------------------------------------------------------------
		normalx=(double **)new long[mesh->GetFaceNumber()];
		normaly=(double **)new long[mesh->GetFaceNumber()];
		normalz=(double **)new long[mesh->GetFaceNumber()];
		i=0;
		for(PosFace=(mesh->GetFaceList()).GetHeadPosition();PosFace!=NULL;i++)	{
			face=(QMeshFace *)((mesh->GetFaceList()).GetNext(PosFace));
			face->SetIndexNo(i+1);
			num=face->GetEdgeNum();
			normalx[i]=new double[num];	normaly[i]=new double[num];	normalz[i]=new double[num];
			for(k=0;k<num;k++) {face->GetPlaneEquation(normalx[i][k],normaly[i][k],normalz[i][k],dd);}
		}
		//--------------------------------------------------------------------------------------
		//	Processing for vertex normals
		for(PosNode=(mesh->GetNodeList()).GetHeadPosition();PosNode!=NULL;) {
			node=(QMeshNode *)((mesh->GetNodeList()).GetNext(PosNode));
			if (node->GetAttribFlag(7)) {
				if (_compLocalNormalForSharpVertex(node,normalx,normaly,normalz)) continue;
			}
			node->CalNormal(normal);
			for(PosFace=node->GetFaceList().GetHeadPosition();PosFace!=NULL;) {
				face=(QMeshFace *)(node->GetFaceList().GetNext(PosFace));
				i=face->GetIndexNo()-1;
				num=face->GetEdgeNum();
				for(k=0;k<num;k++) {
					if (face->GetNodeRecordPtr(k)!=node) continue;
					normalx[i][k]=normal[0];
					normaly[i][k]=normal[1];
					normalz[i][k]=normal[2];
				}
			}
		}

		//--------------------------------------------------------------------------------------
		//	Displaying the polygons
		glBegin(GL_TRIANGLES);
		id=0;
		for(PosFace=(mesh->GetFaceList()).GetHeadPosition();PosFace!=NULL;id++)
		{
			face=(QMeshFace *)((mesh->GetFaceList()).GetNext(PosFace));

			if (meshIndex==0) 
				{rr=gg=bb=0.6f;	}
			else 
				{rr=gg=bb=0.8f;	rr=0.0f;}
			glColor3f(rr,gg,bb);

			num=face->GetEdgeNum();
			for(k=0;k<num-2;k++) {
				for(i=0;i<3;i++) {
					if (i<1) {
						node=face->GetNodeRecordPtr(i);
						xx=normalx[id][i];	yy=normaly[id][i];	zz=normalz[id][i];
					}
					else {
						node=face->GetNodeRecordPtr(i+k);
						xx=normalx[id][i+k];yy=normaly[id][i+k];zz=normalz[id][i+k];
					}
					glNormal3d(xx,yy,zz);
					node->GetCoord3D(xx,yy,zz);
					glVertex3d(xx,yy,zz);
				}
			}
		}
		glEnd();

		//--------------------------------------------------------------------------------------
		//	Free the memory of normal array
		num=mesh->GetFaceNumber();
		for(i=0;i<num;i++) {delete (double*)normalx[i];	delete (double*)normaly[i];	delete (double*)normalz[i];}
		delete []normalx;	delete []normaly;	delete []normalz;
		for(PosNode=(mesh->GetNodeList()).GetHeadPosition();PosNode!=NULL;) {
			node=(QMeshNode *)((mesh->GetNodeList()).GetNext(PosNode));
			node->SetAttribFlag(7,false);
		}
		for(PosEdge=(mesh->GetEdgeList()).GetHeadPosition();PosEdge!=NULL;) {
			edge=(QMeshEdge*)((mesh->GetEdgeList()).GetNext(PosEdge));
			edge->SetAttribFlag(1,false);
		}
	}

	for(Pos=meshList.GetHeadPosition();Pos!=NULL;meshIndex++) {
		mesh=(QMeshPatch *)(meshList.GetNext(Pos));

		for(PosFace=mesh->GetNodeList().GetHeadPosition();PosFace!=NULL;) {
			QMeshNode *node=(QMeshNode *)(mesh->GetNodeList().GetNext(PosFace));
			if (!(node->GetAttribFlag(7))) continue;
			node->GetCoord3D(xx,yy,zz);
			glColor3f(1,0,0);
			drawBox(xx,yy,zz,.1f);
		}
	}

	glEndList();
}

void PMBody::_buildDrawMeshList()
{
	GLKPOSITION Pos;
	GLKPOSITION PosEdge;

	if (meshList.GetCount()==0) return;

	glNewList(m_drawListID+1, GL_COMPILE);
	glDisable(GL_LIGHTING);
	glLineWidth(1.0);

	glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
	glColor3f(0.0,0.0,0.0);
	glEnable(GL_POLYGON_OFFSET_LINE);
	glPolygonOffset(1.0,1.0);

	QMeshEdge *edge;
	QMeshPatch *mesh;
	double xx,yy,zz;

	glBegin(GL_LINES);
	int meshIndex=0;
	for(Pos=meshList.GetHeadPosition();Pos!=NULL;meshIndex++) {
		mesh=(QMeshPatch *)(meshList.GetNext(Pos));
//		if (meshIndex!=0) continue;
		if ((activePatch!=0) && (meshIndex!=activePatch-1)) continue;
		if ((activeMaterialRegion!=0) && (mesh->GetMaterial(false)!=activeMaterialRegion)
			 && (mesh->GetMaterial(true)!=activeMaterialRegion)) continue;

		for(PosEdge=(mesh->GetEdgeList()).GetHeadPosition();PosEdge!=NULL;) {
			edge=(QMeshEdge *)((mesh->GetEdgeList()).GetNext(PosEdge));

			QMeshFace *leftFace=edge->GetLeftFace();
			QMeshFace *rightFace=edge->GetRightFace();

			edge->GetStartPoint()->GetCoord3D(xx,yy,zz);
			glVertex3d(xx,yy,zz);
			edge->GetEndPoint()->GetCoord3D(xx,yy,zz);
			glVertex3d(xx,yy,zz);
		}
	}
	glEnd();

	glEndList();
}

void PMBody::_buildDrawProfileList()
{
	GLKPOSITION Pos;	
	GLKPOSITION PosEdge;
	double xx,yy,zz;

	if (meshList.GetCount()==0) return;

	glNewList(m_drawListID+2, GL_COMPILE);

	glDisable(GL_LIGHTING);
	glLineWidth(2.0);

	glColor3f(0,0,0);
	glBegin(GL_LINES);	
	int meshIndex=0;
	for(Pos=meshList.GetHeadPosition();Pos!=NULL;meshIndex++) {
		QMeshPatch *mesh=(QMeshPatch *)(meshList.GetNext(Pos));
		if ((activePatch!=0) && (meshIndex!=activePatch-1)) continue;
		if ((activeMaterialRegion!=0) && (mesh->GetMaterial(false)!=activeMaterialRegion)
			 && (mesh->GetMaterial(true)!=activeMaterialRegion)) continue;

		for(PosEdge=(mesh->GetEdgeList()).GetHeadPosition();PosEdge!=NULL;) {
			QMeshEdge *edge=(QMeshEdge *)((mesh->GetEdgeList()).GetNext(PosEdge));
			if ((edge->GetLeftFace()==NULL || edge->GetRightFace()==NULL)//) {
				|| (edge->GetLeftFace()->m_nIdentifiedPatchIndex!=edge->GetRightFace()->m_nIdentifiedPatchIndex)) {
				edge->GetStartPoint()->GetCoord3D(xx,yy,zz);
				glVertex3d(xx,yy,zz);
				edge->GetEndPoint()->GetCoord3D(xx,yy,zz);
				glVertex3d(xx,yy,zz);
			}
		}
	}
	glEnd();

	glLineWidth(4.0);
	glColor3f(0,0,1);
	glBegin(GL_LINES);	
	for(Pos=meshList.GetHeadPosition();Pos!=NULL;) {
		QMeshPatch *mesh=(QMeshPatch *)(meshList.GetNext(Pos));

		for(PosEdge=(mesh->GetEdgeList()).GetHeadPosition();PosEdge!=NULL;) {
			QMeshEdge *edge=(QMeshEdge *)((mesh->GetEdgeList()).GetNext(PosEdge));
			if (!(edge->GetAttribFlag(1))) continue;

			edge->GetStartPoint()->GetCoord3D(xx,yy,zz);
			glVertex3d(xx,yy,zz);
			edge->GetEndPoint()->GetCoord3D(xx,yy,zz);
			glVertex3d(xx,yy,zz);
		}
	}
	glEnd();

	glLineWidth(1.0);

	glEndList();
}

void PMBody::_buildDrawPreMeshList()
{
	GLKPOSITION Pos;
	GLKPOSITION PosFace;
	QMeshFace *face;
	QMeshNode *node;
	QMeshPatch *mesh;
	int i,num;	double xx,yy,zz,dd;

	if (meshList.GetCount()==0) return;

	glNewList(m_drawListID+3, GL_COMPILE);

	glColor3f(1.0f,1.0f,1.0f);

	glEnable(GL_NORMALIZE);

	for(Pos=meshList.GetHeadPosition();Pos!=NULL;) {
		mesh=(QMeshPatch *)(meshList.GetNext(Pos));
		if ((activeMaterialRegion!=0) && (mesh->GetMaterial(false)!=activeMaterialRegion)
			 && (mesh->GetMaterial(true)!=activeMaterialRegion)) continue;

		for(PosFace=(mesh->GetFaceList()).GetHeadPosition();PosFace!=NULL;)
		{
			face=(QMeshFace *)((mesh->GetFaceList()).GetNext(PosFace));
//			continue;
//			if (face->GetAttribFlag(4)) glColor3f(1,0,0); else continue;
			num=face->GetEdgeNum();
			glBegin(GL_POLYGON);
			for(i=0;i<num;i++) 
			{
				node=face->GetNodeRecordPtr(i);
				face->GetPlaneEquation(xx,yy,zz,dd);
				glNormal3d(xx,yy,zz);
				node->GetCoord3D(xx,yy,zz);
				glVertex3d(xx,yy,zz);
			}
			glEnd();
		}
	}

	glEndList();
}

void PMBody::drawShade()
{
	if (meshList.IsEmpty()) {glDeleteLists(m_drawListID, 4); m_drawListID=-1; return;}
	glCallList(m_drawListID);
}

void PMBody::drawMesh()
{
	if (meshList.IsEmpty()) {glDeleteLists(m_drawListID, 4); m_drawListID=-1; return;}
	glCallList(m_drawListID+1);
}

void PMBody::drawProfile()
{
	if (meshList.IsEmpty()) {glDeleteLists(m_drawListID, 4); m_drawListID=-1; return;}
	glCallList(m_drawListID+2);
}

void PMBody::drawPreMesh()
{
	return;
	if (meshList.IsEmpty()) {glDeleteLists(m_drawListID, 4); m_drawListID=-1; return;}
	glCallList(m_drawListID+3);
}

void PMBody::drawBox(float xx, float yy, float zz, float r)
{
	glBegin(GL_QUADS);

	glNormal3f(0.0f,0.0f,-1.0f);
	glVertex3f(xx-r,yy-r,zz-r);
	glVertex3f(xx-r,yy+r,zz-r);
	glVertex3f(xx+r,yy+r,zz-r);
	glVertex3f(xx+r,yy-r,zz-r);

	glNormal3f(0.0f,0.0f,1.0f);
	glVertex3f(xx-r,yy-r,zz+r);
	glVertex3f(xx+r,yy-r,zz+r);
	glVertex3f(xx+r,yy+r,zz+r);
	glVertex3f(xx-r,yy+r,zz+r);
		
	glNormal3f(-1.0f,0.0f,0.0f);
	glVertex3f(xx-r,yy-r,zz-r);
	glVertex3f(xx-r,yy-r,zz+r);
	glVertex3f(xx-r,yy+r,zz+r);
	glVertex3f(xx-r,yy+r,zz-r);
		
	glNormal3f(1.0f,0.0f,0.0f);
	glVertex3f(xx+r,yy-r,zz-r);
	glVertex3f(xx+r,yy+r,zz-r);
	glVertex3f(xx+r,yy+r,zz+r);
	glVertex3f(xx+r,yy-r,zz+r);
		
	glNormal3f(0.0f,-1.0f,0.0f);
	glVertex3f(xx-r,yy-r,zz-r);
	glVertex3f(xx+r,yy-r,zz-r);
	glVertex3f(xx+r,yy-r,zz+r);
	glVertex3f(xx-r,yy-r,zz+r);
		
	glNormal3f(0.0f,1.0f,0.0f);
	glVertex3f(xx-r,yy+r,zz-r);
	glVertex3f(xx-r,yy+r,zz+r);
	glVertex3f(xx+r,yy+r,zz+r);
	glVertex3f(xx+r,yy+r,zz-r);

	glEnd();
}

void PMBody::drawHighLight()
{
}

void PMBody::ClearAll()
{
	GLKPOSITION Pos;

	for(Pos=meshList.GetHeadPosition();Pos!=NULL;) {
		QMeshPatch *mesh=(QMeshPatch *)(meshList.GetNext(Pos));
		delete mesh;
	}
	meshList.RemoveAll();
}

void PMBody::computeRange()
{
	double range=0.0,ll,xx,yy,zz;
	GLKPOSITION Pos;
	GLKPOSITION PosNode;

	for(Pos=meshList.GetHeadPosition();Pos!=NULL;) {
		QMeshPatch *mesh=(QMeshPatch *)(meshList.GetNext(Pos));
		for(PosNode=(mesh->GetNodeList()).GetHeadPosition();PosNode!=NULL;) {
			QMeshNode *node=(QMeshNode *)((mesh->GetNodeList()).GetNext(PosNode));

			node->GetCoord3D(xx,yy,zz);
			ll=xx*xx+yy*yy+zz*zz;

			if (ll>range) range=ll;
		}
	}

	m_range=(float)(sqrt(range));
}

bool PMBody::_compLocalNormalForSharpVertex(QMeshNode* node, double **normalx, 
											double **normaly, double **normalz)
{
	GLKPOSITION Pos;
	QMeshEdge* edge; QMeshEdge* currentedge = NULL;
	QMeshEdge* startedge = NULL; QMeshEdge* nextedge = NULL;
	double xx,yy,zz,dd;
	QMeshFace* face, *currentface;  
	double *nx,*ny,*nz;
	int *facenum;
	int num,j,i,edgeNum,id;

	for(Pos = node->GetFaceList().GetHeadPosition(); Pos!=NULL;) {
		face = (QMeshFace*)(node->GetFaceList().GetNext(Pos));
		face->m_nIdentifiedPatchIndex = -1;
	}
	for(Pos = node->GetEdgeList().GetHeadPosition(); Pos!=NULL;) {
		edge = (QMeshEdge*)(node->GetEdgeList().GetNext(Pos));
		if(edge->GetAttribFlag(1)) {currentedge = edge; break;}
	}

	startedge = currentedge;
	num = 0;
	do{
		if(currentedge->GetStartPoint() == node)
			currentface = currentedge->GetLeftFace();
		else
			currentface = currentedge->GetRightFace();

		if(currentface == NULL) return false;
		currentface->m_nIdentifiedPatchIndex = num;

		int edgeNum = currentface->GetEdgeNum();
		nextedge = NULL;
		for(i=0;i<edgeNum;i++) {
			if(currentface->GetEdgeRecordPtr(i)==currentedge) continue;
			if( (currentface->GetEdgeRecordPtr(i)->GetStartPoint()==node) 
				|| (currentface->GetEdgeRecordPtr(i)->GetEndPoint()==node)) 
				{nextedge = currentface->GetEdgeRecordPtr(i); break;}
		}

		if (nextedge == NULL) {
			printf("Warning: The model is wrong in topology!\n");
			return false;
		}
		if(nextedge->GetAttribFlag(1)) num++;
		currentedge = nextedge;
	}while(currentedge!=startedge);


	nx = new double[num];
	ny = new double[num];
	nz = new double[num];
	facenum = new int[num];
	for(i=0;i<num;i++){
		nx[i] = ny[i] = nz[i] = 0;
		facenum[i]=0;
	}
	for(Pos=node->GetFaceList().GetHeadPosition();Pos!=NULL;) {
		face=(QMeshFace*)(node->GetFaceList().GetNext(Pos));
		face->GetPlaneEquation(xx,yy,zz,dd);
		i = face->m_nIdentifiedPatchIndex;
		if (i==-1) {
			printf("Warning:The model is wrong in topology!\n");
			delete []nx; delete []ny; delete []nz; delete []facenum;
			return false;
		}
		nx[i]+=xx;
		ny[i]+=yy;
		nz[i]+=zz;
		facenum[i]=facenum[i]+1;
	}

	for(i=0;i<num;i++) {
		nx[i] = nx[i]/(double)(facenum[i]);
		ny[i] = ny[i]/(double)(facenum[i]);
		nz[i] = nz[i]/(double)(facenum[i]);
	}

	for(Pos = node->GetFaceList().GetHeadPosition(); Pos!=NULL;)
	{
		face = (QMeshFace*)(node->GetFaceList().GetNext(Pos));
		i = face->m_nIdentifiedPatchIndex;
		if (i==-1) continue;
		id=face->GetIndexNo()-1;
		edgeNum=face->GetEdgeNum();
		for(j=0;j<edgeNum;j++) {if (face->GetNodeRecordPtr(j)==node) break;}
		normalx[id][j]=nx[i];
		normaly[id][j]=ny[i];
		normalz[id][j]=nz[i];
	}

	delete []nx; delete []ny; delete []nz; delete []facenum;

	return true;
}

void PMBody::_changeValueToColor(double maxValue, double minValue, double Value, 
								 float & nRed, float & nGreen, float & nBlue)
{
//	Value=fabs(Value);

	if (Value<minValue) 
	{
		nRed=0.0;
		nGreen=0.0;
		nBlue=0.0;
		return;
	}

	if ((maxValue-minValue)<0.000000000001)
	{
		nRed=0.0;
		nGreen=0.0;
		nBlue=1.0;
		return;
	}

	double temp=(Value-minValue)/(maxValue-minValue);

//	nRed=(float)(1.0-temp);	nGreen=(float)(1.0-temp); nBlue=(float)(1.0-temp);	return;

	if (temp>0.75)
	{
		nRed=1;
		nGreen=(float)(1.0-(temp-0.75)/0.25);	
		if (nGreen<0) nGreen=0.0f;
		nBlue=0;
		return;
	}
	if (temp>0.5)
	{
		nRed=(float)((temp-0.5)/0.25);
		nGreen=1;
		nBlue=0;
		return;
	}
	if (temp>0.25)
	{
		nRed=0;
		nGreen=1;
		nBlue=(float)(1.0-(temp-0.25)/0.25);
		return;
	}
	else
	{
		nRed=0;
		nGreen=(float)(temp/0.25);
		nBlue=1;
	}
}
