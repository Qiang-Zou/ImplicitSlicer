#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <fstream>

#include "QMeshPatch.h"
#include "QMeshFace.h"
#include "QMeshEdge.h"
#include "QMeshNode.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

QMeshPatch::QMeshPatch()
{
	indexno=0;
	for(int i=0;i<8;i++) flags[i]=false;
	nodeList.RemoveAll();
	edgeList.RemoveAll();
	faceList.RemoveAll();
	m_materialNegativeDir=1;	m_materialPositiveDir=0;

	int num=faceList.GetCount();
}

QMeshPatch::~QMeshPatch()
{
	ClearAll();
}

//////////////////////////////////////////////////////////////////////
// Implementation
//////////////////////////////////////////////////////////////////////

void QMeshPatch::ClearAll()
{
	GLKPOSITION Pos;

	for(Pos=faceList.GetHeadPosition();Pos!=NULL;) {
		QMeshFace* face=(QMeshFace*)(faceList.GetNext(Pos));
		delete face;
	}
	faceList.RemoveAll();

	for(Pos=edgeList.GetHeadPosition();Pos!=NULL;) {
		QMeshEdge* edge=(QMeshEdge*)(edgeList.GetNext(Pos));
		delete edge;
	}
	edgeList.RemoveAll();

	for(Pos=nodeList.GetHeadPosition();Pos!=NULL;) {
		QMeshNode* node=(QMeshNode*)(nodeList.GetNext(Pos));
		delete node;
	}
	nodeList.RemoveAll();
}

void QMeshPatch::InverseOrientation()
{
	GLKPOSITION Pos;
    QMeshEdge *edgeArray[MAX_EDGE_NUM];
    bool edgeDir[MAX_EDGE_NUM];

	for(Pos=faceList.GetHeadPosition();Pos!=NULL;) {
		QMeshFace *face=(QMeshFace *)(faceList.GetNext(Pos));
		int i,eNum=face->GetEdgeNum();
		for(i=0;i<eNum;i++) {
			edgeArray[eNum-1-i]=face->GetEdgeRecordPtr(i);
			edgeDir[eNum-1-i]=face->IsNormalDirection(i);
		}
		for(i=0;i<eNum;i++) {
			face->SetEdgeRecordPtr(i,edgeArray[i]);
			face->SetDirectionFlag(i,!(edgeDir[i]));
		}
		face->CalPlaneEquation();
	}
	//------------------------------------------------------------------------------
	for(Pos=nodeList.GetHeadPosition();Pos!=NULL;) {
		QMeshNode *node=(QMeshNode *)(nodeList.GetNext(Pos));
		node->GetEdgeList().RemoveAll();
	}
	//------------------------------------------------------------------------------
	for(Pos=edgeList.GetHeadPosition();Pos!=NULL;) {
		QMeshEdge *edge=(QMeshEdge *)(edgeList.GetNext(Pos));
		edge->SetLeftFace(NULL);	edge->SetRightFace(NULL);
		edge->GetStartPoint()->AddEdge(edge);	edge->GetEndPoint()->AddEdge(edge);
	}
	//------------------------------------------------------------------------------
	for(Pos=faceList.GetHeadPosition();Pos!=NULL;) {
		QMeshFace *face=(QMeshFace *)(faceList.GetNext(Pos));
		int i,eNum=face->GetEdgeNum();
		for(i=0;i<eNum;i++) {
			if (face->IsNormalDirection(i))
				face->GetEdgeRecordPtr(i)->SetLeftFace(face);
			else
				face->GetEdgeRecordPtr(i)->SetRightFace(face);
		}
	}
}

bool QMeshPatch::inputOFFFile(char* filename)
{
	FILE *fp;
	char buf[100];
	GLKPOSITION Pos;
	GLKPOSITION PosNode;	
	QMeshNode *node,*startNode,*endNode;
	QMeshEdge *edge;
	QMeshFace *face;
	QMeshNode **nodeArray;
	float xx,yy,zz;
	int faceNum,nodeNum,i;

	fp = fopen(filename, "r");
    if(!fp) {
	    printf("===============================================");
        printf("Can not open the data file, called in OFF File Import!");
	    printf("===============================================");
	    return false;
	}
	ClearAll();

	fscanf(fp,"%s\n",buf);
	fscanf(fp,"%d %d %d\n",&nodeNum,&faceNum,&i);

	for(i=0;i<nodeNum;i++) {
		fscanf(fp,"%f %f %f\n",&xx,&yy,&zz);
		node=new QMeshNode;
		node->SetMeshPatchPtr(this);
		node->SetCoord3D(xx,yy,zz);
		node->SetIndexNo(nodeList.GetCount()+1);
		nodeList.AddTail(node);
	}

	nodeArray=(QMeshNode**)new long[nodeNum];
	i=0;
	for(Pos=nodeList.GetHeadPosition();Pos!=NULL;i++) {
		node=(QMeshNode*)(nodeList.GetNext(Pos));
		nodeArray[i]=node;
	}

	for(i=0;i<faceNum;i++) {
		int num,nodeIndex;
		fscanf(fp,"%d ",&num);
//		num=3;
		if (num>2) {
			face=new QMeshFace;
			face->SetMeshPatchPtr(this);
			face->SetIndexNo(faceList.GetCount()+1);
			faceList.AddTail(face);

			for(int j=0;j<num;j++) {
                if (j==(num-1)) fscanf(fp,"%d\n",&nodeIndex);
                else fscanf(fp,"%d ",&nodeIndex);

				(face->GetAttachedList()).AddTail(nodeArray[nodeIndex]);
			}
		}
	}

	fclose(fp);

	delete [](QMeshNode**)nodeArray;

	//---------------------------------------------------------------------
	//	Build the topology
	//---------------------------------------------------------------------
	//	Step 1: build the edges
	for(Pos=faceList.GetHeadPosition();Pos!=NULL;) {
		face=(QMeshFace*)(faceList.GetNext(Pos));

		int edgeNum=(face->GetAttachedList()).GetCount();
		face->SetEdgeNum(edgeNum);

		nodeArray=(QMeshNode**)new long[edgeNum];
		i=0;
		for(PosNode=(face->GetAttachedList()).GetHeadPosition();PosNode!=NULL;i++) {
			nodeArray[i]=(QMeshNode*)((face->GetAttachedList()).GetNext(PosNode));
			(nodeArray[i]->GetFaceList()).AddTail(face);
		}
		for(i=0;i<edgeNum;i++) {
            edge=NULL;
            startNode=nodeArray[i];
            endNode=nodeArray[(i+1)%edgeNum];
			bool bDir;
			for(PosNode=(startNode->GetEdgeList()).GetHeadPosition();PosNode!=NULL;) {
				QMeshEdge *temp=(QMeshEdge *)((startNode->GetEdgeList()).GetNext(PosNode));
				if ((temp->GetStartPoint()==startNode) && (temp->GetEndPoint()==endNode)) {
					edge=temp;	bDir=true;
				}
				else if ((temp->GetStartPoint()==endNode) && (temp->GetEndPoint()==startNode)) {
					edge=temp;	bDir=false;
				}
			}
			if (edge) {
				face->SetEdgeRecordPtr(i,edge);
				if (bDir) {
					face->SetDirectionFlag(i,true);
					edge->SetLeftFace(face);
				}
				else {
					face->SetDirectionFlag(i,false);
					edge->SetRightFace(face);
				}
			}
			else {
				edge=new QMeshEdge;
				edge->SetMeshPatchPtr(this);
				edge->SetStartPoint(startNode);
				edge->SetEndPoint(endNode);
				edge->SetIndexNo(edgeList.GetCount()+1);
				edgeList.AddTail(edge);

				edge->SetLeftFace(face);
				face->SetEdgeRecordPtr(i,edge);
				face->SetDirectionFlag(i,true);
				(startNode->GetEdgeList()).AddTail(edge);
				(endNode->GetEdgeList()).AddTail(edge);
			}
		}

		delete [](QMeshNode**)nodeArray;

		face->GetAttachedList().RemoveAll();
	}
	//---------------------------------------------------------------------
	//	Step 2: compute the normal
	for(Pos=faceList.GetHeadPosition();Pos!=NULL;) {
		face=(QMeshFace*)(faceList.GetNext(Pos));
		face->CalPlaneEquation();
	}

	return true;
}

bool QMeshPatch::inputMFile(char* filename)
{
	FILE *fp;
	char linebuf[256],buf[100];
	GLKPOSITION Pos;
	GLKPOSITION PosNode;	
	int i,index,index1,index2,index3;
	QMeshNode *node,*startNode,*endNode;
	QMeshEdge *edge;
	QMeshFace *face;
	QMeshNode **nodeArray;
	float xx,yy,zz;
//	float minX,maxX,minY,maxY,minZ,maxZ;

	fp = fopen(filename, "r");
    if(!fp) {
	    printf("===============================================\n");
	    printf("Can not open the data file - OBJ File Import!\n");
	    printf("===============================================\n");
	    return false;
	}

	ClearAll();
	while(!feof(fp)) {
		sprintf(buf,"");
		sprintf(linebuf,"");
		fgets(linebuf, 255, fp);
		sscanf(linebuf,"%s",buf);
	
		if (strcmp(buf,"Vertex")==0 )
		{
			sscanf(linebuf, "%s %d %f %f %f \n", buf, &index, &xx, &yy, &zz);
//			xx=xx*100.0f;	yy=yy*100.0f;	zz=zz*100.0f;

			node=new QMeshNode;
			node->SetMeshPatchPtr(this);
			node->SetCoord3D(xx,yy,zz);
			node->SetIndexNo(nodeList.GetCount()+1);
			nodeList.AddTail(node);
		}
	}
	fclose(fp);

	int nodeNum=nodeList.GetCount();
	nodeArray=(QMeshNode**)new long[nodeNum];
	i=0;
	for(Pos=nodeList.GetHeadPosition();Pos!=NULL;i++) {
		node=(QMeshNode*)(nodeList.GetNext(Pos));
		nodeArray[i]=node;
	}

	fp = fopen(filename, "r");
	while(!feof(fp)) {
		sprintf(buf,"");
		sprintf(linebuf,"");
		fgets(linebuf, 255, fp);
		sscanf(linebuf,"%s",buf);
		
		if (strcmp(buf,"Face")==0 )
		{
			sscanf(linebuf, "%s %d %d %d %d \n", buf, &index, &index1, &index2, &index3);

			face=new QMeshFace;
			face->SetMeshPatchPtr(this);
			face->SetIndexNo(faceList.GetCount()+1);
			faceList.AddTail(face);
			
			(face->GetAttachedList()).AddTail(nodeArray[index1-1]);
			(face->GetAttachedList()).AddTail(nodeArray[index2-1]);
			(face->GetAttachedList()).AddTail(nodeArray[index3-1]);
		}
	}
	fclose(fp);

	delete [](QMeshNode**)nodeArray;

	//---------------------------------------------------------------------
	//	Build the topology
	//---------------------------------------------------------------------
	//	Step 1: build the edges
	for(Pos=faceList.GetHeadPosition();Pos!=NULL;) {
		face=(QMeshFace*)(faceList.GetNext(Pos));

		int edgeNum=(face->GetAttachedList()).GetCount();
		face->SetEdgeNum(edgeNum);

		nodeArray=(QMeshNode**)new long[edgeNum];
		i=0;
		for(PosNode=(face->GetAttachedList()).GetHeadPosition();PosNode!=NULL;i++) {
			nodeArray[i]=(QMeshNode*)((face->GetAttachedList()).GetNext(PosNode));
			(nodeArray[i]->GetFaceList()).AddTail(face);
		}
		for(i=0;i<edgeNum;i++) {
			edge=NULL;	startNode=nodeArray[i];	endNode=nodeArray[(i+1)%edgeNum];
			bool bDir;
			for(PosNode=(startNode->GetEdgeList()).GetHeadPosition();PosNode!=NULL;) {
				QMeshEdge *temp=(QMeshEdge *)((startNode->GetEdgeList()).GetNext(PosNode));
				if ((temp->GetStartPoint()==startNode) && (temp->GetEndPoint()==endNode)) {
					edge=temp;	bDir=true;
				}
				else if ((temp->GetStartPoint()==endNode) && (temp->GetEndPoint()==startNode)) {
					edge=temp;	bDir=false;
				}
			}
			if (edge) {
				face->SetEdgeRecordPtr(i,edge);
				if (bDir) {
					face->SetDirectionFlag(i,true);
					edge->SetLeftFace(face);
				}
				else {
					face->SetDirectionFlag(i,false);
					edge->SetRightFace(face);
				}
			}
			else {
				edge=new QMeshEdge;
				edge->SetMeshPatchPtr(this);
				edge->SetStartPoint(startNode);
				edge->SetEndPoint(endNode);
				edge->SetIndexNo(edgeList.GetCount()+1);
				edgeList.AddTail(edge);

				edge->SetLeftFace(face);
				face->SetEdgeRecordPtr(i,edge);
				face->SetDirectionFlag(i,true);
				(startNode->GetEdgeList()).AddTail(edge);
				(endNode->GetEdgeList()).AddTail(edge);
			}
		}

		delete [](QMeshNode**)nodeArray;

		face->GetAttachedList().RemoveAll();
	}
	//---------------------------------------------------------------------
	//	Step 2: compute the normal
	for(Pos=faceList.GetHeadPosition();Pos!=NULL;) {
		face=(QMeshFace*)(faceList.GetNext(Pos));
		face->CalPlaneEquation();
	}

	return true;
}

bool QMeshPatch::inputPLY2File(char* filename)
{
	FILE *fp;
	GLKPOSITION Pos;
	GLKPOSITION PosNode;	
	int i,j,nodeNum,faceNum,index,edgeNum;
	QMeshNode *node,*startNode,*endNode;
	QMeshEdge *edge;
	QMeshFace *face;
	QMeshNode **nodeArray;
	float xx,yy,zz;
//	float minX,maxX,minY,maxY,minZ,maxZ;

	fp = fopen(filename, "r");
    if(!fp) {
	    printf("===============================================\n");
	    printf("Can not open the data file - PLY2 File Import!\n");
	    printf("===============================================\n");
	    return false;
	}

	ClearAll();
	fscanf(fp,"%d\n",&nodeNum);
	fscanf(fp,"%d\n",&faceNum);

	nodeArray=(QMeshNode**)new long[nodeNum];
	for(i=0;i<nodeNum;i++) {
		fscanf(fp,"%f %f %f \n", &xx, &yy, &zz);

		node=new QMeshNode;
		node->SetMeshPatchPtr(this);
		node->SetCoord3D(xx,yy,zz);
		node->SetIndexNo(nodeList.GetCount()+1);
		nodeList.AddTail(node);
		nodeArray[i]=node;
	}

	for(i=0;i<faceNum;i++) {
		fscanf(fp,"%d ",&edgeNum);

		face=new QMeshFace;
		face->SetMeshPatchPtr(this);
		face->SetIndexNo(faceList.GetCount()+1);
		faceList.AddTail(face);

		for(j=0;j<edgeNum;j++) {
			fscanf(fp,"%d ",&index);
			(face->GetAttachedList()).AddTail(nodeArray[index]);
		}
	}

	fclose(fp);

	delete [](QMeshNode**)nodeArray;

	//---------------------------------------------------------------------
	//	Build the topology
	//---------------------------------------------------------------------
	//	Step 1: build the edges
	for(Pos=faceList.GetHeadPosition();Pos!=NULL;) {
		face=(QMeshFace*)(faceList.GetNext(Pos));

		int edgeNum=(face->GetAttachedList()).GetCount();
		face->SetEdgeNum(edgeNum);

		nodeArray=(QMeshNode**)new long[edgeNum];
		i=0;
		for(PosNode=(face->GetAttachedList()).GetHeadPosition();PosNode!=NULL;i++) {
			nodeArray[i]=(QMeshNode*)((face->GetAttachedList()).GetNext(PosNode));
			(nodeArray[i]->GetFaceList()).AddTail(face);
		}
		for(i=0;i<edgeNum;i++) {
			edge=NULL;	startNode=nodeArray[i];	endNode=nodeArray[(i+1)%edgeNum];
			bool bDir;
			for(PosNode=(startNode->GetEdgeList()).GetHeadPosition();PosNode!=NULL;) {
				QMeshEdge *temp=(QMeshEdge *)((startNode->GetEdgeList()).GetNext(PosNode));
				if ((temp->GetStartPoint()==startNode) && (temp->GetEndPoint()==endNode)) {
					edge=temp;	bDir=true;
				}
				else if ((temp->GetStartPoint()==endNode) && (temp->GetEndPoint()==startNode)) {
					edge=temp;	bDir=false;
				}
			}
			if (edge) {
				face->SetEdgeRecordPtr(i,edge);
				if (bDir) {
					face->SetDirectionFlag(i,true);
					edge->SetLeftFace(face);
				}
				else {
					face->SetDirectionFlag(i,false);
					edge->SetRightFace(face);
				}
			}
			else {
				edge=new QMeshEdge;
				edge->SetMeshPatchPtr(this);
				edge->SetStartPoint(startNode);
				edge->SetEndPoint(endNode);
				edge->SetIndexNo(edgeList.GetCount()+1);
				edgeList.AddTail(edge);

				edge->SetLeftFace(face);
				face->SetEdgeRecordPtr(i,edge);
				face->SetDirectionFlag(i,true);
				(startNode->GetEdgeList()).AddTail(edge);
				(endNode->GetEdgeList()).AddTail(edge);
			}
		}

		delete [](QMeshNode**)nodeArray;

		face->GetAttachedList().RemoveAll();
	}
	//---------------------------------------------------------------------
	//	Step 2: compute the normal
	for(Pos=faceList.GetHeadPosition();Pos!=NULL;) {
		face=(QMeshFace*)(faceList.GetNext(Pos));
		face->CalPlaneEquation();
	}
	for(Pos=edgeList.GetHeadPosition();Pos!=NULL;) {
		edge=(QMeshEdge*)(edgeList.GetNext(Pos));
		if ((edge->GetLeftFace()) && (edge->GetRightFace())) {
			edge->SetAttribFlag(0,false);
		}
		else {
			edge->SetAttribFlag(0,true);
			edge->GetStartPoint()->SetAttribFlag(0,true);
			edge->GetEndPoint()->SetAttribFlag(0,true);
		}
	}

	return true;
}

void QMeshPatch::constructionFromVerFaceTable(int nodeNum, float *nodeTable, int faceNum, unsigned int* faceTable)
{
	QMeshNode *node,*startNode,*endNode;
	QMeshEdge *edge;
	QMeshFace *face;
	QMeshNode **nodeArray;
	GLKPOSITION Pos;
	GLKPOSITION PosNode;
	int i;

	nodeArray=(QMeshNode**)new long[nodeNum];
	for(i=0;i<nodeNum;i++) {
		node=new QMeshNode;		nodeArray[i]=node;
		node->SetMeshPatchPtr(this);
		node->SetCoord3D(nodeTable[i*3],nodeTable[i*3+1],nodeTable[i*3+2]);
		node->SetCoord3D_last(nodeTable[i*3],nodeTable[i*3+1],nodeTable[i*3+2]);
		node->SetIndexNo(nodeList.GetCount()+1);
//		node->SetAttribFlag(4);
		nodeList.AddTail(node);
	}
//	delete [](QMeshNode**)nodeArray;	return;
	//--------------------------------------------------------------------------------------------------------
	for(i=0;i<faceNum;i++) {
		face=new QMeshFace;
		face->SetMeshPatchPtr(this);
		face->SetIndexNo(faceList.GetCount()+1);
		faceList.AddTail(face);

		(face->GetAttachedList()).AddTail(nodeArray[faceTable[i*4+0]-1]);	//printf("%d ",faceTable[i*4+0]-1);
		(face->GetAttachedList()).AddTail(nodeArray[faceTable[i*4+1]-1]);	//printf("%d ",faceTable[i*4+1]-1);
		(face->GetAttachedList()).AddTail(nodeArray[faceTable[i*4+2]-1]);	//printf("%d ",faceTable[i*4+2]-1);
		if (faceTable[i*4+3]==0) continue;
		(face->GetAttachedList()).AddTail(nodeArray[faceTable[i*4+3]-1]);	//printf("%d ",faceTable[i*4+3]-1);
	}
	delete [](QMeshNode**)nodeArray;	

	//--------------------------------------------------------------------------------------------------------
	//	Build the topology
	//--------------------------------------------------------------------------------------------------------
	//	Step 1: build the edges
	int faceIndex=0;
	for(Pos=faceList.GetHeadPosition();Pos!=NULL;faceIndex++) {
		face=(QMeshFace*)(faceList.GetNext(Pos));

		int edgeNum=(face->GetAttachedList()).GetCount();
		face->SetEdgeNum(edgeNum);

		nodeArray=(QMeshNode**)new long[edgeNum];
		i=0;
		for(PosNode=(face->GetAttachedList()).GetHeadPosition();PosNode!=NULL;i++) {
			nodeArray[i]=(QMeshNode*)((face->GetAttachedList()).GetNext(PosNode));
			(nodeArray[i]->GetFaceList()).AddTail(face);
		}

		for(i=0;i<edgeNum;i++) {
			edge=NULL;	startNode=nodeArray[i];	endNode=nodeArray[(i+1)%edgeNum];
			bool bDir;
			for(PosNode=(startNode->GetEdgeList()).GetHeadPosition();PosNode!=NULL;) {
				QMeshEdge *temp=(QMeshEdge *)((startNode->GetEdgeList()).GetNext(PosNode));
				if ((temp->GetStartPoint()==startNode) && (temp->GetEndPoint()==endNode) && (temp->GetLeftFace()==NULL)) {
					edge=temp;	bDir=true;
				}
				else if ((temp->GetStartPoint()==endNode) && (temp->GetEndPoint()==startNode) && (temp->GetRightFace()==NULL)) {
					edge=temp;	bDir=false;
				}
			}
			if (edge && bDir) {
				face->SetEdgeRecordPtr(i,edge);
				face->SetDirectionFlag(i,true);
				edge->SetLeftFace(face);
			}
			else if (edge && (!bDir)) {
				face->SetEdgeRecordPtr(i,edge);
				face->SetDirectionFlag(i,false);
				edge->SetRightFace(face);
			}
			else {
				edge=new QMeshEdge;
				edge->SetMeshPatchPtr(this);
				edge->SetStartPoint(startNode);
				edge->SetEndPoint(endNode);
				edge->SetIndexNo(edgeList.GetCount()+1);
				edgeList.AddTail(edge);

				edge->SetLeftFace(face);
				face->SetEdgeRecordPtr(i,edge);
				face->SetDirectionFlag(i,true);
				(startNode->GetEdgeList()).AddTail(edge);
				(endNode->GetEdgeList()).AddTail(edge);
			}
		}

		delete [](QMeshNode**)nodeArray;
		face->GetAttachedList().RemoveAll();
	}
	//---------------------------------------------------------------------
	//	Step 2: compute the normal
	for(Pos=faceList.GetHeadPosition();Pos!=NULL;) {
		face=(QMeshFace*)(faceList.GetNext(Pos));
		face->CalPlaneEquation();
	}
	for(Pos=edgeList.GetHeadPosition();Pos!=NULL;) {
		edge=(QMeshEdge*)(edgeList.GetNext(Pos));
		if ((edge->GetLeftFace()) && (edge->GetRightFace())) continue;
		edge->SetAttribFlag(0);
		edge->GetStartPoint()->SetAttribFlag(0);
		edge->GetEndPoint()->SetAttribFlag(0);
    }
}

vector<double> QMeshPatch::getBBox()
{
    vector<double> bbox(6, 0);
    if (GetNodeNumber() == 0) return bbox;

    bbox[0] = bbox[1] = bbox[2] = numeric_limits<double>::max();
    bbox[3] = bbox[4] = bbox[5] = numeric_limits<double>::min();
    for (auto iter(this->VertBegin()); iter != this->VertEnd(); iter = iter->Next()) {
        double xx, yy, zz;
        auto vert = GetVertAt(iter);
        vert->GetCoord3D(xx, yy, zz);
        if (xx<bbox[0]) bbox[0] = xx;
        if (xx>bbox[3]) bbox[3] = xx;
        if (yy<bbox[1]) bbox[1] = yy;
        if (yy>bbox[4]) bbox[4] = yy;
        if (zz<bbox[2]) bbox[2] = zz;
        if (zz>bbox[5]) bbox[5] = zz;
    }

    return bbox;
}

double QMeshPatch::CompAverageEdgeLength()
{
    double avgLength = 0;
    size_t div(0);
    for (auto Pos = GetEdgeList().GetHeadPosition(); Pos != NULL;) {
        QMeshEdge *edge = (QMeshEdge *)(GetEdgeList().GetNext(Pos));
        ++div;
        avgLength += edge->CalLength();
    }
    return avgLength / div;
}

bool QMeshPatch::inputOBJFile(char* filename)
{
	FILE *fp;
	char fields[MAX_EDGE_NUM][255];
	char linebuf[256],buf[100];
	GLKPOSITION Pos;
	GLKPOSITION PosNode;	
	int i;
	QMeshNode *node,*startNode,*endNode;
	QMeshEdge *edge;
	QMeshFace *face;
	QMeshNode **nodeArray;
	float xx,yy,zz,ww;
//	float minX,maxX,minY,maxY,minZ,maxZ;

	fp = fopen(filename, "r");
    if(!fp) {
	    printf("===============================================\n");
	    printf("Can not open the data file - OBJ File Import!\n");
	    printf("===============================================\n");
	    return false;
	}

	ClearAll();
	while(!feof(fp)) {
		sprintf(buf,"");
		sprintf(linebuf,"");
		fgets(linebuf, 255, fp);
		sscanf(linebuf,"%s",buf);
	
		if ( (strlen(buf)==1) && (buf[0]=='v') )
		{
            sscanf(linebuf, "%s %f %f %f \n", buf, &xx, &yy, &zz);

			node=new QMeshNode;
			node->SetMeshPatchPtr(this);
			node->SetCoord3D(xx,yy,zz);
			node->SetCoord3D_last(xx,yy,zz);
			node->SetIndexNo(nodeList.GetCount()+1);
			nodeList.AddTail(node);
            node->GetFaceList().RemoveAll();
		}
	}
	fclose(fp);

	int nodeNum=nodeList.GetCount();
	nodeArray=(QMeshNode**)new long[nodeNum];
	i=0;
	for(Pos=nodeList.GetHeadPosition();Pos!=NULL;i++) {
		node=(QMeshNode*)(nodeList.GetNext(Pos));
		nodeArray[i]=node;
	}

	fp = fopen(filename, "r");
	while(!feof(fp)) {
		sprintf(buf,"");
		sprintf(linebuf,"");
		fgets(linebuf, 255, fp);
		sscanf(linebuf,"%s ",buf);
		
		if ( (strlen(buf)==1) && (buf[0]=='f') )
		{
			char seps[]=" \n";
			char seps2[]="/";
			char *token;
			char linebuf2[255];

			strcpy(linebuf2,linebuf);

			int num=0;
			token = strtok( linebuf, seps );
			while(token!=NULL) {token=strtok(NULL,seps); num++;}
			num=num-1;
            //printf("face node-num:%d \n",num);

			if (num>MAX_EDGE_NUM) continue;
			if (num<1) continue;

			face=new QMeshFace;
			face->SetMeshPatchPtr(this);
			face->SetIndexNo(faceList.GetCount()+1);
			faceList.AddTail(face);

			token = strtok( linebuf2, seps );	
			for(i=0;i<num;i++) {
				token = strtok( NULL, seps );
				strcpy(fields[i],token);
			}

			bool bValid=true;
			for(i=0;i<num;i++) {
				token = strtok( fields[i], seps2 );
				int nodeIndex=atoi(token);
                if (nodeIndex<=0) continue; // to protect in case if there is something wrong

				(face->GetAttachedList()).AddTail(nodeArray[nodeIndex-1]);
			}
			if (!bValid) {delete face; faceList.RemoveTail(); continue;}

			bool bDegenerated=false;
			for(Pos=face->GetAttachedList().GetHeadPosition();Pos!=NULL;) {
				QMeshNode *pNode=(QMeshNode *)(face->GetAttachedList().GetNext(Pos));
				GLKPOSITION Pos2=Pos;
				for(;Pos2!=NULL;) {
					QMeshNode *qNode=(QMeshNode *)(face->GetAttachedList().GetNext(Pos2));
					if ((pNode==qNode)) {
						bDegenerated=true;
						break;
					}
				}
				if (bDegenerated) break;
			}
			if (bDegenerated) {
				faceList.RemoveTail();
				delete face;
			}
		}
	}
	fclose(fp);

	delete [](QMeshNode**)nodeArray;

	//---------------------------------------------------------------------
	//	Build the topology
	//---------------------------------------------------------------------
	//	Step 1: build the edges
	for(Pos=faceList.GetHeadPosition();Pos!=NULL;) {
		face=(QMeshFace*)(faceList.GetNext(Pos));

		int edgeNum=(face->GetAttachedList()).GetCount();
		face->SetEdgeNum(edgeNum);

		nodeArray=(QMeshNode**)new long[edgeNum];
		i=0;
		for(PosNode=(face->GetAttachedList()).GetHeadPosition();PosNode!=NULL;i++) {
			nodeArray[i]=(QMeshNode*)((face->GetAttachedList()).GetNext(PosNode));
			(nodeArray[i]->GetFaceList()).AddTail(face);
		}

		for(i=0;i<edgeNum;i++) {
			edge=NULL;	startNode=nodeArray[i];	endNode=nodeArray[(i+1)%edgeNum];
			bool bDir;
			for(PosNode=(startNode->GetEdgeList()).GetHeadPosition();PosNode!=NULL;) {
				QMeshEdge *temp=(QMeshEdge *)((startNode->GetEdgeList()).GetNext(PosNode));
				if ((temp->GetStartPoint()==startNode) && (temp->GetEndPoint()==endNode) && (temp->GetLeftFace()==NULL)) {
					edge=temp;	bDir=true;
				}
				else if ((temp->GetStartPoint()==endNode) && (temp->GetEndPoint()==startNode) && (temp->GetRightFace()==NULL)) {
					edge=temp;	bDir=false;
				}
			}
			if (edge && bDir) {
				face->SetEdgeRecordPtr(i,edge);
				face->SetDirectionFlag(i,true);
				edge->SetLeftFace(face);
			}
			else if (edge && (!bDir)) {
				face->SetEdgeRecordPtr(i,edge);
				face->SetDirectionFlag(i,false);
				edge->SetRightFace(face);
			}
			else {
				edge=new QMeshEdge;
				edge->SetMeshPatchPtr(this);
				edge->SetStartPoint(startNode);
				edge->SetEndPoint(endNode);
				edge->SetIndexNo(edgeList.GetCount()+1);
				edgeList.AddTail(edge);

				edge->SetLeftFace(face);
				face->SetEdgeRecordPtr(i,edge);
				face->SetDirectionFlag(i,true);
				(startNode->GetEdgeList()).AddTail(edge);
				(endNode->GetEdgeList()).AddTail(edge);
			}
		}

		delete [](QMeshNode**)nodeArray;
		face->GetAttachedList().RemoveAll();
	}
	//---------------------------------------------------------------------
	//	Step 2: compute the normal
	for(Pos=faceList.GetHeadPosition();Pos!=NULL;) {
		face=(QMeshFace*)(faceList.GetNext(Pos));
		face->CalPlaneEquation();
	}
	for(Pos=edgeList.GetHeadPosition();Pos!=NULL;) {
		edge=(QMeshEdge*)(edgeList.GetNext(Pos));
		if ((edge->GetLeftFace()) && (edge->GetRightFace())) continue;
		edge->SetAttribFlag(0);
		edge->GetStartPoint()->SetAttribFlag(0);
		edge->GetEndPoint()->SetAttribFlag(0);
	}

	return true;
}

void QMeshPatch::outputTrglOBJFile(char* filename)
{
	FILE *fp;
	GLKPOSITION Pos;
	QMeshNode *node;
	QMeshFace *face;
	double xx,yy,zz;
	int i,num,index;

	fp = fopen(filename, "w");
    if(!fp)
	{
		printf("===============================================\n");
	    printf("Can not open the data file - OBJ File Export!\n");
		printf("===============================================\n");
	    return;
	}

	fprintf(fp,"# The units used in this file are meters.\n");
	
	i=1;
	for(Pos=nodeList.GetHeadPosition();Pos!=NULL;i++) {
		node=(QMeshNode *)(nodeList.GetNext(Pos));
		node->GetCoord3D(xx,yy,zz);
		node->SetIndexNo(i);
//		fprintf(fp,"v %.5f %.5f %.5f\n",(float)yy,(float)zz,(float)xx);
		fprintf(fp,"v %.5f %.5f %.5f\n",(float)xx,(float)yy,(float)zz);
//		fprintf(fp,"v %.12f %.12f %.12f\n",(float)zz,(float)xx,(float)yy);
	}

	fprintf(fp,"\n");
	
	for(Pos=faceList.GetHeadPosition();Pos!=NULL;) {
		face=(QMeshFace *)(faceList.GetNext(Pos));
		num=face->GetEdgeNum();
		
		fprintf(fp,"f ");
		index=face->GetNodeRecordPtr(0)->GetIndexNo();
		fprintf(fp,"%d ",index);
		index=face->GetNodeRecordPtr(1)->GetIndexNo();
		fprintf(fp,"%d ",index);
		index=face->GetNodeRecordPtr(2)->GetIndexNo();
		fprintf(fp,"%d ",index);
		fprintf(fp,"\n");

		for(i=3;i<num;i++) {
			fprintf(fp,"f ");
			index=face->GetNodeRecordPtr(0)->GetIndexNo();
			fprintf(fp,"%d ",index);
			index=face->GetNodeRecordPtr(i-1)->GetIndexNo();
			fprintf(fp,"%d ",index);
			index=face->GetNodeRecordPtr(i)->GetIndexNo();
			fprintf(fp,"%d ",index);
			fprintf(fp,"\n");
		}
	}

    fclose(fp);
}

void QMeshPatch::outputOFFFile(char *filename)
{
    double x, y, z;

    // do writing
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        printf("===============================================\n");
        printf("Can not open the data file - OFF File Export!\n");
        printf("===============================================\n");
        return;
    }
    // 1-2 line: node numbers edge numbers
    outfile << "OFF" << endl;
    outfile << GetNodeNumber() << " " << GetFaceNumber() << " " << 0;
    // vertices
    int i(0);
    for (auto iter = VertBegin(); iter != VertEnd(); iter = iter->Next()) {
        auto node = GetVertAt(iter);
        node->GetCoord3D(x, y, z);
        node->SetIndexNo(i++);
        outfile << endl << x << " " << y << " " << z;
    }

    // faces (indices starting at 0)
    int num;
    auto t = faceList;
    for (auto iter = FaceBegin(); iter != FaceEnd(); iter = iter->Next()) {
        auto face = (QMeshFace*)(iter->data);
        num = face->GetEdgeNum();
        outfile << endl << num;
        for(i=0;i<num;i++) outfile << " " << face->GetNodeRecordPtr(i)->GetIndexNo();
    }
    outfile.close();
}

void QMeshPatch::outputOBJFile(char* filename, bool b2D)
{
	FILE *fp;
	GLKPOSITION Pos;
	QMeshNode *node;
	QMeshFace *face;
	double xx,yy,zz;
	int i,num,index;

	fp = fopen(filename, "w");
    if(!fp)
	{
		printf("===============================================\n");
	    printf("Can not open the data file - OBJ File Export!\n");
		printf("===============================================\n");
	    return;
	}

	fprintf(fp,"# The units used in this file are meters.\n");
	
	i=1;
	for(Pos=nodeList.GetHeadPosition();Pos!=NULL;i++) {
		node=(QMeshNode *)(nodeList.GetNext(Pos));
		node->GetCoord3D(xx,yy,zz);
		node->SetIndexNo(i);
		if (b2D) {
			node->GetCoord2D(xx,yy);
			fprintf(fp,"v %.5f %.5f 0\n",(float)xx,(float)yy); continue;
		}
//		fprintf(fp,"v %.5f %.5f %.5f\n",(float)yy,(float)zz,(float)xx);
		fprintf(fp,"v %.12f %.12f %.12f\n",xx,yy,zz);
	}

	fprintf(fp,"\n");
	
	for(Pos=faceList.GetHeadPosition();Pos!=NULL;) {
		face=(QMeshFace *)(faceList.GetNext(Pos));
		num=face->GetEdgeNum();
		fprintf(fp,"f ");
		for(i=0;i<num;i++) {
			index=face->GetNodeRecordPtr(i)->GetIndexNo();
			fprintf(fp,"%d ",index);
		}
		fprintf(fp,"\n");
	}

	fclose(fp);
}

bool QMeshPatch::GetAttribFlag( const int whichBit )
{
	return flags[whichBit];
}

void QMeshPatch::SetAttribFlag( const int whichBit, const bool toBe )
{
	flags[whichBit]=toBe;
}

int QMeshPatch::GetIndexNo() //from 1 to n
{
	return indexno;
}

void QMeshPatch::SetIndexNo( const int _index )
{
	indexno=_index;
}

int QMeshPatch::GetFaceNumber()
{
	return faceList.GetCount();	//from 1 to n
}

QMeshFace* QMeshPatch::GetFaceRecordPtr(int No)	//from 1 to n
{
	if( (No < 1) || (No > faceList.GetCount()))    return  NULL;
    return (QMeshFace *)faceList.GetAt(faceList.FindIndex(No-1));
}

GLKObList& QMeshPatch::GetFaceList()
{
	return faceList;
}

int QMeshPatch::GetEdgeNumber()	//from 1 to n
{
	return edgeList.GetCount();
}

QMeshEdge* QMeshPatch::GetEdgeRecordPtr(int No)	//from 1 to n
{
	if( (No < 1) || (No > edgeList.GetCount()))    return  NULL;
    return (QMeshEdge *)edgeList.GetAt(edgeList.FindIndex(No-1));
}

GLKObList& QMeshPatch::GetEdgeList() 
{
	return edgeList;
}

int QMeshPatch::GetNodeNumber()	//from 1 to n
{
	return nodeList.GetCount();
}

QMeshNode* QMeshPatch::GetNodeRecordPtr(int No)	//from 1 to n
{
	if( (No < 1) || (No > nodeList.GetCount()))    return  NULL;
    return (QMeshNode *)nodeList.GetAt(nodeList.FindIndex(No-1));
}

GLKObList& QMeshPatch::GetNodeList() 
{
	return nodeList;
}

void QMeshPatch::SetMaterial(bool bDir, int material)
{
	if (bDir)
		m_materialPositiveDir=material;
	else
		m_materialNegativeDir=material;
}

int QMeshPatch::GetMaterial(bool bDir)
{
	if (bDir)
		return m_materialPositiveDir;
	else
		return m_materialNegativeDir;
}

