#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>

#include "QMeshLib/QSurfaceMesh.h"
#include "QMeshLib/QMesh/QMeshPatch.h"
#include "QMeshLib/QMesh/QMeshFace.h"
#include "QMeshLib/QMesh/QMeshEdge.h"
#include "QMeshLib/QMesh/QMeshNode.h"

#include "PMBody.h"
#include "MeshWorksDataBoard.h"

MeshWorksDataBoard::MeshWorksDataBoard(void)
{
	m_polyMeshBody = NULL;
}

MeshWorksDataBoard::~MeshWorksDataBoard(void)
{
	m_bVertexNormalShading = true;
}

void MeshWorksDataBoard::InputOBJFile(char* filename)
{
    QSurfaceMesh *newMesh;

    newMesh = new QSurfaceMesh;
	if (newMesh->inputOBJFile(filename)) {
		m_polyMeshBody->meshList.AddTail(newMesh);
		m_polyMeshBody->computeRange();

		printf("-----------------------------------------------\n");
		printf("Mesh Face Number: %d\n", newMesh->GetFaceNumber());
		printf("Mesh Edge Number: %d\n", newMesh->GetEdgeNumber());
		printf("Mesh Node Number: %d\n", newMesh->GetNodeNumber());
	}
	else {
		delete newMesh;
	}
}

void MeshWorksDataBoard::_mallocIndexArray(int ***&indexArray, int xNum, int yNum, int zNum)
{
	int i, j, k;

	indexArray = (int***)new long[xNum];
	for (i = 0; i<xNum; i++) {
		indexArray[i] = (int**)new long[yNum];
		for (j = 0; j<yNum; j++) {
			indexArray[i][j] = new int[zNum];
			for (k = 0; k<zNum; k++) indexArray[i][j][k] = 0;
		}
	}
}

void MeshWorksDataBoard::_freeIndexArray(int ***&indexArray, int xNum, int yNum, int zNum)
{
	int i, j;

	for (i = 0; i<xNum; i++) {
		for (j = 0; j<yNum; j++) {
			delete[](int*)(indexArray[i][j]);
		}
		delete[](int**)(indexArray[i]);
	}
	delete[](int***)(indexArray);

	indexArray = NULL;
}

void MeshWorksDataBoard::InputOFFFile(char* filename)
{
    QSurfaceMesh *newMesh;

    newMesh = new QSurfaceMesh;
	if (newMesh->inputOFFFile(filename)) {
		m_polyMeshBody->meshList.AddTail(newMesh);
		m_polyMeshBody->computeRange();

		printf("-----------------------------------------------\n");
		printf("Mesh Face Number: %d\n", newMesh->GetFaceNumber());
		printf("Mesh Edge Number: %d\n", newMesh->GetEdgeNumber());
		printf("Mesh Node Number: %d\n", newMesh->GetNodeNumber());
	}
	else {
		delete newMesh;
	}
}

void MeshWorksDataBoard::InputWIOFile(char* filename)
{
	FILE* fp;	char buffer[100];	GLKPOSITION Pos;
	int meshNum, faceNum, nodeNum, edgeNum, i, j;	float xx, yy, zz;

	fp = fopen(filename, "r");
	if (!fp) { printf("Error: Can not open the data file\n"); return; }

	fscanf(fp, "%d\n", &meshNum);
	fscanf(fp, "%s\n", buffer);

	for (i = 0; i<meshNum; i++) {
        QSurfaceMesh *mesh = new QSurfaceMesh;	mesh->SetIndexNo(i);
		m_polyMeshBody->meshList.AddTail(mesh);
		fscanf(fp, "%d %d %d \n", &nodeNum, &edgeNum, &faceNum);

		//-----------------------------------------------------------------------------
		//	loading each QMeshNode
		QMeshNode **nodeArr;
		if (nodeNum>0) nodeArr = (QMeshNode **)new long[nodeNum];
		for (j = 0; j<nodeNum; j++) {
			QMeshNode *node = new QMeshNode;	node->SetMeshPatchPtr(mesh);
			node->SetIndexNo(j);
			nodeArr[j] = node;
			fscanf(fp, "%f %f %f \n", &xx, &yy, &zz);
			node->SetCoord3D(xx, yy, zz);		node->SetAttribFlag(0);
			(mesh->GetNodeList()).AddTail(node);
		}

		//-----------------------------------------------------------------------------
		//	loading each QMeshEdge
		QMeshEdge** edgeArr;
		if (edgeNum>0) edgeArr = (QMeshEdge **)new long[edgeNum];
		for (j = 0; j<edgeNum; j++) {
			QMeshEdge *edge = new QMeshEdge;	edge->SetMeshPatchPtr(mesh);
			int n1, n2;
			fscanf(fp, "%d %d\n", &n1, &n2);
			edgeArr[j] = edge;	edge->SetIndexNo(j);
			edge->SetStartPoint(nodeArr[n1 - 1]);		edge->SetEndPoint(nodeArr[n2 - 1]);
			if (n1>nodeNum || n2>nodeNum) {
				int yu = 0;
				yu++;
			}
			(mesh->GetEdgeList()).AddTail(edge);
		}

		//-----------------------------------------------------------------------------
		//	loading each QMeshFace
		QMeshFace** faceArr;
		if (faceNum>0) faceArr = (QMeshFace **)new long[faceNum];
		for (j = 0; j<faceNum; j++) {
			QMeshFace *face = new QMeshFace;	face->SetMeshPatchPtr(mesh);
			int n1, n2, n3, b1, b2, b3;
			fscanf(fp, "%d %d %d %d %d %d\n", &n1, &n2, &n3, &b1, &b2, &b3);
			faceArr[j] = face;	face->SetIndexNo(j);
			face->SetEdgeNum(3);
			face->SetEdgeRecordPtr(0, edgeArr[n1 - 1]);
			face->SetEdgeRecordPtr(1, edgeArr[n2 - 1]);
			face->SetEdgeRecordPtr(2, edgeArr[n3 - 1]);
			if (b1>0)
			{
				face->SetDirectionFlag(0, true);	edgeArr[n1 - 1]->SetLeftFace(face);
			}
			else
			{
				face->SetDirectionFlag(0, false);	edgeArr[n1 - 1]->SetRightFace(face);
			}
			if (b2>0)
			{
				face->SetDirectionFlag(1, true);	edgeArr[n2 - 1]->SetLeftFace(face);
			}
			else
			{
				face->SetDirectionFlag(1, false);	edgeArr[n2 - 1]->SetRightFace(face);
			}
			if (b3>0)
			{
				face->SetDirectionFlag(2, true);	edgeArr[n3 - 1]->SetLeftFace(face);
			}
			else
			{
				face->SetDirectionFlag(2, false);	edgeArr[n3 - 1]->SetRightFace(face);
			}
			(mesh->GetFaceList()).AddTail(face);

			face->GetNodeRecordPtr(0)->AddFace(face);
			face->GetNodeRecordPtr(1)->AddFace(face);
			face->GetNodeRecordPtr(2)->AddFace(face);
			face->CalPlaneEquation();
		}

		//-----------------------------------------------------------------------------
		//	free the memory
		if (nodeNum>0) delete (QMeshNode **)nodeArr;
		if (edgeNum>0) delete (QMeshEdge **)edgeArr;
		if (faceNum>0) delete (QMeshFace **)faceArr;

		//-----------------------------------------------------------------------------
		//	fill in the information of topology
		for (Pos = mesh->GetEdgeList().GetHeadPosition(); Pos != NULL;) {
			QMeshEdge *edge = (QMeshEdge *)(mesh->GetEdgeList().GetNext(Pos));
			edge->GetStartPoint()->AddEdge(edge);	edge->GetEndPoint()->AddEdge(edge);
			if (edge->GetLeftFace() == NULL || edge->GetRightFace() == NULL) {
				edge->SetAttribFlag(0, true);
				edge->GetStartPoint()->SetAttribFlag(0, true);
				edge->GetEndPoint()->SetAttribFlag(0, true);
			}
		}
	}

	fclose(fp);

	m_polyMeshBody->computeRange();
}

void MeshWorksDataBoard::OutputWIOFile(char* filename)
{
	FILE *fp;	GLKPOSITION PosMesh;	GLKPOSITION Pos;	double xx, yy, zz;
	int i, meshNum, nodeNum, edgeNum, faceNum;

	if (!(fp = fopen(filename, "w"))) return;
	meshNum = m_polyMeshBody->meshList.GetCount();
	fprintf(fp, "%d\n-----------------------------------------------------------\n", meshNum);
	for (PosMesh = m_polyMeshBody->meshList.GetHeadPosition(); PosMesh != NULL;) {
        QSurfaceMesh *mesh = (QSurfaceMesh *)(m_polyMeshBody->meshList.GetNext(PosMesh));
		nodeNum = mesh->GetNodeNumber();
		edgeNum = mesh->GetEdgeNumber();
		faceNum = mesh->GetFaceNumber();
		fprintf(fp, "%d %d %d\n", nodeNum, edgeNum, faceNum);

		i = 0; for (Pos = mesh->GetNodeList().GetHeadPosition(); Pos != NULL; i++) {
			QMeshNode* node = (QMeshNode*)(mesh->GetNodeList().GetNext(Pos));	node->SetIndexNo(i);
			node->GetCoord3D(xx, yy, zz);
			fprintf(fp, "%f %f %f\n", (float)xx, (float)yy, (float)zz);
		}
		i = 0; for (Pos = mesh->GetEdgeList().GetHeadPosition(); Pos != NULL; i++) {
			QMeshEdge* edge = (QMeshEdge*)(mesh->GetEdgeList().GetNext(Pos));	edge->SetIndexNo(i);
			fprintf(fp, "%d %d\n", (edge->GetStartPoint()->GetIndexNo() + 1), (edge->GetEndPoint()->GetIndexNo() + 1));
		}
		i = 0; for (Pos = mesh->GetFaceList().GetHeadPosition(); Pos != NULL; i++) {
			QMeshFace* face = (QMeshFace*)(mesh->GetFaceList().GetNext(Pos));	face->SetIndexNo(i);
			fprintf(fp, "%d %d %d ", (face->GetEdgeRecordPtr(0)->GetIndexNo() + 1),
				(face->GetEdgeRecordPtr(1)->GetIndexNo() + 1), (face->GetEdgeRecordPtr(2)->GetIndexNo() + 1));
			for (int j = 0; j<3; j++)
			if (face->IsNormalDirection(j)) fprintf(fp, "%d ", 1); else fprintf(fp, "%d ", 0);
			fprintf(fp, "\n");
		}
	}
	fclose(fp);
}

void MeshWorksDataBoard::OutputOBJFile(char* filename)
{
	FILE *fp;	GLKPOSITION PosMesh;	GLKPOSITION Pos;	double xx, yy, zz;
	int i, index, num, nodeStIndex;

	for (PosMesh = m_polyMeshBody->meshList.GetHeadPosition(); PosMesh != NULL;) {
        QSurfaceMesh *mesh = (QSurfaceMesh *)(m_polyMeshBody->meshList.GetNext(PosMesh));
		i = 0; for (Pos = mesh->GetNodeList().GetHeadPosition(); Pos != NULL; i++) {
			QMeshNode* node = (QMeshNode*)(mesh->GetNodeList().GetNext(Pos));	node->SetIndexNo(i);
		}
	}

	fp = fopen(filename, "w");
	if (!fp) {
		printf("===============================================\n");
		printf("Can not open the data file - OBJ File Export!\n");
		printf("===============================================\n");
		return;
	}
	fprintf(fp, "# The units used in this file are meters.\n");

	for (PosMesh = m_polyMeshBody->meshList.GetHeadPosition(); PosMesh != NULL;) {
        QSurfaceMesh *mesh = (QSurfaceMesh *)(m_polyMeshBody->meshList.GetNext(PosMesh));
		for (Pos = mesh->GetNodeList().GetHeadPosition(); Pos != NULL;) {
			QMeshNode* node = (QMeshNode*)(mesh->GetNodeList().GetNext(Pos));
			node->GetCoord3D(xx, yy, zz);
			fprintf(fp, "v %.12f %.12f %.12f\n", xx, yy, zz);
		}
	}
	fprintf(fp, "\n");

	nodeStIndex = 1;
	for (PosMesh = m_polyMeshBody->meshList.GetHeadPosition(); PosMesh != NULL;) {
        QSurfaceMesh *mesh = (QSurfaceMesh *)(m_polyMeshBody->meshList.GetNext(PosMesh));

		for (Pos = mesh->GetFaceList().GetHeadPosition(); Pos != NULL;) {
			QMeshFace *face = (QMeshFace *)(mesh->GetFaceList().GetNext(Pos));
			num = face->GetEdgeNum();
			fprintf(fp, "f ");
			for (i = 0; i<num; i++) {
				index = nodeStIndex + face->GetNodeRecordPtr(i)->GetIndexNo();
				fprintf(fp, "%d ", index);
			}
			fprintf(fp, "\n");
		}
		nodeStIndex += mesh->GetNodeNumber();
	}
}
