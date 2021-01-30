// VSA.cpp: implementation of the Variational Shape Approximation class.
//
//////////////////////////////////////////////////////////////////////
#include "VSA.h"

#include <math.h>
#include <iostream>

#include "../GLKLib/GLKObList.h"
#include "../GLKLib/GLKGeometry.h"
#include "../GLKLib/GLKMatrixlib.h"

#include "../QMeshLib/QMesh/QMeshNode.h"
#include "../QMeshLib/QMesh/QMeshEdge.h"
#include "../QMeshLib/QMesh/QMeshFace.h"
#include "../QMeshLib/QMesh/QMeshPatch.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
VSAHeapNode::VSAHeapNode()
{
    index=0;
}
VSAHeapNode::~VSAHeapNode()
{
}

VSANode::VSANode()
{
    ProxyIndicator=0;
    for(int i=0; i<3; i++) { centerpnt[i]=0.0; }
    area=0.0;

    meshobj = NULL;		//this will not release memory when we call deconstructor
}
VSANode::~VSANode()
{
}

VSA::VSA()
{
    m_meshpatch = NULL;		//this will not release memory when we call deconstructor
    m_vsanodes = NULL;
    m_proxies = NULL;
}
VSA::~VSA()
{
    if(m_vsanodes!=NULL) {
        for(int i=0; i<m_facenum; i++)
            delete (VSANode *)m_vsanodes[i];
        delete [] (VSANode **)m_vsanodes;
    }
    if(m_proxies!=NULL) {
        for(int i=0; i<m_regionnum; i++)
            delete (double *)m_proxies[i];
        delete [] (double **)m_proxies;
    }
}

//////////////////////////////////////////////////////////////////////
// Implementation for VSANode
//////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////
// Implementation for VSA
//////////////////////////////////////////////////////////////////////
void VSA::InitializeVSA(QMeshPatch *parapatch, int pararegionnum, short paradim)
{
    m_dim = paradim;
    m_regionnum = pararegionnum;
    if(m_regionnum==0)		//we should have at least one region
        m_regionnum++;
    if(m_dim==2)
        m_facenum = parapatch->GetEdgeList().GetCount();
    else
        m_facenum = parapatch->GetFaceList().GetCount();

    m_meshpatch = parapatch;

    m_vsanodes = (VSANode **)new long[m_facenum];		//Note: we suppose the index for edges and faces are ready for use
    for(int i=0; i<m_facenum; i++)
    {
        m_vsanodes[i] = new VSANode;
    }

    m_proxies = (double **)new long[m_regionnum];
    for(int i=0; i<m_regionnum; i++)
    {
        if(m_dim==2)
            m_proxies[i] = new double[3];
        else
            m_proxies[i] = new double[4];
    }
}
void VSA::PerformVSA2D(int iternum, double desireddistterror, bool needcontourbdcheck)
{
    GLKGeometry geo;
    GLKMatrixLib mal;

    //build relationship between VSA heap Nodes and QMesh entities
    QMeshEdge *tempedge;
    double startpntcoord[3];
    double endpntcoord[3];
    for(GLKPOSITION Pos=m_meshpatch->GetEdgeList().GetHeadPosition(); Pos!=NULL; )
    {
        tempedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetNext(Pos);
        tempedge->regionind = -1;
        tempedge->GetStartPoint()->GetCoord3D(startpntcoord[0], startpntcoord[1], startpntcoord[2]);
        tempedge->GetEndPoint()->GetCoord3D(endpntcoord[0], endpntcoord[1], endpntcoord[2]);

        m_vsanodes[tempedge->GetIndexNo()]->meshobj = tempedge;
        m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator = -1;
        m_vsanodes[tempedge->GetIndexNo()]->area = tempedge->GetLength();                 //Note: here we suppose the length of edge has already been calculated, which may not true for other appl
        m_vsanodes[tempedge->GetIndexNo()]->centerpnt[0] = (startpntcoord[0]+endpntcoord[0])/2.0;			//Note: here, we suppose the y axis is trivial axis
        m_vsanodes[tempedge->GetIndexNo()]->centerpnt[2] = (startpntcoord[2]+endpntcoord[2])/2.0;
    }

    //select m_regionnum starting faces and initialize the proxies
    int steplength = m_meshpatch->GetEdgeList().GetCount()/m_regionnum;
    for(int i=0; i<m_regionnum; i++)
    {
        GLKPOSITION tempos = m_meshpatch->GetEdgeList().FindIndex(i*steplength);
        tempedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetAt(tempos);
        tempedge->GetStartPoint()->GetCoord3D(startpntcoord[0], startpntcoord[1], startpntcoord[2]);
        tempedge->GetEndPoint()->GetCoord3D(endpntcoord[0], endpntcoord[1], endpntcoord[2]);

        m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator = i;
        geo.CalLineEquation(m_proxies[i][0], m_proxies[i][1], m_proxies[i][2], startpntcoord[0], startpntcoord[2], endpntcoord[0], endpntcoord[2]);
    }


    GLKHeap *heap=new GLKHeap(m_facenum*2,true);
    double ***covmatrixarray = (double ***)new long[m_regionnum];
    for(int i=0; i<m_regionnum; i++)
        mal.CreateMatrix(covmatrixarray[i], 2, 2);
    double *totalareaforregions = new double[m_regionnum];
    double *regioncenterx = new double[m_regionnum];
    double *regioncenterz = new double[m_regionnum];
    double *mindisttforregions = new double[m_regionnum];
    double *regiondistortionerrors = new double[m_regionnum];
    VSANode **tempvsanodelist = (VSANode **)new long[m_regionnum];

    //the iteration to clustering, fitting and partition
    double localmaxregiondistterror;
    int maxdisttregion;
    int maxdisttvsanodeind;
    localmaxregiondistterror = this->VSACoreIterations(iternum, heap, covmatrixarray, totalareaforregions, regioncenterx, regioncenterz,
                                                       mindisttforregions, regiondistortionerrors, tempvsanodelist, maxdisttregion, maxdisttvsanodeind, needcontourbdcheck);

    //perform region insertion if necessary
    //while(localmaxregiondistterror>desireddistterror)
    //{
    //	m_regionnum++;		//increase the regions number

    //	//deal with proxies
    //	double **newproxies = (double **)new long[m_regionnum];
    //	for(int i=0; i<m_regionnum; i++)
    //		newproxies[i] = new double[3];
    //	for(int i=0; i<m_regionnum; i++)
    //	{
    //		if(i<=maxdisttregion)
    //		{
    //			newproxies[i][0] = m_proxies[i][0];		newproxies[i][1] = m_proxies[i][1];		newproxies[i][2] = m_proxies[i][2];
    //		}
    //		else if(i==maxdisttregion+1)
    //		{
    //			GLKPOSITION tempos = m_meshpatch->GetEdgeList().FindIndex(maxdisttvsanodeind);
    //			tempedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetAt(tempos);
    //			tempedge->GetStartPoint()->GetCoord3D(startpntcoord[0], startpntcoord[1], startpntcoord[2]);
    //			tempedge->GetEndPoint()->GetCoord3D(endpntcoord[0], endpntcoord[1], endpntcoord[2]);

    //			geo.CalLineEquation(newproxies[i][0], newproxies[i][1], newproxies[i][2], startpntcoord[0], startpntcoord[2], endpntcoord[0], endpntcoord[2]);
    //		}
    //		else
    //		{
    //			newproxies[i][0] = m_proxies[i-1][0];		newproxies[i][1] = m_proxies[i-1][1];		newproxies[i][2] = m_proxies[i-1][2];
    //		}
    //	}
    //	for(int i=0; i<m_regionnum-1; i++)
    //		delete (double *)m_proxies[i];
    //	delete [] (double **)m_proxies;
    //	m_proxies = newproxies;

    //	//deal with the region label
    //	for(GLKPOSITION Pos=m_meshpatch->GetEdgeList().GetHeadPosition(); Pos!=NULL; )
    //	{
    //		tempedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetNext(Pos);
    //		if(m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator>maxdisttregion)
    //		{
    //			m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator = 1+m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator;
    //		}
    //	}
    //	m_vsanodes[maxdisttvsanodeind]->ProxyIndicator = maxdisttregion+1;

    //	//deal with the dimension of dynamic arrays
    //	delete [] totalareaforregions;
    //	delete [] regioncenterx;
    //	delete [] regioncenterz;
    //	delete [] mindisttforregions;
    //	delete [] regiondistortionerrors;
    //	delete [] (VSANode **)tempvsanodelist;
    //	totalareaforregions = new double[m_regionnum];
    //	regioncenterx = new double[m_regionnum];
    //	regioncenterz = new double[m_regionnum];
    //	mindisttforregions = new double[m_regionnum];
    //	regiondistortionerrors = new double[m_regionnum];
    //	tempvsanodelist = (VSANode **)new long[m_regionnum];

    //	//deal with the covmatrixarray
    //	double ***newcovmatrixarray = (double ***)new long[m_regionnum];
    //	for(int i=0; i<m_regionnum-1; i++)
    //		newcovmatrixarray[i] = covmatrixarray[i];
    //	mal.CreateMatrix(newcovmatrixarray[m_regionnum-1], 2, 2);
    //	delete [] (double ***)covmatrixarray;
    //	covmatrixarray = newcovmatrixarray;

    //	//continue to call VSA iterations function for several iterations
    //	localmaxregiondistterror = this->VSACoreIterations(5, heap, covmatrixarray, totalareaforregions, regioncenterx, regioncenterz,
    //			mindisttforregions, regiondistortionerrors, tempvsanodelist, maxdisttregion, maxdisttvsanodeind, needcontourbdcheck);
    //}

    delete heap;
    for(int i=0; i<m_regionnum; i++)
        mal.DeleteMatrix(covmatrixarray[i], 2, 2);
    delete [] (double ***)covmatrixarray;
    delete [] totalareaforregions;
    delete [] regioncenterx;
    delete [] regioncenterz;
    delete [] mindisttforregions;
    delete [] regiondistortionerrors;
    delete [] (VSANode **)tempvsanodelist;

    return ;
}

double VSA::VSACoreIterations(int maxiter, GLKHeap *heap, double ***covmatrixarray, double *totalareaforregions, double *regioncenterx, double *regioncenterz, double *mindisttforregions,
                              double *regiondistortionerrors, VSANode **tempvsanodelist, int &maxdistterrorregion, int &maxdisttvsanodeind, bool needcontourbdcheck)
{
    GLKGeometry geo;
    GLKMatrixLib mal;
    QMeshNode *startpnt;
    QMeshNode *endpnt;
    QMeshEdge *edgewithstartpnt;
    QMeshEdge *edgewithendpnt;
    VSAHeapNode *tempvsaheapnode;
    VSAHeapNode *topheapnode;

    double startpntcoord[3];
    double endpntcoord[3];
    double tempdistterror;

    QMeshEdge *tempedge;
    double p1[2], p2[2];
    double **tempmatrix;
    mal.CreateMatrix(tempmatrix, 2, 2);

    double maxregiondistterror = -1.0e10;
    double maxfacedistterror = -1.0e10;

    // a small function
    auto ComputeInertiaMatrixOfSegment = [](const double *p1, const double *p2, double length, double **result) {
        if(length==0.0)
            length = sqrt((p2[0]-p1[0])*(p2[0]-p1[0])+(p2[1]-p1[1])*(p2[1]-p1[1]));

        double pbar[2] = {(p2[0]-p1[0])/2.0, (p2[1]-p1[1])/2.0};

        double matrixA[2][2] = {{p2[0]-p1[0], 0.0}, {p2[1]-p1[1], 0.0}};

        result[0][0] = length*0.333333333333333333333333*matrixA[0][0]*matrixA[0][0]+length*(p1[0]*p1[0]+p1[0]*pbar[0]+pbar[0]*p1[0]);
        result[0][1] = length*0.333333333333333333333333*matrixA[0][0]*matrixA[1][0]+length*(p1[0]*p1[1]+p1[0]*pbar[1]+pbar[0]*p1[1]);
        result[1][0] = length*0.333333333333333333333333*matrixA[1][0]*matrixA[0][0]+length*(p1[1]*p1[0]+p1[1]*pbar[0]+pbar[1]*p1[0]);
        result[1][1] = length*0.333333333333333333333333*matrixA[1][0]*matrixA[1][0]+length*(p1[1]*p1[1]+p1[1]*pbar[1]+pbar[1]*p1[1]);
    };

    int count = 0;
    while(count<maxiter)
    {
        count++;
        //push the faces adjacent to the seed faces into the queue
        for(GLKPOSITION Pos=m_meshpatch->GetEdgeList().GetHeadPosition(); Pos!=NULL; )
        {
            tempedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetNext(Pos);
            if(m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator!=-1)
            {
                startpnt = tempedge->GetStartPoint();
                endpnt = tempedge->GetEndPoint();
                /*if(startpnt->GetEdgeList().GetCount()==1)
                                        edgewithstartpnt = NULL;
                                else
                                {*/
                edgewithstartpnt = (QMeshEdge *)startpnt->GetEdgeList().GetHead();
                if(edgewithstartpnt!=NULL && edgewithstartpnt->GetIndexNo()==tempedge->GetIndexNo())
                    edgewithstartpnt = (QMeshEdge *)startpnt->GetEdgeList().GetTail();
                if(needcontourbdcheck)
                {
                    if(m_meshpatch->GetEdgeList().Find(edgewithstartpnt)==NULL)             //just check if the edge is contained in current patch
                        edgewithstartpnt = NULL;
                }
                /*}
                                if(endpnt->GetEdgeList().GetCount()==1)
                                        edgewithendpnt = NULL;
                                else
                                {*/
                edgewithendpnt = (QMeshEdge *)endpnt->GetEdgeList().GetTail();
                if(edgewithendpnt!=NULL && edgewithendpnt->GetIndexNo()==tempedge->GetIndexNo())
                    edgewithendpnt = (QMeshEdge *)endpnt->GetEdgeList().GetHead();
                if(needcontourbdcheck)
                {
                    if(m_meshpatch->GetEdgeList().Find(edgewithendpnt)==NULL)
                        edgewithendpnt = NULL;
                }
                /*}*/

                if(edgewithstartpnt!=NULL)
                {
                    if(m_vsanodes[edgewithstartpnt->GetIndexNo()]->ProxyIndicator==-1)
                    {
                        startpnt = edgewithstartpnt->GetStartPoint();
                        endpnt = edgewithstartpnt->GetEndPoint();
                        startpnt->GetCoord3D(startpntcoord[0], startpntcoord[1], startpntcoord[2]);
                        endpnt->GetCoord3D(endpntcoord[0], endpntcoord[1], endpntcoord[2]);
                        tempvsaheapnode = new VSAHeapNode();

                        tempvsaheapnode->attachedObj = m_vsanodes[edgewithstartpnt->GetIndexNo()];
                        tempvsaheapnode->whichproxyagainst = m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator;
                        tempvsaheapnode->SetValue(this->EvaluateDistortionError2D(m_proxies[tempvsaheapnode->whichproxyagainst],
                                                  startpntcoord[0], startpntcoord[2], endpntcoord[0], endpntcoord[2], m_vsanodes[edgewithstartpnt->GetIndexNo()]->area));
                        heap->Insert(tempvsaheapnode);
                    }
                }
                if(edgewithendpnt!=NULL)
                {
                    if(m_vsanodes[edgewithendpnt->GetIndexNo()]->ProxyIndicator==-1)
                    {
                        startpnt = edgewithendpnt->GetStartPoint();
                        endpnt = edgewithendpnt->GetEndPoint();
                        startpnt->GetCoord3D(startpntcoord[0], startpntcoord[1], startpntcoord[2]);
                        endpnt->GetCoord3D(endpntcoord[0], endpntcoord[1], endpntcoord[2]);
                        tempvsaheapnode = new VSAHeapNode();

                        tempvsaheapnode->attachedObj = m_vsanodes[edgewithendpnt->GetIndexNo()];
                        tempvsaheapnode->whichproxyagainst = m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator;
                        tempvsaheapnode->SetValue(this->EvaluateDistortionError2D(m_proxies[tempvsaheapnode->whichproxyagainst],
                                                  startpntcoord[0], startpntcoord[2], endpntcoord[0], endpntcoord[2], m_vsanodes[edgewithendpnt->GetIndexNo()]->area));
                        heap->Insert(tempvsaheapnode);
                    }
                }
            }
        }

        //pop one face each time from queue and label it until the queue is empty
        VSANode *tempvsanode;
        while(!heap->ListEmpty())
        {
            tempvsaheapnode = (VSAHeapNode *)heap->RemoveTop();
            topheapnode = tempvsaheapnode;
            tempvsanode = (VSANode *)tempvsaheapnode->attachedObj;
            if(tempvsanode->ProxyIndicator==-1)
            {
                tempvsanode->ProxyIndicator = tempvsaheapnode->whichproxyagainst;

                tempedge = (QMeshEdge *)tempvsanode->meshobj;

                startpnt = tempedge->GetStartPoint();
                endpnt = tempedge->GetEndPoint();
                /*if(startpnt->GetEdgeList().GetCount()==1)
                                        edgewithstartpnt = NULL;
                                else
                                {*/
                edgewithstartpnt = (QMeshEdge *)startpnt->GetEdgeList().GetHead();
                if(edgewithstartpnt!=NULL && edgewithstartpnt->GetIndexNo()==tempedge->GetIndexNo())
                    edgewithstartpnt = (QMeshEdge *)startpnt->GetEdgeList().GetTail();
                if(needcontourbdcheck)
                {
                    if(m_meshpatch->GetEdgeList().Find(edgewithstartpnt)==NULL)             //just check if the edge is contained in current patch
                        edgewithstartpnt = NULL;
                }
                /*}
                                if(endpnt->GetEdgeList().GetCount()==1)
                                        edgewithendpnt = NULL;
                                else
                                {*/
                edgewithendpnt = (QMeshEdge *)endpnt->GetEdgeList().GetTail();
                if(edgewithendpnt!=NULL && edgewithendpnt->GetIndexNo()==tempedge->GetIndexNo())
                    edgewithendpnt = (QMeshEdge *)endpnt->GetEdgeList().GetHead();
                if(needcontourbdcheck)
                {
                    if(m_meshpatch->GetEdgeList().Find(edgewithendpnt)==NULL)
                        edgewithendpnt = NULL;
                }
                /*}*/

                if(edgewithstartpnt!=NULL)
                {
                    if(m_vsanodes[edgewithstartpnt->GetIndexNo()]->ProxyIndicator==-1)
                    {
                        startpnt = edgewithstartpnt->GetStartPoint();
                        endpnt = edgewithstartpnt->GetEndPoint();
                        startpnt->GetCoord3D(startpntcoord[0], startpntcoord[1], startpntcoord[2]);
                        endpnt->GetCoord3D(endpntcoord[0], endpntcoord[1], endpntcoord[2]);
                        tempvsaheapnode = new VSAHeapNode();

                        tempvsaheapnode->attachedObj = m_vsanodes[edgewithstartpnt->GetIndexNo()];
                        tempvsaheapnode->whichproxyagainst = m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator;
                        tempvsaheapnode->SetValue(this->EvaluateDistortionError2D(m_proxies[tempvsaheapnode->whichproxyagainst],
                                                  startpntcoord[0], startpntcoord[2], endpntcoord[0], endpntcoord[2], m_vsanodes[edgewithstartpnt->GetIndexNo()]->area));
                        heap->Insert(tempvsaheapnode);
                    }
                }
                if(edgewithendpnt!=NULL)
                {
                    if(m_vsanodes[edgewithendpnt->GetIndexNo()]->ProxyIndicator==-1)
                    {
                        startpnt = edgewithendpnt->GetStartPoint();
                        endpnt = edgewithendpnt->GetEndPoint();
                        startpnt->GetCoord3D(startpntcoord[0], startpntcoord[1], startpntcoord[2]);
                        endpnt->GetCoord3D(endpntcoord[0], endpntcoord[1], endpntcoord[2]);
                        tempvsaheapnode = new VSAHeapNode();

                        tempvsaheapnode->attachedObj = m_vsanodes[edgewithendpnt->GetIndexNo()];
                        tempvsaheapnode->whichproxyagainst = m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator;
                        tempvsaheapnode->SetValue(this->EvaluateDistortionError2D(m_proxies[tempvsaheapnode->whichproxyagainst],
                                                  startpntcoord[0], startpntcoord[2], endpntcoord[0], endpntcoord[2], m_vsanodes[edgewithendpnt->GetIndexNo()]->area));
                        heap->Insert(tempvsaheapnode);
                    }
                }
            }

            delete topheapnode;
        }

        //for each region, fit a proxy
        for(int i=0; i<m_regionnum; i++)
        {
            covmatrixarray[i][0][0] = 0.0;	covmatrixarray[i][0][1] = 0.0;	covmatrixarray[i][1][0] = 0.0;	covmatrixarray[i][1][1] = 0.0;
            tempmatrix[0][0] = 0.0;	tempmatrix[0][1] = 0.0;	tempmatrix[1][0] = 0.0;	tempmatrix[1][1] = 0.0;
            totalareaforregions[i] = 0.0;
            regioncenterx[i] = 0.0;
            regioncenterz[i] = 0.0;
        }

        for(GLKPOSITION Pos=m_meshpatch->GetEdgeList().GetHeadPosition(); Pos!=NULL; )
        {
            tempedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetNext(Pos);
            tempedge->GetStartPoint()->GetCoord3D(startpntcoord[0], startpntcoord[1], startpntcoord[2]);
            tempedge->GetEndPoint()->GetCoord3D(endpntcoord[0], endpntcoord[1], endpntcoord[2]);
            p1[0] = startpntcoord[0];		p1[1] = startpntcoord[2];
            p2[0] = endpntcoord[0];		p2[1] = endpntcoord[2];
            ComputeInertiaMatrixOfSegment(p1, p2, tempedge->GetLength(), tempmatrix);

            covmatrixarray[m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator][0][0] += tempmatrix[0][0];
            covmatrixarray[m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator][0][1] += tempmatrix[0][1];
            covmatrixarray[m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator][1][0] += tempmatrix[1][0];
            covmatrixarray[m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator][1][1] += tempmatrix[1][1];

            totalareaforregions[m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator] += m_vsanodes[tempedge->GetIndexNo()]->area;
            regioncenterx[m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator] += m_vsanodes[tempedge->GetIndexNo()]->centerpnt[0]*m_vsanodes[tempedge->GetIndexNo()]->area;
            regioncenterz[m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator] += m_vsanodes[tempedge->GetIndexNo()]->centerpnt[2]*m_vsanodes[tempedge->GetIndexNo()]->area;
        }

        double eigvalues[2];
        for(int i=0; i<m_regionnum; i++)
        {
            regioncenterx[i] = regioncenterx[i]/totalareaforregions[i];
            regioncenterz[i] = regioncenterz[i]/totalareaforregions[i];

            covmatrixarray[i][0][0] -= totalareaforregions[i]*regioncenterx[i]*regioncenterx[i];
            covmatrixarray[i][0][1] -= totalareaforregions[i]*regioncenterx[i]*regioncenterz[i];
            covmatrixarray[i][1][0] -= totalareaforregions[i]*regioncenterz[i]*regioncenterx[i];
            covmatrixarray[i][1][1] -= totalareaforregions[i]*regioncenterz[i]*regioncenterz[i];

            covmatrixarray[i][0][0] *= 1.0e10;
            covmatrixarray[i][0][1] *= 1.0e10;
            covmatrixarray[i][1][0] *= 1.0e10;
            covmatrixarray[i][1][1] *= 1.0e10;

            tempmatrix[0][0] = 0.0;  tempmatrix[0][1] = 0.0;  tempmatrix[1][0] = 0.0;  tempmatrix[1][1] = 0.0;
            eigvalues[0] = 0.0;	eigvalues[1] = 0.0;
            mal.JacobianEigensystemSolver(covmatrixarray[i], 2, tempmatrix, eigvalues, 1.0e-5, 30);
            if(eigvalues[0]<eigvalues[1])
            {
                m_proxies[i][0] = tempmatrix[0][0];		m_proxies[i][1] = tempmatrix[1][0];		m_proxies[i][2] = -(m_proxies[i][0]*regioncenterx[i]+m_proxies[i][1]*regioncenterz[i]);
            }
            else
            {
                m_proxies[i][0] = tempmatrix[0][1];		m_proxies[i][1] = tempmatrix[1][1];		m_proxies[i][2] = -(m_proxies[i][0]*regioncenterx[i]+m_proxies[i][1]*regioncenterz[i]);
            }
        }

        //find a new set of seed faces and clean the label of the faces
        for(int i=0; i<m_regionnum; i++)
        {
            mindisttforregions[i] = 1.0e10;
            tempvsanodelist[i] = NULL;
            regiondistortionerrors[i] = 0.0;
        }
        for(GLKPOSITION Pos=m_meshpatch->GetEdgeList().GetHeadPosition(); Pos!=NULL; )
        {
            tempedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetNext(Pos);
            tempedge->GetStartPoint()->GetCoord3D(startpntcoord[0], startpntcoord[1], startpntcoord[2]);
            tempedge->GetEndPoint()->GetCoord3D(endpntcoord[0], endpntcoord[1], endpntcoord[2]);
            tempdistterror = this->EvaluateDistortionError2D(m_proxies[m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator], startpntcoord[0], startpntcoord[2], endpntcoord[0], endpntcoord[2],
                    m_vsanodes[tempedge->GetIndexNo()]->area);

            if(tempdistterror < mindisttforregions[m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator])
            {
                mindisttforregions[m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator] = tempdistterror;
                tempvsanodelist[m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator] = m_vsanodes[tempedge->GetIndexNo()];
            }

            regiondistortionerrors[m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator] += tempdistterror;

            tempedge->regionind = m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator;		//just used to verify the 2D VSA
        }

        //prepare information for insertion region
        maxregiondistterror = -1.0e10;
        maxfacedistterror = -1.0e10;
        for(int i=0; i<m_regionnum; i++)
        {
            if(regiondistortionerrors[i]>maxregiondistterror)
            {
                maxregiondistterror = regiondistortionerrors[i];
                maxdistterrorregion = i;
            }
        }
        for(GLKPOSITION Pos=m_meshpatch->GetEdgeList().GetHeadPosition(); Pos!=NULL; )
        {
            tempedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetNext(Pos);
            if(m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator==maxdistterrorregion)
            {
                tempedge->GetStartPoint()->GetCoord3D(startpntcoord[0], startpntcoord[1], startpntcoord[2]);
                tempedge->GetEndPoint()->GetCoord3D(endpntcoord[0], endpntcoord[1], endpntcoord[2]);
                tempdistterror = this->EvaluateDistortionError2D(m_proxies[maxdistterrorregion], startpntcoord[0], startpntcoord[2], endpntcoord[0], endpntcoord[2],
                        m_vsanodes[tempedge->GetIndexNo()]->area);
                if(tempdistterror>maxfacedistterror && !(m_vsanodes[tempedge->GetIndexNo()]->centerpnt[0]==tempvsanodelist[maxdistterrorregion]->centerpnt[0]&&
                                                         m_vsanodes[tempedge->GetIndexNo()]->centerpnt[2]==tempvsanodelist[maxdistterrorregion]->centerpnt[2]))
                {
                    maxfacedistterror = tempdistterror;
                    maxdisttvsanodeind = tempedge->GetIndexNo();
                }
            }
        }

        for(GLKPOSITION Pos=m_meshpatch->GetEdgeList().GetHeadPosition(); Pos!=NULL; )
        {
            tempedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetNext(Pos);
            m_vsanodes[tempedge->GetIndexNo()]->ProxyIndicator = -1;
        }
        for(int i=0; i<m_regionnum; i++)
        {
            tempvsanodelist[i]->ProxyIndicator = i;
        }
    }

    mal.DeleteMatrix(tempmatrix, 2, 2);

    return maxregiondistterror;
}

float VSA::EvaluateDistortionError2D(double *proxy, double v1, double v2, double u1, double u2, double length)
{
    double pntonline[2];
    if(abs(proxy[0])>0.1)
    {
        pntonline[1] = 0.0;  pntonline[0] = -proxy[2]/proxy[0];
    }
    else
    {
        pntonline[0] = 0.0;  pntonline[1] = -proxy[2]/proxy[1];
    }

    double d1 = (v1-pntonline[0])*proxy[0]+(v2-pntonline[1])*proxy[1];
    double d2 = (u1-pntonline[0])*proxy[0]+(u2-pntonline[1])*proxy[1];

    float result = 0.33333333333333333*length*(d1*d1+d2*d2+d1*d2);

    return result;
}

void VSA::SimplifyMeshBasedOnVSARegions2D()
{
    QMeshEdge *tempedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetHead();
    QMeshEdge *firstedge = tempedge;
    int lastregionindex = firstedge->regionind;

    QMeshPatch *newpatch = new QMeshPatch();

    double lastpntcoord[3];
    double currpntcoord[3];
    QMeshNode *startingnode;
    QMeshNode *lastnode;
    QMeshNode *currentnode;
    QMeshEdge *newtempedge;
    int count=0;

    if(m_regionnum==1 || m_regionnum==2)
    {
        for(GLKPOSITION Pos=m_meshpatch->GetEdgeList().GetHeadPosition(); Pos!=NULL; )
        {
            tempedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetNext(Pos);
            count++;
            if(count==1)
            {
                tempedge->GetStartPoint()->GetCoord3D(lastpntcoord[0], lastpntcoord[1], lastpntcoord[2]);
                lastnode = new QMeshNode();
                lastnode->SetCoord3D(lastpntcoord[0], lastpntcoord[1], lastpntcoord[2]);

                startingnode = lastnode;
                lastregionindex = tempedge->regionind;
            }
            else
            {
                currentnode = new QMeshNode();
                tempedge->GetStartPoint()->GetCoord3D(currpntcoord[0], currpntcoord[1], currpntcoord[2]);
                currentnode->SetCoord3D(currpntcoord[0], currpntcoord[1], currpntcoord[2]);

                newtempedge = new QMeshEdge();
                newtempedge->SetStartPoint(lastnode);
                newtempedge->SetEndPoint(currentnode);
                newtempedge->regionind = lastregionindex;

                lastnode->GetEdgeList().AddTail(newtempedge);
                currentnode->GetEdgeList().AddTail(newtempedge);

                newpatch->GetEdgeList().AddTail(newtempedge);
                newpatch->GetNodeList().AddTail(lastnode);

                lastnode = currentnode;
                lastregionindex = tempedge->regionind;
            }
        }
        newpatch->GetNodeList().AddTail(lastnode);
    }
    else
    {
        for(GLKPOSITION Pos=m_meshpatch->GetEdgeList().GetHeadPosition(); Pos!=NULL; )
        {
            tempedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetNext(Pos);
            if(tempedge->regionind!=lastregionindex /*true*/)
            {
                count++;
                if(count==1)
                {
                    tempedge->GetStartPoint()->GetCoord3D(lastpntcoord[0], lastpntcoord[1], lastpntcoord[2]);
                    lastnode = new QMeshNode();
                    lastnode->SetCoord3D(lastpntcoord[0], lastpntcoord[1], lastpntcoord[2]);

                    startingnode = lastnode;
                    lastregionindex = tempedge->regionind;
                }
                else
                {
                    currentnode = new QMeshNode();
                    tempedge->GetStartPoint()->GetCoord3D(currpntcoord[0], currpntcoord[1], currpntcoord[2]);
                    currentnode->SetCoord3D(currpntcoord[0], currpntcoord[1], currpntcoord[2]);

                    newtempedge = new QMeshEdge();
                    newtempedge->SetStartPoint(lastnode);
                    newtempedge->SetEndPoint(currentnode);
                    newtempedge->regionind = lastregionindex;

                    lastnode->GetEdgeList().AddTail(newtempedge);
                    currentnode->GetEdgeList().AddTail(newtempedge);

                    newpatch->GetEdgeList().AddTail(newtempedge);
                    newpatch->GetNodeList().AddTail(lastnode);

                    lastnode = currentnode;
                    lastregionindex = tempedge->regionind;
                }
            }
        }

        //construct connection segment(s)
        if(tempedge->regionind==firstedge->regionind)
        {
            newtempedge = new QMeshEdge();
            newtempedge->SetStartPoint(lastnode);
            newtempedge->SetEndPoint(startingnode);
            newtempedge->regionind = lastregionindex;

            lastnode->GetEdgeList().AddTail(newtempedge);
            startingnode->GetEdgeList().AddTail(newtempedge);

            newpatch->GetEdgeList().AddTail(newtempedge);
            newpatch->GetNodeList().AddTail(lastnode);
        }
        else
        {
            QMeshNode *firstnode = new QMeshNode();
            double firstcoord[3];
            firstedge->GetStartPoint()->GetCoord3D(firstcoord[0], firstcoord[1], firstcoord[2]);
            firstnode->SetCoord3D(firstcoord[0], firstcoord[1], firstcoord[2]);

            //first connection segment
            newtempedge = new QMeshEdge();
            newtempedge->SetStartPoint(lastnode);
            newtempedge->SetEndPoint(firstnode);
            newtempedge->regionind = lastregionindex;

            lastnode->GetEdgeList().AddTail(newtempedge);
            firstnode->GetEdgeList().AddTail(newtempedge);

            newpatch->GetEdgeList().AddTail(newtempedge);
            newpatch->GetNodeList().AddTail(lastnode);

            //second connection segment
            newtempedge = new QMeshEdge();
            newtempedge->SetStartPoint(firstnode);
            newtempedge->SetEndPoint(startingnode);
            newtempedge->regionind = firstedge->regionind;

            firstnode->GetEdgeList().AddTail(newtempedge);
            startingnode->GetEdgeList().AddTail(newtempedge);

            newpatch->GetEdgeList().AddTail(newtempedge);
            newpatch->GetNodeList().AddTail(firstnode);
        }
    }

    //set index for entities in new patch
    count=0;
    for(GLKPOSITION Pos=newpatch->GetEdgeList().GetHeadPosition(); Pos!=NULL; )
    {
        tempedge = (QMeshEdge *)newpatch->GetEdgeList().GetNext(Pos);
        tempedge->SetIndexNo(count);
        count++;
    }
    count=0;
    QMeshNode *tempnode;
    QMeshNode *firstnode;
    for(GLKPOSITION Pos=newpatch->GetNodeList().GetHeadPosition(); Pos!=NULL; )
    {
        tempnode = (QMeshNode *)newpatch->GetNodeList().GetNext(Pos);
        if(count==0)
            firstnode = tempnode;
        tempnode->SetIndexNo(count);
        count++;
    }
    firstnode->GetCoord3D(currpntcoord[0], currpntcoord[1], currpntcoord[2]);
    QMeshNode *newnode = new QMeshNode();
    newnode->SetCoord3D(currpntcoord[0], currpntcoord[1], currpntcoord[2]);
    newnode->SetIndexNo(count);
    newpatch->GetNodeList().AddTail(newnode);

    m_meshpatch = newpatch;
}

void VSA::SimplifyMeshBasedOnVSARegions2D_LineIntsct()
{
    GLKGeometry geo;

    QMeshEdge *tempedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetHead();
    QMeshEdge *firstedge = tempedge;
    int lastregionindex = firstedge->regionind;

    QMeshPatch *newpatch = new QMeshPatch();

    double lastpntcoord[3];
    double currpntcoord[3];
    QMeshNode *startingnode;
    QMeshNode *lastnode;
    QMeshNode *currentnode;
    QMeshEdge *newtempedge;
    int count=0;
    for(GLKPOSITION Pos=m_meshpatch->GetEdgeList().GetHeadPosition(); Pos!=NULL; )
    {
        tempedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetNext(Pos);
        if(tempedge->regionind!=lastregionindex)
        {
            count++;
            if(count==1)
            {
                tempedge->GetStartPoint()->GetCoord3D(lastpntcoord[0], lastpntcoord[1], lastpntcoord[2]);
                if(!geo.CalTwoLinesIntersection(m_proxies[lastregionindex][0], m_proxies[lastregionindex][1], m_proxies[lastregionindex][2],
                                                m_proxies[tempedge->regionind][0], m_proxies[tempedge->regionind][1], m_proxies[tempedge->regionind][2], lastpntcoord[0], lastpntcoord[2]))
                {
                    std::cout<<"degenerate case detected!!! Two line have no intersection !!!"<<std::endl;
                }
                lastnode = new QMeshNode();
                lastnode->SetCoord3D(lastpntcoord[0], lastpntcoord[1], lastpntcoord[2]);

                startingnode = lastnode;
                lastregionindex = tempedge->regionind;
            }
            else
            {
                currentnode = new QMeshNode();
                tempedge->GetStartPoint()->GetCoord3D(currpntcoord[0], currpntcoord[1], currpntcoord[2]);
                if(!geo.CalTwoLinesIntersection(m_proxies[lastregionindex][0], m_proxies[lastregionindex][1], m_proxies[lastregionindex][2],
                                                m_proxies[tempedge->regionind][0], m_proxies[tempedge->regionind][1], m_proxies[tempedge->regionind][2], currpntcoord[0], currpntcoord[2]))
                {
                    std::cout<<"degenerate case detected!!! Two line have no intersection !!!"<<std::endl;
                }
                currentnode->SetCoord3D(currpntcoord[0], currpntcoord[1], currpntcoord[2]);

                newtempedge = new QMeshEdge();
                newtempedge->SetStartPoint(lastnode);
                newtempedge->SetEndPoint(currentnode);
                newtempedge->regionind = lastregionindex;

                lastnode->GetEdgeList().AddTail(newtempedge);
                currentnode->GetEdgeList().AddTail(newtempedge);

                newpatch->GetEdgeList().AddTail(newtempedge);
                newpatch->GetNodeList().AddTail(lastnode);

                lastnode = currentnode;
                lastregionindex = tempedge->regionind;
            }
        }
    }

    //construct connection segment(s)
    if(tempedge->regionind==firstedge->regionind)
    {
        newtempedge = new QMeshEdge();
        newtempedge->SetStartPoint(lastnode);
        newtempedge->SetEndPoint(startingnode);
        newtempedge->regionind = lastregionindex;

        lastnode->GetEdgeList().AddTail(newtempedge);
        startingnode->GetEdgeList().AddTail(newtempedge);

        newpatch->GetEdgeList().AddTail(newtempedge);
        newpatch->GetNodeList().AddTail(lastnode);
    }
    else
    {
        QMeshNode *firstnode = new QMeshNode();
        double firstcoord[3];
        firstedge->GetStartPoint()->GetCoord3D(firstcoord[0], firstcoord[1], firstcoord[2]);
        if(!geo.CalTwoLinesIntersection(m_proxies[lastregionindex][0], m_proxies[lastregionindex][1], m_proxies[lastregionindex][2],
                                        m_proxies[firstedge->regionind][0], m_proxies[firstedge->regionind][1], m_proxies[firstedge->regionind][2], firstcoord[0], firstcoord[2]))
        {
            std::cout<<"degenerate case detected!!! Two line have no intersection !!!"<<std::endl;
        }
        firstnode->SetCoord3D(firstcoord[0], firstcoord[1], firstcoord[2]);

        //first connection segment
        newtempedge = new QMeshEdge();
        newtempedge->SetStartPoint(lastnode);
        newtempedge->SetEndPoint(firstnode);
        newtempedge->regionind = lastregionindex;

        lastnode->GetEdgeList().AddTail(newtempedge);
        firstnode->GetEdgeList().AddTail(newtempedge);

        newpatch->GetEdgeList().AddTail(newtempedge);
        newpatch->GetNodeList().AddTail(lastnode);

        //second connection segment
        newtempedge = new QMeshEdge();
        newtempedge->SetStartPoint(firstnode);
        newtempedge->SetEndPoint(startingnode);
        newtempedge->regionind = firstedge->regionind;

        firstnode->GetEdgeList().AddTail(newtempedge);
        startingnode->GetEdgeList().AddTail(newtempedge);

        newpatch->GetEdgeList().AddTail(newtempedge);
        newpatch->GetNodeList().AddTail(firstnode);
    }

    //set index for entities in new patch
    count=0;
    for(GLKPOSITION Pos=newpatch->GetEdgeList().GetHeadPosition(); Pos!=NULL; )
    {
        tempedge = (QMeshEdge *)newpatch->GetEdgeList().GetNext(Pos);
        tempedge->SetIndexNo(count);
        count++;
    }
    count=0;
    QMeshNode *tempnode;
    QMeshNode *firstnode;
    for(GLKPOSITION Pos=newpatch->GetNodeList().GetHeadPosition(); Pos!=NULL; )
    {
        tempnode = (QMeshNode *)newpatch->GetNodeList().GetNext(Pos);
        if(count==0)
            firstnode = tempnode;
        tempnode->SetIndexNo(count);
        count++;
    }
    firstnode->GetCoord3D(currpntcoord[0], currpntcoord[1], currpntcoord[2]);
    QMeshNode *newnode = new QMeshNode();
    newnode->SetCoord3D(currpntcoord[0], currpntcoord[1], currpntcoord[2]);
    newnode->SetIndexNo(count);
    newpatch->GetNodeList().AddTail(newnode);

    m_meshpatch = newpatch;
}

void VSA::BinaryImageInOutCorrection(GLKArray *stickindexx, GLKArray *stickindexz, double *biorigin, double bigridwidth, int baseedgenumber)
{
    if(m_regionnum==1 || m_regionnum==2)
        return;

    QMeshEdge *tempedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetHead();
    QMeshEdge *firstedge = tempedge;
    QMeshNode *firstnode = firstedge->GetStartPoint();
    int lastregionindex = firstedge->regionind;

    double lastpntcoord[3];
    double currpntcoord[3];
    bool isviolate = false;
    double stickstart[2];
    double stickend[2];
    QMeshNode *startingnode;
    QMeshNode *lastnode = NULL;
    QMeshNode *currentnode;
    QMeshEdge *newtempedge;
    fakeregionnum = m_regionnum;
    QMeshEdge *startingedgeoflastregion;
    int count=0;
    GLKPOSITION Pos;
    for(Pos=m_meshpatch->GetEdgeList().GetHeadPosition(); Pos!=NULL; )
    {
        tempedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetNext(Pos);
        if(tempedge->regionind!=lastregionindex)
        {
            count++;
            if(count==1)
            {
                lastnode = tempedge->GetStartPoint();

                startingnode = lastnode;
                lastregionindex = tempedge->regionind;
            }
            else
            {
                currentnode = tempedge->GetStartPoint();
                currentnode->GetCoord3D(currpntcoord[0], currpntcoord[1], currpntcoord[2]);
                lastnode->GetCoord3D(lastpntcoord[0], lastpntcoord[1], lastpntcoord[2]);

                for(int i=lastnode->GetIndexNo()+1; i<currentnode->GetIndexNo(); i++)
                {
                    stickstart[0] = biorigin[0]+stickindexx->GetIntAt(2*i)*bigridwidth;
                    stickstart[1] = biorigin[1]+stickindexz->GetIntAt(2*i)*bigridwidth;
                    stickend[0] = biorigin[0]+stickindexx->GetIntAt(2*i+1)*bigridwidth;
                    stickend[1] = biorigin[1]+stickindexz->GetIntAt(2*i+1)*bigridwidth;
                    if(!this->_IsTwoSegmentIntersect(lastpntcoord[0], lastpntcoord[2], currpntcoord[0], currpntcoord[2], stickstart[0], stickstart[1], stickend[0], stickend[1]))
                    {
                        this->_RecursiveLocalVSA(m_meshpatch, lastnode, currentnode, fakeregionnum,
                                                 stickindexx, stickindexz, biorigin, bigridwidth);
                        break;
                    }
                }

                lastnode = currentnode;
                lastregionindex = tempedge->regionind;
                startingedgeoflastregion = tempedge;
            }
        }
    }

    if(lastnode!=NULL)
    {
        GLKPOSITION startedgepos = m_meshpatch->GetEdgeList().Find(startingedgeoflastregion);
        GLKPOSITION endedgepos = m_meshpatch->GetEdgeList().Find(tempedge);
        for(Pos=startedgepos; Pos!=endedgepos; )
        {
            QMeshEdge *temlocaledge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetNext(Pos);
            temlocaledge->regionind = fakeregionnum;
        }
        fakeregionnum++;

        currentnode = tempedge->GetEndPoint();
        currentnode->GetCoord3D(currpntcoord[0], currpntcoord[1], currpntcoord[2]);
        lastnode->GetCoord3D(lastpntcoord[0], lastpntcoord[1], lastpntcoord[2]);

        int speciallimit;
        if(lastnode->GetIndexNo() > currentnode->GetIndexNo())
            speciallimit = currentnode->GetIndexNo()+m_meshpatch->GetEdgeList().GetCount()-1;
        else
            speciallimit = currentnode->GetIndexNo();
        for(int i=lastnode->GetIndexNo()+1; i<speciallimit; i++)
        {
            stickstart[0] = biorigin[0]+stickindexx->GetIntAt(2*i)*bigridwidth;
            stickstart[1] = biorigin[1]+stickindexz->GetIntAt(2*i)*bigridwidth;
            stickend[0] = biorigin[0]+stickindexx->GetIntAt(2*i+1)*bigridwidth;
            stickend[1] = biorigin[1]+stickindexz->GetIntAt(2*i+1)*bigridwidth;
            if(!this->_IsTwoSegmentIntersect(lastpntcoord[0], lastpntcoord[2], currpntcoord[0], currpntcoord[2], stickstart[0], stickstart[1], stickend[0], stickend[1]))
            {
                this->_RecursiveLocalVSA(m_meshpatch, lastnode, currentnode, fakeregionnum,
                                         stickindexx, stickindexz, biorigin, bigridwidth);
                break;
            }
        }

        lastnode = firstnode;
        currentnode = startingnode;
        currentnode->GetCoord3D(currpntcoord[0], currpntcoord[1], currpntcoord[2]);
        lastnode->GetCoord3D(lastpntcoord[0], lastpntcoord[1], lastpntcoord[2]);

        for(int i=lastnode->GetIndexNo()+1; i<currentnode->GetIndexNo(); i++)
        {
            stickstart[0] = biorigin[0]+stickindexx->GetIntAt(2*i)*bigridwidth;
            stickstart[1] = biorigin[1]+stickindexz->GetIntAt(2*i)*bigridwidth;
            stickend[0] = biorigin[0]+stickindexx->GetIntAt(2*i+1)*bigridwidth;
            stickend[1] = biorigin[1]+stickindexz->GetIntAt(2*i+1)*bigridwidth;
            if(!this->_IsTwoSegmentIntersect(lastpntcoord[0], lastpntcoord[2], currpntcoord[0], currpntcoord[2], stickstart[0], stickstart[1], stickend[0], stickend[1]))
            {
                this->_RecursiveLocalVSA(m_meshpatch, lastnode, currentnode, fakeregionnum,
                                         stickindexx, stickindexz, biorigin, bigridwidth);
                break;
            }
        }
    }
}

bool VSA::_IsTwoSegmentIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
{
    GLKGeometry geo;

    double a1,b1,c1,a2,b2,c2,xx,yy;

    geo.CalLineEquation(a1,b1,c1,x1,y1,x2,y2);
    geo.CalLineEquation(a2,b2,c2,x3,y3,x4,y4);

    if (!(geo.CalTwoLinesIntersection(a1,b1,c1,a2,b2,c2,xx,yy))) return false;

    double u1;
    if (x3==x4)
        u1=(yy-y3)/(y4-y3);
    else
        u1=(xx-x3)/(x4-x3);

    double u2;
    if (x1==x2)
        u2=(yy-y1)/(y2-y1);
    else
        u2=(xx-x1)/(x2-x1);

    if ((u1>=0.0) && (u1<=1.0) && (u2>=0.0) && (u2<=1.0)) return true;

    return false;
}

int VSA::_RecursiveLocalVSA(QMeshPatch *parapatch,  QMeshNode *startnode, QMeshNode *endnode, int &currregionnum, 
                            GLKArray *stickindexx, GLKArray *stickindexz, double *biorigin, double bigridwidth)
{
    QMeshPatch *tempatch = new QMeshPatch();
    QMeshEdge *temedge = (QMeshEdge *)startnode->GetEdgeList().GetTail();
    GLKPOSITION startedgepos = parapatch->GetEdgeList().Find(temedge);
    temedge = (QMeshEdge *)endnode->GetEdgeList().GetHead();
    GLKPOSITION endedgepos = parapatch->GetEdgeList().Find(temedge);

    GLKPOSITION Pos;
    int localcount = 0;
    for(Pos=startedgepos; Pos!=endedgepos; )
    {
        temedge = (QMeshEdge *)parapatch->GetEdgeList().GetNext(Pos);
        temedge->SetIndexNo(localcount);
        localcount++;
        tempatch->GetEdgeList().AddTail(temedge);
    }
    temedge = (QMeshEdge *)parapatch->GetEdgeList().GetNext(Pos);
    temedge->SetIndexNo(localcount);
    tempatch->GetEdgeList().AddTail(temedge);

    VSA newvsaoper;
    newvsaoper.InitializeVSA(tempatch, 2, 2);
    newvsaoper.PerformVSA2D(8, 1.0e8, true);

    //find out the boundary point bewteen two region
    temedge = (QMeshEdge *)tempatch->GetEdgeList().GetTail();
    int secondregionid = temedge->regionind;
    temedge = (QMeshEdge *)tempatch->GetEdgeList().GetHead();
    int firstregionind = temedge->regionind;
    QMeshNode *midpnt;
    for(Pos=tempatch->GetEdgeList().GetHeadPosition(); Pos!=NULL; )
    {
        temedge = (QMeshEdge *)tempatch->GetEdgeList().GetNext(Pos);
        if(temedge->regionind!=firstregionind)
        {
            midpnt = temedge->GetStartPoint();
            break;
        }
    }

    double startpntcoord[3];
    double midpntcoord[3];
    double endpntcoord[3];
    startnode->GetCoord3D(startpntcoord[0], startpntcoord[1], startpntcoord[2]);
    midpnt->GetCoord3D(midpntcoord[0], midpntcoord[1], midpntcoord[2]);
    endnode->GetCoord3D(endpntcoord[0], endpntcoord[1], endpntcoord[2]);
    double stickstart[2];
    double stickend[2];
    int firstregionnum = 1;
    int secondregionnum = 1;


    for(int i=startnode->GetIndexNo()+1; i<midpnt->GetIndexNo(); i++)
    {
        stickstart[0] = biorigin[0]+stickindexx->GetIntAt(2*i)*bigridwidth;
        stickstart[1] = biorigin[1]+stickindexz->GetIntAt(2*i)*bigridwidth;
        stickend[0] = biorigin[0]+stickindexx->GetIntAt(2*i+1)*bigridwidth;
        stickend[1] = biorigin[1]+stickindexz->GetIntAt(2*i+1)*bigridwidth;
        if(!this->_IsTwoSegmentIntersect(startpntcoord[0], startpntcoord[2], midpntcoord[0], midpntcoord[2], stickstart[0], stickstart[1], stickend[0], stickend[1]))
        {
            firstregionnum = this->_RecursiveLocalVSA(m_meshpatch, startnode, midpnt, currregionnum, stickindexx, stickindexz, biorigin, bigridwidth);
        }
    }

    int speciallimit;
    if(midpnt->GetIndexNo() > endnode->GetIndexNo())
        speciallimit = endnode->GetIndexNo()+m_meshpatch->GetEdgeList().GetCount()-1;
    else
        speciallimit = endnode->GetIndexNo();
    for(int i=midpnt->GetIndexNo()+1; i<speciallimit; i++)
    {
        stickstart[0] = biorigin[0]+stickindexx->GetIntAt(2*i)*bigridwidth;
        stickstart[1] = biorigin[1]+stickindexz->GetIntAt(2*i)*bigridwidth;
        stickend[0] = biorigin[0]+stickindexx->GetIntAt(2*i+1)*bigridwidth;
        stickend[1] = biorigin[1]+stickindexz->GetIntAt(2*i+1)*bigridwidth;
        if(!this->_IsTwoSegmentIntersect(midpntcoord[0], midpntcoord[2], endpntcoord[0], endpntcoord[2], stickstart[0], stickstart[1], stickend[0], stickend[1]))
        {
            secondregionnum = this->_RecursiveLocalVSA(m_meshpatch, midpnt, endnode, currregionnum, stickindexx, stickindexz, biorigin, bigridwidth);
        }
    }
    //increase region index to the number bigger than current region number by 1 or 2
    if(firstregionnum==1)
    {
        for(Pos=tempatch->GetEdgeList().GetHeadPosition(); Pos!=NULL; )
        {
            temedge = (QMeshEdge *)tempatch->GetEdgeList().GetNext(Pos);
            if(temedge->regionind==firstregionind)
            {
                temedge->regionind = currregionnum;
            }
        }
    }
    currregionnum++;
    if(secondregionnum==1)
    {
        for(Pos=tempatch->GetEdgeList().GetHeadPosition(); Pos!=NULL; )
        {
            temedge = (QMeshEdge *)tempatch->GetEdgeList().GetNext(Pos);
            if(temedge->regionind==secondregionid)
            {
                temedge->regionind = currregionnum;
            }
        }
    }
    currregionnum++;

    tempatch->GetEdgeList().RemoveAll();
    delete tempatch;

    return firstregionnum+secondregionnum;
}

double VSA::DistortionErrorCorrection(double desireddistterror)
{
    if(m_regionnum==1 || m_regionnum==2)
        return 0.0;

    GLKGeometry geo;

    QMeshEdge *tempedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetHead();
    QMeshEdge *firstedge = tempedge;
    int lastregionindex = firstedge->regionind;

    double lastpntcoord[3];
    double currpntcoord[3];
    double v1[3];
    double v2[3];
    double fakeproxy[3];
    double length;
    double regiondistterror;
    QMeshNode *startingnode;
    QMeshNode *lastnode;
    QMeshNode *currentnode;
    GLKPOSITION startedgepos;
    GLKPOSITION endedgepos;
    int count=0;
    double largerdistterror = 0.0;
    for(GLKPOSITION Pos=m_meshpatch->GetEdgeList().GetHeadPosition(); Pos!=NULL; )
    {
        tempedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetNext(Pos);
        if(tempedge->regionind!=lastregionindex)
        {
            count++;
            if(count==1)
            {
                lastnode = tempedge->GetStartPoint();

                startingnode = lastnode;
                lastregionindex = tempedge->regionind;
            }
            else
            {
                currentnode = tempedge->GetStartPoint();
                currentnode->GetCoord3D(currpntcoord[0], currpntcoord[1], currpntcoord[2]);
                lastnode->GetCoord3D(lastpntcoord[0], lastpntcoord[1], lastpntcoord[2]);
                geo.CalLineEquation(fakeproxy[0], fakeproxy[1], fakeproxy[2], lastpntcoord[0], lastpntcoord[2], currpntcoord[0], currpntcoord[2]);

                tempedge = (QMeshEdge *)lastnode->GetEdgeList().GetTail();
                startedgepos = m_meshpatch->GetEdgeList().Find(tempedge);
                tempedge = (QMeshEdge *)currentnode->GetEdgeList().GetHead();
                endedgepos = m_meshpatch->GetEdgeList().Find(tempedge);

                regiondistterror = 0.0;
                GLKPOSITION Pos;
                for(Pos=startedgepos; Pos!=endedgepos; )
                {
                    tempedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetNext(Pos);
                    tempedge->GetStartPoint()->GetCoord3D(v1[0], v1[1], v1[2]);
                    tempedge->GetEndPoint()->GetCoord3D(v2[0], v2[1], v2[2]);
                    length = sqrt((v1[0]-v2[0])*(v1[0]-v2[0])+(v1[2]-v2[2])*(v1[2]-v2[2]));
                    regiondistterror += this->EvaluateDistortionError2D(fakeproxy, v1[0], v1[2], v2[0], v2[2], length);
                }
                tempedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetNext(Pos);
                tempedge->GetStartPoint()->GetCoord3D(v1[0], v1[1], v1[2]);
                tempedge->GetEndPoint()->GetCoord3D(v2[0], v2[1], v2[2]);
                length = sqrt((v1[0]-v2[0])*(v1[0]-v2[0])+(v1[2]-v2[2])*(v1[2]-v2[2]));
                regiondistterror += this->EvaluateDistortionError2D(fakeproxy, v1[0], v1[2], v2[0], v2[2], length);

                if(regiondistterror>desireddistterror)
                {
                    this->_RecursiveLocalVSAForDistterror(m_meshpatch, lastnode, currentnode, fakeregionnum,desireddistterror, largerdistterror);
                    /*std::cout<<"distortion error correction performed"<<std::endl;*/
                }
                else
                {
                    if(regiondistterror>largerdistterror)
                        largerdistterror = regiondistterror;
                }

                lastnode = currentnode;
                lastregionindex = tempedge->regionind;
            }
        }
    }
    return largerdistterror;
}

int VSA::_RecursiveLocalVSAForDistterror(QMeshPatch *parapatch, QMeshNode *startnode, QMeshNode *endnode, int &currregionnum, double desireddistterror, double &largerdistterror)
{
    GLKGeometry geo;

    QMeshPatch *tempatch = new QMeshPatch();
    QMeshEdge *temedge = (QMeshEdge *)startnode->GetEdgeList().GetTail();
    GLKPOSITION startedgepos = parapatch->GetEdgeList().Find(temedge);
    temedge = (QMeshEdge *)endnode->GetEdgeList().GetHead();
    GLKPOSITION endedgepos = parapatch->GetEdgeList().Find(temedge);

    GLKPOSITION Pos;
    int localcount = 0;
    for(Pos=startedgepos; Pos!=endedgepos; )
    {
        temedge = (QMeshEdge *)parapatch->GetEdgeList().GetNext(Pos);
        temedge->SetIndexNo(localcount);
        localcount++;
        tempatch->GetEdgeList().AddTail(temedge);
    }
    temedge = (QMeshEdge *)parapatch->GetEdgeList().GetNext(Pos);
    temedge->SetIndexNo(localcount);
    tempatch->GetEdgeList().AddTail(temedge);

    VSA newvsaoper;
    newvsaoper.InitializeVSA(tempatch, 2, 2);
    newvsaoper.PerformVSA2D(8, 1.0e8, true);

    //find out the boundary point bewteen two region
    temedge = (QMeshEdge *)tempatch->GetEdgeList().GetHead();
    int firstregionind = temedge->regionind;
    QMeshNode *midpnt;
    for(Pos=tempatch->GetEdgeList().GetHeadPosition(); Pos!=NULL; )
    {
        temedge = (QMeshEdge *)tempatch->GetEdgeList().GetNext(Pos);
        if(temedge->regionind!=firstregionind)
        {
            midpnt = temedge->GetStartPoint();
            break;
        }
    }

    double startpntcoord[3];
    double midpntcoord[3];
    double endpntcoord[3];
    startnode->GetCoord3D(startpntcoord[0], startpntcoord[1], startpntcoord[2]);
    midpnt->GetCoord3D(midpntcoord[0], midpntcoord[1], midpntcoord[2]);
    endnode->GetCoord3D(endpntcoord[0], endpntcoord[1], endpntcoord[2]);
    double v1[3];
    double v2[3];
    double fakeproxy[3];
    double length;
    double regiondistterror;
    int firstregionnum = 1;
    int secondregionnum = 1;

    geo.CalLineEquation(fakeproxy[0], fakeproxy[1], fakeproxy[2], startpntcoord[0], startpntcoord[2], midpntcoord[0], midpntcoord[2]);
    temedge = (QMeshEdge *)startnode->GetEdgeList().GetTail();
    startedgepos = m_meshpatch->GetEdgeList().Find(temedge);
    temedge = (QMeshEdge *)midpnt->GetEdgeList().GetHead();
    endedgepos = m_meshpatch->GetEdgeList().Find(temedge);
    regiondistterror = 0.0;
    for(Pos=startedgepos; Pos!=endedgepos; )
    {
        temedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetNext(Pos);
        temedge->GetStartPoint()->GetCoord3D(v1[0], v1[1], v1[2]);
        temedge->GetEndPoint()->GetCoord3D(v2[0], v2[1], v2[2]);
        length = sqrt((v1[0]-v2[0])*(v1[0]-v2[0])+(v1[2]-v2[2])*(v1[2]-v2[2]));
        regiondistterror += this->EvaluateDistortionError2D(fakeproxy, v1[0], v1[2], v2[0], v2[2], length);
    }
    temedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetNext(Pos);
    temedge->GetStartPoint()->GetCoord3D(v1[0], v1[1], v1[2]);
    temedge->GetEndPoint()->GetCoord3D(v2[0], v2[1], v2[2]);
    length = sqrt((v1[0]-v2[0])*(v1[0]-v2[0])+(v1[2]-v2[2])*(v1[2]-v2[2]));
    regiondistterror += this->EvaluateDistortionError2D(fakeproxy, v1[0], v1[2], v2[0], v2[2], length);
    if(regiondistterror>desireddistterror)
        firstregionnum = this->_RecursiveLocalVSAForDistterror(m_meshpatch, startnode, midpnt, fakeregionnum,desireddistterror, largerdistterror);
    else
    {
        if(regiondistterror>largerdistterror)
            largerdistterror = regiondistterror;
    }

    geo.CalLineEquation(fakeproxy[0], fakeproxy[1], fakeproxy[2], midpntcoord[0], midpntcoord[2], endpntcoord[0], endpntcoord[2]);
    temedge = (QMeshEdge *)midpnt->GetEdgeList().GetTail();
    startedgepos = m_meshpatch->GetEdgeList().Find(temedge);
    temedge = (QMeshEdge *)endnode->GetEdgeList().GetHead();
    endedgepos = m_meshpatch->GetEdgeList().Find(temedge);
    regiondistterror = 0.0;
    for(Pos=startedgepos; Pos!=endedgepos; )
    {
        temedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetNext(Pos);
        temedge->GetStartPoint()->GetCoord3D(v1[0], v1[1], v1[2]);
        temedge->GetEndPoint()->GetCoord3D(v2[0], v2[1], v2[2]);
        length = sqrt((v1[0]-v2[0])*(v1[0]-v2[0])+(v1[2]-v2[2])*(v1[2]-v2[2]));
        regiondistterror += this->EvaluateDistortionError2D(fakeproxy, v1[0], v1[2], v2[0], v2[2], length);
    }
    temedge = (QMeshEdge *)m_meshpatch->GetEdgeList().GetNext(Pos);
    temedge->GetStartPoint()->GetCoord3D(v1[0], v1[1], v1[2]);
    temedge->GetEndPoint()->GetCoord3D(v2[0], v2[1], v2[2]);
    length = sqrt((v1[0]-v2[0])*(v1[0]-v2[0])+(v1[2]-v2[2])*(v1[2]-v2[2]));
    regiondistterror += this->EvaluateDistortionError2D(fakeproxy, v1[0], v1[2], v2[0], v2[2], length);
    if(regiondistterror>desireddistterror)
        secondregionnum = this->_RecursiveLocalVSAForDistterror(m_meshpatch, midpnt, endnode, fakeregionnum,desireddistterror, largerdistterror);
    else
    {
        if(regiondistterror>largerdistterror)
            largerdistterror = regiondistterror;
    }
    //increase region index to the number bigger than current region number by 1 or 2
    if(firstregionnum==1)
    {
        for(Pos=tempatch->GetEdgeList().GetHeadPosition(); Pos!=NULL; )
        {
            temedge = (QMeshEdge *)tempatch->GetEdgeList().GetNext(Pos);
            if(temedge->regionind==firstregionind)
            {
                temedge->regionind = currregionnum;
            }
        }
    }
    currregionnum++;
    if(secondregionnum==1)
    {
        for(Pos=tempatch->GetEdgeList().GetHeadPosition(); Pos!=NULL; )
        {
            temedge = (QMeshEdge *)tempatch->GetEdgeList().GetNext(Pos);
            if(temedge->regionind!=currregionnum-1)
            {
                temedge->regionind = currregionnum;
            }
        }
    }
    currregionnum++;

    tempatch->GetEdgeList().RemoveAll();
    delete tempatch;

    return firstregionnum+secondregionnum;
}
