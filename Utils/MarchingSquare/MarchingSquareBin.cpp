#include "MarchingSquareBin.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <set>
#include <queue>
#include <list>
#include <memory>
#include <time.h>

#include <Eigen/Core>

#include "../GLKLib/GLKGeometry.h"

#include "../QMeshLib/QMesh/QMeshNode.h"
#include "../QMeshLib/QMesh/QMeshEdge.h"
#include "../QMeshLib/QMesh/QMeshFace.h"
#include "../QMeshLib/QMesh/QMeshPatch.h"

#include "VSA.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
MarchingSquareBin::MarchingSquareBin(std::unordered_map<std::string, double> &m_image_, size_t xres_, size_t zres_,
                                     double xzgridwidth_, double height_, double *imageorigin_)
    : m_image(&m_image_),
      xres(xres_),
      zres(zres_),
      xzgridwidth(xzgridwidth_),
      height(height_)
{
    imageorigin[0] = imageorigin_[0];
    imageorigin[1] = imageorigin_[1];

    stickindexx = new GLKArray(50, 50, 1);;
    stickindexz = new GLKArray(50, 50, 1);;
    m_contour = NULL;
}

MarchingSquareBin::~MarchingSquareBin()
{
    if(stickindexx!=NULL) delete stickindexx;
    if(stickindexz!=NULL) delete stickindexz;
    if(m_contour!=NULL) delete m_contour;
}

void MarchingSquareBin::doContouringWithOptim()
{
    double smootherrorcoeff = 0.001;
    int VSAiter = 20;
    int simpratio = 10;
    double movementcoef = 0.65;
    double desireddistortionerror = 0.1 * xzgridwidth * xzgridwidth;

    FindAllSticks();
    basicContouring();
    BiLaplacianSmoothOp(smootherrorcoeff, movementcoef);
    PerformVSAOnContour(VSAiter, desireddistortionerror, simpratio);
    convert2XYPlane();
}

void MarchingSquareBin::doContouring() {
    FindAllSticks();
    basicContouring();
    convert2XYPlane();
}

void MarchingSquareBin::FindAllSticks()
{
    // stick direction: from an outside node to an inside node
    for(int i=0; i<xres; i++) {
        for(int j=0; j<zres-1; j++)  {
            bool v1(m_image->find(std::to_string(i) + " " + std::to_string(j)) == m_image->end());
            bool v2(m_image->find(std::to_string(i) + " " + std::to_string(j + 1)) == m_image->end());
            if(v1 == !v2) {
                stickindexx->Add(i);
                stickindexx->Add(i);
                if(v1) {stickindexz->Add(j); stickindexz->Add(j+1);}
                else { stickindexz->Add(j+1); stickindexz->Add(j);}
            }
        }
    }

    for(int j=0; j<zres; j++) {
        for(int i=0; i<xres-1; i++) {
            bool v1(m_image->find(std::to_string(i) + " " + std::to_string(j)) == m_image->end());
            bool v2(m_image->find(std::to_string(i+1) + " " + std::to_string(j)) == m_image->end());
            if(v1 == !v2)  {
                stickindexz->Add(j);
                stickindexz->Add(j);
                if(v1) {stickindexx->Add(i); stickindexx->Add(i+1);}
                else {stickindexx->Add(i+1);stickindexx->Add(i);}
            }
        }
    }
}

void MarchingSquareBin::basicContouring()
{
    //int ambcount = 0;			//this variable is temp used for count the occurance of ambiguous cases

    int firstx, firstz;
    int secondx, secondz;
    int startx, startz;
    int endx, endz;
    int nextstartx, nextstartz;
    int nextendx, nextendz;
    QMeshNode *firstnode;
    QMeshNode *lastnode;
    QMeshNode *nextnode;
    QMeshEdge *newedge;
    int currsticknum;
    int temstartx, temstartz, temendx, temendz;

    int localcount = 0;

    double startpntsurfacedist;
    double endpntsurfacedist;

    m_contour = new QMeshPatch();

    currsticknum = stickindexx->GetSize()/2;
    int dealtsticknum = 0;
    while(dealtsticknum!=currsticknum) {
        size_t edgeCount(0);
        QMeshFace* curFace = new QMeshFace();
        curFace->SetMeshPatchPtr(m_contour);
        curFace->SetIndexNo(m_contour->GetFaceNumber());
        m_contour->GetFaceList().AddTail(curFace);

        startx = firstx = stickindexx->GetIntAt(dealtsticknum*2);
        startz = firstz = stickindexz->GetIntAt(dealtsticknum*2);
        endx = secondx = stickindexx->GetIntAt(dealtsticknum*2+1);
        endz = secondz = stickindexz->GetIntAt(dealtsticknum*2+1);
        nextstartx = 0;		nextstartz = 0;		nextendx = 0;		nextendz = 0;
        lastnode = new QMeshNode();
        lastnode->SetCoord3D(imageorigin[0]+xzgridwidth*(startx+endx)/2.0, height, imageorigin[1]+xzgridwidth*(startz+endz)/2.0);
        firstnode = lastnode;
        firstnode->AddFace(curFace);

        while(true) {
            // marching squares look-up table
            if(startx==endx) {
                if(startz<endz) {
                    bool v1(m_image->find(std::to_string(endx+1) + " " + std::to_string(endz)) == m_image->end());
                    bool v2(m_image->find(std::to_string(startx+1) + " " + std::to_string(startz)) == m_image->end());
                    if(v1) {
                        nextstartx = endx+1;
                        nextstartz = endz;
                        nextendx = endx;
                        nextendz = endz;
                        /*if(m_image[startx+1][startz]==false)
                                                        ambcount++;*/
                    }
                    else if(v2) {
                        nextstartx = startx+1;
                        nextstartz = startz;
                        nextendx = endx+1;
                        nextendz = endz;
                    }
                    else  {
                        nextstartx = startx;
                        nextstartz = startz;
                        nextendx = startx+1;
                        nextendz = startz;
                    }
                }
                else {
                    bool v1(m_image->find(std::to_string(endx-1) + " " + std::to_string(endz)) == m_image->end());
                    bool v2(m_image->find(std::to_string(startx-1) + " " + std::to_string(startz)) == m_image->end());
                    if(v1)  {
                        nextstartx = endx-1;
                        nextstartz = endz;
                        nextendx = endx;
                        nextendz = endz;
                        /*if(m_image[startx-1][startz]==false)
                                                        ambcount++;*/
                    }
                    else if(v2) {
                        nextstartx = startx-1;
                        nextstartz = startz;
                        nextendx = endx-1;
                        nextendz = endz;
                    }
                    else {
                        nextstartx = startx;
                        nextstartz = startz;
                        nextendx = startx-1;
                        nextendz = startz;
                    }
                }
            }
            else {
                if(startx<endx)  {
                    bool v1(m_image->find(std::to_string(endx) + " " + std::to_string(endz-1)) == m_image->end());
                    bool v2(m_image->find(std::to_string(startx) + " " + std::to_string(startz-1)) == m_image->end());
                    if(v1)  {
                        nextstartx = endx;
                        nextstartz = endz-1;
                        nextendx = endx;
                        nextendz = endz;
                        /*if(m_image[startx][startz-1]==false)
                                                        ambcount++;*/
                    }
                    else if(v2) {
                        nextstartx = startx;
                        nextstartz = startz-1;
                        nextendx = endx;
                        nextendz = endz-1;
                    }
                    else {
                        nextstartx = startx;
                        nextstartz = startz;
                        nextendx = startx;
                        nextendz = startz-1;
                    }
                }
                else {
                    bool v1(m_image->find(std::to_string(endx) + " " + std::to_string(endz+1)) == m_image->end());
                    bool v2(m_image->find(std::to_string(startx) + " " + std::to_string(startz+1)) == m_image->end());
                    if(v1) {
                        nextstartx = endx;
                        nextstartz = endz+1;
                        nextendx = endx;
                        nextendz = endz;
                        /*if(m_image[startx][startz+1]==false)
                                                        ambcount++;*/
                    }
                    else if(v2) {
                        nextstartx = startx;
                        nextstartz = startz+1;
                        nextendx = endx;
                        nextendz = endz+1;
                    }
                    else {
                        nextstartx = startx;
                        nextstartz = startz;
                        nextendx = startx;
                        nextendz = startz+1;
                    }
                }
            }
            // end of the look-up table

            for(int i=dealtsticknum; i<currsticknum; i++) {
                temstartx = stickindexx->GetIntAt(i*2);
                temstartz = stickindexz->GetIntAt(i*2);
                temendx = stickindexx->GetIntAt(i*2+1);
                temendz = stickindexz->GetIntAt(i*2+1);
                if(startx==temstartx && startz==temstartz && endx==temendx && endz==temendz) {
                    temstartx = stickindexx->GetIntAt(dealtsticknum*2);
                    temstartz = stickindexz->GetIntAt(dealtsticknum*2);
                    temendx = stickindexx->GetIntAt(dealtsticknum*2+1);
                    temendz = stickindexz->GetIntAt(dealtsticknum*2+1);
                    stickindexx->SetAt(dealtsticknum*2, startx);
                    stickindexz->SetAt(dealtsticknum*2, startz);
                    stickindexx->SetAt(dealtsticknum*2+1, endx);
                    stickindexz->SetAt(dealtsticknum*2+1, endz);
                    stickindexx->SetAt(i*2, temstartx);
                    stickindexz->SetAt(i*2, temstartz);
                    stickindexx->SetAt(i*2+1, temendx);
                    stickindexz->SetAt(i*2+1, temendz);
                    break;
                }
            }

            // construct the edge
            nextnode = new QMeshNode();
            nextnode->AddFace(curFace);
            //middle point
            nextnode->SetCoord3D(imageorigin[0]+xzgridwidth*(nextstartx+nextendx)/2.0, height, imageorigin[1]+xzgridwidth*(nextstartz+nextendz)/2.0);

            newedge = new QMeshEdge();
            newedge->SetLeftFace(curFace);
            newedge->SetStartPoint(lastnode);
            lastnode->AddEdge(newedge);
            m_contour->GetEdgeList().AddTail(newedge);
            m_contour->GetNodeList().AddTail(lastnode);
            lastnode->SetIndexNo(localcount);
            newedge->SetIndexNo(localcount);
            curFace->SetEdgeRecordPtr(edgeCount, newedge);
            curFace->SetDirectionFlag(edgeCount, true);
            edgeCount++;
            localcount++;
            dealtsticknum++;

            startx = nextstartx;
            startz = nextstartz;
            endx = nextendx;
            endz = nextendz;
            lastnode = nextnode;

            if(nextstartx!=firstx || nextstartz!=firstz || nextendx!=secondx || nextendz!=secondz) {
                newedge->SetEndPoint(nextnode);
                nextnode->AddEdge(newedge);
            }
            else {
                newedge->SetEndPoint(firstnode);
                firstnode->GetEdgeList().AddHead(newedge);
                delete nextnode;
                break;
            }
        }
    }

    //if(ambcount>0) std::cout<<ambcount<<std::endl;
}

void MarchingSquareBin::BiLaplacianSmoothOp(double smootherrorcoeff, double movementcoeffict)
{
    GLKGeometry geo;

    double pnta[3];
    double pntb[3];
    double pntc[3];
    double pntd[3];
    double pnte[3];
    double movevct[3];
    QMeshNode *currnode;
    QMeshNode *nodea;
    QMeshNode *nodeb;
    QMeshNode *noded;
    QMeshNode *nodee;
    QMeshEdge *firstleveledgeb;
    QMeshEdge *firstleveledged;
    QMeshEdge *secondleveledgea;
    QMeshEdge *secondleveledgee;
    QMeshEdge *temedge;
    int startx, startz, endx, endz;
    int localcounter;
    int biggerind;

    double accummovement = 0.0;
    double tempcoord[3];

    double lengthedgeb, lengthedged, targetpnt[3];
    int nodenum = m_contour->GetNodeList().GetCount();
    double **temnodescoord = (double **)new long[nodenum];
    for(int i=0; i<nodenum; i++)
        temnodescoord[i] = new double[3];
    for(GLKPOSITION Pos=m_contour->GetEdgeList().GetHeadPosition(); Pos!=NULL;) {
        temedge = (QMeshEdge *)(m_contour->GetEdgeList().GetNext(Pos));
        temedge->CalLength();
    }

    int itercounter = 0;
    while(itercounter<=30) { // maximum steps
        ++itercounter;
        localcounter = 0;
        for(GLKPOSITION Pos=m_contour->GetNodeList().GetHeadPosition(); Pos!=NULL; ) {
            //find the point a, b, c, d, e.
            currnode = (QMeshNode *)(m_contour->GetNodeList().GetNext(Pos));
            firstleveledgeb = (QMeshEdge *)(currnode->GetEdgeList().GetHead());
            firstleveledged = (QMeshEdge *)(currnode->GetEdgeList().GetTail());
            if(firstleveledgeb->GetStartPoint()->GetIndexNo()!=currnode->GetIndexNo())
                nodeb = firstleveledgeb->GetStartPoint();
            else
                nodeb = firstleveledgeb->GetEndPoint();

            if(firstleveledged->GetStartPoint()->GetIndexNo()!=currnode->GetIndexNo())
                noded = firstleveledged->GetStartPoint();
            else
                noded = firstleveledged->GetEndPoint();

            temedge = (QMeshEdge *)(nodeb->GetEdgeList().GetHead());
            if(temedge->GetIndexNo()!=firstleveledgeb->GetIndexNo())
                secondleveledgea = temedge;
            else
                secondleveledgea = (QMeshEdge *)(nodeb->GetEdgeList().GetTail());

            temedge = (QMeshEdge *)(noded->GetEdgeList().GetHead());
            if(temedge->GetIndexNo()!=firstleveledged->GetIndexNo())
                secondleveledgee = temedge;
            else
                secondleveledgee = (QMeshEdge *)(noded->GetEdgeList().GetTail());

            if(secondleveledgea->GetStartPoint()->GetIndexNo()!=nodeb->GetIndexNo())
                nodea = secondleveledgea->GetStartPoint();
            else
                nodea = secondleveledgea->GetEndPoint();

            if(secondleveledgee->GetStartPoint()->GetIndexNo()!=noded->GetIndexNo())
                nodee = secondleveledgee->GetStartPoint();
            else
                nodee = secondleveledgee->GetEndPoint();

            currnode->GetCoord3D(pntc[0], pntc[1], pntc[2]);
            nodea->GetCoord3D(pnta[0], pnta[1], pnta[2]);
            nodeb->GetCoord3D(pntb[0], pntb[1], pntb[2]);
            noded->GetCoord3D(pntd[0], pntd[1], pntd[2]);
            nodee->GetCoord3D(pnte[0], pnte[1], pnte[2]);

            //improved smoothing technique by only sliding vertices on sticks
            startx = stickindexx->GetIntAt(localcounter*2);		startz = stickindexz->GetIntAt(localcounter*2);
            endx = stickindexx->GetIntAt(localcounter*2+1);		endz = stickindexz->GetIntAt(localcounter*2+1);
            pnta[0] = imageorigin[0]+startx*xzgridwidth;		pnta[2] = imageorigin[1]+startz*xzgridwidth;
            pnte[0] = imageorigin[0]+endx*xzgridwidth;			pnte[2] = imageorigin[1]+endz*xzgridwidth;
            bool intersect;
            intersect = geo.CalTwoLineSegmentsIntersection(pnta[0], pnta[2], pnte[0], pnte[2], pntb[0], pntb[2], pntd[0], pntd[2], targetpnt[0], targetpnt[2]);
            if(!intersect) {
                if(startx==endx) {
                    targetpnt[0] = pnta[0];
                    targetpnt[2] = pntb[2];
                }
                else {
                    targetpnt[0] = pntb[0];
                    targetpnt[2] = pnta[2];
                }
            }
            temnodescoord[localcounter][0] = pntc[0]+movementcoeffict*movementcoeffict*(targetpnt[0]-pntc[0]);
            temnodescoord[localcounter][2] = pntc[2]+movementcoeffict*movementcoeffict*(targetpnt[2]-pntc[2]);
            temnodescoord[localcounter][1] = height;

            ++localcounter;
        }

        //update the coords for all nodes
        localcounter = 0;
        accummovement = 0.0;
        for(GLKPOSITION Pos=m_contour->GetNodeList().GetHeadPosition(); Pos!=NULL; ) {
            currnode = (QMeshNode *)(m_contour->GetNodeList().GetNext(Pos));
            currnode->GetCoord3D(tempcoord[0], tempcoord[1], tempcoord[2]);
            accummovement += sqrt((tempcoord[0]-temnodescoord[localcounter][0])*(tempcoord[0]-temnodescoord[localcounter][0])+
                    (tempcoord[1]-temnodescoord[localcounter][1])*(tempcoord[1]-temnodescoord[localcounter][1])+
                    (tempcoord[2]-temnodescoord[localcounter][2])*(tempcoord[2]-temnodescoord[localcounter][2]));
            currnode->SetCoord3D(temnodescoord[localcounter][0], temnodescoord[localcounter][1], temnodescoord[localcounter][2]);
            localcounter++;
        }
        if(smootherrorcoeff < (accummovement/m_contour->GetNodeList().GetCount())/xzgridwidth<smootherrorcoeff) break;
    }

    for(int i=0; i<nodenum; i++) delete [] (double *)temnodescoord[i];
    delete [] (double **)temnodescoord;
}

double MarchingSquareBin::PerformVSAOnContour(int iter, double paradistterror, int simpratio)
{
    GLKObList temppatchlist;
    temppatchlist.RemoveAll();
    double maxdistterror = 0.0;
    double tempdistterror;

    GLKGeometry geo;
    QMeshNode *firstnode;
    QMeshNode *startnode;
    QMeshNode *endnode;
    double startnodecoord[3];
    double endnodecoord[3];
    double edgevec[3];
    double resultnormal[3];
    double tempnormal[3];
    int localcount = 0;
    double dir[3] = {0.0, 1.0, 0.0};

    int simpcontourcount = 0;
    int contouredgecount = 0;
    QMeshEdge *curredge;
    QMeshNode *currnode;
    double startpntcoord[3];
    double currentendpntcoord[3];
    QMeshPatch *temppatch = new QMeshPatch();
    for(GLKPOSITION Pos=m_contour->GetEdgeList().GetHeadPosition(); Pos!=NULL; ) {
        curredge = (QMeshEdge *)(m_contour->GetEdgeList().GetNext(Pos));
        curredge->SetIndexNo(contouredgecount);
        //curredge->regionind = contouredgecount;       //this is used only when you want to cancle VSA segmentation
        contouredgecount++;

        if(contouredgecount==1)
            curredge->GetStartPoint()->GetCoord3D(startpntcoord[0], startpntcoord[1],startpntcoord[2]);

        temppatch->GetEdgeList().AddTail(curredge);
        curredge->GetEndPoint()->GetCoord3D(currentendpntcoord[0], currentendpntcoord[1], currentendpntcoord[2]);
        if(currentendpntcoord[0]==startpntcoord[0] && currentendpntcoord[1]==startpntcoord[1] && currentendpntcoord[2]==startpntcoord[2]){
            temppatchlist.AddTail(temppatch);
            temppatch = new QMeshPatch();
            contouredgecount = 0;
            simpcontourcount++;
        }
    }
    /*std::cout<<simpcontourcount<<std::endl;*/
    delete temppatch;

    //perform VSA and construct a simplified mesh
    QMeshPatch *simplifiedpatch = new QMeshPatch();
    int totaledgecount = 0;			//this variable count the number of total edge in order to make use of the global stick index stored in stickindexarrays
    for(GLKPOSITION Pos=temppatchlist.GetHeadPosition(); Pos!=NULL; ) {
        temppatch = (QMeshPatch *)(temppatchlist.GetNext(Pos));

        VSA vsaoper;
        vsaoper.InitializeVSA(temppatch, temppatch->GetEdgeList().GetCount()/simpratio, 2);
        vsaoper.PerformVSA2D(iter, paradistterror, false);

        vsaoper.BinaryImageInOutCorrection(stickindexx, stickindexz, imageorigin, xzgridwidth, totaledgecount);
        tempdistterror = vsaoper.DistortionErrorCorrection(paradistterror);
        if(tempdistterror>maxdistterror)
            maxdistterror = tempdistterror;
        totaledgecount += temppatch->GetEdgeList().GetCount();

        //this peice of function simplify the original patch and put the simplification result into the new "simplifiedpatch"
        vsaoper.SimplifyMeshBasedOnVSARegions2D();
        for(GLKPOSITION Pos=vsaoper.m_meshpatch->GetEdgeList().GetHeadPosition(); Pos!=NULL; ) {
            curredge = (QMeshEdge *)(vsaoper.m_meshpatch->GetEdgeList().GetNext(Pos));
            simplifiedpatch->GetEdgeList().AddTail(curredge);
        }

        localcount = 0;
        for(GLKPOSITION Pos=vsaoper.m_meshpatch->GetNodeList().GetHeadPosition(); Pos!=NULL; ) {
            currnode = (QMeshNode *)(vsaoper.m_meshpatch->GetNodeList().GetNext(Pos));
            simplifiedpatch->GetNodeList().AddHead(currnode);

            //calculate normal for each vertex
            if(localcount==0)  {
                firstnode = currnode;
                startnode = firstnode;
                startnode->SetNormal(0.0, 0.0, 0.0);
            }
            if(localcount!=0) {
                endnode = currnode;

                startnode->GetCoord3D(startnodecoord[0], startnodecoord[1], startnodecoord[2]);
                endnode->GetCoord3D(endnodecoord[0], endnodecoord[1], endnodecoord[2]);
                edgevec[0] = endnodecoord[0]-startnodecoord[0];		edgevec[1] = endnodecoord[1]-startnodecoord[1];		edgevec[2] = endnodecoord[2]-startnodecoord[2];
                geo.VectorProduct(edgevec, dir, resultnormal);

                endnode->SetNormal(resultnormal[0], resultnormal[1], resultnormal[2]);
                startnode->GetNormal(tempnormal[0], tempnormal[1], tempnormal[2]);
                tempnormal[0] += resultnormal[0];  tempnormal[1] += resultnormal[1];  resultnormal[2] += resultnormal[2];
                startnode->SetNormal(tempnormal[0], tempnormal[1], tempnormal[2]);

                startnode = endnode;
            }
            localcount++;
        }
        if(vsaoper.m_meshpatch->GetNodeList().GetCount()>0)  {
            firstnode->GetNormal(tempnormal[0], tempnormal[1], tempnormal[2]);
            tempnormal[0] += resultnormal[0];  tempnormal[1] += resultnormal[1];  resultnormal[2] += resultnormal[2];
            firstnode->SetNormal(tempnormal[0], tempnormal[1], tempnormal[2]);
            startnode->SetNormal(tempnormal[0], tempnormal[1], tempnormal[2]);
        }

        vsaoper.m_meshpatch->GetEdgeList().RemoveAll();
        vsaoper.m_meshpatch->GetNodeList().RemoveAll();
        delete vsaoper.m_meshpatch;

        temppatch->GetEdgeList().RemoveAll();
        delete temppatch;
    }

    delete m_contour;
    m_contour = simplifiedpatch;
    localcount = 0;
    for(GLKPOSITION Pos=m_contour->GetEdgeList().GetHeadPosition(); Pos!=NULL; ) {
        curredge = (QMeshEdge *)(m_contour->GetEdgeList().GetNext(Pos));
        curredge->SetIndexNo(localcount);
        localcount++;
    }
    localcount = 0;
    for(GLKPOSITION Pos=m_contour->GetNodeList().GetHeadPosition(); Pos!=NULL; ) {
        currnode = (QMeshNode *)(m_contour->GetNodeList().GetNext(Pos));
        currnode->SetIndexNo(localcount);
        localcount++;
    }

    return maxdistterror;
}

void MarchingSquareBin::convert2XYPlane()
{
    double x, y, z;
    for(auto Pos=m_contour->GetNodeList().GetHeadPosition(); Pos!=NULL; ) {
        auto currnode = (QMeshNode *)(m_contour->GetNodeList().GetNext(Pos));
        currnode->GetCoord3D(x, y, z);
        std::swap(y, z);
        currnode->SetCoord3D(x, y, z);
    }
}

void MarchingSquareBin::writeContours(const std::string& path)
{
    std::shared_ptr<char> cstr(new char[path.length() + 1]);
    strcpy(cstr.get(), path.c_str());
    m_contour->outputOFFFile(cstr.get());
}

QMeshPatch *MarchingSquareBin::getContours(bool isDeleteData) {
    if (!isDeleteData) return m_contour;
    else {
        auto temp = m_contour;
        m_contour = NULL;
        return temp;
    }
}
