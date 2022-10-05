#define _CRT_SECURE_NO_DEPRECATE

#if defined (__APPLE__)
#include <GLUT/glut.h>
#include <sys/uio.h>
#include <dirent.h>
#else
#include <GL/glew.h>
#include <GL/GLAux.h>
#include <GL/glut.h>
#include <io.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <iostream>

#include "GLKLib/GLK.h"
#include "GLKLib/GLKCameraTool.h"

#include "QMeshLib/QMesh/QMeshPatch.h"
#include "QMeshLib/QSurfaceMesh.h"

#include "PMBody.h"
#include "MeshWorksDataBoard.h"

#include "LatticeModeler.h"
#include "call_cuda.h"
using namespace  std;

#define _MENU_QUIT						10001
#define _MENU_FILE_OPEN					10002
#define _MENU_FILE_SAVE					10003

#define _MENU_VIEW_ISOMETRIC			10101
#define _MENU_VIEW_FRONT				10102
#define _MENU_VIEW_BACK					10103
#define _MENU_VIEW_TOP					10104
#define _MENU_VIEW_BOTTOM				10105
#define _MENU_VIEW_LEFT					10106
#define _MENU_VIEW_RIGHT				10107
#define _MENU_VIEW_ORBITPAN				10108
#define _MENU_VIEW_ZOOMWINDOW			10109
#define _MENU_VIEW_ZOOMIN				10110
#define _MENU_VIEW_ZOOMOUT				10111
#define _MENU_VIEW_ZOOMALL				10112
#define _MENU_VIEW_PROFILE				10113
#define _MENU_VIEW_SHADE				10114
#define _MENU_VIEW_MESH					10115
#define _MENU_VIEW_AXIS					10116
#define _MENU_VIEW_COORD				10117
#define _MENU_VIEW_MESHSMOOTHSHADING	10118
#define _MENU_VIEW_SNAPSHOT				10121

#define _MENU_MESH_MAKECENTER			10201
#define _MENU_MESH_THICKENING			10202
#define _MENU_MESH_HOLEFILLING			10203
#define _MENU_MESH_UNIREMESHING			10204
#define _MENU_MESH_LATTICE              10205
#define _MENU_MESH_TPMS                 10206
#define _MENU_MESH_CLEARALL				10299

GLK _pGLK;
MeshWorksDataBoard _pDataBoard;
int _pMainWnd;

extern void menuEvent(int idCommand);

//#if defined (__APPLE__)
//#define CCL_DEFAULT_FOLDER_LOCATION     "Data/"
//#else
//#define CCL_DEFAULT_FOLDER_LOCATION     "D:\\data\\lattice\\Data\\"
/////////////////////////////////////////////////////////////////////////////////////////////////
////
//#ifdef _DEBUG
//#define _CRTDBG_MAP_ALLOC	// for memory-leak detection
//#include <stdlib.h>
//#include <crtdbg.h>
//#endif
////
/////////////////////////////////////////////////////////////////////////////////////////////////
//#endif

#define DEFAULT_FOLDER_LOCATION     "Data/"

#if defined (__APPLE__)
#define GET_REAL_PATH(X) realpath((X))     
#else
#include <cstdlib>
#define GET_REAL_PATH(X) _fullpath((char*)malloc(1024*sizeof(char)) ,(X), 1024)     
#endif 
#define CCL_DEFAULT_FOLDER_LOCATION  GET_REAL_PATH(DEFAULT_FOLDER_LOCATION)

void displayCoordinate(int x, int y)
{
    double wx,wy,wz;

    _pGLK.screen_to_wcl(x, y, wx, wy, wz);
    _pGLK.m_currentCoord[0]=(float)wx;
    _pGLK.m_currentCoord[1]=(float)wy;
    _pGLK.m_currentCoord[2]=(float)wz;

    //	printf("(%.2f, %.2f, %.2f)\n",(float)wx,(float)wy,(float)wz);

    _pGLK.refresh();
}

void specialKeyboardFunc(int key, int x, int y)
{
    pick_event pe;
    switch(_pGLK.m_mouseState) {
    case 1:{pe.nFlags=GLUT_LEFT_BUTTON;
    }break;
    case 2:{pe.nFlags=GLUT_MIDDLE_BUTTON;
    }break;
    case 3:{pe.nFlags=GLUT_RIGHT_BUTTON;
    }break;
    }
    pe.x=(double)x;	pe.y=(double)y;
    pe.nChar=-key;

    _pGLK.m_nModifier=0;
    switch(glutGetModifiers()) {
    case GLUT_ACTIVE_SHIFT:{_pGLK.m_nModifier=1;}break;
    case GLUT_ACTIVE_CTRL:{_pGLK.m_nModifier=2;	}break;
    case GLUT_ACTIVE_ALT:{_pGLK.m_nModifier=3;	}break;
    }

    if (_pGLK.GetCurrentTool()) _pGLK.GetCurrentTool()->process_event(KEY_PRESS,pe);
}

void keyboardFunc(unsigned char key, int x, int y)
{
    //------------------------------------------------------------------
    //	Hot Key Processing
    switch(key) {
    case 1:{	// ctrl+a
        menuEvent(_MENU_VIEW_ZOOMALL); return;
    }break;
    case 4:{	// ctrl+d
        //			menuEvent(_MENU_VIEW_GPUCPUPNTSDISP); return;
    }break;
    case 12:{	// ctrl+l
        menuEvent(_MENU_VIEW_MESHSMOOTHSHADING); return;
    }break;
    case 15:{	// ctrl+o
        menuEvent(_MENU_FILE_OPEN); return;
    }break;
    case 16:{	// ctrl+p
        //			menuEvent(_MENU_SPHS_STEPSIM); return;
    }break;
    case 18:{	// ctrl+r
        menuEvent(_MENU_VIEW_ORBITPAN); return;
    }break;
    case 19:{	// ctrl+s
        menuEvent(_MENU_FILE_SAVE); return;
    }break;
    case 23:{	// ctrl+w
        menuEvent(_MENU_VIEW_ZOOMWINDOW); return;
    }break;
    case 26:{	// ctrl+z
        menuEvent(_MENU_VIEW_SNAPSHOT); return;
    }break;
    }

    pick_event pe;
    switch(_pGLK.m_mouseState) {
    case 1:{pe.nFlags=GLUT_LEFT_BUTTON;
    }break;
    case 2:{pe.nFlags=GLUT_MIDDLE_BUTTON;
    }break;
    case 3:{pe.nFlags=GLUT_RIGHT_BUTTON;
    }break;
    }
    pe.x=(double)x;	pe.y=(double)y;
    pe.nChar=key;

    _pGLK.m_nModifier=0;
    switch(glutGetModifiers()) {
    case GLUT_ACTIVE_SHIFT:{_pGLK.m_nModifier=1;}break;
    case GLUT_ACTIVE_CTRL:{_pGLK.m_nModifier=2;	}break;
    case GLUT_ACTIVE_ALT:{_pGLK.m_nModifier=3;	}break;
    }

    if (_pGLK.GetCurrentTool()) _pGLK.GetCurrentTool()->process_event(KEY_PRESS,pe);
}

void motionFunc(int x, int y)
{
    if (_pGLK.m_mouseState==0) return;

    pick_event pe;
    switch(_pGLK.m_mouseState) {
    case 1:{pe.nFlags=GLUT_LEFT_BUTTON;
    }break;
    case 2:{pe.nFlags=GLUT_MIDDLE_BUTTON;
    }break;
    case 3:{pe.nFlags=GLUT_RIGHT_BUTTON;
    }break;
    }
    pe.x=(double)x;
    pe.y=(double)y;

    if (_pGLK.m_bCoordDisp) displayCoordinate(x,y);
    if (_pGLK.GetCurrentTool()) _pGLK.GetCurrentTool()->process_event(MOUSE_MOVE,pe);
}

void passiveMotionFunc(int x, int y)
{
    pick_event pe;
    pe.nFlags=-1;
    pe.x=(double)x;
    pe.y=(double)y;
    if (_pGLK.m_bCoordDisp) displayCoordinate(x,y);
    if (_pGLK.GetCurrentTool()) _pGLK.GetCurrentTool()->process_event(MOUSE_MOVE,pe);
}

void mouseFunc(int button, int state, int x, int y)
{
    if (state==GLUT_DOWN) {
        pick_event pe;
        _pGLK.m_nModifier=0;
        switch(glutGetModifiers()) {
        case GLUT_ACTIVE_SHIFT:{_pGLK.m_nModifier=1;}break;
        case GLUT_ACTIVE_CTRL:{_pGLK.m_nModifier=2;	}break;
        case GLUT_ACTIVE_ALT:{_pGLK.m_nModifier=3;	}break;
        }
        if (button==GLUT_LEFT_BUTTON) {pe.nFlags=GLUT_LEFT_BUTTON;_pGLK.m_mouseState=1;}
        if (button==GLUT_MIDDLE_BUTTON) {pe.nFlags=GLUT_MIDDLE_BUTTON;_pGLK.m_mouseState=2;}
        if (button==GLUT_RIGHT_BUTTON) {pe.nFlags=GLUT_RIGHT_BUTTON;_pGLK.m_mouseState=3;}
        pe.x=(double)x;
        pe.y=(double)y;
        if (_pGLK.GetCurrentTool()) _pGLK.GetCurrentTool()->process_event(MOUSE_BUTTON_DOWN,pe);
    }
    else if (state==GLUT_UP) {
        pick_event pe;
        _pGLK.m_nModifier=0;
        switch(glutGetModifiers()) {
        case GLUT_ACTIVE_SHIFT:{_pGLK.m_nModifier=1;}break;
        case GLUT_ACTIVE_CTRL:{_pGLK.m_nModifier=2;	}break;
        case GLUT_ACTIVE_ALT:{_pGLK.m_nModifier=3;	}break;
        }
        if (button==GLUT_LEFT_BUTTON) pe.nFlags=GLUT_LEFT_BUTTON;
        if (button==GLUT_MIDDLE_BUTTON) pe.nFlags=GLUT_MIDDLE_BUTTON;
        if (button==GLUT_RIGHT_BUTTON) pe.nFlags=GLUT_RIGHT_BUTTON;
        pe.x=(double)x;
        pe.y=(double)y;
        if (_pGLK.GetCurrentTool()) _pGLK.GetCurrentTool()->process_event(MOUSE_BUTTON_UP,pe);

        _pGLK.m_mouseState=0;
    }
}

void animationFunc()
{
    /*	if (_pDataBoard.m_vdFieldCudaBody) {
        int activeSlide=_pDataBoard.m_vdFieldCudaBody->GetActiveSlice();
        activeSlide++;
        _pDataBoard.m_vdFieldCudaBody->SetActiveSlice(activeSlide);
        _pDataBoard.m_vdFieldCudaBody->BuildGLList();
        _pGLK.refresh();
    }*/
}

void visibleFunc(int visible)
{
    if (visible==GLUT_VISIBLE)
        glutIdleFunc(animationFunc);
    //		glutIdleFunc(NULL);
    else
        glutIdleFunc(NULL);
}

void displayFunc(void)
{
    _pGLK.refresh();
}

void reshapeFunc(int w, int h) 
{
    _pGLK.Reshape(w,h);
}

void initFunc()
{
    _pGLK.Initialization();

    //	_pGLK.SetAxisDisplay(false);
    _pGLK.SetMesh(false);
    //	_pGLK.SetClearColor(1.0,1.0,1.0);

    GLKCameraTool *myTool=new GLKCameraTool(&_pGLK,ORBITPAN);
    _pGLK.clear_tools();
    _pGLK.set_tool(myTool);
}

#if defined (__APPLE__)
bool fileChosenByList(char directorty[], char selectedFileName[])
{
    DIR *dirp;      struct dirent *dp;      int fileNum=0;
    int colNum=3,colsize;
    
    //--------------------------------------------------------------------------------------
    //  The following lines list out all the files in the folder
    colsize=80/colNum-4;
    if ((dirp = opendir(directorty)) == NULL) {
        printf("Error: couldn't open '%s'\n",directorty);   return false;
    }
    do{
        if ((dp = readdir(dirp)) != NULL) {
            printf( "%*d: %s %*s", 2, fileNum++, dp->d_name, colsize-strlen(dp->d_name), " ");
            if ((fileNum%colNum)==0) printf("\n");
        }
    }while(dp!=NULL);
    closedir(dirp);
    
    //--------------------------------------------------------------------------------------
    //  The following lines select the file according to user input
    int inputNum;
    char inputStr[200];
    printf("\nPlease select the file name for import: ");
    //    scanf("%s",inputStr);	printf("\n"); sscanf(inputStr,"%d",&inputNum);
    inputNum = 5;
    if (inputNum<0 || inputNum>=fileNum) {
        printf("Incorrect Input!!!\n");
        return false;
    }
    fileNum=0;
    dirp = opendir(directorty);
    while((dp = readdir(dirp)) != NULL) {
        if (fileNum==inputNum) {
            strcpy(selectedFileName,dp->d_name);
            break;
        }
        fileNum++;
    }
    closedir(dirp);
    
    printf("----------------------------------------\nSelected File: %s\n",selectedFileName);
    return true;
}

bool isFileExist(char directorty[], char filename[])
{
    DIR *dirp;      struct dirent *dp;
    
    if ((dirp = opendir(directorty)) == NULL) {
        printf("Error: couldn't open '.'\n");
        return false;
    }
    while((dp = readdir(dirp)) != NULL) {
        if (strcmp(dp->d_name, filename) == 0) {
            closedir(dirp);
            return true;
        }
    }
    closedir(dirp);
    
    return false;
}

#else

bool fileChosenByList(char directorty[], char selectedFileName[])
{
    struct _finddata_t c_file;
	intptr_t hFile;
    long fileNum=0;
    char filespef[200];
    int colNum=3,colsize;
    
    colsize=80/colNum-4;
    strcpy(filespef,directorty);
    strcat(filespef,"*.*");

    hFile = _findfirst( filespef, &c_file );
    if( (hFile) == -1 ) {
        // printf( "No file is found!\n");
        printf("Error: couldn't open '%s'\n",directorty);
        return false;
    }

	std::cout << "dir: " << directorty << std::endl;

	printf( "%*d: %s %*s", 2, fileNum++, c_file.name, colsize-strlen(c_file.name), " ");
    while(_findnext( hFile, &c_file )!=-1) {
        printf( "%*d: %s %*s", 2, fileNum++, c_file.name, colsize-strlen(c_file.name), " ");
        if ((fileNum%colNum)==0) printf("\n");
    }
    _findclose(hFile);

    int inputNum;	char inputStr[200];
    printf("\nPlease select the file name for import: ");
    scanf("%s",inputStr);	printf("\n");	sscanf(inputStr,"%d",&inputNum);
    if (inputNum<0 || inputNum>=fileNum) {printf("Incorrect Input!!!\n"); return false;}
    
    fileNum=0;
    if( (hFile = _findfirst( filespef, &c_file )) == -1 ) {return false;}
    if (inputNum!=0) {
        fileNum++;
        while(_findnext( hFile, &c_file )!=-1) {
            if (fileNum==inputNum) break;
            fileNum++;
        }
    }
    _findclose(hFile);
    strcpy(selectedFileName,c_file.name);
    
    printf("----------------------------------------\nSelected File: %s\n",selectedFileName);
    return true;
}

bool isFileExist(char dir[], char filename[])
{
    char fullfilename[1024];
    sprintf(fullfilename,"%s%s",dir,filename);
    
    struct _finddata_t c_file;
    long hFile;
    if( (hFile = _findfirst( fullfilename, &c_file )) == -1L ) {
        return false;
    }
    return true;
}
#endif


//---------------------------------------------------------------------------------
//	The following functions are for menu processing
void menuFuncFileSave()
{
    char filename[1024],exstr[4],name[256],directory[256],answer[10];
    
    strcpy(directory,CCL_DEFAULT_FOLDER_LOCATION);
    printf("\nPlease specify the file name for export: ");
    scanf("%s",name);    printf("\n");
    
    if (isFileExist(directory,name)) {
        printf( "The file - %s has been found, do you want to overwite it? (y/n)\n", name);
        scanf("%s",answer);
        if (answer[0]!='y' && answer[0]!='Y') return;
    }
    strcpy(filename,directory);	strcat(filename,name);
    
    int length=(int)(strlen(filename));
    exstr[0]=filename[length-3];
    exstr[1]=filename[length-2];
    exstr[2]=filename[length-1];
    exstr[3]='\0';
    
    if (strcmp(exstr, "obj") == 0) {	//.obj file
        _pDataBoard.OutputOBJFile(filename);
    }
    else if (strcmp(exstr, "wio") == 0) {	//.wio file
        _pDataBoard.OutputWIOFile(filename);
    }
    else {
        printf("Warning: incorrect file extension, no file is saved!\n");
    }
}

void menuFuncFileImageSnapShot()
{
    //	int sx,sy;	_pGLK.GetSize(sx,sy);
    //	GLKAVIGenerator::SnapShot(sx, sy, "Data/Snapshot.bmp");
}

void menuFuncFileOpen()
{
    char filename[1024],exstr[4],name[256],directory[256];
    
    strcpy(directory,CCL_DEFAULT_FOLDER_LOCATION);
    if (!fileChosenByList(directory,name)) return;
    if (!isFileExist(directory, name)) {printf("The file - %s is not found!\n", name); return;}
    strcpy(filename,directory);	strcat(filename,name);
    
    int length=(int)(strlen(filename));
    exstr[0]=filename[length-3];
    exstr[1]=filename[length-2];
    exstr[2]=filename[length-1];
    exstr[3]='\0';

    if (strcmp(exstr, "obj") == 0) {	//.obj file
        if (!(_pDataBoard.m_polyMeshBody))
            _pDataBoard.m_polyMeshBody = new PMBody;
        else
            _pGLK.DelDisplayObj2(_pDataBoard.m_polyMeshBody);
        _pDataBoard.InputOBJFile(filename);
        long time = clock();
        _pDataBoard.m_polyMeshBody->BuildGLList(_pDataBoard.m_bVertexNormalShading);
        printf("Build GL List Time (ms): %ld\n", clock() - time); time = clock();
        printf("--------------------------------------------\n");
        _pGLK.AddDisplayObj(_pDataBoard.m_polyMeshBody, true);
    }
    else if (strcmp(exstr, "wio") == 0) {	//.wio file
        if (!(_pDataBoard.m_polyMeshBody))
            _pDataBoard.m_polyMeshBody = new PMBody;
        else
            _pGLK.DelDisplayObj2(_pDataBoard.m_polyMeshBody);
        _pDataBoard.InputWIOFile(filename);
        long time = clock();
        _pDataBoard.m_polyMeshBody->BuildGLList(_pDataBoard.m_bVertexNormalShading);
        printf("Build GL List Time (ms): %ld\n", clock() - time); time = clock();
        printf("--------------------------------------------\n");
        _pGLK.AddDisplayObj(_pDataBoard.m_polyMeshBody, true);
    }
    else if (strcmp(exstr, "off") == 0) {	//.off file
        if (!(_pDataBoard.m_polyMeshBody))
            _pDataBoard.m_polyMeshBody = new PMBody;
        else
            _pGLK.DelDisplayObj2(_pDataBoard.m_polyMeshBody);
        _pDataBoard.InputOFFFile(filename);
        long time = clock();
        _pDataBoard.m_polyMeshBody->BuildGLList(_pDataBoard.m_bVertexNormalShading);
        printf("Build GL List Time (ms): %ld\n", clock() - time); time = clock();
        printf("--------------------------------------------\n");
        _pGLK.AddDisplayObj(_pDataBoard.m_polyMeshBody, true);
    }

    _pGLK.zoom_all_in_view();
}

void menuFuncQuit()
{
    exit(0);
}

void menuFuncMeshMakeCenter()
{
    if (!(_pDataBoard.m_polyMeshBody)) { printf("None mesh is found!\n");	return; }

    double bndBox[6];	long time = clock();
    _pDataBoard.m_polyMeshBody->CompBoundingBox(bndBox);
    _pDataBoard.m_polyMeshBody->Transformation(-(bndBox[0]+bndBox[1])*0.5,-(bndBox[2]+bndBox[3])*0.5, -(bndBox[4]+bndBox[5])*0.5);
    //printf("Total Processing Time (ms): %ld\n\n", clock() - time);
    _pDataBoard.m_polyMeshBody->BuildGLList(_pDataBoard.m_bVertexNormalShading);
    _pGLK.refresh();
}

void menuFuncMeshThickening()
{
    //	if (!(_pDataBoard.m_polyMeshBody)) { printf("None mesh is found!\n");	return; }
    //	char buf[20];	float offset;	double avgLength;	int nRes;
    //    QSurfaceMesh *mesh = (QSurfaceMesh *)(_pDataBoard.m_polyMeshBody->GetMeshList().GetTail());
    //    if (mesh == NULL) { printf("None mesh is found!\n");	return; }

    //	avgLength = _pDataBoard.m_polyMeshBody->CompAverageEdgeLength();
    //	printf("Offset (in terms of average triangular edge length - %lf) = ?", avgLength);
    //	scanf("%s", buf);
    //	printf("%s\n", buf);
    //	sscanf(buf, "%f", &offset);
    //	//--------------------------------------------------------------------------------------------
    //	printf("\nPlease specify the resolution for sampling: ");
    //	scanf("%s", buf);	printf("\n");	sscanf(buf, "%d", &nRes);
    //	if (nRes <= 0) { printf("Incorrect InputL: %d!!!\n", nRes); return; }
    //	printf("Sampling Resolution: %d\n", nRes);

    //	ThickeningMeshPatch meshThickening;
    //    QSurfaceMesh *offsetMesh;
    //	long time = clock();
    //    meshThickening.DoPatchThickening(mesh, (double)offset*avgLength, nRes, offsetMesh);
    //	printf("Total Processing Time (ms): %ld\n\n", clock() - time);

    //	_pGLK.DelDisplayObj2(_pDataBoard.m_polyMeshBody);
    //	_pDataBoard.m_polyMeshBody->computeRange();
    //	_pDataBoard.m_polyMeshBody->BuildGLList(_pDataBoard.m_bVertexNormalShading);
    //	_pGLK.AddDisplayObj(_pDataBoard.m_polyMeshBody, true);
}

void menuFuncLatticeModeling()
{
    char filename[1024],name[256],directory[256];

    strcpy(directory,CCL_DEFAULT_FOLDER_LOCATION);
    if (!fileChosenByList(directory,name)) return;
    if (!isFileExist(directory, name)) {printf("The file - %s is not found!\n", name); return;}
    strcpy(filename,directory);	strcat(filename,name);

    LatticeModeler modeler;
    modeler.sliceModelImplicitStream(std::string(filename));
}


void menuFuncMeshHoleFilling()
{
    //    if (!(_pDataBoard.m_polyMeshBody)) { printf("None mesh is found!\n");	return; }
    //    char buf[20];	float offset;	double avgLength;	int nRes;
    //    QSurfaceMesh *mesh = (QSurfaceMesh *)(_pDataBoard.m_polyMeshBody->GetMeshList().GetTail());
    //    if (mesh == NULL) { printf("None mesh is found!\n");	return; }
    
    //    ThickeningMeshPatch meshThickening;
    //    long time = clock();
    //    meshThickening.DoHoleTriangulation(mesh);
    //    printf("Total Processing Time (ms): %ld\n\n", clock() - time);
    
    //    _pGLK.DelDisplayObj2(_pDataBoard.m_polyMeshBody);
    //    _pDataBoard.m_polyMeshBody->computeRange();
    //    _pDataBoard.m_polyMeshBody->BuildGLList(_pDataBoard.m_bVertexNormalShading);
    //    _pGLK.AddDisplayObj(_pDataBoard.m_polyMeshBody, true);
}

void menuEvent(int idCommand)
{
    switch (idCommand) {
    case _MENU_QUIT:menuFuncQuit();
        break;

        //--------------------------------------------------------------------
        //	File related
    case _MENU_FILE_OPEN:menuFuncFileOpen();
        break;
    case _MENU_FILE_SAVE:menuFuncFileSave();
        break;

        //--------------------------------------------------------------------
        //	View related
    case _MENU_VIEW_ISOMETRIC:_pGLK.SetViewDirection(VD_ISOMETRICVIEW);
        break;
    case _MENU_VIEW_FRONT:_pGLK.SetViewDirection(VD_FRONTVIEW);
        break;
    case _MENU_VIEW_BACK:_pGLK.SetViewDirection(VD_BACKVIEW);
        break;
    case _MENU_VIEW_TOP:_pGLK.SetViewDirection(VD_TOPVIEW);
        break;
    case _MENU_VIEW_BOTTOM:_pGLK.SetViewDirection(VD_BOTTOMVIEW);
        break;
    case _MENU_VIEW_LEFT:_pGLK.SetViewDirection(VD_LEFTVIEW);
        break;
    case _MENU_VIEW_RIGHT:_pGLK.SetViewDirection(VD_RIGHTVIEW);
        break;
    case _MENU_VIEW_ORBITPAN:{
        GLKCameraTool *myTool = new GLKCameraTool(&_pGLK, ORBITPAN);
        _pGLK.clear_tools();
        _pGLK.set_tool(myTool);
    }break;
    case _MENU_VIEW_ZOOMWINDOW:{
        GLKCameraTool *myTool = new GLKCameraTool(&_pGLK, ZOOMWINDOW);
        _pGLK.clear_tools();
        _pGLK.set_tool(myTool);
    }break;
    case _MENU_VIEW_ZOOMIN:_pGLK.zoom(1.5);
        break;
    case _MENU_VIEW_ZOOMOUT:_pGLK.zoom(0.75);
        break;
    case _MENU_VIEW_ZOOMALL:_pGLK.zoom_all_in_view();
        break;
    case _MENU_VIEW_PROFILE:{_pGLK.SetProfile(!(_pGLK.GetProfile())); _pGLK.refresh();
    }break;
    case _MENU_VIEW_SHADE:{
        _pGLK.SetShading(!(_pGLK.GetShading()));
        if (_pGLK.GetShading()) {
            //			if (_pDataBoard.m_nurbsSurfBody!=NULL) _pDataBoard.m_nurbsSurfBody->BuildGLList(true);
            //			if (_pDataBoard.m_polyMeshBody!=NULL) _pDataBoard.m_polyMeshBody->BuildGLList(true);
        }
        _pGLK.refresh();
    }break;
    case _MENU_VIEW_MESH:{
        _pGLK.SetMesh(!(_pGLK.GetMesh()));
        if (_pGLK.GetMesh()) {
            //			if (_pDataBoard.m_nurbsSurfBody!=NULL) _pDataBoard.m_nurbsSurfBody->BuildGLList(false);
            //			if (_pDataBoard.m_polyMeshBody!=NULL) _pDataBoard.m_polyMeshBody->BuildGLList(false);
        }
        _pGLK.refresh();
    }break;
    case _MENU_VIEW_AXIS:{_pGLK.SetAxisDisplay(!(_pGLK.GetAxisDisplay())); _pGLK.refresh();
    }break;
    case _MENU_VIEW_COORD:{_pGLK.m_bCoordDisp = !(_pGLK.m_bCoordDisp); _pGLK.refresh();
    }break;
    case _MENU_VIEW_MESHSMOOTHSHADING:{
        _pDataBoard.m_bVertexNormalShading = !(_pDataBoard.m_bVertexNormalShading);
        if (_pDataBoard.m_polyMeshBody) {
            _pDataBoard.m_polyMeshBody->BuildGLList(_pDataBoard.m_bVertexNormalShading);
            _pGLK.refresh();
        }
    }break;
    case _MENU_VIEW_SNAPSHOT:{menuFuncFileImageSnapShot();
    }break;

        //--------------------------------------------------------------------
        //	Mesh-Processing Function related
    case _MENU_MESH_MAKECENTER:{menuFuncMeshMakeCenter();
    }break;
    case _MENU_MESH_THICKENING:{menuFuncMeshThickening();
    }break;
    case _MENU_MESH_LATTICE:{menuFuncLatticeModeling();
    }break;
    case _MENU_MESH_HOLEFILLING:{menuFuncMeshHoleFilling();
    }break;
    case _MENU_MESH_CLEARALL:{
        char answer[20];
        printf("Warning: are you sure that you are going to clear all mesh objects? (y/n)\n");	scanf("%s", answer);
        if (answer[0] == 'y' || answer[0] == 'Y') { _pGLK.DelDisplayObj(_pDataBoard.m_polyMeshBody);	_pDataBoard.m_polyMeshBody = NULL; }
    }break;

        //--------------------------------------------------------------------
        //	Other Functions

    }
}

int buildPopupMenu (void)
{
    int mainMenu,fileSubMenu,viewSubMenu,meshSubMenu;

    fileSubMenu = glutCreateMenu(menuEvent);
    glutAddMenuEntry("Open\tCtrl+O", _MENU_FILE_OPEN);
    glutAddMenuEntry("Save\tCtrl+S", _MENU_FILE_SAVE);

    viewSubMenu = glutCreateMenu(menuEvent);
    glutAddMenuEntry("Isometric", _MENU_VIEW_ISOMETRIC);
    glutAddMenuEntry("Front", _MENU_VIEW_FRONT);
    glutAddMenuEntry("Back", _MENU_VIEW_BACK);
    glutAddMenuEntry("Top", _MENU_VIEW_TOP);
    glutAddMenuEntry("Bottom", _MENU_VIEW_BOTTOM);
    glutAddMenuEntry("Left", _MENU_VIEW_LEFT);
    glutAddMenuEntry("Right", _MENU_VIEW_RIGHT);
    glutAddMenuEntry("----",-1);
    glutAddMenuEntry("Orbot and Pan\tCtrl+R",_MENU_VIEW_ORBITPAN);
    glutAddMenuEntry("Zoom Window\tCtrl+W",_MENU_VIEW_ZOOMWINDOW);
    glutAddMenuEntry("Zoom In",_MENU_VIEW_ZOOMIN);
    glutAddMenuEntry("Zoom Out",_MENU_VIEW_ZOOMOUT);
    glutAddMenuEntry("Zoom All\tCtrl+A",_MENU_VIEW_ZOOMALL);
    glutAddMenuEntry("----",-1);
    glutAddMenuEntry("Profile",_MENU_VIEW_PROFILE);
    glutAddMenuEntry("Shade",_MENU_VIEW_SHADE);
    glutAddMenuEntry("Mesh",_MENU_VIEW_MESH);
    glutAddMenuEntry("----",-1);
    glutAddMenuEntry("Axis Frame",_MENU_VIEW_AXIS);
    glutAddMenuEntry("Coordinate",_MENU_VIEW_COORD);
    glutAddMenuEntry("----",-1);
    glutAddMenuEntry("Smooth Shading\tCtrl+L", _MENU_VIEW_MESHSMOOTHSHADING);
    glutAddMenuEntry("Image Snap Shot\tCtrl+Z", _MENU_VIEW_SNAPSHOT);

    meshSubMenu = glutCreateMenu(menuEvent);
    glutAddMenuEntry("Lattice Modeling", _MENU_MESH_LATTICE);
    glutAddMenuEntry("----", -1);
    glutAddMenuEntry("Thickening", _MENU_MESH_THICKENING);
    glutAddMenuEntry("----", -1);
    glutAddMenuEntry("Hole Filling", _MENU_MESH_HOLEFILLING);
    glutAddMenuEntry("Uniform Remeshing", _MENU_MESH_UNIREMESHING);
    glutAddMenuEntry("----", -1);
    glutAddMenuEntry("Make Centralized", _MENU_MESH_MAKECENTER);
    glutAddMenuEntry("----", -1);
    glutAddMenuEntry("Clear All Meshes", _MENU_MESH_CLEARALL);

    mainMenu = glutCreateMenu(menuEvent);
    glutAddSubMenu("File", fileSubMenu);
    glutAddSubMenu("View", viewSubMenu);
    glutAddSubMenu("Mesh", meshSubMenu);
    glutAddMenuEntry("----",-1);
    glutAddMenuEntry("Quit", _MENU_QUIT);

    return mainMenu;
}

//---------------------------------------------------------------------------------
//	The major function of a program
int main(int argc, char *argv[])
{
    glutInit(&argc, argv);
    //    glutInitDisplayMode(GLUT_DEPTH | GLUT_RGB | GLUT_DOUBLE | GLUT_MULTISAMPLE | GLUT_STENCIL);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_RGBA | GLUT_ALPHA | GLUT_DOUBLE | GLUT_STENCIL);

    _pMainWnd=glutCreateWindow("MeshWorks ver 0.1");
    glutDisplayFunc(displayFunc);
    glutReshapeWindow(1000, 750);
    glutMouseFunc(mouseFunc);
    glutMotionFunc(motionFunc);
    glutPassiveMotionFunc(passiveMotionFunc);
    glutKeyboardFunc(keyboardFunc);
    glutSpecialFunc(specialKeyboardFunc);
    glutReshapeFunc(reshapeFunc);
    glutVisibilityFunc(visibleFunc);

    initFunc();
    _pGLK.SetClearColor(0.35f,0.35f,0.35f);
    //	_pGLK.SetClearColor(1.0f,1.0f,1.0f);
    _pGLK.SetForegroundColor(1.0f,1.0f,1.0f);
    _pGLK.m_bCoordDisp=false;

    //	_pGLK.SetProfile(false);

    displayFunc();
    
#if defined (__APPLE__)
#else
    if(glewInit() != GLEW_OK) {
        printf("glewInit failed. Exiting...\n");
        return false;
    }
    if (glewIsSupported("GL_VERSION_2_0")) {
        printf("\nReady for OpenGL 2.0\n");
        printf("-------------------------------------------------\n");
        printf("GLSL will be used to speed up sampling\n");
    }
    else {
        printf("OpenGL 2.0 not supported\n");
        return false;
    }
#endif
    
    printf("PntWorks Started\n");
    printf("--------------------------------------------------\n");
    printf("Please select the following functions by hot-keys:\n\n");
    printf("Ctrl - O      Open\n");
    printf("Ctrl - R      Orbit and Pan\n");
    printf("Ctrl - W      Zoom Window\n");
    printf("--------------------------------------------------\n");

    buildPopupMenu();
    glutAttachMenu(GLUT_RIGHT_BUTTON);

    glutSwapBuffers();
    glutMainLoop();

    return 0;             /* ANSI C requires main to return int. */
}
