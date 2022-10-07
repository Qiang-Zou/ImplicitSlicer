#ifndef CALL_CUDA_H
#define CALL_CUDA_H
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include "malloc.h"
extern "C" void call_krFDMContouring_SubtractSolidRegion(bool *gridNodes, bool *suptNodes, bool *seedNodes, int nodeNum, int3 imageRes, int iy);
extern "C" void call_krFDMContouring_CopyNodesrom3Dto2D(bool *m_2DNodes, bool *m_3DNodes, int nodeNum, int3 imageRes, int iy);
extern "C" void call_krFDMContouring_Filter4(bool* outNodes, bool *inNodesA, bool *inNodesB, int nodeNum, int3 imageRes, int iy);
extern "C" void call_krFDMContouring_integrateImageintoGrossImage(bool* outNodes, bool *inNodes, int nodeNum, int3 imageRes, int iy);
extern "C" void call_krFDMContouring_CopyNodesrom2Dto3D(bool *m_2DNodes, bool *m_3DNodes, int nodeNum, int3 imageRes, int iy);
extern "C" void call_krFDMContouring_Filter5(bool* outNodes, bool *inNodesA, bool *gridNodes, int nodeNum, int iy, int3 imageRes);

extern "C" void call_krSLAContouring_Initialization(bool *tempImg, bool *targetImg, bool *gridNodes, int nodeNum, int3 imageRes, int iy);
extern "C" void call_krSLAContouring_Filter1(bool *assistImg, bool *tempImg, bool *targetImg, int nodeNum, int2 suptRes, int suptRadius, int iy);
extern "C" void call_krSLAContouring_OrthoSearchRemainAnchorZ(bool *assistImg, bool *tempImg, bool *targetImg, int nodeNum, int2 indRes, int2 imgRes, int iy);
extern "C" void call_krSLAContouring_OrthoSearchRemainAnchorX(bool *assistImg, bool *tempImg, bool *targetImg, int nodeNum, int2 indRes, int2 imgRes, int iy);
extern "C" void call_krSLAContouring_Filter5(bool *gridNodes, bool *tempImg, bool *suptNodes, unsigned int *linkIndex, int nodeNum, int3 imageRes, int iy);
extern "C" void call_krSLAContouring_FindAllLinks(unsigned int *linkIndex, unsigned int *linkLayerC, short *linkLayerD, short2 *linkID, bool *tempImg, unsigned int *count, int nodeNum, int3 imageRes);
extern "C" void call_krSLAContouring_RelateAllLinksBetweenLayers(unsigned int *linkIndex, unsigned int *linkLayerC, bool *gridNodes, bool *suptNodes, int nodeNum, int3 imageRes);
extern "C" void call_krFDMContouring_Dilation(bool *seedNodes, bool* output, int nodeNum, int3 imageRes, double realThreshold, int gridRadius, int iy);
extern "C" void call_krFDMContouring_VerticalSpptPxlProp(bool *gridNodes, bool *suptNodes, int nodeNum, int3 imageRes, int iy);

class call_func
{
public:
	call_func(void) {};
	~call_func(void) {};

	static void setdev_ptr(unsigned int* LinkIndex, int nodeNum, unsigned int& LinkNum);
	static void setdev_ptr2(unsigned int* linkLayerC, unsigned int LinkNum, int3 imgRes);

};

#endif // CALL_CUDA_H
