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
#endif // CALL_CUDA_H
