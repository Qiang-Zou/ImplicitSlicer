#include "LDNIcudaOperation.h"
#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include <sys/stat.h>

#include "cuda.h"
#include "cutil.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/remove.h>
#include <thrust/sequence.h>
#include <thrust/fill.h>
#include <thrust/count.h>
#include <thrust/functional.h>
#include <thrust/sequence.h>
#include <thrust/iterator/reverse_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/unique.h>
#include <cstdlib>

#include "GLKLib/GLKGeometry.h"


#include "LDNIcudaSolid.h"
#include "call_cuda.h"


extern __global__ void krFDMContouring_SubtractSolidRegion(bool *gridNodes, bool *suptNodes, bool *seedNodes, int nodeNum, int3 imageRes, int iy);
extern __global__ void krFDMContouring_CopyNodesrom3Dto2D(bool *m_2DNodes, bool *m_3DNodes, int nodeNum, int3 imageRes, int iy);
extern __global__ void krFDMContouring_Filter4(bool* outNodes, bool *inNodesA, bool *inNodesB, int nodeNum, int3 imageRes, int iy);
extern __device__ float interpoint(int x1, int y1, int x2, int y2, int x0);
extern __global__ void krFDMContouring_integrateImageintoGrossImage(bool* outNodes, bool *inNodes, int nodeNum, int3 imageRes, int iy);
extern __global__ void krFDMContouring_CopyNodesrom2Dto3D(bool *m_2DNodes, bool *m_3DNodes, int nodeNum, int3 imageRes, int iy);
extern __global__ void krFDMContouring_Filter5(bool* outNodes, bool *inNodesA, bool *gridNodes, int nodeNum, int iy, int3 imageRes);
extern __global__ void krFDMContouring_Dilation(bool *seedNodes, bool* output, int nodeNum, int3 imageRes, double realThreshold, int gridRadius, int iy);
extern __global__ void krFDMContouring_Erosion(bool *gridNodes, bool* output, int nodeNum, int3 imageRes, double realThreshold, int gridRadius);
extern __global__ void krSLAContouring_Initialization(bool *tempImg, bool *targetImg, bool *gridNodes, int nodeNum, int3 imageRes, int iy);
extern __global__ void krSLAContouring_Filter1(bool *assistImg, bool *tempImg, bool *targetImg, int nodeNum, int2 suptRes, int suptRadius, int iy);
extern __global__ void krSLAContouring_OrthoSearchRemainAnchorZ(bool *assistImg, bool *tempImg, bool *targetImg, int nodeNum, int2 indRes, int2 imgRes, int iy);
extern __global__ void krSLAContouring_OrthoSearchRemainAnchorX(bool *assistImg, bool *tempImg, bool *targetImg, int nodeNum, int2 indRes, int2 imgRes, int iy);
extern __global__ void krSLAContouring_FillAnchorValue(bool* seedNodes, bool* inNodes, unsigned int* _value, unsigned int marker, int nodeNum, int iy);
extern __global__ void krSLAContouring_GetAnchorPoint(bool* targetImg, unsigned int* _value, unsigned int* anchorPt, int2 imgRes, int nodeNum, int iy);
extern __global__ void krSLAContouring_FillImageValue(bool* tempImg, bool* suptImg, unsigned int* _value, int2 imgRes, int nodeNum, unsigned int init, int iy);
extern __global__ void krSLAContouring_Filter5(bool *gridNodes, bool *tempImg, bool *suptNodes, unsigned int *linkIndex, int nodeNum, int3 imageRes, int iy);
extern __global__ void krSLAContouring_InitializedDistanceMap(bool *seedNodes, short2 *dMap, int nodeNum, int2 imageRes, int iy, int texwidth);
extern __global__ void krSLAContouring_FindAllLinks(unsigned int *linkIndex, unsigned int *linkLayerC, short *linkLayerD, short2 *linkID, bool *tempImg, unsigned int *count, int nodeNum, int3 imageRes);
extern __global__ void krSLAContouring_RelateAllLinksBetweenLayers(unsigned int *linkIndex, unsigned int *linkLayerC, bool *gridNodes, bool *suptNodes, int nodeNum, int3 imageRes);
extern __global__ void krFDMContouring_VerticalSpptPxlProp(bool *gridNodes, bool *suptNodes, int nodeNum, int3 imageRes, int iy);
#define MARKER      -32768
#define BLOCKSIZE	64
#define TILE_DIM	32
#define BLOCK_ROWS	8
texture<short2> disTexColor;
texture<short2> disTexLinks;
#define TOID(x, y, size)    (__mul24((y), (size)) + (x))
struct transformOp : public thrust::unary_function<unsigned int, unsigned int>
{
	unsigned int w, h;

	__host__ __device__
		transformOp(unsigned int _w, unsigned int _h) : w(_w), h(_h) {}


	__host__ __device__
		unsigned int operator()(unsigned int index)
	{
		unsigned int m, n;

		n = index % w;
		m = index / w;


		return (n*h + m);
	}
};

__global__ void krFDMContouring_InitializedValue(short2 *dMap, int nodeNum, int value)
{
    int index=threadIdx.x+blockIdx.x*blockDim.x;

    while (index<nodeNum) {

        dMap[index] = make_short2(value, value);

        index += blockDim.x * gridDim.x;
    }
}
__global__ void krFDMContouring_InitializedDistanceMap(bool *gridNodes, short2 *dMap, int nodeNum, int realNum, int3 imageRes, int iy, int texwidth)
{
    int index=threadIdx.x+blockIdx.x*blockDim.x;
    unsigned int ix,iz, id;
    while (index<nodeNum) {

        //if (index < realNum)
        //{

            ix = index%imageRes.x;	iz = index/imageRes.x;
            id = iz*texwidth + ix;
            //if (index == 6911)
            //	printf("ix iz %d %d %d\n", ix, iz, id);
            if (gridNodes[iz*imageRes.x*imageRes.y+iy*imageRes.x+ix])
            {

                dMap[id].x = ix;
                dMap[id].y = iz;
                //dMap[index].x = MARKER;
                //dMap[index].y = MARKER;
                //if (iy == 137)
                //	printf("ix iy iz : %d %d %d\n", ix, iy, iz);

            }
            /*else
            {
                dMap[index].x = MARKER;
                dMap[index].y = MARKER;

            }*/
        /*}
        else
        {
            dMap[index].x = MARKER;
            dMap[index].y = MARKER;
        }*/


        index += blockDim.x * gridDim.x;
    }
}
__global__ void krFDMContouring_kernelFloodDown(short2 *output, int size, int bandSize)
{
    int tx = blockIdx.x * blockDim.x + threadIdx.x;
    int ty = blockIdx.y * bandSize;
    int id = TOID(tx, ty, size);

    short2 pixel1, pixel2;

    pixel1 = make_short2(MARKER, MARKER);

    for (int i = 0; i < bandSize; i++, id += size) {

        //if (id >= imageSize) break;
        pixel2 = tex1Dfetch(disTexColor, id);




        if (pixel2.x != MARKER)
            pixel1 = pixel2;

        //if (TOID(tx, ty, size) == 4351)
        //if (id == 6911)
        //	printf("flood down %d %d %d %d\n", pixel1.x, pixel1.y,  pixel2.x, pixel2.y);

        output[id] = pixel1;


    }
}
__global__ void krFDMContouring_kernelFloodUp(short2 *output, int size, int bandSize)
{
    int tx = blockIdx.x * blockDim.x + threadIdx.x;
    int ty = (blockIdx.y+1) * bandSize - 1;
    int id = TOID(tx, ty, size);

    short2 pixel1, pixel2;
    int dist1, dist2;

    pixel1 = make_short2(MARKER, MARKER);

    for (int i = 0; i < bandSize; i++, id -= size) {
        //if (id >= imageSize) continue;
        //if (TOID(tx, ty, size) == 8191)
        //	printf("111 %d %d \n", pixel1.x, pixel1.y);

        dist1 = abs(pixel1.y - ty + i);

        //if (TOID(tx, ty, size) == 8191)
        //	printf("222 %d \n", dist1);

        pixel2 = tex1Dfetch(disTexColor, id);
        dist2 = abs(pixel2.y - ty + i);

        if (dist2 < dist1)
            pixel1 = pixel2;


        //if (TOID(tx, ty, size) == 8191)
        //		printf("aaa %d %d %d %d %d %d %d %d %d\n", pixel2.x, pixel2.y, dist1, dist2, pixel1.x, pixel1.y, ty, i, id);


        output[id] = pixel1;
    }
}
__global__ void krFDMContouring_kernelPropagateInterband(short2 *output, int size, int bandSize)
{
    int tx = blockIdx.x * blockDim.x + threadIdx.x;
    int inc = __mul24(bandSize, size);
    int ny, nid, nDist;
    short2 pixel;

    // Top row, look backward
    int ty = __mul24(blockIdx.y, bandSize);
    int topId = TOID(tx, ty, size);
    int bottomId = TOID(tx, ty + bandSize - 1, size);

    //if (topId >= imageSize)
    //	pixel = make_short2(MARKER, MARKER);
    //else
    pixel = tex1Dfetch(disTexColor, topId);

    int myDist = abs(pixel.y - ty);

    for (nid = bottomId - inc; nid >= 0; nid -= inc) {
        //if (nid >= imageSize) continue;
        pixel = tex1Dfetch(disTexColor, nid);



        if (pixel.x != MARKER) {
            nDist = abs(pixel.y - ty);

            if (nDist < myDist)
                output[topId] = pixel;

            //if (topId == 4095) printf("?? %d %d \n", pixel.x, pixel.y);
            break;
        }
    }

    // Last row, look downward
    ty = ty + bandSize - 1;

    //if (bottomId >= imageSize) return;

    pixel = tex1Dfetch(disTexColor, bottomId);
    myDist = abs(pixel.y - ty);

    for (ny = ty + 1, nid = topId + inc; ny < size; ny += bandSize, nid += inc) {
        //if (nid >= imageSize) break;
        pixel = tex1Dfetch(disTexColor, nid);

        //if (nid == 4351)
        //	printf("bbb %d %d\n", pixel.x, pixel.y);


        if (pixel.x != MARKER) {
            nDist = abs(pixel.y - ty);

            if (nDist < myDist)
                output[bottomId] = pixel;
            //if (bottomId == 4095) printf("!? %d %d %d\n", pixel.x, pixel.y, nid);
            break;
        }
    }
}
__global__ void krFDMContouring_kernelUpdateVertical(short2 *output, int size, int band, int bandSize)
{
    int tx = blockIdx.x * blockDim.x + threadIdx.x;
    int ty = blockIdx.y * bandSize;

    short2 top = tex1Dfetch(disTexLinks, TOID(tx, ty, size));
    short2 bottom = tex1Dfetch(disTexLinks, TOID(tx, ty + bandSize - 1, size));
    short2 pixel;

    int dist, myDist;

    int id = TOID(tx, ty, size);

    for (int i = 0; i < bandSize; i++, id += size) {
        pixel = tex1Dfetch(disTexColor, id);



        myDist = abs(pixel.y - (ty + i));

        dist = abs(top.y - (ty + i));
        if (dist < myDist) { myDist = dist; pixel = top; }


        dist = abs(bottom.y - (ty + i));
        if (dist < myDist) pixel = bottom;

        output[id] = pixel;
    }
}
__global__ void krFDMContouring_kernelTranspose(short2 *data, int size)
{
    __shared__ short2 block1[TILE_DIM][TILE_DIM + 1];
    __shared__ short2 block2[TILE_DIM][TILE_DIM + 1];

    int blockIdx_y = blockIdx.x;
    int blockIdx_x = blockIdx.x+blockIdx.y;

    if (blockIdx_x >= gridDim.x)
        return ;

    int blkX, blkY, x, y, id1, id2;
    short2 pixel;

    blkX = __mul24(blockIdx_x, TILE_DIM);
    blkY = __mul24(blockIdx_y, TILE_DIM);

    x = blkX + threadIdx.x;
    y = blkY + threadIdx.y;
    id1 = __mul24(y, size) + x;

    x = blkY + threadIdx.x;
    y = blkX + threadIdx.y;
    id2 = __mul24(y, size) + x;

    // read the matrix tile into shared memory
    for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
        block1[threadIdx.y + i][threadIdx.x] = tex1Dfetch(disTexColor, id1 + __mul24(i, size));
        block2[threadIdx.y + i][threadIdx.x] = tex1Dfetch(disTexColor, id2 + __mul24(i, size));


    }

    __syncthreads();

    // write the transposed matrix tile to global memory
    for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
        pixel = block1[threadIdx.x][threadIdx.y + i];
        data[id2 + __mul24(i, size)] = make_short2(pixel.y, pixel.x);


        pixel = block2[threadIdx.x][threadIdx.y + i];
        data[id1 + __mul24(i, size)] = make_short2(pixel.y, pixel.x);

    }
}
__global__ void krFDMContouring_kernelProximatePoints(short2 *stack, int size, int bandSize)
{
    int tx = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;
    int ty = __mul24(blockIdx.y, bandSize);
    int id = TOID(tx, ty, size);
    int lasty = -1;
    short2 last1, last2, current;
    float i1, i2;


    last1.y = -1; last2.y = -1;

    for (int i = 0; i < bandSize; i++, id += size) {
        current = tex1Dfetch(disTexColor, id);


        if (current.x != MARKER) {

            while (last2.y >= 0) {

                i1 = interpoint(last1.x, last2.y, last2.x, lasty, tx);
                i2 = interpoint(last2.x, lasty, current.x, current.y, tx);



                if (i1 < i2)
                    break;

                lasty = last2.y; last2 = last1;

                if (last1.y >= 0)
                    last1 = stack[TOID(tx, last1.y, size)];
            }

            last1 = last2; last2 = make_short2(current.x, lasty); lasty = current.y;

            stack[id] = last2;
        }
    }

    // Store the pointer to the tail at the last pixel of this band
    if (lasty != ty + bandSize - 1)
        stack[TOID(tx, ty + bandSize - 1, size)] = make_short2(MARKER, lasty);
}
__global__ void krFDMContouring_kernelCreateForwardPointers(short2 *output, int size, int bandSize)
{
    int tx = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;
    int ty = __mul24(blockIdx.y+1, bandSize) - 1;
    int id = TOID(tx, ty, size);
    int lasty = -1, nexty;
    short2 current;

    // Get the tail pointer

    current = tex1Dfetch(disTexLinks, id);

    if (current.x == MARKER)
        nexty = current.y;
    else
        nexty = ty;

    for (int i = 0; i < bandSize; i++, id -= size)
    {

        if (ty - i == nexty) {
            current = make_short2(lasty, tex1Dfetch(disTexLinks, id).y);


            output[id] = current;

            lasty = nexty;
            nexty = current.y;
        }
    }

        // Store the pointer to the head at the first pixel of this band
    if (lasty != ty - bandSize + 1)
            output[id + size] = make_short2(lasty, MARKER);
}
__global__ void krFDMContouring_kernelMergeBands(short2 *output, int size, int bandSize)
{
    int tx = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;
    int band1 = blockIdx.y * 2;
    int band2 = band1 + 1;
    int firsty, lasty;
    short2 last1, last2, current;
    // last1 and last2: x component store the x coordinate of the site,
    // y component store the backward pointer
    // current: y component store the x coordinate of the site,
    // x component store the forward pointer

    // Get the two last items of the first list
    lasty = __mul24(band2, bandSize) - 1;


    last2 = make_short2(tex1Dfetch(disTexColor, TOID(tx, lasty, size)).x,
        tex1Dfetch(disTexLinks, TOID(tx, lasty, size)).y);



    if (last2.x == MARKER) {
        lasty = last2.y;

        if (lasty >= 0)
            last2 = make_short2(tex1Dfetch(disTexColor, TOID(tx, lasty, size)).x,
            tex1Dfetch(disTexLinks, TOID(tx, lasty, size)).y);
        else
            last2 = make_short2(MARKER, MARKER);
    }

    if (last2.y >= 0) {
        // Second item at the top of the stack
        last1 = make_short2(tex1Dfetch(disTexColor, TOID(tx, last2.y, size)).x,
            tex1Dfetch(disTexLinks, TOID(tx, last2.y, size)).y);
    }

    // Get the first item of the second band
    firsty = __mul24(band2, bandSize);
    current = make_short2(tex1Dfetch(disTexLinks, TOID(tx, firsty, size)).x,
        tex1Dfetch(disTexColor, TOID(tx, firsty, size)).x);

    if (current.y == MARKER) {
        firsty = current.x;

        if (firsty >= 0)
            current = make_short2(tex1Dfetch(disTexLinks, TOID(tx, firsty, size)).x,
            tex1Dfetch(disTexColor, TOID(tx, firsty, size)).x);
        else
            current = make_short2(MARKER, MARKER);
    }

    float i1, i2;

    // Count the number of item in the second band that survive so far.
    // Once it reaches 2, we can stop.
    int top = 0;

    while (top < 2 && current.y >= 0) {
        // While there's still something on the left
        while (last2.y >= 0) {
            i1 = interpoint(last1.x, last2.y, last2.x, lasty, tx);
            i2 = interpoint(last2.x, lasty, current.y, firsty, tx);

            if (i1 < i2)
                break;

            lasty = last2.y; last2 = last1;
            top--;

            if (last1.y >= 0)
                last1 = make_short2(tex1Dfetch(disTexColor, TOID(tx, last1.y, size)).x,
                output[TOID(tx, last1.y, size)].y);
        }

        // Update the current pointer

        output[TOID(tx, firsty, size)] = make_short2(current.x, lasty);

        if (lasty >= 0)
        {

            output[TOID(tx, lasty, size)] = make_short2(firsty, last2.y);
        }

        last1 = last2; last2 = make_short2(current.y, lasty); lasty = firsty;
        firsty = current.x;

        top = max(1, top + 1);

        // Advance the current pointer to the next one
        if (firsty >= 0)
            current = make_short2(tex1Dfetch(disTexLinks, TOID(tx, firsty, size)).x,
            tex1Dfetch(disTexColor, TOID(tx, firsty, size)).x);
        else
            current = make_short2(MARKER, MARKER);
    }

    // Update the head and tail pointer.
    firsty = __mul24(band1, bandSize);
    lasty = __mul24(band2, bandSize);
    current = tex1Dfetch(disTexLinks, TOID(tx, firsty, size));

    if (current.y == MARKER && current.x < 0) {	// No head?
        last1 = tex1Dfetch(disTexLinks, TOID(tx, lasty, size));

        if (last1.y == MARKER)
            current.x = last1.x;
        else
            current.x = lasty;


        output[TOID(tx, firsty, size)] = current;
    }

    firsty = __mul24(band1, bandSize) + bandSize - 1;
    lasty = __mul24(band2, bandSize) + bandSize - 1;
    current = tex1Dfetch(disTexLinks, TOID(tx, lasty, size));

    if (current.x == MARKER && current.y < 0) {	// No tail?
        last1 = tex1Dfetch(disTexLinks, TOID(tx, firsty, size));

        if (last1.x == MARKER)
            current.y = last1.y;
        else
            current.y = firsty;


        output[TOID(tx, lasty, size)] = current;
    }
}
__global__ void krFDMContouring_kernelDoubleToSingleList(short2 *output, int size)
{
    int tx = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;
    int ty = blockIdx.y;
    int id = TOID(tx, ty, size);

    output[id] = make_short2(tex1Dfetch(disTexColor, id).x, tex1Dfetch(disTexLinks, id).y);

}
__global__ void krFDMContouring_kernelColor(short2 *output, int size)
{
    __shared__ short2 s_last1[BLOCKSIZE], s_last2[BLOCKSIZE];
    __shared__ int s_lasty[BLOCKSIZE];

    int col = threadIdx.x;
    int tid = threadIdx.y;
    int tx = __mul24(blockIdx.x, blockDim.x) + col;
    int dx, dy, lasty;
    unsigned int best, dist;
    short2 last1, last2;
    int count = 0;
    if (tid == blockDim.y - 1) {
        lasty = size - 1;


        last2 = tex1Dfetch(disTexColor, __mul24(lasty, size) + tx);

        if (last2.x == MARKER) {
            lasty = last2.y;
            last2 = tex1Dfetch(disTexColor, __mul24(lasty, size) + tx);
        }

        if (last2.y >= 0)
            last1 = tex1Dfetch(disTexColor, __mul24(last2.y, size) + tx);

        s_last1[col] = last1; s_last2[col] = last2; s_lasty[col] = lasty;
    }

    __syncthreads();

    count = 0;
    for (int ty = size - 1 - tid; ty >= 0; ty -= blockDim.y) {
        last1 = s_last1[col]; last2 = s_last2[col]; lasty = s_lasty[col];


        dx = last2.x - tx; dy = lasty - ty;
        best = dist = __mul24(dx, dx) + __mul24(dy, dy);

        while (last2.y >= 0) {
            dx = last1.x - tx; dy = last2.y - ty;
            dist = __mul24(dx, dx) + __mul24(dy, dy);

            if (dist > best)
                break;
            count++;
            best = dist; lasty = last2.y; last2 = last1;

            if (last2.y >= 0)
                last1 = tex1Dfetch(disTexColor, __mul24(last2.y, size) + tx);
        }

        __syncthreads();


        output[TOID(tx, ty, size)] = make_short2(last2.x, lasty);

        if (tid == blockDim.y - 1) {
            s_last1[col] = last1; s_last2[col] = last2; s_lasty[col] = lasty;
        }

        __syncthreads();
    }
}


void LDNIcudaOperation::LDNIFDMContouring_BuildDistanceMap(bool *gridNodes, int layerID, int imageSize[], short2 *&disTexturesA, short2 *&disTexturesB,  int disTexSize)
{


    int phase1Band  = 16;
    int phase2Band  = 16;
    int phase3Band  = 16;
    long time = clock();

    krFDMContouring_InitializedValue<<<BLOCKS_PER_GRID,THREADS_PER_BLOCK>>>( disTexturesA, disTexSize*disTexSize, MARKER);

    //Initialize the Distance Map
    krFDMContouring_InitializedDistanceMap<<<BLOCKS_PER_GRID,THREADS_PER_BLOCK>>>(gridNodes, disTexturesA, imageSize[0]*imageSize[2], imageSize[0]*imageSize[2],
                                                                    make_int3(imageSize[0],imageSize[1],imageSize[2]), layerID, disTexSize);

    //printf("Start computing distance field of layer (%d) ......\n",layerID);
    //Phase 1:  Flood vertically in their own bands
    dim3 block = dim3(BLOCKSIZE);
    dim3 grid = dim3(disTexSize / block.x, phase1Band);


    cudaBindTexture(0, disTexColor, disTexturesA);
    krFDMContouring_kernelFloodDown<<< grid, block >>>(disTexturesB, disTexSize, disTexSize / phase1Band);
    cudaUnbindTexture(disTexColor);


    cudaBindTexture(0, disTexColor, disTexturesB);
    krFDMContouring_kernelFloodUp<<< grid, block >>>(disTexturesB, disTexSize, disTexSize / phase1Band);


    //Phase 1:  Passing information between bands
    grid = dim3(disTexSize / block.x, phase1Band);
    krFDMContouring_kernelPropagateInterband<<< grid, block >>>(disTexturesA, disTexSize, disTexSize / phase1Band);


    cudaBindTexture(0, disTexLinks, disTexturesA);
    krFDMContouring_kernelUpdateVertical<<< grid, block >>>(disTexturesB, disTexSize, phase1Band, disTexSize / phase1Band);
    cudaUnbindTexture(disTexLinks);
    cudaUnbindTexture(disTexColor);


    //Phase 1: Transpose
    block = dim3(TILE_DIM, BLOCK_ROWS);
    grid = dim3(disTexSize / TILE_DIM, disTexSize / TILE_DIM);

    cudaBindTexture(0, disTexColor, disTexturesB);
    krFDMContouring_kernelTranspose<<< grid, block >>>(disTexturesB, disTexSize);
    cudaUnbindTexture(disTexColor);



    //Phase 2: Compute proximate points locally in each band
    block = dim3(BLOCKSIZE);
    grid = dim3(disTexSize / block.x, phase2Band);
    cudaBindTexture(0, disTexColor, disTexturesB);
    krFDMContouring_kernelProximatePoints<<< grid, block >>>(disTexturesA, disTexSize, disTexSize / phase2Band);
    cudaBindTexture(0, disTexLinks, disTexturesA);
    krFDMContouring_kernelCreateForwardPointers<<< grid, block >>>(disTexturesA, disTexSize, disTexSize / phase2Band);

    //Phase 2:  Repeatly merging two bands into one
    for (int noBand = phase2Band; noBand > 1; noBand /= 2) {
        grid = dim3(disTexSize / block.x, noBand / 2);
        krFDMContouring_kernelMergeBands<<< grid, block >>>(disTexturesA, disTexSize, disTexSize / noBand);
    }

    //Phase 2:  Replace the forward link with the X coordinate of the seed to remove the need of looking at the other texture. We need it for coloring.
    grid = dim3(disTexSize / block.x, disTexSize);
    krFDMContouring_kernelDoubleToSingleList<<< grid, block >>>(disTexturesA, disTexSize);
    cudaUnbindTexture(disTexLinks);
    cudaUnbindTexture(disTexColor);


    //Phase 3:
    block = dim3(BLOCKSIZE / phase2Band, phase2Band);
    grid = dim3(disTexSize / block.x);
    cudaBindTexture(0, disTexColor, disTexturesA);
    krFDMContouring_kernelColor<<< grid, block >>>(disTexturesB, disTexSize);
    cudaUnbindTexture(disTexColor);


    //Phase 3: Transpose
    block = dim3(TILE_DIM, BLOCK_ROWS);
    grid = dim3(disTexSize / TILE_DIM, disTexSize / TILE_DIM);

    cudaBindTexture(0, disTexColor, disTexturesB);
    krFDMContouring_kernelTranspose<<< grid, block >>>(disTexturesB, disTexSize);
    cudaUnbindTexture(disTexColor);

    //--------------------------------------------------------------------------------------------------------------------------------------------------
    //printf("Finish computing distance field of layer. %ld(ms)\n", (clock()-time)/(CLOCKS_PER_SEC/1000));



}

extern "C" void call_krFDMContouring_Dilation(bool *seedNodes, bool* output, int nodeNum, int3 imageRes, double realThreshold, int gridRadius, int iy)
{
	krFDMContouring_Dilation << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (seedNodes, output, nodeNum, imageRes, realThreshold, gridRadius, iy);
}
__global__ void krFDMContouring_Dilation(bool *seedNodes, bool* output, int nodeNum, int3 imageRes, double realThreshold, int gridRadius, int iy)
{
    int index=threadIdx.x+blockIdx.x*blockDim.x;
    unsigned int ix, iz;
    int m, n;

    while (index<nodeNum) {
        ix = index%imageRes.x;		iz = index/imageRes.x;

        if (seedNodes[index])
        {
            //if (iy == 87) printf("--%d %d %d %d\n", iy, index, ix, iz);
            for(m=-gridRadius; m<=gridRadius; m++)
            {
                for(n=-gridRadius; n<=gridRadius; n++)
                {
                    if(ix+m<0 || ix+m>=imageRes.x || iz+n<0 || iz+n>=imageRes.z || (m==0 && n==0))
                        continue;
                    else if(m*m+n*n > realThreshold*realThreshold)
                        continue;

                    output[(iz+n)*imageRes.x+(ix+m)] = true;
                    //if (iy == 137 && ix == 74 && iz == 36) printf("%d %d %d\n", index, ix, iz);
                }
            }
        }

        index += blockDim.x * gridDim.x;

    }
}

__global__ void krFDMContouring_Filter1(bool *gridNodes, bool *outNodes, int nodeNum, int3 imageRes, int iy)
{
    int index=threadIdx.x+blockIdx.x*blockDim.x;
    unsigned int ix, iz;
    bool node;

    while (index<nodeNum) {
        ix = index%imageRes.x;	iz = index/imageRes.x;

        node = outNodes[index];
        //outNodes[index] = node && (!gridNodes[iz*imageRes.x*imageRes.y+iy*imageRes.x+ix]);
        outNodes[index] = node && (!gridNodes[index]);
        //Exclude self-support region
        //if ((node && (!gridNodes[index])) && iy == 125)
        //	printf("%d %d %d \n", ix, iy, iz);

        index += blockDim.x * gridDim.x;

    }

}

__global__ void krFDMContouring_Filter2(short2 *disMap, bool *outNodes, bool *suptNodes, int nodeNum, int3 imageRes, bool *RegionFlag, double t, int disTexSize, int iy)
{
    int index=threadIdx.x+blockIdx.x*blockDim.x;
    unsigned int ix, iz;
    bool node;
    short2 d;
    double dist;

    while (index<nodeNum) {
        ix = index%imageRes.x;	iz = index/imageRes.x;

        node = outNodes[index];
        if (node)
        {
            //d = disMap[index];
            d = disMap[iz*disTexSize+ix];
            d.x = ix - d.x;
            d.y = iz - d.y;
            dist = d.x*d.x + d.y*d.y;

            //if (iy == 137 && ix == 73) printf("%d %d %d %d %f\n", ix, iz, d.x, d.y, dist);

            if(suptNodes[index] && dist <= t)
            {
                RegionFlag[0] = false;
            }
            else
                outNodes[index] = false;

        }

        index += blockDim.x * gridDim.x;

    }

}

__global__ void krFDMContouring_Filter3(bool* outNodes, bool *seedNodes, int nodeNum)
{
    int index=threadIdx.x+blockIdx.x*blockDim.x;
    unsigned int ix, iz;
    bool node;

    while (index<nodeNum) {


        node = outNodes[index];
        if (node)
        {
            seedNodes[index] = true;
            node = true;
        }
        else
        {
            node = seedNodes[index];
        }

        outNodes[index] = node;

        index += blockDim.x * gridDim.x;

    }

}

extern "C" void call_krFDMContouring_SubtractSolidRegion(bool *gridNodes, bool *suptNodes, bool *seedNodes, int nodeNum, int3 imageRes, int iy)
{
    krFDMContouring_SubtractSolidRegion<<<BLOCKS_PER_GRID,THREADS_PER_BLOCK>>>(gridNodes, suptNodes, seedNodes, nodeNum, imageRes, iy);
}

__global__ void krFDMContouring_SubtractSolidRegion(bool *gridNodes, bool *suptNodes, bool *seedNodes, int nodeNum, int3 imageRes, int iy)
{
    int index=threadIdx.x+blockIdx.x*blockDim.x;
    unsigned int ix, iz;

    while (index<nodeNum) {
        ix = index%imageRes.x;	//iy = (index/imageRes.x)%imageRes.y+1;
        iz = index/imageRes.x;
        //iz = index/(imageRes.x*imageRes.y);

		if (gridNodes[iz*imageRes.x*imageRes.y + iy * imageRes.x + ix])
		{
			suptNodes[iz*imageRes.x*(imageRes.y - 1) + iy * imageRes.x + ix] = false;
			seedNodes[index] = true;
		}
		else if (gridNodes[iz*imageRes.x*imageRes.y + (iy + 1)*imageRes.x + ix])
		{
			suptNodes[iz*imageRes.x*(imageRes.y - 1) + iy * imageRes.x + ix] = true;
		}
		else
		{
			suptNodes[iz*imageRes.x*(imageRes.y - 1) + iy * imageRes.x + ix] = false;
		}
        /*if (gridNodes[iy*imageRes.x*imageRes.z+iz*imageRes.x+ix])
        {
            suptNodes[iy*imageRes.x*imageRes.z+iz*imageRes.x+ix] = false;
            seedNodes[index] = true;
        }
        else if (gridNodes[(iy+1)*imageRes.x*imageRes.z+iz*imageRes.x+ix])
        {
            suptNodes[iy*imageRes.x*imageRes.z+iz*imageRes.x+ix] = true;
        }
        else
        {
            suptNodes[iy*imageRes.x*imageRes.z+iz*imageRes.x+ix] = false;
        }*/

        index += blockDim.x * gridDim.x;
    }
}

extern "C" void call_krFDMContouring_CopyNodesrom3Dto2D(bool *m_2DNodes, bool *m_3DNodes, int nodeNum, int3 imageRes, int iy)
{
    krFDMContouring_CopyNodesrom3Dto2D<<<BLOCKS_PER_GRID,THREADS_PER_BLOCK>>>(m_2DNodes, m_3DNodes, nodeNum, imageRes, iy);
}
__global__ void krFDMContouring_CopyNodesrom3Dto2D(bool *m_2DNodes, bool *m_3DNodes, int nodeNum, int3 imageRes, int iy)
{
    int index=threadIdx.x+blockIdx.x*blockDim.x;
    unsigned int ix, iz;

    while (index<nodeNum) {
        ix = index%imageRes.x;
        iz = index/imageRes.x;

        m_2DNodes[index] = m_3DNodes[iz*imageRes.x*imageRes.y+iy*imageRes.x+ix];

        index += blockDim.x * gridDim.x;
    }
}

void LDNIcudaOperation::LDNIFDMContouring_GrowthAndSwallow(bool *gridNodes, int i, int imageSize[], double t, double nSampleWidth, bool *&suptRegion, bool *&solidRegion, short2 *distanceMapA, short2 *distanceMapB , int disTexSize, bool *tempNodes, bool *tempSeeds)
{
    bool* m_RegionFinish;

//	bool *tempNodes;
//	bool *tempSeeds;

    double realThreshold = (2.5*nSampleWidth-nSampleWidth)/nSampleWidth;
    int gridRadius = (int)floor(realThreshold);
    int3 imgRes = make_int3(imageSize[0], imageSize[1], imageSize[2]);
    double realT = (t-nSampleWidth)/nSampleWidth;
    bool* bRegionFinish = (bool*)malloc(sizeof(bool));
    //long time = clock();

    LDNIFDMContouring_BuildDistanceMap(gridNodes, i, imageSize, distanceMapA, distanceMapB, disTexSize);
    //cudaThreadSynchronize();
//	printf("1 %ld(ms)\n", clock()-time);

    //CUDA_SAFE_CALL( cudaMalloc( (void**)&(tempNodes), imageSize[0]*imageSize[2]*sizeof(bool) ) );
    CUDA_SAFE_CALL( cudaMemset( (void*)tempNodes, false, imageSize[0]*imageSize[2]*sizeof(bool) ) );
    //CUDA_SAFE_CALL( cudaMalloc( (void**)&(tempSeeds), imageSize[0]*imageSize[2]*sizeof(bool) ) );
    CUDA_SAFE_CALL( cudaMemset( (void*)tempSeeds, false, imageSize[0]*imageSize[2]*sizeof(bool) ) );
    CUDA_SAFE_CALL( cudaMalloc( (void**)&(m_RegionFinish), sizeof(bool) ) );
    CUDA_SAFE_CALL( cudaMemset( (void*)m_RegionFinish, false, sizeof(bool) ) );


    CUDA_SAFE_CALL( cudaMemcpy( tempNodes, solidRegion, imageSize[0]*imageSize[2]*sizeof(bool), cudaMemcpyDeviceToDevice ) );

    bRegionFinish[0] = false;

    //time = clock();

    while(!bRegionFinish[0])
    {
        bRegionFinish[0] = true;
        CUDA_SAFE_CALL( cudaMemset( (void*)m_RegionFinish, true, sizeof(bool) ) );

        //dilation
        CUDA_SAFE_CALL( cudaMemcpy( tempSeeds, tempNodes, imageSize[0]*imageSize[2]*sizeof(bool), cudaMemcpyDeviceToDevice ) );
        krFDMContouring_Dilation<<<BLOCKS_PER_GRID,THREADS_PER_BLOCK>>>(tempSeeds, tempNodes, imageSize[0]*imageSize[2], imgRes, realThreshold, gridRadius, i);
        krFDMContouring_Filter1<<<BLOCKS_PER_GRID,THREADS_PER_BLOCK>>>(solidRegion, tempNodes, imageSize[0]*imageSize[2], imgRes, i);


        //intersection with original input image
        krFDMContouring_Filter2<<<BLOCKS_PER_GRID,THREADS_PER_BLOCK>>>(distanceMapB, tempNodes, suptRegion, imageSize[0]*imageSize[2], imgRes, m_RegionFinish, realT*realT,  disTexSize, i);
        CUDA_SAFE_CALL( cudaMemcpy(bRegionFinish, m_RegionFinish, sizeof(bool), cudaMemcpyDeviceToHost) );
        krFDMContouring_Filter3<<<BLOCKS_PER_GRID,THREADS_PER_BLOCK>>>(tempNodes, solidRegion, imageSize[0]*imageSize[2]);


    }

    krFDMContouring_Filter1<<<BLOCKS_PER_GRID,THREADS_PER_BLOCK>>>(solidRegion, suptRegion, imageSize[0]*imageSize[2], imgRes, i);
    //cudaThreadSynchronize();
    //printf("2 %ld(ms)\n", (clock()-time)/(CLOCKS_PER_SEC/1000));



    //cudaFree(tempNodes);
    //cudaFree(tempSeeds);
    cudaFree(m_RegionFinish);


}

extern "C" void call_krFDMContouring_Filter4(bool* outNodes, bool *inNodesA, bool *inNodesB, int nodeNum, int3 imageRes, int iy)
{
    krFDMContouring_Filter4<<<BLOCKS_PER_GRID,THREADS_PER_BLOCK>>>(outNodes, inNodesA, inNodesB, nodeNum, imageRes, iy);
}
__global__ void krFDMContouring_Filter4(bool* outNodes, bool *inNodesA, bool *inNodesB, int nodeNum, int3 imageRes, int iy)
{
    int index=threadIdx.x+blockIdx.x*blockDim.x;
    unsigned int ix, iz;
    bool node;

    while (index<nodeNum) {
        ix = index%imageRes.x;	iz = index/imageRes.x;

        node = inNodesA[iz*imageRes.x*imageRes.y+iy*imageRes.x+ix] && (!inNodesB[index]);
        outNodes[index] = node;
        inNodesA[iz*imageRes.x*imageRes.y+iy*imageRes.x+ix] = inNodesB[index];

        index += blockDim.x * gridDim.x;
    }

}
extern "C" void call_krFDMContouring_integrateImageintoGrossImage(bool* outNodes, bool *inNodes, int nodeNum, int3 imageRes, int iy)
{
    krFDMContouring_integrateImageintoGrossImage<<<BLOCKS_PER_GRID,THREADS_PER_BLOCK>>>(outNodes, inNodes, nodeNum, imageRes, iy);
}
__global__ void krFDMContouring_integrateImageintoGrossImage(bool* outNodes, bool *inNodes, int nodeNum, int3 imageRes, int iy)
{
    int index=threadIdx.x+blockIdx.x*blockDim.x;
    unsigned int ix, iz;
    bool node;

    while (index<nodeNum) {

        ix = index%imageRes.x;	iz = index/imageRes.x;
        node = outNodes[index];

        if (inNodes[iz*imageRes.x*imageRes.y+iy*imageRes.x+ix] && !node)
            outNodes[index] = true;

        index += blockDim.x * gridDim.x;
    }


}
extern "C" void call_krFDMContouring_CopyNodesrom2Dto3D(bool *m_2DNodes, bool *m_3DNodes, int nodeNum, int3 imageRes, int iy)
{
    krFDMContouring_CopyNodesrom2Dto3D<<<BLOCKS_PER_GRID,THREADS_PER_BLOCK>>>(m_2DNodes, m_3DNodes, nodeNum, imageRes, iy);
}
__global__ void krFDMContouring_CopyNodesrom2Dto3D(bool *m_2DNodes, bool *m_3DNodes, int nodeNum, int3 imageRes, int iy)
{
    int index=threadIdx.x+blockIdx.x*blockDim.x;
    unsigned int ix, iz;

    while (index<nodeNum) {
        ix = index%imageRes.x;
        iz = index/imageRes.x;

        m_3DNodes[iz*imageRes.x*imageRes.y+iy*imageRes.x+ix] = m_2DNodes[index];

        //if (iy == 0 && m_2DNodes[index])
        //	printf("copy %d %d %d %d %d %d %d \n", ix, iy, iz, iz*imageRes.x*imageRes.y+iy*imageRes.x+ix, index, imageRes.x, imageRes.y);
        /*if (m_3DNodes[iz*imageRes.x*imageRes.y+iy*imageRes.x+ix])
                printf("warning!!! %d %d %d \n", ix, iy, iz);
        if (m_2DNodes[index])
            printf("2warning!!! %d %d %d \n", ix, iy, iz);*/
        /*if (m_3DNodes[iz*imageRes.x*imageRes.y+27*imageRes.x+ix]) printf("warning!!! %d %d %d \n", ix, iy, iz);
        if (m_3DNodes[iz*imageRes.x*imageRes.y+26*imageRes.x+ix]) printf("warning!!! %d %d %d \n", ix, iy, iz);
        if (m_3DNodes[iz*imageRes.x*imageRes.y+29*imageRes.x+ix]) printf("warning!!! %d %d %d \n", ix, iy, iz);
        if (m_3DNodes[iz*imageRes.x*imageRes.y+30*imageRes.x+ix]) printf("warning!!! %d %d %d \n", ix, iy, iz);*/

        index += blockDim.x * gridDim.x;
    }
}

void LDNIcudaOperation::LDNIFDMContouring_Closing(int imageSize[], double t, double nSampleWidth, bool *&inputNodes, int i, bool *tempNodes, bool *tempSeeds)
{
	double realThreshold = (t - nSampleWidth) / nSampleWidth;
	int gridRadius = (int)floor(realThreshold);
	int3 imgRes = make_int3(imageSize[0], imageSize[1], imageSize[2]);
	int3 suptimgRes = make_int3(imageSize[0], imageSize[1] - 1, imageSize[2]);

	//long time = clock();
	//bool *tempNodes;
	//bool *tempSeeds;
	//CUDA_SAFE_CALL( cudaMalloc( (void**)&(tempNodes), imageSize[0]*imageSize[2]*sizeof(bool) ) );
	CUDA_SAFE_CALL(cudaMemset((void*)tempNodes, false, imageSize[0] * imageSize[2] * sizeof(bool)));
	//CUDA_SAFE_CALL( cudaMalloc( (void**)&(tempSeeds), imageSize[0]*imageSize[2]*sizeof(bool) ) );
	CUDA_SAFE_CALL(cudaMemset((void*)tempSeeds, false, imageSize[0] * imageSize[2] * sizeof(bool)));
	cudaThreadSynchronize();
	//printf("1 %ld(ms) \n", (clock()-time)/(CLOCKS_PER_SEC/1000));	time = clock();
	//CUDA_SAFE_CALL( cudaMemcpy( tempNodes, inputNodes, imageSize[0]*imageSize[2]*sizeof(bool), cudaMemcpyDeviceToDevice ) );
	krFDMContouring_CopyNodesrom3Dto2D << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (tempNodes, inputNodes, imageSize[0] * imageSize[2], suptimgRes, i);
	CUDA_SAFE_CALL(cudaMemcpy(tempSeeds, tempNodes, imageSize[0] * imageSize[2] * sizeof(bool), cudaMemcpyDeviceToDevice));


	krFDMContouring_Dilation << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (tempNodes, tempSeeds, imageSize[0] * imageSize[2], imgRes, realThreshold, gridRadius, i);
	cudaThreadSynchronize();
	//printf("2 %ld(ms) \n", (clock()-time)/(CLOCKS_PER_SEC/1000));	time = clock();

	CUDA_SAFE_CALL(cudaMemcpy(tempNodes, tempSeeds, imageSize[0] * imageSize[2] * sizeof(bool), cudaMemcpyDeviceToDevice));
	krFDMContouring_Erosion << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (tempNodes, tempSeeds, imageSize[0] * imageSize[2], imgRes, realThreshold, gridRadius);

	krFDMContouring_CopyNodesrom2Dto3D << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (tempSeeds, inputNodes, imageSize[0] * imageSize[2], suptimgRes, i);
	cudaThreadSynchronize();
	//printf("3 %ld(ms) \n", (clock()-time)/(CLOCKS_PER_SEC/1000));	time = clock();


	//cudaFree(tempNodes);
}

extern "C" void call_krFDMContouring_Filter5(bool* outNodes, bool *inNodesA, bool *gridNodes, int nodeNum, int iy, int3 imageRes)
{
	krFDMContouring_Filter5<<<BLOCKS_PER_GRID,THREADS_PER_BLOCK>>>(outNodes, inNodesA, gridNodes, nodeNum, iy, imageRes);
}
__global__ void krFDMContouring_Filter5(bool* outNodes, bool *inNodesA, bool *gridNodes, int nodeNum, int iy, int3 imageRes)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int ix, iz;
	bool node;

	while (index < nodeNum) {

		ix = index % imageRes.x;	iz = index / imageRes.x;

		node = outNodes[iz*imageRes.x*(imageRes.y - 1) + iy * imageRes.x + ix];
		outNodes[iz*imageRes.x*(imageRes.y - 1) + iy * imageRes.x + ix] = node && (!inNodesA[index]) && (!gridNodes[iz*imageRes.x*imageRes.y + iy * imageRes.x + ix]);

		index += blockDim.x * gridDim.x;
	}

}


__device__ float interpoint(int x1, int y1, int x2, int y2, int x0)
{
    float xM = float(x1 + x2) / 2.0f;
    float yM = float(y1 + y2) / 2.0f;
    float nx = x2 - x1;
    float ny = y2 - y1;

    return yM + nx * (xM - x0) / ny;
}

__global__ void krFDMContouring_Erosion(bool *gridNodes, bool* output, int nodeNum, int3 imageRes, double realThreshold, int gridRadius)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int ix, iz;
	int m, n;
	bool bflag;

	while (index < nodeNum) {
		ix = index % imageRes.x;		iz = index / imageRes.x;

		if (gridNodes[iz*imageRes.x + ix])
		{
			bflag = true;
			for (m = -gridRadius; m <= gridRadius; m++)
			{
				for (n = -gridRadius; n <= gridRadius; n++)
				{
					if (ix + m < 0 || ix + m >= imageRes.x || iz + n < 0 || iz + n >= imageRes.z || (m == 0 && n == 0))
						continue;
					else if (m*m + n * n > realThreshold*realThreshold)
						continue;

					if (!gridNodes[(iz + n)*imageRes.x + (ix + m)])
					{
						output[index] = false;
						bflag = false;
						break;
					}
				}
				if (!bflag)
					break;
			}
		}
		index += blockDim.x * gridDim.x;

	}
}

extern "C" void call_krSLAContouring_Initialization(bool *tempImg, bool *targetImg, bool *gridNodes, int nodeNum, int3 imageRes, int iy)
{
	krSLAContouring_Initialization << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (tempImg, targetImg, gridNodes, nodeNum, imageRes, iy);
}
__global__ void krSLAContouring_Initialization(bool *tempImg, bool *targetImg, bool *gridNodes, int nodeNum, int3 imageRes, int iy)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int ix, iz;
	bool node;

	while (index < nodeNum) {
		ix = index % imageRes.x;	//iy = (index/imageRes.x)%imageRes.y+1;
		iz = index / imageRes.x;
		//iz = index/(imageRes.x*imageRes.y);

		node = gridNodes[iz*imageRes.x*imageRes.y + iy * imageRes.x + ix];

		tempImg[index] = node;
		targetImg[index] = gridNodes[iz*imageRes.x*imageRes.y + (iy + 1)*imageRes.x + ix] && (!node);


		//if ((gridNodes[iz*imageRes.x*imageRes.y+(iy+1)*imageRes.x+ix] && (!node)) && iy == 134)
		//if (iy == 0 && targetImg[index])
		//	printf(" %d %d %d \n", ix, iy, iz);

		index += blockDim.x * gridDim.x;
	}
}

void LDNIcudaOperation::LDNISLAContouring_GrowthAndSwallow(double t, bool *&suptNodes, bool *&seedNodes, int i, int imageSize[], double nSampleWidth, short2 *disTextureA, short2 *disTextureB, int disTexSize)
{
	bool* m_RegionFinish;
	bool *tempNodes;
	bool *tempSeeds;

	double realThreshold = (2.5*nSampleWidth - nSampleWidth) / nSampleWidth;
	int gridRadius = (int)floor(realThreshold);
	int3 imgRes = make_int3(imageSize[0], imageSize[1], imageSize[2]);
	double realT = (t - nSampleWidth) / nSampleWidth;
	bool* bRegionFinish = (bool*)malloc(sizeof(bool));

	LDNISLAContouring_BuildDistanceMap(seedNodes, i, imageSize, disTextureA, disTextureB, disTexSize);

	CUDA_SAFE_CALL(cudaMalloc((void**)&(tempNodes), imageSize[0] * imageSize[2] * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)tempNodes, false, imageSize[0] * imageSize[2] * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(tempSeeds), imageSize[0] * imageSize[2] * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)tempSeeds, false, imageSize[0] * imageSize[2] * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(m_RegionFinish), sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)m_RegionFinish, false, sizeof(bool)));


	CUDA_SAFE_CALL(cudaMemcpy(tempNodes, seedNodes, imageSize[0] * imageSize[2] * sizeof(bool), cudaMemcpyDeviceToDevice));

	bRegionFinish[0] = false;



	while (!bRegionFinish[0])
	{
		bRegionFinish[0] = true;
		CUDA_SAFE_CALL(cudaMemset((void*)m_RegionFinish, true, sizeof(bool)));

		//dilation
		CUDA_SAFE_CALL(cudaMemcpy(tempSeeds, tempNodes, imageSize[0] * imageSize[2] * sizeof(bool), cudaMemcpyDeviceToDevice));
		krFDMContouring_Dilation << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (tempSeeds, tempNodes, imageSize[0] * imageSize[2], imgRes, realThreshold, gridRadius, i);
		krFDMContouring_Filter1 << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (seedNodes, tempNodes, imageSize[0] * imageSize[2], imgRes, i);


		//intersection with original input image
		krFDMContouring_Filter2 << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (disTextureB, tempNodes, suptNodes, imageSize[0] * imageSize[2], imgRes, m_RegionFinish, realT*realT, disTexSize, i);
		CUDA_SAFE_CALL(cudaMemcpy(bRegionFinish, m_RegionFinish, sizeof(bool), cudaMemcpyDeviceToHost));
		krFDMContouring_Filter3 << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (tempNodes, seedNodes, imageSize[0] * imageSize[2]);


	}

	//printf("test \n");
	krFDMContouring_Filter1 << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (seedNodes, suptNodes, imageSize[0] * imageSize[2], imgRes, i);
	//printf("test 2\n");



	cudaFree(tempNodes);
	cudaFree(tempSeeds);
	cudaFree(m_RegionFinish);
}

extern "C" void call_krSLAContouring_Filter1(bool *assistImg, bool *tempImg, bool *targetImg, int nodeNum, int2 suptRes, int suptRadius, int iy)
{
	krSLAContouring_Filter1 << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (assistImg, tempImg, targetImg, nodeNum, suptRes, suptRadius, iy);
}
__global__ void krSLAContouring_Filter1(bool *assistImg, bool *tempImg, bool *targetImg, int nodeNum, int2 suptRes, int suptRadius, int iy)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int ix, iz;
	bool node;
	int idx, idz;

	while (index < nodeNum) {

		ix = index % suptRes.x;	iz = index / suptRes.x;

		idx = ix * suptRadius;	idz = iz * suptRadius;
		if (targetImg[idz*suptRes.x + idx])
		{
			assistImg[idz*suptRes.x + idx] = true;
			tempImg[idz*suptRes.x + idx] = true;
		}


		index += blockDim.x * gridDim.x;

	}

}

extern "C" void call_krSLAContouring_OrthoSearchRemainAnchorZ(bool *assistImg, bool *tempImg, bool *targetImg, int nodeNum, int2 indRes, int2 imgRes, int iy)
{
	krSLAContouring_OrthoSearchRemainAnchorZ << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (assistImg, tempImg, targetImg, nodeNum,
		indRes, imgRes, iy);
}
__global__ void krSLAContouring_OrthoSearchRemainAnchorZ(bool *assistImg, bool *tempImg, bool *targetImg, int nodeNum, int2 indRes, int2 imgRes, int iy)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int ix, i, accumx, accumz, localcount;
	bool bstart;

	while (index < nodeNum) {
		ix = index * indRes.x;
		bstart = false;
		accumx = 0;
		accumz = 0;
		localcount = 0;
		for (i = 0; i < indRes.y; i++)
		{
			if (targetImg[i*imgRes.x + ix])
			{
				bstart = true;
				localcount++;
				accumx += ix;
				accumz += i;
			}
			if (!targetImg[i*imgRes.x + ix] && bstart)
			{
				bstart = false;
				accumx /= localcount;
				accumz /= localcount;
				assistImg[accumz*imgRes.x + accumx] = true;
				tempImg[accumz*imgRes.x + accumx] = true;
				localcount = 0;
				accumx = 0;
				accumz = 0;
			}
		}



		index += blockDim.x * gridDim.x;
	}
}

extern "C" void call_krSLAContouring_OrthoSearchRemainAnchorX(bool *assistImg, bool *tempImg, bool *targetImg, int nodeNum, int2 indRes, int2 imgRes, int iy)
{
	krSLAContouring_OrthoSearchRemainAnchorX << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (assistImg, tempImg, targetImg, nodeNum,
		indRes, imgRes, iy);
}
__global__ void krSLAContouring_OrthoSearchRemainAnchorX(bool *assistImg, bool *tempImg, bool *targetImg, int nodeNum, int2 indRes, int2 imgRes, int iy)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int iz, i, accumx, accumz, localcount;
	bool bstart;

	while (index < nodeNum) {
		iz = index * indRes.x;
		bstart = false;
		accumx = 0;
		accumz = 0;
		localcount = 0;
		for (i = 0; i < indRes.y; i++)
		{
			if (targetImg[iz*imgRes.x + i])
			{
				bstart = true;
				localcount++;
				accumx += i;
				accumz += iz;
			}
			if (!targetImg[iz*imgRes.x + i] && bstart)
			{
				bstart = false;
				accumx /= localcount;
				accumz /= localcount;
				assistImg[accumz*imgRes.x + accumx] = true;
				tempImg[accumz*imgRes.x + accumx] = true;
				localcount = 0;
				accumx = 0;
				accumz = 0;
			}
		}



		index += blockDim.x * gridDim.x;
	}
}

void LDNIcudaOperation::LDNISLAContouring_ThirdClassCylinder(double threshold, bool *&targetImg, bool *&assistImg, bool *&tempImg, int2 imgRes, double nSampleWidth, int i, short2 *disTextureA, short2 *disTextureB, int disTexSize)
{


	int imageSize[3] = { imgRes.x, 0, imgRes.y };
	int nodeNum = imgRes.x*imgRes.y;
	int regionNum;
	unsigned int MAX_MARKER = max(imgRes.x, imgRes.y)*max(imgRes.x, imgRes.y) + 10;//nodeNum+10;//((1 << 31)-1)*2 +1; // 4294967295

	thrust::device_ptr<bool> target_ptr(targetImg);
	int c_count = thrust::count(target_ptr, target_ptr + nodeNum, 1);



	if (c_count > 0)
	{
		//printf("start...\n");
		//printf("count %d %d \n", c_count, i);
		double realThreshold = (2.5*nSampleWidth - nSampleWidth) / nSampleWidth;
		int gridRadius = (int)floor(realThreshold);
		krFDMContouring_Dilation << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (targetImg, assistImg, nodeNum, make_int3(imgRes.x, i, imgRes.y), realThreshold, gridRadius, i);



		unsigned int* anchorPtIndex;
		CUDA_SAFE_CALL(cudaMalloc((void**)&(anchorPtIndex), nodeNum * sizeof(unsigned int)));
		CUDA_SAFE_CALL(cudaMemset((void*)anchorPtIndex, 0, nodeNum * sizeof(unsigned int)));

		unsigned int* anchorPtValue;
		CUDA_SAFE_CALL(cudaMalloc((void**)&(anchorPtValue), nodeNum * sizeof(unsigned int)));
		CUDA_SAFE_CALL(cudaMemset((void*)anchorPtValue, 0, nodeNum * sizeof(unsigned int)));


		// Initialize flooding index
		thrust::device_ptr<unsigned int> index_ptr(anchorPtIndex);
		thrust::sequence(index_ptr, index_ptr + nodeNum, 0, 1);

		// Initialize flooding value
		krSLAContouring_FillAnchorValue << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (targetImg, assistImg, anchorPtValue, MAX_MARKER, nodeNum, i);

		// region flooding
		//if (i==87)
		regionNum = LDNIFlooding_Color2DFlooding(anchorPtIndex, assistImg, anchorPtValue, nodeNum, MAX_MARKER, imgRes);

		// get the anchor point

		thrust::fill(index_ptr, index_ptr + nodeNum, MAX_MARKER);
		krSLAContouring_GetAnchorPoint << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (targetImg, anchorPtValue, anchorPtIndex, imgRes, nodeNum, i);
		// fill image
		CUDA_SAFE_CALL(cudaMemset((void*)assistImg, false, nodeNum * sizeof(bool)));
		//if (i==68)
		krSLAContouring_FillImageValue << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (assistImg, tempImg, anchorPtIndex, imgRes, nodeNum, MAX_MARKER, i);
		CUDA_SAFE_CALL(cudaMemset((void*)targetImg, false, nodeNum * sizeof(bool)));

		//if (i==104)
		//	krFDMContouring_Test2D<<<BLOCKS_PER_GRID,THREADS_PER_BLOCK>>>(tempImg, nodeNum, imgRes);

		//
		LDNISLAContouring_GrowthAndSwallow(threshold, targetImg, assistImg, i, imageSize, nSampleWidth, disTextureA, disTextureB, disTexSize);





		cudaFree(anchorPtIndex);
		cudaFree(anchorPtValue);
		//printf("End...\n");

	}








}

extern "C" void call_krSLAContouring_Filter5(bool *gridNodes, bool *tempImg, bool *suptNodes, unsigned int *linkIndex, int nodeNum, int3 imageRes, int iy)
{
	krSLAContouring_Filter5 << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (gridNodes, tempImg, suptNodes, linkIndex, nodeNum, imageRes, iy);
}
__global__ void krSLAContouring_Filter5(bool *gridNodes, bool *tempImg, bool *suptNodes, unsigned int *linkIndex, int nodeNum, int3 imageRes, int iy)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int ix, iz, id, s_id, prevs_id;
	bool node, node2, node3;

	while (index < nodeNum) {

		ix = index % imageRes.x;	iz = index / imageRes.x;

		id = iz * imageRes.x*imageRes.y + iy * imageRes.x + ix;
		s_id = iz * imageRes.x*(imageRes.y - 1) + iy * imageRes.x + ix;
		prevs_id = iz * imageRes.x*(imageRes.y - 1) + (iy + 1)*imageRes.x + ix;

		node = tempImg[index];
		if (node)
		{
			atomicAdd(&linkIndex[index], 1);
		}

		if (iy + 1 >= imageRes.y - 1) node3 = false;
		else node3 = suptNodes[prevs_id];

		node2 = (node3 && !gridNodes[id]);

		suptNodes[s_id] = node || node2;

		index += blockDim.x * gridDim.x;

	}
}




__global__ void krSLAContouring_FillAnchorValue(bool* seedNodes, bool* inNodes, unsigned int* _value, unsigned int marker, int nodeNum, int iy)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int ind;

	while (index < nodeNum) {

		if (inNodes[index] || seedNodes[index])
		{
			//if (iy == 68)
			//	printf("fill anchor %d %d %d \n", index, index%213, index/213);
			_value[index] = index;

			//if (iy >= 330)
			//	printf("fill anchor %d %d %d \n", index, index%151, index/151);
		}
		else
		{
			_value[index] = marker;
		}


		index += blockDim.x * gridDim.x;
	}

}

__global__ void krSLAContouring_GetAnchorPoint(bool* targetImg, unsigned int* _value, unsigned int* anchorPt, int2 imgRes, int nodeNum, int iy)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int ind, ix, iz;

	while (index < nodeNum) {

		if (targetImg[index])
		{
			ix = index % imgRes.x; iz = index / imgRes.x;

			//if (iy == 68)
			//	printf("%d %d %d %d %d\n", ix, iz, index,_value[index], ix*imgRes.x+iz);

			//atomicMin( &anchorPt[_value[index]], ix*imgRes.x+iz);
			atomicMin(&anchorPt[_value[index]], ix*max(imgRes.x, imgRes.y) + iz);
		}

		index += blockDim.x * gridDim.x;
	}

}

__global__ void krSLAContouring_FillImageValue(bool* tempImg, bool* suptImg, unsigned int* _value, int2 imgRes, int nodeNum, unsigned int init, int iy)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int ix, iz, value;

	while (index < nodeNum) {

		value = _value[index];


		if (value < init)
		{
			ix = value % max(imgRes.x, imgRes.y);	iz = value / max(imgRes.x, imgRes.y);

			tempImg[ix*imgRes.x + iz] = true;
			suptImg[ix*imgRes.x + iz] = true;

		}


		index += blockDim.x * gridDim.x;
	}

}

void LDNIcudaOperation::LDNISLAContouring_BuildDistanceMap(bool *seedNodes, int i, int imageSize[], short2 *&disTexturesA, short2 *&disTexturesB, int disTexSize)
{
	int phase1Band = 16;
	int phase2Band = 16;
	int phase3Band = 16;

	krFDMContouring_InitializedValue << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (disTexturesA, disTexSize*disTexSize, MARKER);

	krSLAContouring_InitializedDistanceMap << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (seedNodes, disTexturesA, imageSize[0] * imageSize[2], make_int2(imageSize[0], imageSize[2]), i, disTexSize);

	dim3 block = dim3(BLOCKSIZE);
	dim3 grid = dim3(disTexSize / block.x, phase1Band);


	cudaBindTexture(0, disTexColor, disTexturesA);
	krFDMContouring_kernelFloodDown << < grid, block >> > (disTexturesB, disTexSize, disTexSize / phase1Band);
	cudaUnbindTexture(disTexColor);


	cudaBindTexture(0, disTexColor, disTexturesB);
	krFDMContouring_kernelFloodUp << < grid, block >> > (disTexturesB, disTexSize, disTexSize / phase1Band);


	//Phase 1:  Passing information between bands
	grid = dim3(disTexSize / block.x, phase1Band);
	krFDMContouring_kernelPropagateInterband << < grid, block >> > (disTexturesA, disTexSize, disTexSize / phase1Band);


	cudaBindTexture(0, disTexLinks, disTexturesA);
	krFDMContouring_kernelUpdateVertical << < grid, block >> > (disTexturesB, disTexSize, phase1Band, disTexSize / phase1Band);
	cudaUnbindTexture(disTexLinks);
	cudaUnbindTexture(disTexColor);


	//Phase 1: Transpose
	block = dim3(TILE_DIM, BLOCK_ROWS);
	grid = dim3(disTexSize / TILE_DIM, disTexSize / TILE_DIM);

	cudaBindTexture(0, disTexColor, disTexturesB);
	krFDMContouring_kernelTranspose << < grid, block >> > (disTexturesB, disTexSize);
	cudaUnbindTexture(disTexColor);



	//Phase 2: Compute proximate points locally in each band
	block = dim3(BLOCKSIZE);
	grid = dim3(disTexSize / block.x, phase2Band);
	cudaBindTexture(0, disTexColor, disTexturesB);
	krFDMContouring_kernelProximatePoints << < grid, block >> > (disTexturesA, disTexSize, disTexSize / phase2Band);
	cudaBindTexture(0, disTexLinks, disTexturesA);
	krFDMContouring_kernelCreateForwardPointers << < grid, block >> > (disTexturesA, disTexSize, disTexSize / phase2Band);

	//Phase 2:  Repeatly merging two bands into one
	for (int noBand = phase2Band; noBand > 1; noBand /= 2) {
		grid = dim3(disTexSize / block.x, noBand / 2);
		krFDMContouring_kernelMergeBands << < grid, block >> > (disTexturesA, disTexSize, disTexSize / noBand);
	}

	//Phase 2:  Replace the forward link with the X coordinate of the seed to remove the need of looking at the other texture. We need it for coloring.
	grid = dim3(disTexSize / block.x, disTexSize);
	krFDMContouring_kernelDoubleToSingleList << < grid, block >> > (disTexturesA, disTexSize);
	cudaUnbindTexture(disTexLinks);
	cudaUnbindTexture(disTexColor);


	//Phase 3: 
	block = dim3(BLOCKSIZE / phase2Band, phase2Band);
	grid = dim3(disTexSize / block.x);
	cudaBindTexture(0, disTexColor, disTexturesA);
	krFDMContouring_kernelColor << < grid, block >> > (disTexturesB, disTexSize);
	cudaUnbindTexture(disTexColor);


	//Phase 3: Transpose
	block = dim3(TILE_DIM, BLOCK_ROWS);
	grid = dim3(disTexSize / TILE_DIM, disTexSize / TILE_DIM);

	cudaBindTexture(0, disTexColor, disTexturesB);
	krFDMContouring_kernelTranspose << < grid, block >> > (disTexturesB, disTexSize);
	cudaUnbindTexture(disTexColor);

}

unsigned int LDNIcudaOperation::LDNIFlooding_Color2DFlooding(unsigned int *&index, bool *&key, unsigned int *&value, unsigned int arrsize, unsigned int init, int2 imgRes)
{
	bool bflag[4] = { true };
	bool result = true;
	unsigned int* verify;
	CUDA_SAFE_CALL(cudaMalloc((void**)&(verify), arrsize * sizeof(unsigned int)));
	CUDA_SAFE_CALL(cudaMemset((void*)verify, 0, arrsize * sizeof(unsigned int)));

	thrust::device_ptr<bool> key_ptr(key);
	thrust::device_ptr<unsigned int> index_ptr(index);
	thrust::device_ptr<unsigned int> value_ptr(value);
	thrust::device_ptr<unsigned int> verify_ptr(verify);

	thrust::equal_to<bool> binary_pred;
	thrust::minimum<unsigned int>  binary_op;




	while (true)
	{
		//Step 1: forward flooding
		{
			//Sorting value and key base on index
			thrust::copy(index_ptr, index_ptr + arrsize, verify_ptr);
			thrust::sort_by_key(index_ptr, index_ptr + arrsize, value_ptr);
			thrust::sort_by_key(verify_ptr, verify_ptr + arrsize, key_ptr);


			//Copying the value before scan and compare with the result after scan, to see any changes
			thrust::copy(value_ptr, value_ptr + arrsize, verify_ptr);
			thrust::inclusive_scan_by_key(key_ptr, key_ptr + arrsize, value_ptr, value_ptr, binary_pred, binary_op); // in-place scan
			bflag[0] = thrust::equal(value_ptr, value_ptr + arrsize, verify_ptr);
		}

		//Step 2: downward flooding
		{
			//transform the index
			thrust::transform(index_ptr, index_ptr + arrsize, index_ptr, transformOp(imgRes.x, imgRes.y));

			//Sorting value and key base on index
			thrust::copy(index_ptr, index_ptr + arrsize, verify_ptr);
			thrust::sort_by_key(index_ptr, index_ptr + arrsize, value_ptr);
			thrust::sort_by_key(verify_ptr, verify_ptr + arrsize, key_ptr);


			//Copying the value before scan and compare with the result after scan, to see any changes
			thrust::copy(value_ptr, value_ptr + arrsize, verify_ptr);
			thrust::inclusive_scan_by_key(key_ptr, key_ptr + arrsize, value_ptr, value_ptr, binary_pred, binary_op); // in-place scan
			bflag[1] = thrust::equal(value_ptr, value_ptr + arrsize, verify_ptr);

			thrust::transform(index_ptr, index_ptr + arrsize, index_ptr, transformOp(imgRes.y, imgRes.x));
		}


		//Step 3: backware flooding
		{

			thrust::reverse(value_ptr, value_ptr + arrsize);
			thrust::reverse(key_ptr, key_ptr + arrsize);

			thrust::copy(index_ptr, index_ptr + arrsize, verify_ptr);
			thrust::sort_by_key(index_ptr, index_ptr + arrsize, value_ptr);
			thrust::sort_by_key(verify_ptr, verify_ptr + arrsize, key_ptr);

			//Copying the value before scan and compare with the result after scan, to see any changes
			thrust::copy(value_ptr, value_ptr + arrsize, verify_ptr);
			thrust::inclusive_scan_by_key(key_ptr, key_ptr + arrsize, value_ptr, value_ptr, binary_pred, binary_op); // in-place scan
			bflag[2] = thrust::equal(value_ptr, value_ptr + arrsize, verify_ptr);

			//reverse back, same as index
			thrust::reverse(value_ptr, value_ptr + arrsize);
			thrust::reverse(key_ptr, key_ptr + arrsize);
		}



		//Step 3: upward flooding
		{
			thrust::transform(index_ptr, index_ptr + arrsize, index_ptr, transformOp(imgRes.x, imgRes.y));

			thrust::reverse(value_ptr, value_ptr + arrsize);
			thrust::reverse(key_ptr, key_ptr + arrsize);

			thrust::copy(index_ptr, index_ptr + arrsize, verify_ptr);
			thrust::sort_by_key(index_ptr, index_ptr + arrsize, value_ptr);
			thrust::sort_by_key(verify_ptr, verify_ptr + arrsize, key_ptr);

			//Copying the value before scan and compare with the result after scan, to see any changes
			thrust::copy(value_ptr, value_ptr + arrsize, verify_ptr);
			thrust::inclusive_scan_by_key(key_ptr, key_ptr + arrsize, value_ptr, value_ptr, binary_pred, binary_op); // in-place scan
			bflag[3] = thrust::equal(value_ptr, value_ptr + arrsize, verify_ptr);

			//reverse back, same as index
			thrust::reverse(value_ptr, value_ptr + arrsize);
			thrust::reverse(key_ptr, key_ptr + arrsize);
			//thrust::reverse(index_ptr, index_ptr + arrsize);

			//transform the index
			thrust::transform(index_ptr, index_ptr + arrsize, index_ptr, transformOp(imgRes.y, imgRes.x));
		}



		//printf("flood 1 %d %d %d %d\n", bflag[0], bflag[1], bflag[2], bflag[3]);
		//break;
		//result = bflag[0] || bflag[1] || bflag[2] || bflag[3]; // false = no more flooding area
		if (bflag[0] && bflag[1] && bflag[2] && bflag[3]) break;

	}

	thrust::copy(index_ptr, index_ptr + arrsize, verify_ptr);
	thrust::sort_by_key(index_ptr, index_ptr + arrsize, value_ptr);
	thrust::sort_by_key(verify_ptr, verify_ptr + arrsize, key_ptr);


	return 0;





}

__global__ void krSLAContouring_InitializedDistanceMap(bool *seedNodes, short2 *dMap, int nodeNum, int2 imageRes, int iy, int texwidth)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int ix, iz, id;
	while (index < nodeNum) {

		ix = index % imageRes.x;	iz = index / imageRes.x;
		id = iz * texwidth + ix;

		if (seedNodes[index])
		{

			dMap[id].x = ix;
			dMap[id].y = iz;
		}

		index += blockDim.x * gridDim.x;
	}
}

extern "C" void call_krSLAContouring_FindAllLinks(unsigned int *linkIndex, unsigned int *linkLayerC, short *linkLayerD, short2 *linkID, bool *tempImg, unsigned int *count, int nodeNum, int3 imageRes)
{
	krSLAContouring_FindAllLinks << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (linkIndex, linkLayerC, linkLayerD, linkID, tempImg, count, nodeNum, imageRes);
}
__global__ void krSLAContouring_FindAllLinks(unsigned int *linkIndex, unsigned int *linkLayerC, short *linkLayerD, short2 *linkID, bool *tempImg, unsigned int *count, int nodeNum, int3 imageRes)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int ix, iy, iz, id, st;
	bool node;

	while (index < nodeNum) {

		ix = index % imageRes.x;	iy = (index / imageRes.x) % imageRes.y;
		iz = index / (imageRes.x*imageRes.y);


		if (tempImg[index])
		{
			id = atomicAdd(&count[iz*imageRes.x + ix], 1);
			st = linkIndex[iz*imageRes.x + ix];

			linkLayerC[st + id] = iy;
			linkLayerD[st + id] = iy;
			linkID[st + id] = make_short2(ix, iz);

			//printf("find all %d %d %d %d %d\n", index, ix, iz, iy, st+id);
		}


		index += blockDim.x * gridDim.x;

	}
}

extern "C" void call_krSLAContouring_RelateAllLinksBetweenLayers(unsigned int *linkIndex, unsigned int *linkLayerC, bool *gridNodes, bool *suptNodes, int nodeNum, int3 imageRes)
{
	krSLAContouring_RelateAllLinksBetweenLayers << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (linkIndex, linkLayerC, gridNodes, suptNodes, nodeNum, imageRes);
}
__global__ void krSLAContouring_RelateAllLinksBetweenLayers(unsigned int *linkIndex, unsigned int *linkLayerC, bool *gridNodes, bool *suptNodes, int nodeNum, int3 imageRes)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int ix, iy, iz, id, s_id, st, num, i;
	bool node, node2;

	while (index < nodeNum) {

		ix = index % imageRes.x;	iy = (index / imageRes.x) % imageRes.y;
		iz = index / (imageRes.x*imageRes.y);


		s_id = iz * imageRes.x*(imageRes.y - 1) + (iy + 1)*imageRes.x + ix; // Mark: could be error when iy+1 >= imageRes.y-1

		if ((iy + 1) >= imageRes.y - 1) node2 = false;
		else node2 = suptNodes[s_id];

		if (node2 && !gridNodes[index])
		{
			st = linkIndex[iz*imageRes.x + ix];
			num = linkIndex[iz*imageRes.x + ix + 1] - st;

			for (i = 0; i < num; i++)
			{
				id = atomicMin(&linkLayerC[st + i], iy);
				//if (ix == 27 && iz == 48)
				//	printf("?? %d %d %d %d \n", ix, iz, iy, id);
			}

		}

		index += blockDim.x * gridDim.x;

	}
}

extern "C" void call_krFDMContouring_VerticalSpptPxlProp(bool *gridNodes, bool *suptNodes, int nodeNum, int3 imageRes, int iy)
{
	krFDMContouring_VerticalSpptPxlProp << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (gridNodes, suptNodes, nodeNum, imageRes, iy);
}
__global__ void krFDMContouring_VerticalSpptPxlProp(bool *gridNodes, bool *suptNodes, int nodeNum, int3 imageRes, int iy)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int ix, iz, id, prev_id, suptid;
	bool node;

	while (index < nodeNum) {
		ix = index % imageRes.x;
		iz = index / imageRes.x;

		id = iz * imageRes.x*imageRes.y + iy * imageRes.x + ix;
		suptid = iz * imageRes.x*(imageRes.y - 1) + iy * imageRes.x + ix;
		prev_id = iz * imageRes.x*(imageRes.y - 1) + (iy + 1)*imageRes.x + ix; // Mark: could be error when iy+1 >= imageRes.y-1

		if (suptNodes[prev_id] && !suptNodes[suptid] && !gridNodes[id])
			suptNodes[suptid] = true;

		index += blockDim.x * gridDim.x;
	}
}

__global__ void krSLAContouring_CounterLink1(short2 *linkID, bool *linkcount, int* count, int nodeNum, int linkNum, int linkThreshold)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int ix, iy, iz, id, st;
	short a, b, c, d;
	bool node;

	while (index < nodeNum) {

		ix = index % linkNum;	iy = index / linkNum;

		if (ix <= iy) {
			index += blockDim.x * gridDim.x;
			continue;
		}

		//printf("%d %d %d %d %d %d %d\n", ix, iy, linkNum, linkID[ix].x, linkID[ix].y, linkID[iy].x, linkID[iy].y);
		a = linkID[ix].x;
		b = linkID[ix].y;
		c = linkID[iy].x;
		d = linkID[iy].y;



		if ((a - c)*(a - c) + (b - d)*(b - d) <= linkThreshold * linkThreshold)
		{
			linkcount[index] = true;
			atomicAdd(count, 1);
			//printf("-- %d %d %d %d %d %d\n", a, b, c, d, ix, iy);
		}

		index += blockDim.x * gridDim.x;

	}
}

__global__ void krSLAContouring_FilterLink1(short2 *linkID, unsigned int *linklayerC, short *linklayerD, int* count, bool *linkcount, int nodeNum, int linkNum, int linkThreshold, int lengthofLayer, int furtherStepLength)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int ix, iy;
	int a, b, c, d, overlapStartinglayerID, overlapEndinglayerID, diff, pixelNum;
	short2 e, f;


	while (index < nodeNum) {
		ix = index % linkNum;	iy = index / linkNum;

		if (linkcount[index])
		{

			a = linklayerC[ix];
			b = linklayerC[iy];



			if (a >= b)
				overlapStartinglayerID = a;
			else
				overlapStartinglayerID = b;


			c = linklayerD[ix];
			d = linklayerD[iy];



			if (c <= d)
				overlapEndinglayerID = c;
			else
				overlapEndinglayerID = d;

			diff = overlapEndinglayerID - overlapStartinglayerID;

			e = linkID[ix];
			f = linkID[iy];




			if (abs(e.x - f.x) != abs(e.y - f.y))
				pixelNum = abs(e.x - f.x) + abs(e.y - f.y) - 1;
			else
				pixelNum = abs(e.x - f.x) - 1;


			//if (f.x == 59 || e.x == 59)
			//{
				//printf("!! %d %d %d %d %d %d %d %d %d %d %d %d\n",  linkID[ix].x,  linkID[ix].y,  linkID[iy].x,  linkID[iy].y , a, b, c, d, overlapEndinglayerID, overlapStartinglayerID, diff, pixelNum);


			//}

			if (diff <= (pixelNum - (lengthofLayer - 1) + 1) / furtherStepLength + 1)
			{
				linkcount[index] = false;

			}
			else
			{

				atomicAdd(count, 1);
			}
		}

		index += blockDim.x * gridDim.x;
	}

}

#define C_EPS 1.0e-8
__device__ bool _calTwoLinesIntersection(double3 l1, double3 l2, double pt[])
{
	double d = l1.x*l2.y - l2.x*l1.y;

	if (fabs(d) < C_EPS) return false;

	pt[0] = -(l1.z*l2.y - l2.z*l1.y) / d;
	pt[1] = (l1.z*l2.x - l2.z*l1.x) / d;

	return true;
}
__device__ double3 _calEquation(float2 v1, float2 v2)
{
	double3 p;
	p.x = v2.y - v1.y;
	p.y = v1.x - v2.x;
	double sqroot = sqrt(p.x*p.x + p.y*p.y);
	p.x = p.x / sqroot;
	p.y = p.y / sqroot;

	if (fabs(p.y) < C_EPS)
	{
		p.x = 1.0;
		p.y = 0.0;
		p.z = -v1.x;
		return p;
	}

	p.z = -(p.y*v1.y + p.x*v1.x);
	return p;


}

__device__ bool _calTwoLineSegmentsIntersection(float2 vMinus1, float2 v_a, float2 v_b, float2 vPlus1, double pt[])
{
	double3 L1, L2;


	L1 = _calEquation(v_a, v_b);
	L2 = _calEquation(vPlus1, vMinus1);


	if (!(_calTwoLinesIntersection(L1, L2, pt))) return false;


	double u1;
	if (fabs(vPlus1.x - vMinus1.x) <= C_EPS)
		u1 = (pt[1] - vPlus1.y) / (vMinus1.y - vPlus1.y);
	else
		u1 = (pt[0] - vPlus1.x) / (vMinus1.x - vPlus1.x);

	double u2;

	if (fabs(v_a.x - v_b.x) <= C_EPS)
		u2 = (pt[1] - v_a.y) / (v_b.y - v_a.y);
	else
		u2 = (pt[0] - v_a.x) / (v_b.x - v_a.x);

	if ((u1 >= 0.0) && (u1 <= 1.0) && (u2 >= 0.0) && (u2 <= 1.0)) return true;

	return false;

}

__global__ void krSLAContouring_CalConnectionMap(short2 *linkID, int i, int j, bool *bflag, int *pixelNum, int nodeNum, int2 imgRes,
	double nSampleWidth)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	short2 st, ed;
	unsigned int ix, iy;

	st = linkID[i];
	ed = linkID[j];

	float2 edgeSt, edgeEd;
	edgeSt.x = (st.x + 0.5)*nSampleWidth;		edgeSt.y = (st.y + 0.5)*nSampleWidth;
	edgeEd.x = (ed.x + 0.5)*nSampleWidth;		edgeEd.y = (ed.y + 0.5)*nSampleWidth;

	float2 edgeT1, edgeT2;
	double pt[2];

	bool intersect;
	//printf("st ed nnn %d %d %d %d \n", st.x, st.y, ed.x, ed.y);

	while (index < nodeNum) {
		ix = index % imgRes.x;	iy = index / imgRes.x;

		if (ix <= max(st.x, ed.x) && ix >= min(st.x, ed.x))
		{
			if (iy <= max(st.y, ed.y) && iy >= min(st.y, ed.y))
			{

				if ((ix == st.x && iy == st.y) || (ix == ed.x && iy == ed.y))
				{
					bflag[index] = true;
					//printf("1 nnn %d %d %d %d %d %d %d \n", st.x, st.y, ed.x, ed.y, ix, iy, index);
					index += blockDim.x * gridDim.x;
					continue;
				}

				edgeT1.x = (ix)*nSampleWidth; edgeT1.y = (iy)*nSampleWidth;
				edgeT2.x = (ix + 1)*nSampleWidth;	edgeT2.y = (iy)*nSampleWidth;

				intersect = _calTwoLineSegmentsIntersection(edgeT1, edgeSt, edgeEd, edgeT2, pt);

				if (intersect)
				{
					bflag[index] = true;
					//printf("2 nnn %d %d %d %d %d %d %d \n", st.x, st.y, ed.x, ed.y, ix, iy, index);
					atomicAdd(pixelNum, 1);
					index += blockDim.x * gridDim.x;
					continue;
				}

				edgeT1.x = (ix)*nSampleWidth; edgeT1.y = (iy)*nSampleWidth;
				edgeT2.x = (ix)*nSampleWidth;	edgeT2.y = (iy + 1)*nSampleWidth;

				intersect = _calTwoLineSegmentsIntersection(edgeT1, edgeSt, edgeEd, edgeT2, pt);

				if (intersect)
				{
					bflag[index] = true;
					//printf("3 nnn %d %d %d %d %d %d %d \n", st.x, st.y, ed.x, ed.y, ix, iy, index);
					atomicAdd(pixelNum, 1);
					index += blockDim.x * gridDim.x;
					continue;
				}

				edgeT1.x = (ix + 1)*nSampleWidth; edgeT1.y = (iy)*nSampleWidth;
				edgeT2.x = (ix + 1)*nSampleWidth;	edgeT2.y = (iy + 1)*nSampleWidth;

				intersect = _calTwoLineSegmentsIntersection(edgeT1, edgeSt, edgeEd, edgeT2, pt);

				if (intersect)
				{
					bflag[index] = true;
					//printf("4 nnn %d %d %d %d %d %d %d \n", st.x, st.y, ed.x, ed.y, ix, iy, index);
					atomicAdd(pixelNum, 1);
					index += blockDim.x * gridDim.x;
					continue;
				}

				edgeT1.x = (ix)*nSampleWidth; edgeT1.y = (iy + 1)*nSampleWidth;
				edgeT2.x = (ix + 1)*nSampleWidth;	edgeT2.y = (iy + 1)*nSampleWidth;

				intersect = _calTwoLineSegmentsIntersection(edgeT1, edgeSt, edgeEd, edgeT2, pt);

				if (intersect)
				{
					bflag[index] = true;
					//printf("5 nnn %d %d %d %d %d %d %d \n", st.x, st.y, ed.x, ed.y, ix, iy, index);
					atomicAdd(pixelNum, 1);
					index += blockDim.x * gridDim.x;
					continue;
				}

				edgeT1.x = (ix)*nSampleWidth; edgeT1.y = (iy + 0.5)*nSampleWidth;
				edgeT2.x = (ix + 1)*nSampleWidth;	edgeT2.y = (iy + 0.5)*nSampleWidth;

				intersect = _calTwoLineSegmentsIntersection(edgeT1, edgeSt, edgeEd, edgeT2, pt);

				if (intersect)
				{
					bflag[index] = true;
					// printf("6 nnn %d %d %d %d %d %d %d \n", st.x, st.y, ed.x, ed.y, ix, iy, index);
					atomicAdd(pixelNum, 1);
					index += blockDim.x * gridDim.x;
					continue;

				}
			}
		}



		index += blockDim.x * gridDim.x;
	}
}

__global__ void krSLAContouring_CheckConnectionMap(short2 *linkID, int i, int j, bool *bflag, int *pixelNum, int nodeNum, int2 imgRes)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	short2 st, ed;
	int ix, iy;
	bool btest1, btest2;

	st = linkID[i];
	ed = linkID[j];

	while (index < nodeNum) {
		ix = index % imgRes.x;	iy = index / imgRes.x;
		btest1 = false;
		btest2 = false;

		if (bflag[index])
		{
			if (ix <= max(st.x, ed.x) && ix >= min(st.x, ed.x))
			{
				if (iy <= max(st.y, ed.y) && iy >= min(st.y, ed.y))
				{
					if ((ix == st.x && iy == st.y) || (ix == ed.x && iy == ed.y))
					{
						index += blockDim.x * gridDim.x;
						continue;
					}

					if (ed.x > st.x)
					{
						if (ed.y > st.y)
						{
							if (bflag[(iy + 1)*imgRes.x + ix])
								btest1 = true;
							if (bflag[(iy)*imgRes.x + ix + 1])
								btest2 = true;

							if (!btest1 && !btest2)
							{
								bflag[(iy)*imgRes.x + ix + 1] = true;
								atomicAdd(pixelNum, 1);
							}
						}
						else if (ed.y < st.y)
						{
							if (bflag[(iy - 1)*imgRes.x + ix])
								btest1 = true;
							if (bflag[(iy)*imgRes.x + ix + 1])
								btest2 = true;

							if (!btest1 && !btest2)
							{
								bflag[(iy)*imgRes.x + ix + 1] = true;
								atomicAdd(pixelNum, 1);
							}
						}
					}
					else if (ed.x < st.x)
					{
						if (ed.y > st.y)
						{
							if (bflag[(iy + 1)*imgRes.x + ix])
								btest1 = true;
							if (bflag[(iy)*imgRes.x + ix - 1])
								btest2 = true;

							if (!btest1 && !btest2)
							{
								bflag[(iy)*imgRes.x + ix - 1] = true;
								atomicAdd(pixelNum, 1);
							}
						}
						else if (ed.y < st.y)
						{
							if (bflag[(iy - 1)*imgRes.x + ix])
								btest1 = true;
							if (bflag[(iy)*imgRes.x + ix - 1])
								btest2 = true;

							if (!btest1 && !btest2)
							{
								bflag[(iy)*imgRes.x + ix - 1] = true;
								atomicAdd(pixelNum, 1);
							}
						}
					}
				}
			}
		}



		index += blockDim.x * gridDim.x;
	}

}

__global__ void krSLAContouring_ConnectionMapOnLayers(bool* gridNodes, bool* tempImg, bool* linkMap, bool* bflagDelete, short2* linkID, int i, int j, int pixelNum, int lengthofLayer, int furtherStepLength, int layerNumforOneCircle, int endlayer, int startlayer, int3 imgRes, int nodeNum)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	int ix, iy, iz, m;
	short2 st, ed;
	int localID;
	int layerID;
	int norm;
	st = linkID[i];
	ed = linkID[j];
	bool bflag = false;

	while (index < nodeNum) {
		ix = index % imgRes.x;	iz = (index / imgRes.x) % imgRes.y;
		iy = index / (imgRes.x*imgRes.y);


		layerID = startlayer + 2 + iy;
		localID = iy % layerNumforOneCircle;
		if (layerNumforOneCircle - localID > endlayer - layerID)
		{
			index += blockDim.x * gridDim.x;
			//if (ix == 0 && iz == 0) printf("@@ %d %d %d %d %d %d\n", layerNumforOneCircle, localID, startlayer, endlayer, layerID, iy);
			continue;
		}

		bflag = false;
		if (linkMap[iz*imgRes.x + ix])//&& iy==4)
		{
			//norm = abs(iz - st.y) + abs(ix - st.x) - 1;
			norm = abs(iz - ed.y) + abs(ix - ed.x) - 1;
			//printf("3: %d %d %d %d %d -- %d %d %d\n", ix, iz, ed.x, ed.y, norm, localID, layerNumforOneCircle-1, pixelNum-1);

			if (norm >= 0 && !((localID == layerNumforOneCircle - 1) && (norm > pixelNum - 1)))
			{

				if (norm >= localID * furtherStepLength && norm < localID*furtherStepLength + lengthofLayer)
				{
					tempImg[index] = true;
					//printf("1: %d %d %d %d %d %d %d %d\n", norm, ix, iz, iy, localID, layerID, startlayer, endlayer);
					bflag = true;
				}
				if (norm <= pixelNum - 1 - (localID*furtherStepLength) && norm > pixelNum - 1 - (localID*furtherStepLength + lengthofLayer))
				{
					tempImg[index] = true;
					//printf("2: %d %d %d %d %d %d %d %d\n", norm, ix, iz, iy, st.x, st.y, ed.x, ed.y);
					bflag = true;
				}


				if (bflag)
				{
					if (gridNodes[iz*imgRes.x*imgRes.z + (layerID)*imgRes.x + ix])
					{
						for (m = 1; m <= layerNumforOneCircle - localID - 1; m++)
						{
							bflagDelete[iy + m] = true;
							//printf("delete %d %d %d %d %d \n", ix, iz, layerID, iy, iy+m);
							//printf("delete %d \n", layerID);
						}

						for (m = iy; m >= iy - localID; m--)
						{
							bflagDelete[m] = true;
							//printf("delete %d \n", layerID);
						}
					}
				}

			}
		}




		index += blockDim.x * gridDim.x;
	}

}

__global__ void krSLAContouring_DeleteMapOnLayers(bool* tempImg, bool* bflagDelete, int layerNumforOneCircle, int2 imgRes, int nodeNum)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	int ix, iy, iz;
	int localID;

	while (index < nodeNum) {
		ix = index % imgRes.x;	iz = (index / imgRes.x) % imgRes.y;
		iy = index / (imgRes.x*imgRes.y);

		localID = iy % layerNumforOneCircle;
		if (bflagDelete[iy] && tempImg[index])
		{
			tempImg[index] = false;
		}

		index += blockDim.x * gridDim.x;
	}


}

__global__ void krSLAContouring_FilterLink2(bool* suptNodes, bool* tempImgs, int startlayer, int endlayer, int3 imgRes, int nodeNum)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	int ix, iy, iz;
	int layerID;

	while (index < nodeNum) {
		ix = index % imgRes.x;	iz = (index / imgRes.x) % imgRes.y;
		iy = index / (imgRes.x*imgRes.y);

		layerID = iy + startlayer + 2;

		if (!suptNodes[iz*imgRes.x*imgRes.z + layerID * imgRes.x + ix] && tempImgs[index])
		{
			suptNodes[iz*imgRes.x*imgRes.z + layerID * imgRes.x + ix] = true;
			//if (iz == 191) printf("filter %d %d %d \n", ix, iy, iz);
		}



		index += blockDim.x * gridDim.x;
	}

}

void LDNIcudaOperation::LDNISLAContouring_GenerateConnectionforCylinders(unsigned int *linkLayerC, short *linkLayerD,
	short2 *linkID, bool *gridNodes, bool *&suptNodes, int imageSize[], int linkThreshold,
	int lengthofLayer, int furtherStepLength, int linkNum, double nSampleWidth)
{
	int *linkCount;
	bool *linkfilter;


	CUDA_SAFE_CALL(cudaMalloc((void**)&(linkCount), sizeof(unsigned int)));
	CUDA_SAFE_CALL(cudaMemset((void*)linkCount, 0, sizeof(unsigned int)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&(linkfilter), linkNum*linkNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)linkfilter, false, linkNum*linkNum * sizeof(bool)));

	krSLAContouring_CounterLink1 << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (linkID, linkfilter, linkCount, linkNum*linkNum, linkNum, linkThreshold);

	/*int *numofLink = (int*)malloc(sizeof(int));
	CUDA_SAFE_CALL( cudaMemcpy( numofLink, linkCount, sizeof(int), cudaMemcpyDeviceToHost ) );

	printf("num of link %d \n", numofLink[0]);*/

	CUDA_SAFE_CALL(cudaMemset((void*)linkCount, 0, sizeof(unsigned int)));
	krSLAContouring_FilterLink1 << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (linkID, linkLayerC, linkLayerD, linkCount, linkfilter, linkNum*linkNum, linkNum, linkThreshold, lengthofLayer, furtherStepLength);


	int *numofLink = (int*)malloc(sizeof(int));
	CUDA_SAFE_CALL(cudaMemcpy(numofLink, linkCount, sizeof(int), cudaMemcpyDeviceToHost));

	printf("num of link %d \n", numofLink[0]);


	bool *blink = (bool*)malloc(linkNum*linkNum * sizeof(bool));
	CUDA_SAFE_CALL(cudaMemcpy(blink, linkfilter, linkNum*linkNum * sizeof(bool), cudaMemcpyDeviceToHost));
	unsigned int *cpu_layerC = (unsigned int*)malloc(linkNum * sizeof(unsigned int));
	CUDA_SAFE_CALL(cudaMemcpy(cpu_layerC, linkLayerC, linkNum * sizeof(unsigned int), cudaMemcpyDeviceToHost));
	short *cpu_layerD = (short*)malloc(linkNum * sizeof(short));
	CUDA_SAFE_CALL(cudaMemcpy(cpu_layerD, linkLayerD, linkNum * sizeof(short), cudaMemcpyDeviceToHost));


	bool *linkMap;
	bool *tempLayer;
	bool *bflagLayer;
	CUDA_SAFE_CALL(cudaMalloc((void**)&(linkMap), imageSize[0] * imageSize[2] * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)linkMap, false, imageSize[0] * imageSize[2] * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)linkCount, 0, sizeof(int)));

	int m, n;
	int *numofPixel = (int*)malloc(sizeof(int));
	int st_layer, ed_layer, numofLayer;
	int layerNumforOneCircle;

	for (int i = 0; i < linkNum*linkNum; i++)
	{
		if (blink[i])
		{
			m = i % linkNum;	n = i / linkNum;
			numofPixel[0] = 0;


			krSLAContouring_CalConnectionMap << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (linkID, m, n, linkMap, linkCount, imageSize[0] * imageSize[2], make_int2(imageSize[0], imageSize[2]), nSampleWidth);


			krSLAContouring_CheckConnectionMap << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (linkID, m, n, linkMap, linkCount, imageSize[0] * imageSize[2], make_int2(imageSize[0], imageSize[2]));

			thrust::device_ptr<bool> target_ptr(linkMap);
			numofPixel[0] = thrust::count(target_ptr, target_ptr + (imageSize[0] * imageSize[2]), 1);




			numofPixel[0] = numofPixel[0] - 2 - 1;



			if (numofPixel[0] > lengthofLayer)
			{

				st_layer = cpu_layerC[m] >= cpu_layerC[n] ? cpu_layerC[m] : cpu_layerC[n];
				ed_layer = cpu_layerD[m] <= cpu_layerD[n] ? cpu_layerD[m] : cpu_layerD[n];
				numofLayer = abs(st_layer - ed_layer) + 1 - 2 - 1; // if(this->layerind<startinglayer+2 || this->layerind>= endinglayer)



				CUDA_SAFE_CALL(cudaMalloc((void**)&(tempLayer), numofLayer*imageSize[0] * imageSize[2] * sizeof(bool)));
				CUDA_SAFE_CALL(cudaMemset((void*)tempLayer, false, numofLayer*imageSize[0] * imageSize[2] * sizeof(bool)));
				CUDA_SAFE_CALL(cudaMalloc((void**)&(bflagLayer), numofLayer * sizeof(bool)));
				CUDA_SAFE_CALL(cudaMemset((void*)bflagLayer, false, numofLayer * sizeof(bool)));


				layerNumforOneCircle = (numofPixel[0] - (lengthofLayer - 1) + 1) / furtherStepLength + 1;

				int ccc = thrust::count(target_ptr, target_ptr + (imageSize[0] * imageSize[2]), 1);
				//printf("num of pixel %d %d %d %d \n", numofPixel[0],numofLayer, layerNumforOneCircle, ccc);

				krSLAContouring_ConnectionMapOnLayers << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (gridNodes, tempLayer, linkMap, bflagLayer, linkID, m, n, numofPixel[0], lengthofLayer, furtherStepLength, layerNumforOneCircle, ed_layer, st_layer, make_int3(imageSize[0], imageSize[2], imageSize[1]), imageSize[0] * imageSize[2] * numofLayer);
				krSLAContouring_DeleteMapOnLayers << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (tempLayer, bflagLayer, layerNumforOneCircle, make_int2(imageSize[0], imageSize[2]), imageSize[0] * imageSize[2] * numofLayer);
				krSLAContouring_FilterLink2 << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (suptNodes, tempLayer, st_layer, ed_layer, make_int3(imageSize[0], imageSize[2], imageSize[1] - 1), imageSize[0] * imageSize[2] * numofLayer);

				/**/
				cudaFree(bflagLayer);
				cudaFree(tempLayer);

				//break;
			}

			CUDA_SAFE_CALL(cudaMemset((void*)linkMap, false, imageSize[0] * imageSize[2] * sizeof(bool)));
			CUDA_SAFE_CALL(cudaMemset((void*)linkCount, 0, sizeof(int)));
		}

	}




	cudaFree(linkMap);
	cudaFree(linkfilter);
	cudaFree(linkCount);

	free(numofLink);
	free(blink);
	free(cpu_layerC);
	free(cpu_layerD);
	free(numofPixel);


}


void call_func::setdev_ptr(unsigned int* LinkIndex, int nodeNum, unsigned int& LinkNum)
{
	thrust::device_ptr<unsigned int> dev_ptr(LinkIndex); //	Wrap raw pointers with dev_ptr
	thrust::exclusive_scan(dev_ptr, dev_ptr + (nodeNum + 1), dev_ptr); //	in-place scan
	LinkNum = dev_ptr[nodeNum];
	printf("max links ----- %d\n", LinkNum);
}


void call_func::setdev_ptr2(unsigned int* linkLayerC, unsigned int LinkNum, int3 imgRes)
{
	thrust::device_ptr<unsigned int> c_ptr(linkLayerC);
	thrust::fill(c_ptr, c_ptr + LinkNum, imgRes.y + 1);
}