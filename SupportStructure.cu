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
extern __global__ void krFDMContouring_Erosion(bool *gridNodes, bool* output, int nodeNum, int3 imageRes, double realThreshold, int gridRadius);

#define MARKER      -32768
#define BLOCKSIZE	64
#define TILE_DIM	32
#define BLOCK_ROWS	8
texture<short2> disTexColor;
texture<short2> disTexLinks;
#define TOID(x, y, size)    (__mul24((y), (size)) + (x))

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
