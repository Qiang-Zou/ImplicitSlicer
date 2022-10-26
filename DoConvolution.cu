#include "cuda.h"
#include "cutil.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "LDNIcudaSolid.h"
#include "call_cuda.h"
#include <Eigen/Dense>
extern __global__ void DoConvolution(SCALAR* digitImage, int gridIdxMinX, int gridIdxMinY, double cellSize, int nodeNum, int2 XY,int3 imgRes, double left, double right, Eigen::Vector3d p0, Eigen::Vector3d p1, double convolRadius, double radius, double layerLevel);

extern "C" void call_DoConvolution(SCALAR* digitImage, int gridIdxMinX, int gridIdxMinY, double cellSize, int nodeNum, int2 XY,int3 imgRes, double left, double right, Eigen::Vector3d p0, Eigen::Vector3d p1, double convolRadius, double radius, double layerLevel)
{	
	DoConvolution << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (digitImage,gridIdxMinX, gridIdxMinY, cellSize, nodeNum, XY,imgRes, left, right, p0,p1, convolRadius, radius, layerLevel);
}
__global__ void DoConvolution(SCALAR* digitImage,int gridIdxMinX,int gridIdxMinY,double cellSize,int nodeNum, int2 XY,int3 imgRes,double left,double right, Eigen::Vector3d p0, Eigen::Vector3d p1,double convolRadius,double radius,double layerLevel)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int rowIdx, colIdx;

	while (index < nodeNum) {
		rowIdx = index % XY.x+ gridIdxMinX;
		colIdx = index / XY.x+ gridIdxMinY;

		// compute convolution value for this point with regards to edge e
		Eigen::Vector3d p(left + cellSize * rowIdx, right + cellSize * colIdx, layerLevel);
		double l((p1 - p0).norm()), a((p0 - p).dot(p1 - p0)), b((p0 - p).norm());
		// compute line-sphere intersection
		l *= l; a *= 2; b = b * b - convolRadius * convolRadius;
		double delta = a * a - 4 * l * b;
		if (delta >= 1e-6) {
			double intersect1 = (-a - sqrt(delta)) / (2 * l);
			double intersect2 = (-a + sqrt(delta)) / (2 * l);
			double start = max(0., min(intersect1, intersect2));
			double end = min(1., max(intersect1, intersect2));
			if (start <= 1 && end >= 0) {
				l = sqrt(l); a /= -2; b += convolRadius * convolRadius; b = sqrt(b);
				SCALAR f = radius / (15 * pow(convolRadius, 4)) *
					(3 * pow(l, 4) * (pow(end, 5) - pow(start, 5)) -
						15 * a * l * l * (pow(end, 4) - pow(start, 4)) +
						20 * a * a * (pow(end, 3) - pow(start, 3)));
				digitImage[rowIdx+colIdx*imgRes.x] += f;
			}
		}
		// end of convolution

		index += blockDim.x * gridDim.x;
	}
}

