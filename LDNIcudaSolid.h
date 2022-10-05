/*
 *  Copyright (C) 2014, Geometric Design and Manufacturing Lab in THE CHINESE UNIVERSITY OF HONG KONG
 *  All rights reserved.
 *   
 *		 http://ldnibasedsolidmodeling.sourceforge.net/
 *  
 *   
 *  Redistribution and use in source and binary forms, with or without modification, 
 *  are permitted provided that the following conditions are met:
 *  
 *  1. Redistributions of source code must retain the above copyright notice, 
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice, 
 *     this list of conditions and the following disclaimer in the documentation 
 *	   and/or other materials provided with the distribution.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 *   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
 *   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
 *   IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
 *   INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
 *   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
 *   OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
 *   WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 *   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY 
 *   OF SUCH DAMAGE.
 */

#ifndef	_CCL_LDNI_CUDA_SOLID
#define	_CCL_LDNI_CUDA_SOLID

#if defined(_WIN32) || defined(WIN32)
#include "GLLib\incldue\GL/glut.h"
#else
#include <GL/glut.h>
#endif


#define THREADS_PER_BLOCK			256
#define BLOCKS_PER_GRID				32

#define	MAX_NUM_OF_SAMPLES_ON_RAY	256

#define MAX(a,b)		(((a)>(b))?(a):(b))
#define MIN(a,b)		(((a)<(b))?(a):(b))
#define MAX3(a,b,c)		(MAX(MAX(a,b),MAX(a,c)))

#define LOGIC_UNION(insideA, insideB)	(insideA || insideB)
#define LOGIC_INTER(insideA, insideB)	(insideA && insideB)
#define LOGIC_SUBTR(insideA, insideB)	(insideA && (!insideB))

#define START_DEPTH(a)	((a > 0)? a - 1: a + 1)
#define END_DEPTH(a) ((a > 0)? a - 2: a + 2)

#define DEGREE_TO_ROTATE(x)		0.0174532922222*x
#define ROTATE_TO_DEGREE(x)		57.295780490443*x





class LDNIcudaSolid
{
public:
	LDNIcudaSolid() {m_res=m_xSampleNum=m_ySampleNum=m_zSampleNum=0;};
	virtual ~LDNIcudaSolid() {FreeMemory();};

	void MallocMemory(int res);
	void FreeMemory();
	void MallocSampleMemory(short nAxis, int sampleNum);
	void CleanUpSamples();

	bool FileSave(char *filename);
	bool FileRead(const char *filename);
	bool ExportLDNFile(char *filename);

	int GetSampleNumber() {return (m_xSampleNum+m_ySampleNum+m_zSampleNum);};
	int GetSampleNumber(short nDir) {int num;
		switch(nDir) {
		case 0:num=m_xSampleNum;break;
		case 1:num=m_ySampleNum;break;
		case 2:num=m_zSampleNum;break;
		}return num;
	};
	void SetSampleNumber(short nDir, int num) {
		switch(nDir) {
		case 0:m_xSampleNum=num;break;
		case 1:m_ySampleNum=num;break;
		case 2:m_zSampleNum=num;break;
		}
	};

	void SetOrigin(float ox, float oy, float oz) {m_origin[0]=ox; m_origin[1]=oy; m_origin[2]=oz;};
	void GetOrigin(float &ox, float &oy, float &oz) {ox=m_origin[0]; oy=m_origin[1]; oz=m_origin[2];};
	int GetResolution() {return m_res;};
	void SetResolution(int res) {m_res=res;};

	void SetSampleWidth(float width) {m_sampleWidth=width;};
	float GetSampleWidth() {return m_sampleWidth;};

	unsigned int* GetIndexArrayPtr(short nAxis) {return dev_indexArray[nAxis];};
	float* GetSampleDepthArrayPtr(short nAxis) {return dev_sampleDepthArray[nAxis];};
	float* GetSampleNxArrayPtr(short nAxis) {return dev_sampleNxArray[nAxis];};
	float* GetSampleNyArrayPtr(short nAxis) {return dev_sampleNyArray[nAxis];};
	
	void SetIndexArrayPtr(short nAxis, unsigned int* ptr) {dev_indexArray[nAxis]=ptr;};
	void SetSampleDepthArrayPtr(short nAxis, float* ptr) {dev_sampleDepthArray[nAxis]=ptr;};
	void SetSampleNxArrayPtr(short nAxis, float* ptr) {dev_sampleNxArray[nAxis]=ptr;};
	void SetSampleNyArrayPtr(short nAxis, float* ptr) {dev_sampleNyArray[nAxis]=ptr;};

	void CopyIndexArrayToHost(short nAxis, unsigned int* &hostIndexArray);
	void CopySampleArrayToHost(short nAxis, float* &hostSampleDepthArray, float* &hostSampleNxArray, float* &hostSampleNyArray);

	void BuildVBOforRendering(GLuint &m_vboPosition, GLuint &m_vboNormal, int &m_vertexNum, bool &m_cudaRegistered);

	void BuildSampleBasedIndexArray(short nAxis, unsigned int* &sampleBasedIndexArray);

	void SetBoundingBox(float box[]) { bbox[0] = box[0]; bbox[1] = box[1]; bbox[2] = box[2]; bbox[3] = box[3]; bbox[4] = box[4]; bbox[5] = box[5]; };
	void GetBoundingBox(float box[]) { box[0] = bbox[0]; box[1] = bbox[1]; box[2] = bbox[2]; box[3] = bbox[3]; box[4] = bbox[4]; box[5] = bbox[5]; };

private:
	unsigned int *dev_indexArray[3];	// dev_indexArray[nAxis][j*m_res+i] specify the starting index of ray(i,j) in the array of samples
										// dev_indexArray[nAxis][m_res*m_res] specifies the number of samples in the nAxis direction
	float *dev_sampleNxArray[3];
	float *dev_sampleNyArray[3];
	float *dev_sampleDepthArray[3];		// NOTE that: the sign of z-component of normal vector is stored in the depth value
	
	float m_origin[3],m_sampleWidth;	
	int m_res,m_xSampleNum,m_ySampleNum,m_zSampleNum;
	float bbox[6]; //Debbie 17/9
};

#endif
