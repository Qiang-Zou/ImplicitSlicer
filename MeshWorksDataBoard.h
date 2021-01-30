#ifndef _CCL_MESHWORKS_DATA_BOARD
#define _CCL_MESHWORKS_DATA_BOARD

class PMBody;

class MeshWorksDataBoard
{
public:
	MeshWorksDataBoard(void);
	virtual ~MeshWorksDataBoard(void);

	void InputOBJFile(char* filename);
	void OutputOBJFile(char* filename);

	void InputOFFFile(char* filename);

	void InputWIOFile(char* filename);
	void OutputWIOFile(char* filename);


public:
	PMBody *m_polyMeshBody;
	bool m_bVertexNormalShading;

private:
	void _mallocIndexArray(int ***&indexArray, int xNum, int yNum, int zNum);
	void _freeIndexArray(int ***&indexArray, int xNum, int yNum, int zNum);
};

#endif
