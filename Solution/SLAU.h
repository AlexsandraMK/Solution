#pragma once
#include "TimeScheme.h"
#include <set>


class NoBlockSM
{

public:
	vector<int> ig, jg;
	vector<double> l, u, d;
	NoBlockSM(InitialData* data);
	int ReadSparseMatrix(string pathFile);

	void WriteSparseMatrix(string pathFile);

	void AddElement(int i, int j, double elem);

	vector<double> MultMatrByVect(vector<double> b);

	
};

class Block_3_SM
{

public:
	vector<int> ig, jg;
	vector<double**> l, u, d;
	int ReadSparseMatrix(string pathFile);
	void WriteSparseMatrix(string pathFile);
	void AddToBlock(double** to, double into[3][3]);
	void AddToBlock(double** to, double** into);
	void AddElement(int i, int j, double elem[3][3]);
	void AddElement(int i, int j, double** elem);
	vector<double> MultMatrByVect(vector<double> b);
	Block_3_SM(InitialData* data);
	void Clear();
};


class FEM
{
public:
	void SolveSLAU(InitialData* data, TimeScheme* scheme);

	FEM(InitialData* data);

	void WriteResultForSolution(vector<double> q, double time);
	void WriteResultForTest(vector<double> q, double time);
	void SolveInAreaForTest(InitialData* data, double time);
	void SolveInArea(InitialData* data, double time);
	double SolveInPoint(InitialData* data, Knot knot);
	Block_3_SM* A;
	NoBlockSM* M;
	Block_3_SM* G;
	vector<double> q;
	vector<double> qx;
	vector<double> qy;
	vector<double> qz;
	vector<double> d;
	vector<double> u;
	


protected:
	set<double> xArea, yArea, zArea;
	vector<int> KnotsIKEs;
	void LOC();
	void AssemblingGlobalM(NoBlockSM* globalMatrix, InitialData* data);
	void AssemblingGlobalBlockG(Block_3_SM* globalMatrix, InitialData* data);
	void AssemblingGlobalG_aa(NoBlockSM* globalMatrix, InitialData* data, axis a1, axis a2, int iCoeff);
	void CreateArea(InitialData* data);
	void FindIKeForArea(InitialData* data);
	int FindIKe(InitialData* data, Knot* knot);
	void CalcA(InitialData* data, TimeScheme* scheme);
	void CalcD(InitialData* data, TimeScheme* scheme);
	void CalcU(InitialData* data, double time);
	int slauSize;
	vector<Knot> knots;
	vector<double> AssemblingGlobalF(InitialData* data, double time);
	void CalcFirstBoundaryConditions(InitialData* data, double time);
	void Symmetrization(int i);
};
