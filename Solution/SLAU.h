#pragma once
#include "TimeScheme.h"

class SLAU
{
public:
	void SolveSLAU(InitialData* data, TimeScheme* scheme);
	SLAU(InitialData* data);


	
	void WriteResultForSolution(vector<double> q, double time);
	void WriteResultForTest(vector<double> q, double time);
	void SolveInAreaForTest(InitialData* data, double time);
	void SolveInArea(InitialData* data, double time);
	int FindIKe(InitialData* data, Knot* knot);
	double SolveInPoint(InitialData* data, Knot knot);
	vector<vector<double>> A;
	vector<vector<double>> M;
	vector<vector<double>> G_xx_sigma;
	vector<vector<double>> G_yy_sigma;
	vector<vector<double>> G_zz_sigma;
	vector<vector<double>> G_xx;
	vector<vector<double>> G_yy;
	vector<vector<double>> G_zz;
	vector<vector<double>> G_xy;
	vector<vector<double>> G_xz;
	vector<vector<double>> G_yz;
	vector<double> q;
	vector<double> qx;
	vector<double> qy;
	vector<double> qz;
	vector<double> d;
	vector<double> u;


protected:

	void LOC();
	vector<vector<double>> AssemblingGlobalM(InitialData* data);
	vector<vector<double>> AssemblingGlobalG(InitialData* data);
	vector<vector<double>> AssemblingGlobalG_aa(InitialData* data, axis a1, axis a2, int iCoeff);
	void CalcA(InitialData* data, TimeScheme* scheme);
	void CalcD(InitialData* data, TimeScheme* scheme);
	void CalcU(InitialData* data, double time);
	int slauSize;
	vector<Knot> knots;
	vector<double> AssemblingGlobalF(InitialData* data, double time);
	void CalcFirstBoundaryConditions(InitialData* data, double time);
};
