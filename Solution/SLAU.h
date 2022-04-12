#pragma once
#include "TimeScheme.h"

class SLAU
{
public:
	void SolveSLAU(InitialData* data, TimeScheme* scheme);
	SLAU(InitialData* data);
	void WriteResultForSolution(vector<double> q, double time);
	void WriteResultForTest(vector<double> q, double time);
	void SolveInArea(InitialData* data, double time);
	double SolveInPoint(InitialData* data, Knot knot);
	vector<vector<double>> A;
	vector<vector<double>> M;
	vector<vector<double>> G;
	vector<double> q;
	vector<double> d;
	vector<double> u;


protected:

	void LOC();
	vector<vector<double>> AssemblingGlobalM(InitialData* data);
	vector<vector<double>> AssemblingGlobalG(InitialData* data);
	void CalcA(InitialData* data, TimeScheme* scheme);
	void CalcD(InitialData* data, TimeScheme* scheme);
	void CalcU(InitialData* data, double time);
	vector<Knot> knots;
	vector<double> AssemblingGlobalF(InitialData* data, double time);
	void CalcFirstBoundaryConditions(InitialData* data, double time);
};
