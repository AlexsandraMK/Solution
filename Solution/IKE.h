#pragma once
#include "vectorInteractions.h"
#include "IntegrationFunc.h"
#include "Knot.h"
#include <string>
using namespace std;



class IKE
{
public:
	vector<int> globalNumsKnots;
	double lambda, sigma, hi;
	
	virtual void SetGlobalKnotNum(int numKnot, Knot coordinatesKnot) = 0;
	/*void SetGlobalKnotNum(int numKnot, Knot coordinatesKnot)
	{
		globalNumsKnots[iterKnots] = numKnot;
		knots[iterKnots] = coordinatesKnot;
		iterKnots++;
	}*/

	/// <summary>
	/// Вычисляет локальную матрицу Г без влияния параметров среды
	/// </summary>
	virtual vector<vector<double>> CalcLocalG() = 0;

	/// <summary>
	/// Вычисляет локальную матрицу М без влияния параметров среды
	/// </summary>
	virtual vector<vector<double>> CalcLocalM() = 0;

	vector<double> CalcLocalF(vector<double> fInPoints)
	{
		vector<double> f;
		f.resize(COUNT_KNOTS);
		vector<vector<double>> M = CalcLocalM();
		f = MultMatrByVect(M, fInPoints);

		return f;
	};

	/// <summary>
	/// Количество необходимых узлов
	/// </summary>
	virtual int GetCountKnots() = 0;

	virtual double SolveInPoint(Knot knot, vector<double> q) = 0;
	virtual bool IsIn(Knot knot) = 0;

protected:
	int COUNT_KNOTS;
	vector<Knot> knots;
	int iterKnots = 0;
};

class Triangle : public IKE
{
private:
	int COUNT_KNOTS = 3;
	vector<vector<double>> CalcAlfaMatrix();
	double CalcArea();
	double CalcDetD();
	double CalcDetD(Knot knotsTriangle[3]);

public:
	Triangle();
	int GetCountKnots();
	vector<vector<double>> CalcLocalG();
	vector<vector<double>> CalcLocalM();

	double CountBasis(int ind, double x, double y);

	double SolveInPoint(Knot knot, vector<double> q);


	bool IsIn(Knot knot);
	void SetGlobalKnotNum(int numKnot, Knot coordinatesKnot);
	static string ToString()
	{
		return "Triangle";
	}
};

class TriangularPrism : public IKE
{

private:
	int COUNT_KNOTS = 6;
	int CalcMu(int ind);
	int CalcNu(int ind);
	Triangle* base;

	vector<vector<double>> CalcMz();

	vector<vector<double>> CalcGz();

public:
	void CreateBase();
	TriangularPrism();
	int GetCountKnots();
	vector<vector<double>> CalcLocalG();
	double SolveInPoint(Knot knot, vector<double> q);
	vector<vector<double>> CalcLocalM();
	bool IsIn(Knot knot);
	void SetGlobalKnotNum(int numKnot, Knot coordinatesKnot);
	static string ToString()
	{
		return "TriangularPrism";
	}
};

class Hexagon : public IKE
{

private:
	int COUNT_KNOTS = 8;
	int CalcMu(int ind);
	int CalcNu(int ind);
	int CalcZeta(int ind);

	double CalcW(int ind, double alpha);

	double CalcPhi(int ind, vector<double> integrationVar);

	double DifferentiationPhi(int ind, NewAxis axis, vector<double> integrationVar);

	function<double(vector<double> integrationVar, int, int)> Mij = [this](vector<double> integrationVar, int i, int j)
	{
		double detJacobian = CalcDetMatrix(CalcJacobian(integrationVar));

		double res = CalcPhi(i, integrationVar) * CalcPhi(j, integrationVar) * detJacobian;
		return res;
	};

	function<double(vector<double> integrationVar, int, int)> Gij = [this](vector<double> integrationVar, int i, int j)
	{
		vector<vector<double>> Jacobian = CalcJacobian(integrationVar);
		double detJacobian = CalcDetMatrix(Jacobian);
		vector<vector<double>> reversed_Jacobian = CalcReverseMatrixWithSize3(Jacobian);

		vector<double> J_grad_i = MultMatrByVect(reversed_Jacobian, CalcGrad(i, integrationVar));
		vector<double> J_grad_j = MultMatrByVect(reversed_Jacobian, CalcGrad(j, integrationVar));

		double res = CalcScalar(J_grad_i, J_grad_j) *detJacobian;
		return res; // Исправлено
	};

	vector<vector<double>> CalcJacobian(vector<double> integrationVar);

	vector<double> CalcGrad(int ind, vector<double> integrationVar);

	vector<double> CalcF(vector<double> integrationVar, Knot knot);

	

public:
	int GetCountKnots();
	Hexagon();
	vector<vector<double>> CalcLocalG();
	vector<vector<double>> CalcLocalM();
	double SolveInPoint(Knot knot, vector<double> q);
	bool IsIn(Knot knot);
	void SetGlobalKnotNum(int numKnot, Knot coordinatesKnot);
	static string ToString()
	{
		return "Hexagon";
	}
};






