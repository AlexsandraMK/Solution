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
	

	void SetGlobalKnotNum(int numKnot, Knot coordinatesKnot)
	{
		globalNumsKnots[iterKnots] = numKnot;
		knots[iterKnots] = coordinatesKnot;
		iterKnots++;
	}

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

public:
	Triangle();
	int GetCountKnots();
	vector<vector<double>> CalcLocalG();
	vector<vector<double>> CalcLocalM();

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

	Triangle CreateBase();

	vector<vector<double>> CalcMz();

	vector<vector<double>> CalcGz();

public:
	TriangularPrism();
	int GetCountKnots();
	vector<vector<double>> CalcLocalG();
	vector<vector<double>> CalcLocalM();
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

		return CalcPhi(i, integrationVar) * CalcPhi(j, integrationVar) * detJacobian;
	};

	function<double(vector<double> integrationVar, int, int)> Gij = [this](vector<double> integrationVar, int i, int j)
	{
		vector<vector<double>> reversed_Jacobian;

		vector<vector<double>> Jacobian = CalcJacobian(integrationVar);
		int sizeJacobian = Jacobian.size();
		reversed_Jacobian.resize(sizeJacobian);
		for (int i = 0; i < sizeJacobian; i++)	reversed_Jacobian[i].resize(sizeJacobian);

		for (int i = 0; i < sizeJacobian; i++)
		{
			for (int j = 0; j < sizeJacobian; j++)
			{
				double min[4]{};
				int k = 0;
				for (int im = 0; im < sizeJacobian; im++)
				{
					for (int jm = 0; jm < sizeJacobian; jm++)
					{
						if (im != i && jm != j)
							min[k++] = Jacobian[im][jm];
					}
				}

				reversed_Jacobian[j][i] = pow(-1, i + j + 2) * (min[0] * min[3] - min[1] * min[2]);
			}
		}



		double detJacobian = CalcDetMatrix(Jacobian);

		vector<double> J_grad_i = MultMatrByVect(reversed_Jacobian, CalcGrad(i, integrationVar));

		vector<double> J_grad_j = MultMatrByVect(reversed_Jacobian, CalcGrad(j, integrationVar));

		return CalcScalar(J_grad_i, J_grad_j) * detJacobian;
	};

	vector<vector<double>> CalcJacobian(vector<double> integrationVar);

	vector<double> CalcGrad(int ind, vector<double> integrationVar);

public:
	int GetCountKnots();
	Hexagon();
	vector<vector<double>> CalcLocalG();
	vector<vector<double>> CalcLocalM();

	static string ToString()
	{
		return "Hexagon";
	}
};






