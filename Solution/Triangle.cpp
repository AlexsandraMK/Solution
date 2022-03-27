#include "IKE.h"

Triangle::Triangle()
{
	globalNumsKnots.resize(COUNT_KNOTS);
	knots.resize(COUNT_KNOTS);
	lambda = sigma = hi = 0;
	iterKnots = 0;
}

int Triangle::GetCountKnots()
{
	return COUNT_KNOTS;
}

double Triangle::CalcDetD()
{
	double detD = (knots[1].x - knots[0].x) * (knots[2].y - knots[0].y) - (knots[2].x - knots[0].x) * (knots[1].y - knots[0].y);
	return detD;
}

vector<vector<double>> Triangle::CalcAlfaMatrix()
{
	Knot top1 = knots[0];
	Knot top2 = knots[1];
	Knot top3 = knots[2];
	double detD = CalcDetD();
	vector<vector<double>> alfaMatrix =
	{
		{(top2.x * top3.y - top3.x * top2.y) / detD,	(top2.y - top3.y) / detD,	(top3.x - top2.x) / detD},
		{(top3.x * top1.y - top1.x * top3.y) / detD,	(top3.y - top1.y) / detD,	(top1.x - top3.x) / detD},
		{(top1.x * top2.y - top2.x * top1.y) / detD,	(top1.y - top2.y) / detD,	(top2.x - top1.x) / detD}
	};

	return alfaMatrix;
}

double Triangle::CalcArea()
{
    double detD = CalcDetD();
    return 1. / 2. * abs(detD);
}

vector<vector<double>> Triangle::CalcLocalG()
{
	double area = CalcArea();
	vector<vector<double>> alfaMatrix = CalcAlfaMatrix();
	vector<vector<double>> G;
	G.resize(COUNT_KNOTS);

	for (int i = 0; i < COUNT_KNOTS; i++)
	{
		// Заполнение строки значением = площади треугольника
		G[i].resize(COUNT_KNOTS, area);
		for (int j = 0; j < COUNT_KNOTS; j++)
		{
			G[i][j] *= alfaMatrix[i][1] * alfaMatrix[j][1] + alfaMatrix[i][2] * alfaMatrix[j][2];
		}
	}

	return G;
}

vector<vector<double>> Triangle::CalcLocalM()
{
	double area = CalcArea();
	double mult = area / 12.;


	vector<vector<double>> M =
	{
		{mult * 2., mult * 1., mult * 1.},
		{mult * 1., mult * 2., mult * 1.},
		{mult * 1., mult * 1., mult * 2.}
	};

	return M;
}