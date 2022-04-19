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

double Triangle::CalcDetD(Knot knotsTriangle[3])
{
	double detD = (knotsTriangle[1].x - knotsTriangle[0].x) * (knotsTriangle[2].y - knotsTriangle[0].y) - (knotsTriangle[2].x - knotsTriangle[0].x) * (knotsTriangle[1].y - knotsTriangle[0].y);
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

	double sum = 0;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			sum += G[i][j];

	if (sum > 1e-14) cout << "Ошибка в Gxy" << globalNumsKnots[0] << "\t" << globalNumsKnots[1] << "\t" << globalNumsKnots[2];

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

double Triangle::CountBasis(int ind, double x, double y)
{
	vector<vector<double>> alpha = CalcAlfaMatrix();

	return alpha[ind][0] + alpha[ind][1] * x + alpha[ind][2] * y;
}

double Triangle::SolveInPoint(Knot knot, vector<double> q)
{
	double bases[] = {	CountBasis(0,knot.x, knot.y),
						CountBasis(1,knot.x, knot.y),
						CountBasis(2,knot.x, knot.y) };

	return 0;
}

bool Triangle::IsIn(Knot knot)
{
	double sum = 0;
	for (int i = 0; i < COUNT_KNOTS; i++)
	{
		Knot* knotsTriangle = new Knot[3]{	knots[i % COUNT_KNOTS],
											knots[(i + 1) % COUNT_KNOTS],
											knot };
		sum += abs(CalcDetD(knotsTriangle));
		delete[] knotsTriangle;
	}

	double area = CalcArea();
	if (fabs(area - sum/2.) <= 1e-10) return true;

	return false;
}

void Triangle::SetGlobalKnotNum(int numKnot, Knot coordinatesKnot)
{
	globalNumsKnots[iterKnots] = numKnot;
	knots[iterKnots] = coordinatesKnot;
	iterKnots++;
}