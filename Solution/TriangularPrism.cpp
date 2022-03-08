#include "IKE.h"

	TriangularPrism::TriangularPrism()
	{
		globalNumsKnots.resize(COUNT_KNOTS);
		lambda = sigma = hi = 0;
		iterKnots = 0;
	}

    vector<vector<double>> TriangularPrism::CalcMz()
    {
        double hz = knots[COUNT_KNOTS - 1].z - knots[0].z;
        double mult = hz / 6;

        vector<vector<double>> Mz = 
        {
            {mult * 2, mult * 1},
            {mult * 1, mult * 2},
        };
      
        return Mz;
    }

    vector<vector<double>> TriangularPrism::CalcGz()
    {
        double hz = knots[COUNT_KNOTS - 1].z - knots[0].z;
        double mult = 1 / hz;

        vector<vector<double>> Gz =
        {
            {mult * 1, mult * (-1)},
            {mult * (-1), mult * 1},
        };

        return Gz;
    }

    int TriangularPrism::CalcMu(int ind) { return (ind - 1) % 3; }
    int TriangularPrism::CalcNu(int ind) { return (ind - 1) / 3; }

    Triangle TriangularPrism::CreateBase()
    {
        Triangle base;
        int baseCountKnots = COUNT_KNOTS / 2;

        for (int i = 0; i < baseCountKnots; i++)
            base.SetGlobalKnotNum(globalNumsKnots[i], knots[i]);

        return base;
    }

    vector<vector<double>> TriangularPrism::CalcLocalM()
    {
        Triangle base = CreateBase();

        vector<vector<double>> Mxy = base.CalcLocalM();

        vector<vector<double>> Mz = CalcMz();

        vector<vector<double>> M;
        M.resize(COUNT_KNOTS);
        int muI = 0, muJ = 0, nuI = 0, nuJ = 0;

        for (int i = 0; i < COUNT_KNOTS; i++)
        {
            M[i].resize(COUNT_KNOTS);
            muI = CalcMu(i);
            nuI = CalcNu(i);

            for (int j = 0; j < COUNT_KNOTS; j++)
            {
                muJ = CalcMu(j);
                nuJ = CalcNu(j);
                M[i][j] = Mxy[muI][muJ] * Mz[nuI][nuJ];
            }
        }
    }

    vector<vector<double>> TriangularPrism::CalcLocalG()
    {
        Triangle base = CreateBase();

        vector<vector<double>> Mxy = base.CalcLocalM();
        vector<vector<double>> Gxy = base.CalcLocalG();

        vector<vector<double>> Mz = CalcMz();
        vector<vector<double>> Gz = CalcGz();

        vector<vector<double>> G;
        G.resize(COUNT_KNOTS);
        int muI = 0, muJ = 0, nuI = 0, nuJ = 0;

        for (int i = 0; i < COUNT_KNOTS; i++)
        {
            G[i].resize(COUNT_KNOTS);
            muI = CalcMu(i);
            nuI = CalcNu(i);

            for (int j = 0; j < COUNT_KNOTS; j++)
            {
                muJ = CalcMu(j);
                nuJ = CalcNu(j);
                G[i][j] = Gxy[muI][muJ] * Mz[nuI][nuJ] + Mxy[muI][muJ] * Gz[nuI][nuJ];
            }
        }
    }





 