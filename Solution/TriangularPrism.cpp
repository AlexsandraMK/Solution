#include "IKE.h"

	TriangularPrism::TriangularPrism()
	{
		globalNumsKnots.resize(COUNT_KNOTS);
        knots.resize(COUNT_KNOTS);
		lambda = sigma = hi = 0;
		iterKnots = 0;

        base = NULL;
	}

    int TriangularPrism::GetCountKnots()
    {
        return COUNT_KNOTS;
    }

    vector<vector<double>> TriangularPrism::CalcMz()
    {
        double hz = knots[COUNT_KNOTS - 1].z - knots[0].z;
        double mult = hz / 6.;

        vector<vector<double>> Mz = 
        {
            {mult * 2., mult * 1.},
            {mult * 1., mult * 2.},
        };
      
        return Mz;
    }

    vector<vector<double>> TriangularPrism::CalcGz()
    {
        double hz = knots[COUNT_KNOTS - 1].z - knots[0].z;
        double mult = 1. / hz;

        vector<vector<double>> Gz =
        {
            {mult * 1., mult * (-1.)},
            {mult * (-1.), mult * 1.},
        };

        double sum = 0;
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                sum += Gz[i][j];

        if (sum != 0) cout << "Îøèáêà â Gz" << globalNumsKnots[0] << "\t" << globalNumsKnots[1] << "\t" << globalNumsKnots[2];

        return Gz;
    }

    int TriangularPrism::CalcMu(int ind) { return ind % 3; }
    int TriangularPrism::CalcNu(int ind) { return ind / 3; }

    void TriangularPrism::CreateBase()
    {
        base = new Triangle();
        int baseCountKnots = COUNT_KNOTS / 2;

        for (int i = 0; i < baseCountKnots; i++)
            base->SetGlobalKnotNum(globalNumsKnots[i], knots[i]);
    }

    vector<vector<double>> TriangularPrism::CalcLocalM()
    {
        if (base == NULL) CreateBase();

        vector<vector<double>> Mxy = base->CalcLocalM();

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

        return M;
    }

    vector<vector<double>> TriangularPrism::CalcLocalG()
    {
        if (base == NULL) CreateBase();
        vector<vector<double>> Mxy = base->CalcLocalM();
        vector<vector<double>> Gxy = base->CalcLocalG();

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

        double sum = 0;
        for (int i = 0; i < COUNT_KNOTS; i++)
            for (int j = 0; j < COUNT_KNOTS; j++)
                sum += G[i][j];

        if (sum > 1e-14) cout << "Îøèáêà â G" << globalNumsKnots[0] << "\t" << globalNumsKnots[1] << "\t" << globalNumsKnots[2];

        return G;
    }


    double TriangularPrism::SolveInPoint(Knot knot, vector<double> q)
    {
        double basesBase[] = { base->CountBasis(0,knot.x, knot.y),
                               base->CountBasis(1,knot.x, knot.y),
                               base->CountBasis(2,knot.x, knot.y) };

        double hz = knots[COUNT_KNOTS - 1].z - knots[0].z;

        double basesZ[] = { (knots[COUNT_KNOTS-1].z - knot.z) / hz,
                            (knot.z - knots[0].z) / hz };

        double u =  q[globalNumsKnots[0]] * basesBase[0] * basesZ[0] +
                    q[globalNumsKnots[1]] * basesBase[1] * basesZ[0] + 
                    q[globalNumsKnots[2]] * basesBase[2] * basesZ[0] + 
                    q[globalNumsKnots[3]] * basesBase[0] * basesZ[1] +
                    q[globalNumsKnots[4]] * basesBase[1] * basesZ[1] +
                    q[globalNumsKnots[5]] * basesBase[2] * basesZ[1];

        return u;
    }




 