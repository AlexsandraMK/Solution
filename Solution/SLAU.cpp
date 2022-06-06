#include "SLAU.h"
#include <iostream>
#include <fstream>
#include "InputFuncs.h"








FEM::FEM(InitialData* data)
{
	cout << "Создаем СЛАУ:\n";
	slauSize = data->knots.size() * 3;
	knots.resize(data->knots.size());
	for (int i = 0; i < data->knots.size(); i++)
		knots[i] = data->knots[i];

	A = new Block_3_SM(data);
	q.resize(slauSize);
	qx.resize(data->knots.size());
	qy.resize(data->knots.size());
	qz.resize(data->knots.size());
	u.resize(slauSize);
	d.resize(slauSize);


	M = new NoBlockSM(data);
	if (M->ReadSparseMatrix("./matrix/Mhi.txt"))
	{
		cout << "\tПодсчет матрицы Mhi:\n";
		AssemblingGlobalM(M, data);
		M->WriteSparseMatrix("./matrix/Mhi.txt");
	}


	//G_xx = new NoBlockSM(data);
	//if (G_xx->ReadSparseMatrix( "./matrix/G_xx.txt"))
	//{
	//	cout << "\tПодсчет матрицы G_xx:\n";
	//	AssemblingGlobalG_aa(G_xx, data, x, x, 0);
	//	G_xx->WriteSparseMatrix("./matrix/G_xx.txt");
	//}

	//G_yy = new NoBlockSM(data);
	//if (G_yy->ReadSparseMatrix("./matrix/G_yy.txt"))
	//{
	//	cout << "\tПодсчет матрицы G_yy:\n";
	//	AssemblingGlobalG_aa(G_yy, data, y, y, 0);
	//	G_yy->WriteSparseMatrix("./matrix/G_yy.txt");
	//}

	//G_zz = new NoBlockSM(data);
	//if (G_zz->ReadSparseMatrix("./matrix/G_zz.txt"))
	//{
	//	cout << "\tПодсчет матрицы G_zz:\n";
	//	AssemblingGlobalG_aa(G_zz, data, z, z, 0);
	//	G_zz->WriteSparseMatrix("./matrix/G_zz.txt");
	//}
	//
	//G_xx_sigma = new NoBlockSM(data);
	//if (G_xx_sigma->ReadSparseMatrix("./matrix/G_xx_sigma.txt"))
	//{

	//	cout << "\tПодсчет матрицы G_xx_sigma:\n";
	//	AssemblingGlobalG_aa(G_xx_sigma, data, x, x, 1);
	//	G_xx_sigma->WriteSparseMatrix("./matrix/G_xx_sigma.txt");
	//}

	//G_yy_sigma = new NoBlockSM(data);
	//if (G_yy_sigma->ReadSparseMatrix("./matrix/G_yy_sigma.txt"))
	//{
	//	cout << "\tПодсчет матрицы G_yy_sigma:\n";
	//	AssemblingGlobalG_aa(G_yy_sigma, data, y, y, 1);
	//	G_yy_sigma->WriteSparseMatrix("./matrix/G_yy_sigma.txt");
	//}

	//G_zz_sigma = new NoBlockSM(data);
	//if (G_zz_sigma->ReadSparseMatrix("./matrix/G_zz_sigma.txt"))
	//{
	//	cout << "\tПодсчет матрицы G_zz_sigma:\n";
	//	AssemblingGlobalG_aa(G_zz_sigma, data, z, z, 1);
	//	G_zz_sigma->WriteSparseMatrix("./matrix/G_zz_sigma.txt");
	//}
	//
	//G_xy = new NoBlockSM(data);
	//if (G_xy->ReadSparseMatrix("./matrix/G_xy.txt"))
	//{
	//	cout << "\tПодсчет матрицы G_xy:\n";
	//	AssemblingGlobalG_aa(G_xy, data, x, y, 0);
	//	G_xy->WriteSparseMatrix("./matrix/G_xy.txt");
	//}

	//G_xz = new NoBlockSM(data);
	//if (G_xz->ReadSparseMatrix("./matrix/G_xz.txt"))
	//{
	//	cout << "\tПодсчет матрицы G_xz:\n";
	//	AssemblingGlobalG_aa(G_xz, data, x, z, 0);
	//	G_xz->WriteSparseMatrix("./matrix/G_xz.txt");
	//}

	//G_yz = new NoBlockSM(data);
	//if (G_yz->ReadSparseMatrix("./matrix/G_yz.txt"))
	//{
	//	cout << "\tПодсчет матрицы G_yz:\n";
	//	AssemblingGlobalG_aa(G_yz, data, y, z, 0);
	//	G_yz->WriteSparseMatrix("./matrix/G_yz.txt");
	//}
	G = new Block_3_SM(data);
	if (G->ReadSparseMatrix("./matrix/G.txt"))
	{
		cout << "\tПодсчет матрицы G:\n";
		AssemblingGlobalBlockG(G, data);
		G->WriteSparseMatrix("./matrix/G.txt");
	}
}

void FEM::AssemblingGlobalBlockG(Block_3_SM* globalMatrix, InitialData* data)
{
	double block[3][3]{};
	int nKEs = data->KEs.size();
	for (int iKE = 0; iKE < nKEs; iKE++)
	{
		IKE* ke = data->KEs[iKE];
		vector<vector<double>> G_xx = ke->CalcLocalG_aa(x, x);
		vector<vector<double>> G_xy = ke->CalcLocalG_aa(x, y);
		vector<vector<double>> G_xz = ke->CalcLocalG_aa(x, z);
		vector<vector<double>> G_yx = ke->CalcLocalG_aa(y, x);
		vector<vector<double>> G_yy = ke->CalcLocalG_aa(y, y);
		vector<vector<double>> G_yz = ke->CalcLocalG_aa(y, z);
		vector<vector<double>> G_zx = ke->CalcLocalG_aa(z, x);
		vector<vector<double>> G_zy = ke->CalcLocalG_aa(z, y);
		vector<vector<double>> G_zz = ke->CalcLocalG_aa(z, z);

		double sigma = data->coeffs[ke->iCoeff].Yung / (2. * (1. + data->coeffs[ke->iCoeff].Poisson));
		double lambda = data->coeffs[ke->iCoeff].Poisson * data->coeffs[ke->iCoeff].Yung / (1. + data->coeffs[ke->iCoeff].Poisson) / (1. - 2. * data->coeffs[ke->iCoeff].Poisson);
		int countKnotsInKe = ke->GetCountKnots();

		for (int i = 0; i < countKnotsInKe; i++)
		{
			int globalI = ke->globalNumsKnots[i];
			for (int j = 0; j < countKnotsInKe; j++)
			{
				int globalJ = ke->globalNumsKnots[j];
				block[0][0] = (lambda + 2 * sigma) * G_xx[i][j] + sigma * (G_yy[i][j] + G_zz[i][j]);
				block[0][1] = lambda * G_xy[i][j] + sigma * G_yx[i][j];
				block[0][2] = lambda * G_xz[i][j] + sigma * G_zx[i][j];

				block[1][0] = lambda * G_yx[i][j] + sigma * G_xy[i][j];
				block[1][1] = (lambda + 2 * sigma) * G_yy[i][j] + sigma * (G_xx[i][j] + G_zz[i][j]);
				block[1][2] = lambda * G_yz[i][j] + sigma * G_zy[i][j];
				block[2][0] = lambda * G_zx[i][j] + sigma * G_xz[i][j];
				block[2][1] = lambda * G_zy[i][j] + sigma * G_yz[i][j];
				block[2][2] = (lambda + 2 * sigma) * G_zz[i][j] + sigma * (G_xx[i][j] + G_yy[i][j]);

				globalMatrix->AddElement(globalI, globalJ, block);
			}
		}
	}
	
}

void  FEM::AssemblingGlobalG_aa(NoBlockSM* globalMatrix, InitialData* data, axis a1, axis a2, int iCoeff)
{
	if (a1 == a2)
	{
		int nKEs = data->KEs.size();
		for (int iKE = 0; iKE < nKEs; iKE++)
		{
			IKE* ke = data->KEs[iKE];
			vector<vector<double>> localMatrix = ke->CalcLocalG_aa(a1, a2);
			double sigma = data->coeffs[ke->iCoeff].Yung / (2. * (1. + data->coeffs[ke->iCoeff].Poisson));
			double lambda = data->coeffs[ke->iCoeff].Poisson * data->coeffs[ke->iCoeff].Yung / (1. + data->coeffs[ke->iCoeff].Poisson) / (1. - 2. * data->coeffs[ke->iCoeff].Poisson);

			double coeff = iCoeff == 0 ? lambda + sigma : sigma;
			string file = "./locals/localG/localG" + std::to_string(iKE) + "_" + std::to_string(a1) + "_" + std::to_string(a2) + ".txt";
			WriteMatrix(localMatrix, file);
			int countKnotsInKe = ke->GetCountKnots();

			for (int j = 0; j < countKnotsInKe; j++)
			{
				int globalJ = ke->globalNumsKnots[j];
				for (int k = 0; k < countKnotsInKe; k++)
				{
					int globalK = ke->globalNumsKnots[k];
					globalMatrix->AddElement(globalJ, globalK, coeff * localMatrix[j][k]);
				}
			}
		}
	}
	else
	{
		int nKEs = data->KEs.size();
		for (int iKE = 0; iKE < nKEs; iKE++)
		{
			IKE* ke = data->KEs[iKE];
			vector<vector<double>> localMatrix = ke->CalcLocalG_aa(a1, a2);
			vector<vector<double>> localMatrix2 = ke->CalcLocalG_aa(a2, a1);
			double sigma = data->coeffs[ke->iCoeff].Yung / (2. * (1. + data->coeffs[ke->iCoeff].Poisson));
			double lambda = data->coeffs[ke->iCoeff].Poisson * data->coeffs[ke->iCoeff].Yung / (1. + data->coeffs[ke->iCoeff].Poisson) / (1. - 2. * data->coeffs[ke->iCoeff].Poisson);

			string file = "./locals/localG/localG" + std::to_string(iKE) + "_" + std::to_string(a1) + "_" + std::to_string(a2) + ".txt";
			WriteMatrix(localMatrix, file);
			int countKnotsInKe = ke->GetCountKnots();

			for (int j = 0; j < countKnotsInKe; j++)
			{
				int globalJ = ke->globalNumsKnots[j];
				for (int k = 0; k < countKnotsInKe; k++)
				{
					int globalK = ke->globalNumsKnots[k];
					globalMatrix->AddElement(globalJ, globalK, lambda * localMatrix[j][k] + sigma * localMatrix2[j][k]);
				}
			}
		}
	}
	

}


void FEM::LOC()
{
	cout << "\n\nВычисления СЛАУ с помощью ЛОС:\n";

	int slauSize = q.size();
	int i;
	double nvzk = 0., alfa = 0., beta = 0., skp = 0., eps = 9.999999682655226e-016;

	vector<double> z, r, p, f;
	z.resize(slauSize);
	r.resize(slauSize);
	p.resize(slauSize);
	f.resize(slauSize);

	double lastnvzk;
	for (i = 0; i < slauSize; i++)
		q[i] = 0.;
	f = A->MultMatrByVect(q);
	f = A->MultMatrByVect(q);

	for (i = 0; i < slauSize; i++)
		z[i] = r[i] = d[i] - f[i];

	p = A->MultMatrByVect(z);
	nvzk = sqrt(CalcScalar(r, r)) / sqrt(CalcScalar(d, d));
	unsigned int k;
	for (k = 1; k < 10000000 && nvzk > eps; k++)
	{
		lastnvzk = nvzk;
		skp = CalcScalar(p, p);
		alfa = CalcScalar(p, r) / skp;

		for (i = 0; i < slauSize; i++)
		{
			q[i] += alfa * z[i];
			r[i] -= alfa * p[i];
		}

		f = A->MultMatrByVect(r);
		beta = -CalcScalar(p, f) / skp;

		for (i = 0; i < slauSize; i++)
		{
			z[i] = r[i] + beta * z[i];
			p[i] = f[i] + beta * p[i];
		}

		nvzk = sqrt(CalcScalar(r, r)) / sqrt(CalcScalar(d, d));
		if (k % 10000 == 0) cout << "iterations: " << k << "\tНевязка:" << nvzk << endl;
	}
	cout << "iterations: " << k << "\tНевязка:" << nvzk << endl;
	cout << "Вычисление СЛАУ закончено.\n\n\n";

}


/// <summary>
/// Вычисление глобальных матриц массы пространственной сетки
/// </summary>
/// <param name="data">Входные данные</param>
void FEM::AssemblingGlobalM(NoBlockSM* globalMatrix, InitialData* data)
{
	int sizeGlobalMatrix = data->knots.size();
	int nKEs = data->KEs.size();

	for (int iKE = 0; iKE < nKEs; iKE++)
	{
		IKE* ke = data->KEs[iKE];

		double hi = /*0;*/data->coeffs[ke->iCoeff].density;

		vector<vector<double>> localMatrix = ke->CalcLocalM();
		//string file = "./locals/localM" + std::to_string(iKE) + ".txt";
		//WriteMatrix(localMatrix, file);
		int countKnotsInKe = ke->GetCountKnots();

		for (int j = 0; j < countKnotsInKe; j++)
		{
			int globalJ = ke->globalNumsKnots[j];
			for (int k = 0; k < countKnotsInKe; k++)
			{
				int globalK = ke->globalNumsKnots[k];
				globalMatrix->AddElement(globalJ, globalK, hi * localMatrix[j][k]);
			}
		}
	}
}


double CalcTimeCoeff3ForHi(TimeScheme* scheme, int nTime)
{
	vector<double> timeToCalc = scheme->time;

	double hiMult = 0;
	double hiDel = 1.;

	for (int i = 1; i < timeToCalc.size();i++)
	{
		if (i != nTime)
		{
			//hiMult += timeToCalc[3] - timeToCalc[i];
			hiDel *= timeToCalc[nTime] - timeToCalc[i];
		}
	}

	return 2. /** hiMult */ / hiDel;
}

double CalcTimeCoeff4ForHi(TimeScheme* scheme, int nTime)
{
	vector<double> timeToCalc = scheme->time;

	double hiMult = 0;
	double hiDel = 1.;

	for (int i = 0; i < timeToCalc.size();i++)
	{
		if (i != nTime)
		{
			hiMult += timeToCalc[3] - timeToCalc[i];
			hiDel *= timeToCalc[nTime] - timeToCalc[i];
		}
	}

	return 2. * hiMult / hiDel;
}

void FEM::CalcA(InitialData* data, TimeScheme* scheme) // Вычисление глобальной A
{
	vector<double> timeToCalc = scheme->time;

	for (int i = 0; i < data->knots.size(); i++)
	{
		A->AddElement(i, i, G->d[i]);
		for (int j = 0; j < 3; j++)
			A->d[i][j][j] += M->d[i] * CalcTimeCoeff4ForHi(scheme, 3);
		/*A->d[i][0][0] = **G_xx->d[i] + **G_xx_sigma->d[i] + **G_yy_sigma->d[i] + **G_zz_sigma->d[i]
				+ **M->d[i] * CalcTimeCoeff4ForHi(scheme, 3);

		A->d[i][1][1] = **G_yy->d[i] + **G_xx_sigma->d[i] + **G_yy_sigma->d[i] + **G_zz_sigma->d[i]
			+ **M->d[i] * CalcTimeCoeff4ForHi(scheme, 3);

		A->d[i][2][2] = **G_zz->d[i] + **G_xx_sigma->d[i] + **G_yy_sigma->d[i] + **G_zz_sigma->d[i]
			+ **M->d[i] * CalcTimeCoeff4ForHi(scheme, 3);

		A->d[i][0][1] = A->d[i][1][0] = **G_xy->d[i];
		A->d[i][0][2] = A->d[i][2][0] = **G_xz->d[i];
		A->d[i][1][2] = A->d[i][2][1] = **G_yz->d[i];*/
	}

	for (int i = 0; i < A->l.size(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3;k++)
			{
				A->l[i][j][k] = G->l[i][j][k];
				A->u[i][j][k] = G->u[i][j][k];
			}

			A->l[i][j][j] += M->l[i] * CalcTimeCoeff4ForHi(scheme, 3);
			A->u[i][j][j] += M->l[i] * CalcTimeCoeff4ForHi(scheme, 3);
		}

		
		/*A->l[i][0][0] 
		A->l[i][0][0] = **G_xx->l[i] + **G_xx_sigma->l[i] + **G_yy_sigma->l[i] + **G_zz_sigma->l[i]
			+ **M->l[i] * CalcTimeCoeff4ForHi(scheme, 3);

		A->l[i][1][1] = **G_yy->l[i] + **G_xx_sigma->l[i] + **G_yy_sigma->l[i] + **G_zz_sigma->l[i]
			+ **M->l[i] * CalcTimeCoeff4ForHi(scheme, 3);

		A->l[i][2][2] = **G_zz->l[i] + **G_xx_sigma->l[i] + **G_yy_sigma->l[i] + **G_zz_sigma->l[i]
			+ **M->l[i] * CalcTimeCoeff4ForHi(scheme, 3);

		A->l[i][0][1] = A->l[i][1][0] = **G_xy->l[i];
		A->l[i][0][2] = A->l[i][2][0] = **G_xz->l[i];
		A->l[i][1][2] = A->l[i][2][1] = **G_yz->l[i];

		A->u[i][0][0] = **G_xx->u[i] + **G_xx_sigma->u[i] + **G_yy_sigma->u[i] + **G_zz_sigma->u[i]
			+ **M->l[i] * CalcTimeCoeff4ForHi(scheme, 3);

		A->u[i][1][1] = **G_yy->u[i] + **G_xx_sigma->u[i] + **G_yy_sigma->u[i] + **G_zz_sigma->u[i]
			+ **M->u[i] * CalcTimeCoeff4ForHi(scheme, 3);

		A->u[i][2][2] = **G_zz->u[i] + **G_xx_sigma->u[i] + **G_yy_sigma->u[i] + **G_zz_sigma->u[i]
			+ **M->u[i] * CalcTimeCoeff4ForHi(scheme, 3);

		A->u[i][0][1] = A->u[i][1][0] = **G_xy->u[i];
		A->u[i][0][2] = A->u[i][2][0] = **G_xz->u[i];
		A->u[i][1][2] = A->u[i][2][1] = **G_yz->u[i];*/
	}

}





void FEM::CalcD(InitialData* data, TimeScheme* scheme) // Вычисление глобальной A
{
	vector<double> timeToCalc = scheme->time;


	vector<double> b;/*AssemblingGlobalF(data, timeToCalc[3]);*/
	b.resize(data->knots.size()*3, 0);
	vector<double>Mq_time0_x, Mq_time1_x, Mq_time2_x, Mq_time0_y, Mq_time1_y, Mq_time2_y, Mq_time0_z, Mq_time1_z, Mq_time2_z;

	Mq_time0_x = M->MultMatrByVect(scheme->qx[0]);
	Mq_time1_x = M->MultMatrByVect(scheme->qx[1]);
	Mq_time2_x = M->MultMatrByVect(scheme->qx[2]);

	Mq_time0_y = M->MultMatrByVect(scheme->qy[0]);
	Mq_time1_y = M->MultMatrByVect(scheme->qy[1]);
	Mq_time2_y = M->MultMatrByVect(scheme->qy[2]);

	Mq_time0_z = M->MultMatrByVect(scheme->qz[0]);
	Mq_time1_z = M->MultMatrByVect(scheme->qz[1]);
	Mq_time2_z = M->MultMatrByVect(scheme->qz[2]);

	for (int i = 0; i < data->knots.size(); i++)
	{

		d[i * 3] = b[i * 3]
			- Mq_time0_x[i] * CalcTimeCoeff4ForHi(scheme, 0)
			- Mq_time1_x[i] * CalcTimeCoeff4ForHi(scheme, 1)
			- Mq_time2_x[i] * CalcTimeCoeff4ForHi(scheme, 2);

		d[i * 3 + 1] = b[i * 3 + 1]
			- Mq_time0_y[i] * CalcTimeCoeff4ForHi(scheme, 0)
			- Mq_time1_y[i] * CalcTimeCoeff4ForHi(scheme, 1)
			- Mq_time2_y[i] * CalcTimeCoeff4ForHi(scheme, 2);

		d[i * 3 + 2] = b[i * 3 + 2]
			- Mq_time0_z[i] * CalcTimeCoeff4ForHi(scheme, 0)
			- Mq_time1_z[i] * CalcTimeCoeff4ForHi(scheme, 1)
			- Mq_time2_z[i] * CalcTimeCoeff4ForHi(scheme, 2);




		//d[i * 3] = b[i * 3]
		//	- Mq_time2_x[i] * CalcTimeCoeff3ForHi(scheme, 2)
		//	- Mq_time1_x[i] * CalcTimeCoeff3ForHi(scheme, 1);

		//d[i * 3 + 1] = b[i * 3 + 1]
		//	- Mq_time2_y[i] * CalcTimeCoeff3ForHi(scheme, 2)
		//	- Mq_time1_y[i] * CalcTimeCoeff3ForHi(scheme, 1);

		//d[i * 3 + 2] = b[i * 3 + 2]
		//	- Mq_time2_z[i] * CalcTimeCoeff3ForHi(scheme, 2)
		//	- Mq_time1_z[i] * CalcTimeCoeff3ForHi(scheme, 1);


		// Проверка на эллиптическую
		/*d[i * 3] = b[i*3];

		d[i * 3 + 1] = b[i * 3 + 1];

		d[i * 3 + 2] = b[i * 3 + 2];*/
	}


}

vector<double> FEM::AssemblingGlobalF(InitialData* data, double time)
{
	vector<double> f;
	f.resize(data->knots.size() * 3);
	int nKEs = data->KEs.size();
	string file;
	for (int iKE = 0; iKE < nKEs; iKE++)
	{
		IKE* ke = data->KEs[iKE];
		int countKnotsInKe = ke->GetCountKnots();
		vector<double> fInKnotsX, fInKnotsY, fInKnotsZ;
		fInKnotsX.resize(ke->GetCountKnots());
		fInKnotsY.resize(ke->GetCountKnots());
		fInKnotsZ.resize(ke->GetCountKnots());
		for (int j = 0; j < countKnotsInKe; j++)
		{
			int globalJ = ke->globalNumsKnots[j];
			fInKnotsX[j] = GetFx(data->knots[globalJ], time);
			fInKnotsY[j] = GetFy(data->knots[globalJ], time);
			fInKnotsZ[j] = GetFz(data->knots[globalJ], time);
		}

		vector<double> localFx = ke->CalcLocalF(fInKnotsX);
		vector<double> localFy = ke->CalcLocalF(fInKnotsY);
		vector<double> localFz = ke->CalcLocalF(fInKnotsZ);
		file = "./locals/localFx" + std::to_string(iKE) + ".txt";
		WriteVector(localFx, file);
		file = "./locals/localFy" + std::to_string(iKE) + ".txt";
		WriteVector(localFy, file);
		file = "./locals/localFz" + std::to_string(iKE) + ".txt";
		WriteVector(localFz, file);
		for (int j = 0; j < countKnotsInKe; j++)
		{
			int globalJ = ke->globalNumsKnots[j] * 3;
			f[globalJ] += localFx[j];
			f[globalJ + 1] += localFy[j];
			f[globalJ + 2] += localFz[j];

		}
	}
	file = "./vectors/b.txt";
	WriteVector(f, file);
	return f;
}

void FEM::CalcFirstBoundaryConditions(InitialData* data, double time)
{
	for (int i = 0; i < data->bounds.size(); i++)
	{
		for (int iKnot = 0; iKnot < 4; iKnot++)
		{
			int global_num_coord = data->bounds[i].globalNum[iKnot];
			for (int r1 = 0; r1 < 3; r1++)
			{
				for (int r2 = 0; r2 < r1; r2++)
					A->d[global_num_coord][r2][r1] = A->d[global_num_coord][r1][r2] = 0;
				A->d[global_num_coord][r1][r1] = 1;
				
			}

			
			for (int j = A->ig[global_num_coord]; j < A->ig[global_num_coord + 1]; j++)
				for (int r1 = 0; r1 < 3; r1++)
				{
					for (int r2 = 0; r2 < 3; r2++)
						A->l[j][r1][r2] = 0;
				}
				
			for (int j = 0; j < A->ig[data->knots.size()]; j++)
				if (A->jg[j] == global_num_coord)
					for (int r1 = 0; r1 < 3; r1++)
					{
						for (int r2 = 0; r2 < 3; r2++)
							A->u[j][r1][r2] = 0;
					}

			d[global_num_coord*3] = u[global_num_coord*3];
			d[global_num_coord*3 + 1] = u[global_num_coord*3 + 1];
			d[global_num_coord*3 + 2] = u[global_num_coord*3 + 2];

			Symmetrization(global_num_coord);
			//d[global_num_coord] = 0.;     // Для решения


			// БОЛШИМ ЧИСЛОМ
			//A[global_num_coord][global_num_coord] = 1e+13;
			//d[global_num_coord] = u[global_num_coord] * 1e+13;
			//A[global_num_coord+1][global_num_coord+1] = 1e+13;
			//d[global_num_coord+1] = u[global_num_coord+1] * 1e+13;
			//A[global_num_coord+2][global_num_coord+2] = 1e+13;
			//d[global_num_coord+2] = u[global_num_coord+2] * 1e+13;
			//d[global_num_coord] = 0.* 10e+14;     // Для решения
		}
	}
}

void FEM::Symmetrization(int i)
{
	for (int j = A->ig[i]; j < A->ig[i + 1]; j++)
	{
		for (int r1 = 0; r1 < 3; r1++)
		{
			for (int r2 = 0; r2 < 3; r2++)
			{
				d[A->jg[j] * 3 + r1] -= d[i * 3+r1] * A->u[j][r1][r2];
				A->u[j][r1][r2] = 0;
			}
		}
	}


	int numKnots = A->d.size();
	for (int iIG = 0; iIG < numKnots; iIG++)
	{
		for (int j = A->ig[iIG]; j < A->ig[iIG+1]; j++)
			if (A->jg[j] == i)
				for (int r1 = 0; r1 < 3; r1++)
				{
					for (int r2 = 0; r2 < 3; r2++)
					{
						d[iIG * 3 + r1] -= d[i * 3] * A->l[j][r1][r2];
						A->l[j][r1][r2] = 0;
					}
				}
	}


	//for (int j = 0; j < i; j++)
	//{
	//	d[j] -= d[i] * A[j][i];
	//	A[j][i] = 0;
	//}

	//for (int j = i+1; j < A.size(); j++)
	//{
	//	d[j] -= d[i] * A[j][i];
	//	A[j][i] = 0;
	//}
}

void FEM::SolveSLAU(InitialData* data, TimeScheme* scheme)
{
	u.resize(slauSize, 0.);
	d.resize(slauSize, 0.);
	
	A->Clear();

	CalcU(data, scheme->time[scheme->time.size() - 1]);
	cout << "Сборка глобальной матрицы А.\n";
	CalcA(data, scheme);
	//A->WriteSparseMatrix("./matrix/A.txt");
	cout << "Сборка глобального вектора d.\n";
	CalcD(data, scheme);
	
	cout << "Учет первых краевых.\n";
	CalcFirstBoundaryConditions(data, scheme->time[scheme->time.size() - 1]);
	//string file = "./matrix/A" + std::to_string(scheme->time[3]) + ".txt";
	//A->WriteSparseMatrix(file);

	//file = "./vectors/d" + std::to_string(scheme->time[3]) + ".txt";
	//WriteVector(d, file);
	cout << "Решение СЛАУ\n";
	LOC();

	for (int i = 0; i < data->knots.size(); i++)
	{
		qx[i] = q[i * 3];
		qy[i] = q[i * 3 + 1];
		qz[i] = q[i * 3 + 2];
	}

}

void FEM::CalcU(InitialData* data, double time)
{
	int slauSize = u.size() / 3;
	for (int j = 0; j < slauSize; j++)
	{
		u[j * 3] = GetUx(data->knots[j], time, 0, data->knotToGo, data->uToGo);
		u[j * 3 + 1] = GetUy(data->knots[j], time, 0, data->knotToGo, data->uToGo);
		u[j * 3 + 2] = GetUz(data->knots[j], time, 0, data->knotToGo, data->uToGo);
	}
}

void FEM::WriteResultForSolution(vector<double> q, double time) //функция вывода в консоль
{
	ofstream out("./res/Result2.txt", ios_base::out | ios_base::app);

	out << endl << "ВРЕМЯ: " << time << endl;
	out << endl << "Результат в узлах (веса):" << endl;
	out << " ___________________________________________________________________ " << endl;

	out.setf(ios::left);
	out.width(15);
	out << "| № элемента " << "  | ";
	out.width(15);
	out << "x" << "| ";
	out.width(15);
	out << "y" << "| ";
	out.width(15);
	out << "z" << "| ";
	out.width(15);
	out << "u*" << "| ";
	//out.width(15);
	//out << "u" << "| ";
	//out.width(15);
	//out << "|u-u*|" << "|";
	out << endl;
	out << "|----------------|----------------|";
	out << "----------------|";
	out << "----------------|----------------|";
	out << endl;

	int slauSize = q.size();
	for (int i = 0; i < slauSize; i++)
	{
		out.setf(ios::left);
		out << "| ";
		out.width(15);
		out << i + 1 << "| ";
		out.width(15);
		out << knots[i].x << "| ";
		out.width(15);
		out << knots[i].y << "| ";
		out.width(15);
		out << knots[i].z << "| ";
		out.width(15);
		out << q[i] << "| ";
		//cout.width(15);
		//cout << u[i] << "| ";
		//cout.width(15);
		//cout << fabs(q[i] - u[i]) << "| ";
		out << endl;
	}

	out.close();
}

void FEM::WriteResultForTest(vector<double> q, double time) //функция вывода в консоль
{
	ofstream out("./res/ResultForTest.txt", ios_base::out | ios_base::app);

	out << endl << "ВРЕМЯ: " << time << endl;
	//out << "Результат в узлах (веса):" << endl;
	//out.setf(ios::left);
	//out.width(15);
	//out << "| № элемента " << "  | ";
	//out.width(15);
	//out << "x" << "| ";
	//out.width(15);
	//out << "y" << "| ";
	//out.width(15);
	//out << "z" << "| ";
	//out.width(15);
	//out << "ux*" << "| ";
	//out.width(15);
	//out << "uy*" << "| ";
	//out.width(15);
	//out << "uz*" << "| ";
	//out.width(15);
	//out << "|ux-ux*| / |ux|" << "|";
	//out.width(15);
	//out << "|uy-uy*| / |uy|" << "|";
	//out.width(15);
	//out << "|uz-uz*| / |uz|" << "|";
	out << endl;

	int slauSize = q.size();
	for (int i = 0; i < slauSize / 3; i++)
	{
		out.setf(ios::left);
		out << "| ";
		out.width(15);
		out << i + 1 << "| ";
		out.width(15);
		out << knots[i].x << "| ";
		out.width(15);
		out << knots[i].y << "| ";
		out.width(15);
		out << knots[i].z << "| ";
		out.width(15);
		out << q[i * 3] << "| ";
		out.width(15);
		out << q[i * 3 + 1] << "| ";
		out.width(15);
		out << q[i * 3 + 2] << "| ";
		//out.width(15);
		//out << fabs(q[i * 3] - u[i * 3]) / fabs(u[i * 3]) << "| ";
		//out.width(15);
		//out << fabs(q[i * 3 + 1] - u[i * 3 + 1]) / fabs(u[i * 3 + 1]) << "| ";
		//out.width(15);
		//out << fabs(q[i * 3 + 2] - u[i * 3 + 2]) / fabs(u[i * 3 + 2]) << "| ";
		out << endl;
	}
}

void FEM::CreateArea(InitialData* data)
{
	for (int i = 0; i < knots.size(); i++)
	{
		xArea.insert(knots[i].x);
		yArea.insert(knots[i].y);
		zArea.insert(knots[i].z);
	}

	double xbeg = *(xArea.begin());
	double ybeg = *(yArea.begin());
	double zbeg = *(zArea.begin());

	double hx = *(--xArea.end()) - xbeg;
	double hy = *(--yArea.end()) - ybeg;
	double hz = *(--zArea.end()) - zbeg;

	//xArea.clear();
	//yArea.clear();
	zArea.clear();

	int nStepsInArea = 50;

	for (int i = 1; i <= nStepsInArea; i++)
	{
		xArea.insert(xbeg + i * hx / nStepsInArea);
		yArea.insert(ybeg + i * hy / nStepsInArea);
		zArea.insert(2.5/*zbeg + i * hz / nStepsInArea*/);
	}


}

void FEM::FindIKeForArea(InitialData* data)
{
	KnotsIKEs.resize(xArea.size() * yArea.size() * zArea.size());

	int i = 0;
	for (set<double> ::iterator iz = zArea.begin(); iz != zArea.end(); iz++)
	{
		for (set<double> ::iterator ix = xArea.begin(); ix != xArea.end(); ix++)
		{
			for (set<double> ::iterator iy = yArea.begin(); iy != yArea.end(); iy++, i++)
			{
				Knot* knot = new Knot(*ix, *iy, *iz);
				KnotsIKEs[i] = FindIKe(data, knot);
			}
		}
	}
}

void FEM::SolveInAreaForTest(InitialData* data, double time)
{

	if (xArea.size() == 0 || yArea.size() == 0 || zArea.size() == 0)
	{
		CreateArea(data);
		FindIKeForArea(data);
	}
	

	string str = "./res/X/ResultAreaX" + std::to_string(time) + ".txt";
	string str1 = "./res/Y/ResultAreaY" + std::to_string(time) + ".txt";
	string str2 = "./res/Z/ResultAreaZ" + std::to_string(time) + ".txt";

	ofstream out(str);
	ofstream out1(str1);
	ofstream out2(str2);
	double result_x, result_y, result_z;
	double trueResult_x, trueResult_y, trueResult_z;
	double diffResult_x, diffResult_y, diffResult_z;

	int i = 0;
	Knot* knot = new Knot();
	set<double> ::iterator iz = zArea.begin();
	for (; iz != zArea.end(); iz++)
	{
		set<double> ::iterator ix = xArea.begin();
		for (; ix != xArea.end(); ix++)
		{
			set<double> ::iterator iy = yArea.begin();
			for (; iy != yArea.end(); iy++, i++)
			{
				knot->x = *ix;
				knot->y = *iy;
				knot->z = *iz;

				if (KnotsIKEs[i] >= 0)
				{
					result_x = data->KEs[KnotsIKEs[i]]->SolveInPoint(*knot, qx);
					result_y = data->KEs[KnotsIKEs[i]]->SolveInPoint(*knot, qy);
					result_z = data->KEs[KnotsIKEs[i]]->SolveInPoint(*knot, qz);



					out.setf(ios::left);
					out1.setf(ios::left);
					out2.setf(ios::left);

					out.width(15);
					out1.width(15);
					out2.width(15);
					out << knot->x << " ";
					out1 << knot->x << " ";
					out2 << knot->x << " ";

					out.width(15);
					out1.width(15);
					out2.width(15);
					out << knot->y << " ";
					out1 << knot->y << " ";
					out2 << knot->y << " ";

					out.precision(10);
					out.width(15);
					out << result_x << " ";

					out1.width(15);
					out1 << result_y << " ";

					out2.width(15);
					out2 << result_z << " ";


					trueResult_x = GetUx(*knot, time,0, data->knotToGo, data->uToGo);
					trueResult_y = GetUy(*knot, time,0, data->knotToGo, data->uToGo);
					trueResult_z = GetUz(*knot, time,0, data->knotToGo, data->uToGo);
					diffResult_x = abs(trueResult_x - result_x);
					diffResult_y = abs(trueResult_y - result_y);
					diffResult_z = abs(trueResult_z - result_z);

					out.width(15);
					out << trueResult_x << " ";

					out1.width(15);
					out1 << trueResult_y << " ";

					out2.width(15);
					out2 << trueResult_z << " ";



					out.width(15);
					out << diffResult_x << " ";

					out1.width(15);
					out1 << diffResult_y << " ";

					out2.width(15);
					out2 << diffResult_z << " ";

					out << endl;
					out1 << endl;
					out2 << endl;
				}
			}
			
		}

	}

	if (out.is_open())  out.close();
	if (out1.is_open())  out1.close();
	if (out2.is_open())  out2.close();
}


void FEM::SolveInArea(InitialData* data, double time) //функция вывода в консоль
{
	set<double> x, y, z;

	for (int i = 0; i < knots.size(); i++)
	{
		x.insert(knots[i].x);
		y.insert(knots[i].y);
		z.insert(knots[i].z);
	}

	double hx = *(--x.end()) - *(x.begin());
	double hy = *(--y.end()) - *(y.begin());
	double hz = *(--z.end()) - *(z.begin());
	double xbeg = *(x.begin());
	double ybeg = *(y.begin());
	double zbeg = *(z.begin());

	//x.clear();
	//y.clear();
	y.clear();

	for (int i = 1; i <= 25; i++)
	{
		z.insert(zbeg + i * hz / 25);
		x.insert(xbeg + i * hx / 25);
		y.insert(0.);
	}

	/*vector<Knot*> areaKnots;
	for (set<double> ::iterator iz = z.begin(); iz != z.end(); iz++ )
	{
		for (set<double> ::iterator ix = x.begin(); ix != x.end(); ix++ )
		{
			for (set<double> ::iterator iy = y.begin(); iy != y.end(); iy++)
			{
				Knot* knot = new Knot(*ix, *iy, *iz);
				areaKnots.push_back(knot);
			}
		}
	}*/

	string str = "ResultArea" + std::to_string(time) + ".txt";

	ofstream out(str);

	/*out << endl << "ВРЕМЯ: " << time << endl;
	out << endl << "Результат в узлах (веса):" << endl;
	out << " _____________________________________________________________________________________ " << endl;

	out.setf(ios::left);
	out.width(15);
	out << "| № элемента " << "  | ";
	out.width(15);
	out << "x" << "| ";
	out.width(15);
	out << "y" << "| ";
	out.width(15);
	out << "z" << "| ";
	out.width(15);
	out << "u*" << "| ";
	out.width(15);
	out << endl;
	out << "|----------------|";
	out << "----------------|----------------|----------------|";
	out << "----------------|";
	out << endl;*/
	/*out.close();*/

	//for (int i = 0; i < areaKnots.size(); i++)
	//{
	int i = 0;
	for (set<double> ::iterator iy = y.begin(); iy != y.end(); iy++)
	{
		for (set<double> ::iterator ix = x.begin(); ix != x.end(); ix++)
		{
			for (set<double> ::iterator iz = z.begin(); iz != z.end(); iz++)
			{
				Knot* knot = new Knot(*ix, *iy, *iz);
				//int iKe = FindIKe(data, areaKnots[i]);
				int iKe = FindIKe(data, knot);


				/*out.open(str, ios_base::out | ios_base::app);*/
				out.setf(ios::left);
				//out << "| ";
				//out.width(15);
				//out << i + 1 << "| ";
				out.width(15);
				out << knot->x << " ";
				out.width(15);
				out << knot->z << " ";
				//out.width(15);
				//out << knot->z << "| ";
				//out.width(15);

				double result;
				if (iKe < 0) result = 0.;
				else  result = data->KEs[iKe]->SolveInPoint(*knot, q);
				//if (fabs(result) < 9e-5) result = 0;
				out << result;

				//if (areaKnots[i]->x == 1) out << " " << data->KEs[16]->SolveInPoint(*areaKnots[i], q);
				//out << "| ";
				out << endl;
				delete knot;
				/*out.close()*/;
			}
		}

	}

	if (out.is_open())out.close();
}

int FEM::FindIKe(InitialData* data, Knot* knot)
{
	int iKe = 0;
	for (; iKe < data->KEs.size(); iKe++)
	{
		if (data->KEs[iKe]->IsIn(*knot))
		{
			return iKe;
		}
	}

	cout << "Точка - " << knot->x << " " << knot->y << " " << knot->z << " - находится вне расчетной области\n";

	return -1;
}

double FEM::SolveInPoint(InitialData* data, Knot knot)
{
	int iKe = 0;
	for (; iKe < data->KEs.size(); iKe++)
	{
		if (data->KEs[iKe]->IsIn(knot))
		{
			return data->KEs[iKe]->SolveInPoint(knot, q);
		}
	}

	cout << "Точка - " << knot.x << " " << knot.y << " " << knot.z << " - находится вне расчетной области";
	return 0;
}




