#include "SLAU.h"
#include <iostream>
#include <fstream>
#include "InputFuncs.h"
#include <set>

void CreateMatrix(vector<vector<double>> &matrix, int size)
{
	matrix.resize(size);
	for (int i = 0; i < size; i++)
		matrix[i].resize(size);
};

SLAU::SLAU(InitialData* data)
{
	cout << "Создаем СЛАУ:\n";
	slauSize = data->knots.size() * 3;
	knots.resize(data->knots.size());
	for (int i = 0; i < data->knots.size(); i++)
		knots[i] = data->knots[i];


	A.resize(slauSize);
	for (int i = 0; i < slauSize; i++)
		A[i].resize(slauSize);
	q.resize(slauSize);
	qx.resize(data->knots.size());
	qy.resize(data->knots.size());
	qz.resize(data->knots.size());
	u.resize(slauSize);
	d.resize(slauSize);

	CreateMatrix(M, knots.size());
	if (ReadMatrix(M, "Mhi.txt"))
	{
		cout << "\tПодсчет матрицы Mhi:\n";
		M = AssemblingGlobalM(data);
		WriteMatrix(M, "Mhi.txt");
	}

	CreateMatrix(G_xx, knots.size());
	if (ReadMatrix(G_xx, "G_xx.txt"))
	{
		cout << "\tПодсчет матрицы G_xx:\n";
		G_xx = AssemblingGlobalG_aa(data, x, x, 0);
		WriteMatrix(G_xx, "G_xx.txt");
	}

	CreateMatrix(G_yy, knots.size());
	if (ReadMatrix(G_yy, "G_yy.txt"))
	{
		cout << "\tПодсчет матрицы G_yy:\n";
		G_yy = AssemblingGlobalG_aa(data, y, y, 0);
		WriteMatrix(G_yy, "G_yy.txt");
	}

	CreateMatrix(G_zz, knots.size());
	if (ReadMatrix(G_zz, "G_zz.txt"))
	{
		cout << "\tПодсчет матрицы G_zz:\n";
		G_zz = AssemblingGlobalG_aa(data, z, z, 0);
		WriteMatrix(G_zz, "G_zz.txt");
	}
	
	CreateMatrix(G_xx_sigma, knots.size());
	if (ReadMatrix(G_xx_sigma, "G_xx_sigma.txt"))
	{

		cout << "\tПодсчет матрицы G_xx_sigma:\n";
		G_xx_sigma = AssemblingGlobalG_aa(data, x, x, 1);
		WriteMatrix(G_xx_sigma, "G_xx_sigma.txt");
	}

	CreateMatrix(G_yy_sigma, knots.size());
	if (ReadMatrix(G_yy_sigma, "G_yy_sigma.txt"))
	{
		cout << "\tПодсчет матрицы G_yy_sigma:\n";
		G_yy_sigma = AssemblingGlobalG_aa(data, y, y, 1);
		WriteMatrix(G_yy_sigma, "G_yy_sigma.txt");
	}


	CreateMatrix(G_zz_sigma, knots.size());
	if (ReadMatrix(G_zz_sigma, "G_zz_sigma.txt"))
	{
		cout << "\tПодсчет матрицы G_zz_sigma:\n";
		G_zz_sigma = AssemblingGlobalG_aa(data, z, z, 1);
		WriteMatrix(G_zz_sigma, "G_zz_sigma.txt");
	}
	

	CreateMatrix(G_xy, knots.size());
	if (ReadMatrix(G_xy, "G_xy.txt"))
	{
		cout << "\tПодсчет матрицы G_xy:\n";
		G_xy = AssemblingGlobalG_aa(data, x, y, 0);
		WriteMatrix(G_xy, "G_xy.txt");
	}

	CreateMatrix(G_xz, knots.size());
	if (ReadMatrix(G_xz, "G_xz.txt"))
	{
		cout << "\tПодсчет матрицы G_xz:\n";
		G_xz = AssemblingGlobalG_aa(data, x, z, 0);
		WriteMatrix(G_xz, "G_xz.txt");
	}

	CreateMatrix(G_yz, knots.size());
	if (ReadMatrix(G_yz, "G_yz.txt"))
	{
		cout << "\tПодсчет матрицы G_yz:\n";
		G_yz = AssemblingGlobalG_aa(data, y, z, 0);
		WriteMatrix(G_yz, "G_yz.txt");
	}


	/*WriteMatrix(M);*/
	/*G = AssemblingGlobalG(data);*/
   /* WriteMatrix(G);*/
}

vector<vector<double>>  SLAU::AssemblingGlobalG_aa(InitialData* data, axis a1, axis a2, int iCoeff)
{
	vector<vector<double>> globalMatrix;
	int sizeGlobalMatrix = data->knots.size();
	globalMatrix.resize(sizeGlobalMatrix);
	for (int i = 0; i < sizeGlobalMatrix; i++)
		globalMatrix[i].resize(sizeGlobalMatrix);

	int nKEs = data->KEs.size();
	for (int iKE = 0; iKE < nKEs; iKE++)
	{
		IKE* ke = data->KEs[iKE];
		vector<vector<double>> localMatrix = ke->CalcLocalG();
		double sigma = data->coeffs[ke->iCoeff].Yung / (2. * (1. + data->coeffs[ke->iCoeff].Poisson));
		double lambda = data->coeffs[ke->iCoeff].Poisson * data->coeffs[ke->iCoeff].Yung / (1. + data->coeffs[ke->iCoeff].Poisson) / (1. - 2. * data->coeffs[ke->iCoeff].Poisson);

		double coeff = iCoeff == 0 ? lambda + sigma : sigma;

		int countKnotsInKe = ke->GetCountKnots();

		for (int j = 0; j < countKnotsInKe; j++)
		{
			int globalJ = ke->globalNumsKnots[j];
			for (int k = 0; k < countKnotsInKe; k++)
			{
				int globalK = ke->globalNumsKnots[k];
				globalMatrix[globalJ][globalK] += coeff * localMatrix[j][k];
			}
		}
	}

	return globalMatrix;
}


void SLAU::LOC()
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
	f = MultMatrByVect(A, q);

	for (i = 0; i < slauSize; i++)
		z[i] = r[i] = d[i] - f[i];

	p = MultMatrByVect(A, z);
	nvzk = sqrt(CalcScalar(r, r)) / sqrt(CalcScalar(d, d));
	unsigned int k;
	for (k = 1; k < 100000 && nvzk > eps; k++)
	{
		lastnvzk = nvzk;
		skp = CalcScalar(p, p);
		alfa = CalcScalar(p, r) / skp;

		for (i = 0; i < slauSize; i++)
		{
			q[i] += alfa * z[i];
			r[i] -= alfa * p[i];
		}

		f = MultMatrByVect(A, r);
		beta = -CalcScalar(p, f) / skp;

		for (i = 0; i < slauSize; i++)
		{
			z[i] = r[i] + beta * z[i];
			p[i] = f[i] + beta * p[i];
		}

		nvzk = sqrt(CalcScalar(r, r)) / sqrt(CalcScalar(d, d));
		cout << "iteration: " << k << "\tНевязка:" << nvzk << endl;
	}

	cout << "Вычисление СЛАУ закончено.\n\n\n";

}

/// <summary>
/// Вычисление глобальных матриц массы пространственной сетки
/// </summary>
/// <param name="data">Входные данные</param>
vector<vector<double>> SLAU::AssemblingGlobalM(InitialData* data)
{
	vector<vector<double>> globalMatrix;
	int sizeGlobalMatrix = data->knots.size();
	globalMatrix.resize(sizeGlobalMatrix);
	for (int i = 0; i < sizeGlobalMatrix; i++)
		globalMatrix[i].resize(sizeGlobalMatrix);

	int nKEs = data->KEs.size();

	for (int iKE = 0; iKE < nKEs; iKE++)
	{
		IKE* ke = data->KEs[iKE];

		double hi = data->coeffs[ke->iCoeff].density;

		vector<vector<double>> localMatrix = ke->CalcLocalM();
		int countKnotsInKe = ke->GetCountKnots();

		for (int j = 0; j < countKnotsInKe; j++)
		{
			int globalJ = ke->globalNumsKnots[j];
			for (int k = 0; k < countKnotsInKe; k++)
			{
				int globalK = ke->globalNumsKnots[k];
				globalMatrix[globalJ][globalK] += hi * localMatrix[j][k];
			}
		}
	}

	return globalMatrix;
}

/// <summary>
/// Вычисление глобальных матриц жесткости пространственной сетки
/// </summary>
/// <param name="data">Входные данные</param>
vector<vector<double>> SLAU::AssemblingGlobalG(InitialData* data)
{
	vector<vector<double>> globalMatrix;
	int sizeGlobalMatrix = data->knots.size();
	globalMatrix.resize(sizeGlobalMatrix);
	for (int i = 0; i < sizeGlobalMatrix; i++)
		globalMatrix[i].resize(sizeGlobalMatrix);

	int nKEs = data->KEs.size();
	for (int iKE = 0; iKE < nKEs; iKE++)
	{
		IKE* ke = data->KEs[iKE];
		vector<vector<double>> localMatrix = ke->CalcLocalG();

		int countKnotsInKe = ke->GetCountKnots();

		for (int j = 0; j < countKnotsInKe; j++)
		{
			int globalJ = ke->globalNumsKnots[j];
			for (int k = 0; k < countKnotsInKe; k++)
			{
				int globalK = ke->globalNumsKnots[k];
				globalMatrix[globalJ][globalK] += localMatrix[j][k];
			}
		}
	}

	return globalMatrix;
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

void SLAU::CalcA(InitialData* data, TimeScheme* scheme) // Вычисление глобальной A
{
	vector<double> timeToCalc = scheme->time;

	for (int i = 0; i < data->knots.size(); i++)
	{
		for (int j = 0; j < data->knots.size(); j++)
		{
			A[i * 3][j * 3] = G_xx[i][j] + G_xx_sigma[i][j] + G_yy_sigma[i][j] + G_zz_sigma[i][j]
				+ M[i][j] * CalcTimeCoeff4ForHi(scheme, 3);

			A[i * 3 + 1][j * 3 + 1] = G_yy[i][j] + G_xx_sigma[i][j] + G_yy_sigma[i][j] + G_zz_sigma[i][j]
				+ M[i][j] * CalcTimeCoeff4ForHi(scheme, 3);

			A[i * 3 + 2][j * 3 + 2] = G_zz[i][j] + G_xx_sigma[i][j] + G_yy_sigma[i][j] + G_zz_sigma[i][j]
				+ M[i][j] * CalcTimeCoeff4ForHi(scheme, 3);

			//A[i * 3][j * 3] = G_xx[i][j] + G_xx_sigma[i][j] + G_yy_sigma[i][j] + G_zz_sigma[i][j]
			//	+ M[i][j] * CalcTimeCoeff3ForHi(scheme, 3);

			//A[i * 3 + 1][j * 3 + 1] = G_yy[i][j] + G_xx_sigma[i][j] + G_yy_sigma[i][j] + G_zz_sigma[i][j]
			//	+ M[i][j] * CalcTimeCoeff3ForHi(scheme, 3);

			//A[i * 3 + 2][j * 3 + 2] = G_zz[i][j] + G_xx_sigma[i][j] + G_yy_sigma[i][j] + G_zz_sigma[i][j]
			//	+ M[i][j] * CalcTimeCoeff4ForHi(scheme, 3); M[i][j] * CalcTimeCoeff3ForHi(scheme, 3);

			A[i * 3][j * 3 + 1] = A[i * 3 + 1][j * 3] = G_xy[i][j];

			A[i * 3][j * 3 + 2] = A[i * 3 + 2][j * 3] = G_xz[i][j];

			A[i * 3 + 1][j * 3 + 2] = A[i * 3 + 2][j * 3 + 1] = G_yz[i][j];
		}
	}

	//cout << endl;
	//for (int i = 0; i < A.size(); i++)
	//    cout << A[13][i] << "\t";
	//cout << endl;
}





void SLAU::CalcD(InitialData* data, TimeScheme* scheme) // Вычисление глобальной A
{
	vector<double> timeToCalc = scheme->time;


	vector<double> b = AssemblingGlobalF(data, timeToCalc[3]);

	vector<double>Mq_time0_x, Mq_time1_x, Mq_time2_x, Mq_time0_y, Mq_time1_y, Mq_time2_y, Mq_time0_z, Mq_time1_z, Mq_time2_z;

	Mq_time0_x = MultMatrByVect(M, scheme->qx[0]);
	Mq_time1_x = MultMatrByVect(M, scheme->qx[1]);
	Mq_time2_x = MultMatrByVect(M, scheme->qx[2]);

	Mq_time0_y = MultMatrByVect(M, scheme->qy[0]);
	Mq_time1_y = MultMatrByVect(M, scheme->qy[1]);
	Mq_time2_y = MultMatrByVect(M, scheme->qy[2]);

	Mq_time0_z = MultMatrByVect(M, scheme->qz[0]);
	Mq_time1_z = MultMatrByVect(M, scheme->qz[1]);
	Mq_time2_z = MultMatrByVect(M, scheme->qz[2]);

	for (int i = 0; i < data->knots.size(); i++)
	{

		/*d[i] = b[i];*/
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
	}


}

vector<double> SLAU::AssemblingGlobalF(InitialData* data, double time)
{
	vector<double> f;
	f.resize(data->knots.size() * 3);
	int nKEs = data->KEs.size();
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

		for (int j = 0; j < countKnotsInKe; j++)
		{
			int globalJ = ke->globalNumsKnots[j] * 3;
			f[globalJ] += localFx[j];
			f[globalJ + 1] += localFy[j];
			f[globalJ + 2] += localFz[j];
		}
	}

	return f;
}

void SLAU::CalcFirstBoundaryConditions(InitialData* data, double time)
{
	for (int i = 0; i < data->bounds.size(); i++)
	{
		for (int j = 0; j < 4; j++)
		{
			int global_num_coord = data->bounds[i].globalNum[j] * 3;
			for (int k = 0; k < data->knots.size(); k++)
			{
				for (int r1 = 0; r1 < 3; r1++)
					for (int r2 = 0; r2 < 3; r2++)
						A[global_num_coord + r1][k * 3 + r2] = 0.;
			}

			A[global_num_coord][global_num_coord] = 1.;
			A[global_num_coord + 1][global_num_coord + 1] = 1.;
			A[global_num_coord + 2][global_num_coord + 2] = 1.;
			d[global_num_coord] = u[global_num_coord];
			d[global_num_coord + 1] = u[global_num_coord + 1];
			d[global_num_coord + 2] = u[global_num_coord + 2];

			//d[global_num_coord] = 0.;     // Для решения


			// БОЛШИМ ЧИСЛОМ
			//A[global_num_coord][global_num_coord] = 10e+14;
			//d[global_num_coord] = u[global_num_coord] * 10e+14;
			//A[global_num_coord+1][global_num_coord+1] = 10e+14;
			//d[global_num_coord+1] = u[global_num_coord+1] * 10e+14;
			//A[global_num_coord+2][global_num_coord+2] = 10e+14;
			//d[global_num_coord+2] = u[global_num_coord+2] * 10e+14;
			//d[global_num_coord] = 0.* 10e+14;     // Для решения
		}
	}
}

void SLAU::SolveSLAU(InitialData* data, TimeScheme* scheme)
{
	u.resize(slauSize, 0.);
	d.resize(slauSize, 0.);
	for (int i = 0; i < slauSize; i++) A[i].resize(slauSize, 0.);


	CalcU(data, scheme->time[scheme->time.size() - 1]);
	cout << "Сборка глобальной матрицы А.\n";
	CalcA(data, scheme);
	cout << "Сборка глобального вектора d.\n";
	CalcD(data, scheme);

	cout << "Учет первых краевых.\n";
	CalcFirstBoundaryConditions(data, scheme->time[scheme->time.size() - 1]);
	//WriteMatrix(A);

	cout << "Решение СЛАУ\n";
	LOC();

	for (int i = 0; i < data->knots.size(); i++)
	{
		qx[i] = q[i * 3];
		qy[i] = q[i * 3 + 1];
		qz[i] = q[i * 3 + 2];
	}

}

void SLAU::CalcU(InitialData* data, double time)
{
	int slauSize = u.size() / 3;
	for (int j = 0; j < slauSize; j++)
	{
		u[j * 3] = GetUx(data->knots[j], time);
		u[j * 3 + 1] = GetUy(data->knots[j], time);
		u[j * 3 + 2] = GetUz(data->knots[j], time);
	}
}

void SLAU::WriteResultForSolution(vector<double> q, double time) //функция вывода в консоль
{
	ofstream out("Result.txt", ios_base::out | ios_base::app);

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

void SLAU::WriteResultForTest(vector<double> q, double time) //функция вывода в консоль
{
	ofstream out("ResultForTest.txt", ios_base::out | ios_base::app);

	out << endl << "ВРЕМЯ: " << time << endl;
	out << "Результат в узлах (веса):" << endl;
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
	out << "ux*" << "| ";
	out.width(15);
	out << "uy*" << "| ";
	out.width(15);
	out << "uz*" << "| ";
	//out.width(15);
	//out << "u" << "| ";
	out.width(15);
	out << "|ux-ux*|" << "|";
	out.width(15);
	out << "|uy-uy*|" << "|";
	out.width(15);
	out << "|uz-uz*|" << "|";
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
		out.width(15);
		out << fabs(q[i * 3] - u[i * 3]) << "| ";
		out.width(15);
		out << fabs(q[i * 3 + 1] - u[i * 3 + 1]) << "| ";
		out.width(15);
		out << fabs(q[i * 3 + 2] - u[i * 3 + 2]) << "| ";
		out << endl;
	}
}

void SLAU::SolveInAreaForTest(InitialData* data, double time)
{
	set<double> x, y, z;

	for (int i = 0; i < knots.size(); i++)
	{
		x.insert(knots[i].x);
		y.insert(knots[i].y);
		z.insert(knots[i].z);
	}

	double xbeg = *(x.begin());
	double ybeg = *(y.begin());
	double zbeg = *(z.begin());

	double hx = *(--x.end()) - xbeg;
	double hy = *(--y.end()) - ybeg;
	double hz = *(--z.end()) - zbeg;

	//x.clear();
	//y.clear();
	z.clear();

	int nStepsInArea = 25;

	for (int i = 1; i <= nStepsInArea; i++)
	{
		x.insert(xbeg + i * hx / nStepsInArea);
		y.insert(ybeg + i * hy / nStepsInArea);
		z.insert(2.5/*zbeg + i * hz / nStepsInArea*/);
	}

	string str = "ResultAreaX" + std::to_string(time) + ".txt";
	string str1 = "ResultAreaY" + std::to_string(time) + ".txt";
	string str2 = "ResultAreaZ" + std::to_string(time) + ".txt";

	ofstream out(str);
	ofstream out1(str1);
	ofstream out2(str2);
	double result_x, result_y, result_z;
	double trueResult_x, trueResult_y, trueResult_z;
	double diffResult_x, diffResult_y, diffResult_z;
	int i = 0;
	for (set<double> ::iterator iz = z.begin(); iz != z.end(); iz++)
	{
		for (set<double> ::iterator ix = x.begin(); ix != x.end(); ix++)
		{
			for (set<double> ::iterator iy = y.begin(); iy != y.end(); iy++, i++)
			{
				Knot* knot = new Knot(*ix, *iy, *iz);
				int iKe = FindIKe(data, knot);

				if (iKe >= 0)
				{
					result_x = data->KEs[iKe]->SolveInPoint(*knot, qx);
					result_y = data->KEs[iKe]->SolveInPoint(*knot, qy);
					result_z = data->KEs[iKe]->SolveInPoint(*knot, qz);

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

					out.width(15);

					out << result_x << " ";

					out1.width(15);

					out1 << result_y << " ";

					out2.width(15);

					out2 << result_z << " ";

					out << endl;
					out1 << endl;
					out2 << endl;
				}

				delete knot;
			}
		}

	}

	if (out.is_open())  out.close();
}


void SLAU::SolveInArea(InitialData* data, double time) //функция вывода в консоль
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

int SLAU::FindIKe(InitialData* data, Knot* knot)
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

double SLAU::SolveInPoint(InitialData* data, Knot knot)
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




