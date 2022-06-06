#include "SLAU.h"

int Block_3_SM::ReadSparseMatrix(string pathFile)
{
	int num_of_knots = d.size();
	double** mat = new double* [num_of_knots * 3]{};
	for (int i = 0; i < num_of_knots * 3; i++)
	{
		mat[i] = new double[num_of_knots * 3]{};
	}

	std::ifstream in(pathFile);
	if (!in.is_open()) return 1;
	for (int i = 0; i < num_of_knots * 3; i++)
	{
		for (int j = 0; j < num_of_knots * 3; j++)
		{
			in >> mat[i][j];
		}
	}

	for (int i = 0; i < num_of_knots; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3;k++)
			{
				d[i][j][k] = mat[i * 3 + j][i * 3 + k];
			}
		}

		for (int j = ig[i]; j < ig[i + 1]; j++)
		{
			for (int r = 0; r < 3; r++)
			{
				for (int k = 0; k < 3;k++)
				{
					l[j][r][k] = mat[i * 3 + r][jg[j] * 3 + k];
					u[j][r][k] = mat[jg[j] * 3 + r][i * 3 + k];
				}
			}
		}
	}

	return 0;
}

void Block_3_SM::WriteSparseMatrix(string pathFile)
{
	int num_of_knots = d.size();
	double** mat = new double* [num_of_knots*3] {};
	for (int i = 0; i < num_of_knots*3; i++)
	{
		mat[i] = new double[num_of_knots*3] {};
	}

	for (int i = 0; i < num_of_knots; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3;k++)
			{
				mat[i * 3 + j][i * 3+k] = d[i][j][k];
			}
		}

		for (int j = ig[i]; j < ig[i + 1]; j++)
		{
			for (int r = 0; r < 3; r++)
			{
				for (int k = 0; k < 3;k++)
				{
					mat[i * 3 + r][jg[j] * 3+ k] = l[j][r][k];
					mat[jg[j]*3 + r][i*3+k] = u[j][r][k];
				}
			}
		}
	}

	std::ofstream out(pathFile);
	out.precision(20);
	out.setf(std::ios::fixed);
	for (int i = 0; i < num_of_knots*3; i++)
	{
		for (int j = 0; j < num_of_knots*3; j++)
		{
			out << mat[i][j] << " ";
		}
		out << "\n";
	}
}

void Block_3_SM::AddToBlock(double** to, double into[3][3])
{
	for (int iBlock = 0; iBlock < 3; iBlock++)
	{
		for (int jBlock = 0; jBlock < 3; jBlock++)
		{
			to[iBlock][jBlock] += into[iBlock][jBlock];
		}
	}
}

void Block_3_SM::AddToBlock(double** to, double** into)
{
	for (int iBlock = 0; iBlock < 3; iBlock++)
	{
		for (int jBlock = 0; jBlock < 3; jBlock++)
		{
			to[iBlock][jBlock] += into[iBlock][jBlock];
		}
	}
}

void Block_3_SM::AddElement(int i, int j, double elem[3][3])
{
	bool found = false;
	if (i == j)
	{
		AddToBlock(d[i], elem);
	}
		
	else if (i < j)
	{
		int m;
		for (m = ig[j]; m < ig[j + 1]; m++)
			if (jg[m] == i) { found = true; break; }
		if (found) AddToBlock(u[m], elem);
	}
	else
	{
		int n;
		for (n = ig[i]; n < ig[i + 1]; n++)
			if (jg[n] == j) { found = true; break; }
		if (found) AddToBlock(l[n], elem);
	}

}

void Block_3_SM::AddElement(int i, int j, double** elem)
{
	bool found = false;
	if (i == j)
	{
		AddToBlock(d[i], elem);
	}

	else if (i < j)
	{
		int m;
		for (m = ig[j]; m < ig[j + 1]; m++)
			if (jg[m] == i) { found = true; break; }
		if (found) AddToBlock(u[m], elem);
	}
	else
	{
		int n;
		for (n = ig[i]; n < ig[i + 1]; n++)
			if (jg[n] == j) { found = true; break; }
		if (found) AddToBlock(l[n], elem);
	}

}

vector<double> Block_3_SM::MultMatrByVect(vector<double> b)
{
	vector<double> v;
	int num_of_knots = d.size();
	v.resize(num_of_knots*3, 0);

	for (int i = 0; i < num_of_knots; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				v[i * 3 + j] += d[i][j][k] * b[i * 3 + k];
			}
		}
	}
		

	for (int i = 0; i < num_of_knots; i++)
		for (int j = ig[i]; j < ig[i + 1]; j++)
		{
			for (int h = 0; h < 3; h++)
			{
				for (int k = 0; k < 3; k++)
				{
					v[i * 3 + h] += l[j][h][k] * b[jg[j] * 3 + k];
					v[jg[j]*3+h] += u[j][h][k] * b[i*3 + k];
				}
			}
		}

	return v;
}

Block_3_SM::Block_3_SM(InitialData* data)
{
	int* list1, * list2;
	int num_of_knots = data->knots.size();
	int* listbeg = new int[num_of_knots];

	for (int i = 0; i < num_of_knots; i++)
		listbeg[i] = -1;

	list1 = new int[num_of_knots * num_of_knots]{};
	list2 = new int[num_of_knots * num_of_knots]{};
	int listsize = -1, iaddr, ind1, ind2, k;
	int num_of_KE = data->KEs.size();
	for (int iel = 0; iel < num_of_KE; iel++) // ï
	{
		int nKnotsInKe = data->KEs[iel]->GetCountKnots();
		for (int i = 0; i < nKnotsInKe; i++) // 
		{
			k = data->KEs[iel]->globalNumsKnots[i]; //
			for (int j = i + 1; j < nKnotsInKe; j++) // need to set N = ?
			{
				ind1 = k;
				ind2 = data->KEs[iel]->globalNumsKnots[j];  //
				if (ind2 < ind1) //
				{
					ind1 = ind2;
					ind2 = k;
				}
				iaddr = listbeg[ind2];
				if (iaddr == -1) // 
				{
					listsize++;
					listbeg[ind2] = listsize;
					list1[listsize] = ind1;
					list2[listsize] = -1;
				}
				else // 
				{
					while (list1[iaddr] < ind1 && list2[iaddr] >= 0)
						iaddr = list2[iaddr];
					if (list1[iaddr] > ind1)  // 
					{                         // 
						listsize++;
						list1[listsize] = list1[iaddr];
						list2[listsize] = list2[iaddr];
						list1[iaddr] = ind1;
						list2[iaddr] = listsize;
					}
					else if (list1[iaddr] < ind1) // 
					{
						listsize++;
						list2[iaddr] = listsize;
						list1[listsize] = ind1;
						list2[listsize] = -1;
					}
				}
			}
		}
	}

	ig.resize(num_of_knots + 1, 0);
	jg.resize(listsize + 1, 0);  // +1???

	for (int i = 0; i < num_of_knots; i++)
	{
		ig[i + 1] = ig[i];

		for (iaddr = listbeg[i]; iaddr != -1; )
		{
			jg[ig[i + 1]] = list1[iaddr];
			ig[i + 1]++;
			iaddr = list2[iaddr];
		}
	}

	delete[] listbeg;
	delete[] list1;
	delete[] list2;

	l.resize(listsize + 1);
	u.resize(listsize + 1);
	d.resize(num_of_knots);

	for (int i = 0; i < l.size(); i++)
	{
		l[i] = new double* [3];
		u[i] = new double* [3];

		for (int j = 0; j < 3; j++)
		{
			l[i][j] = new double[3]{};
			u[i][j] = new double[3]{};
		}
	}

	for (int i = 0; i < d.size(); i++)
	{
		d[i] = new double* [3];
		
		for (int j = 0; j < 3; j++)
		{
			d[i][j] = new double[3]{};
		}
	}
}

void Block_3_SM::Clear()
{
	for (int i = 0; i < l.size(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0;k < 3;k++)
			{
				l[i][j][k] = 0;
				u[i][j][k] = 0;
			}
		}
	}

	for (int i = 0; i < d.size(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0;k < 3;k++)
			{
				d[i][j][k] = 0;
			}
		}
	}

}




NoBlockSM::NoBlockSM(InitialData* data)
{
	int* list1, * list2;
	int num_of_knots = data->knots.size();
	int* listbeg = new int[num_of_knots];

	for (int i = 0; i < num_of_knots; i++)
		listbeg[i] = -1;

	list1 = new int[num_of_knots * num_of_knots]{};
	list2 = new int[num_of_knots * num_of_knots]{};
	int listsize = -1, iaddr, ind1, ind2, k;
	int num_of_KE = data->KEs.size();
	for (int iel = 0; iel < num_of_KE; iel++) // ï
	{
		int nKnotsInKe = data->KEs[iel]->GetCountKnots();
		for (int i = 0; i < nKnotsInKe; i++) // 
		{
			k = data->KEs[iel]->globalNumsKnots[i]; //
			for (int j = i + 1; j < nKnotsInKe; j++) // need to set N = ?
			{
				ind1 = k;
				ind2 = data->KEs[iel]->globalNumsKnots[j];  //
				if (ind2 < ind1) //
				{
					ind1 = ind2;
					ind2 = k;
				}
				iaddr = listbeg[ind2];
				if (iaddr == -1) // 
				{
					listsize++;
					listbeg[ind2] = listsize;
					list1[listsize] = ind1;
					list2[listsize] = -1;
				}
				else // 
				{
					while (list1[iaddr] < ind1 && list2[iaddr] >= 0)
						iaddr = list2[iaddr];
					if (list1[iaddr] > ind1)  // 
					{                         // 
						listsize++;
						list1[listsize] = list1[iaddr];
						list2[listsize] = list2[iaddr];
						list1[iaddr] = ind1;
						list2[iaddr] = listsize;
					}
					else if (list1[iaddr] < ind1) // 
					{
						listsize++;
						list2[iaddr] = listsize;
						list1[listsize] = ind1;
						list2[listsize] = -1;
					}
				}
			}
		}
	}

	ig.resize(num_of_knots + 1, 0);
	jg.resize(listsize + 1, 0);  // +1???

	for (int i = 0; i < num_of_knots; i++)
	{
		ig[i + 1] = ig[i];

		for (iaddr = listbeg[i]; iaddr != -1; )
		{
			jg[ig[i + 1]] = list1[iaddr];
			ig[i + 1]++;
			iaddr = list2[iaddr];
		}
	}

	delete[] listbeg;
	delete[] list1;
	delete[] list2;

	l.resize(listsize + 1,0);
	u.resize(listsize + 1,0);
	d.resize(num_of_knots,0);

}

int NoBlockSM::ReadSparseMatrix(string pathFile)
{
	int num_of_knots = d.size();
	double** mat = new double* [num_of_knots] {};
	for (int i = 0; i < num_of_knots; i++)
	{
		mat[i] = new double[num_of_knots] {};
	}

	std::ifstream in(pathFile);
	if (!in.is_open()) return 1;

	for (int i = 0; i < num_of_knots; i++)
	{
		for (int j = 0; j < num_of_knots; j++)
			in >> mat[i][j];
	}

	in.close();

	for (int i = 0; i < num_of_knots; i++)
	{
		d[i] = mat[i][i];
		for (int j =ig[i]; j < ig[i + 1]; j++)
		{
			l[j] = mat[i][jg[j]];
			u[j] = mat[jg[j]][i];
		}
	}


	return 0;
}

void NoBlockSM::WriteSparseMatrix(string pathFile)
{
	int num_of_knots = d.size();
	double** mat = new double* [num_of_knots] {};
	for (int i = 0; i < num_of_knots; i++)
	{
		mat[i] = new double[num_of_knots] {};
	}

	for (int i = 0; i < num_of_knots; i++)
	{
		mat[i][i] = d[i];
		for (int j = ig[i]; j < ig[i + 1]; j++)
		{
			mat[i][jg[j]] = l[j];
			mat[jg[j]][i] = u[j];
		}
	}

	std::ofstream out(pathFile);
	out.precision(20);
	out.setf(std::ios::fixed);
	for (int i = 0; i < num_of_knots; i++)
	{
		for (int j = 0; j < num_of_knots; j++)
		{
			out << mat[i][j] << " ";
		}
		out << "\n";
	}
}

void NoBlockSM::AddElement(int i, int j, double elem)
{
	bool found = false;
	if (i == j)
		d[i] += elem;
	else if (i < j)
	{
		int m;
		for (m = ig[j]; m < ig[j + 1]; m++)
			if (jg[m] == i) { found = true; break; }
		if (found)
			u[m] += elem; // i-1?
	}
	else
	{
		int n;
		for (n = ig[i]; n < ig[i + 1]; n++)
			if (jg[n] == j) { found = true; break; }
		if (found)
			l[n] += elem; // i-1??
	}

}

vector<double> NoBlockSM::MultMatrByVect(vector<double> b)
{
	vector<double> v;
	int num_of_knots = b.size();
	v.resize(num_of_knots, 0);

	for (int i = 0; i < num_of_knots; i++)
		v[i] = d[i] * b[i];

	for (int i = 0; i < num_of_knots; i++)
		for (int j = ig[i]; j < ig[i + 1]; j++)
		{
			v[i] += l[j] * b[jg[j]];
			v[jg[j]] += u[j] * b[i];
		}

	return v;
}