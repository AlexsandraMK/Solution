#include "InitialData.h"
#include <fstream>
#include <iostream>

/// <summary>
/// ������������� ������ (���������������� � ��������� �����)
/// </summary>
/// <param name="pathFile">��� �����</param>
/// <return name="timeGrid">����� �� �������</return>
InitialData::InitialData()
{
	cout << "������ ������:\n";
	ReadCoordinates("cross.txt");
	ReadCoefficients("coeffs.txt");
	ReadKEs("TriangularPrism.txt");

	ReadKEs("Hexagon.txt");

	ReadBounds("FirstBounds.txt");

	timeGrid = new TimeGrid();
	ReadTime("grid_time.txt");

	cout << "������ ������ ���������.\n";
}

/// <summary>
/// ������ �������������
/// </summary>
/// <param name="pathFile">��� �����</param>
/// <return name="knots">��������� �����</return>
void InitialData::ReadCoefficients(string pathFile)
{
	cout << "\t������ �������������: ";
	ifstream in(pathFile);
	if (!in.is_open())
	{
		cout << "\t���� �� ������\n";
		return;
	}
	int nCoeffs;
	in >> nCoeffs;
	coeffs.resize(nCoeffs);
	for (int i = 0; i < nCoeffs; i++)  in >> coeffs[i].density >> coeffs[i].Yung >> coeffs[i].Poisson;
	in.close();
	cout << "\t�������\n";
}


/// <summary>
/// ������ ��������� �����
/// </summary>
/// <param name="pathFile">��� �����</param>
/// <return name="knots">��������� �����</return>
void InitialData::ReadCoordinates(string pathFile)
{
	cout << "\t������ ���������: ";
	ifstream in(pathFile);
	if (!in.is_open())
	{
		cout << "\t���� �� ������\n";
		return;
	}
	int nKnots;
	in >> nKnots;
	knots.resize(nKnots);
	for (int i = 0; i < nKnots; i++)  in >> knots[i].x >> knots[i].y >> knots[i].z;
	in.close();
	cout << "\t�������\n";
}

/// <summary>
/// ������ �������� ���������
/// </summary>
/// <param name="pathFile">��� �����</param>
/// <param name="knots">���������� �����</param>
/// <return name="KEs">�������� ��������</return>
void InitialData::ReadKEs(string pathFile) // ����                
{
	cout << "\t������ �������� ��������� (���� - " << pathFile << " ): ";
	int indFirstKE = KEs.size();

	ifstream in;
	in.open(pathFile);
	if (!in.is_open())
	{
		cout << "\t���� �� ������\n";
		return;
	}
	int nKEsInFile;
	in >> nKEsInFile;     // ��������� ���������� �������� ���������
	KEs.resize(nKEsInFile + KEs.size()); // �������������� ����������� �������, ������� ������ ��������� �������

	IKE* ke;
	pathFile.erase(pathFile.find_last_of('.'));



	for (int iKE = indFirstKE; iKE < KEs.size(); iKE++)
	{
		if (pathFile == TriangularPrism::ToString())
			ke = new TriangularPrism();
		else
			if (pathFile == Hexagon::ToString())
				ke = new Hexagon();
			else
			{
				cout << "Error";
				exit(0);
			}
		int globalNum;
		int numKnots = ke->GetCountKnots();
		for (int j = 0; j < numKnots; j++)
		{
			in >> globalNum;
			ke->SetGlobalKnotNum(globalNum, knots[globalNum]);
		}

		//in >> ke->lambda;    // ������
		//in >> ke->sigma;     // �����
		//in >> ke->hi;        // ��
		in >> ke->iCoeff;


		KEs[iKE] = ke;
	}
	in.close();
	cout << "\t�������\n";
}


/// <summary>
/// ������ ������
/// </summary>
/// <param name="pathFile">��� �����</param>
/// <return name="bounds">�������</return>
void InitialData::ReadBounds(string pathFile)
{
	cout << "\t������ ������: ";
	ifstream in(pathFile);
	if (!in.is_open())
	{
		cout << "\t���� �� ������\n";
		return;
	}
	int nBounds;
	in >> nBounds;
	bounds.resize(nBounds);
	int boundSize = 4;
	for (int i = 0; i < nBounds; i++)
	{
		for (int j = 0; j < boundSize; j++)
			in >> bounds[i].globalNum[j];
	}
	in.close();
	cout << "\t�������\n";
}

/// <summary>
/// ������ ����� �� �������
/// </summary>
/// <param name="pathFile">��� �����</param>
/// <return name="timeGrid">����� �� �������</return>
void InitialData::ReadTime(string pathFile)
{
	cout << "\t������ �������: ";
	ifstream in(pathFile);
	if (!in.is_open())
	{
		cout << "\t���� �� ������\n";
		return;
	}
	in >> timeGrid->start >> timeGrid->startAfterTime3 >> timeGrid->end >> timeGrid->kAfterTime3 >> timeGrid->nStepsAfterTime3;
	in.close();
	cout << "\t�������\n";
}

double InitialData::ChangeCoordinates(vector<double> q)
{
	double max = 0.;

	for (int i = 0; i < knots.size(); i++)
	{
		knots[i].x += q[i * 3];
		knots[i].y += q[i * 3 + 1];
		knots[i].z += q[i * 3 + 2];
		if (abs(q[i * 3]) > max)
			max = abs(q[i * 3]);
	}

	return max;
}


