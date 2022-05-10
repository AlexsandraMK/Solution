#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <fstream>

static std::vector<double> MultMatrByVect(std::vector< std::vector<double>> matrix, std::vector<double> _vector)
{
	int size = _vector.size();

	std::vector<double> result;
	result.resize(size);

	for (int i = 0; i < size; i++)
	{
		double sum = 0;
		for (int j = 0; j < size; j++)
			sum += matrix[i][j] * _vector[j];
		result[i] = sum;
	}

	return result;
}

static double CalcScalar(std::vector<double> v1, std::vector<double> v2)
{
	int size = v1.size();

	double scalar = 0.;

	for (int i = 0; i < size; i++)
	{
		scalar += v1[i] * v2[i];
	}

	return scalar;
}

static double CalcDetMatrix(std::vector< std::vector<double>> matrix)
{
	if (matrix.size() == 3)
		return matrix[0][0] * matrix[1][1] * matrix[2][2]
				+ matrix[2][0] * matrix[0][1] * matrix[1][2]
				+ matrix[1][0] * matrix[2][1] * matrix[0][2]
				- matrix[2][0] * matrix[1][1] * matrix[0][2]
				- matrix[1][0] * matrix[0][1] * matrix[2][2]
				- matrix[0][0] * matrix[1][2] * matrix[2][1];
}


static std::vector< std::vector<double>> CalcReverseMatrixWithSize3(std::vector< std::vector<double>> matrix)
{
	std::vector< std::vector<double>> reverseMatrix;
	reverseMatrix.resize(3);
	for (int i = 0; i < reverseMatrix.size(); i++)	reverseMatrix[i].resize(reverseMatrix.size(),0);

	double detJacobian = CalcDetMatrix(matrix);

	if (fabs(detJacobian) < 1e-14) 
		std::cout << "ѕроизошла ошибка в подсчете обратного якобиана";

	for (int i = 0; i < reverseMatrix.size(); i++)
	{
		for (int j = 0; j < reverseMatrix.size(); j++)
		{
			double min[4]{};
			int k = 0;
			for (int im = 0; im < reverseMatrix.size(); im++)
			{
				for (int jm = 0; jm < reverseMatrix.size(); jm++)
				{
					if (im != i && jm != j)
						min[k++] = matrix[im][jm];
				}
			}

			reverseMatrix[j][i] = pow(-1, i + j + 2) * (min[0] * min[3] - min[1] * min[2]) / detJacobian;
		}
	}

	return reverseMatrix;
}

static void WriteMatrix(std::vector< std::vector<double>> matrix)
{
	for (int i = 0; i < matrix.size(); i++)
	{
		for (int j = 0; j < matrix.size(); j++)
			std::cout << matrix[i][j] << " ";
		std::cout << "\n";
	}
}

static void SwapLinesInMatrix(std::vector< std::vector<double>> &matrix, int l1, int l2)
{
	double temp;
	for (int i = 0; i < matrix.size(); i++)
	{
		temp = matrix[l1][i];
		matrix[l1][i] = matrix[l2][i];
		matrix[l2][i] = temp;
	}
}

static void SwapLinesInVector(std::vector<double> &vector, int l1, int l2)
{
	double temp = vector[l1];
	vector[l1] = vector[l2];
	vector[l2] = temp;
}


static void WriteMatrix(std::vector< std::vector<double>> matrix, std::string pathFile)
{
	std::ofstream out(pathFile);

	for (int i = 0; i < matrix.size(); i++)
	{
		for (int j = 0; j < matrix.size(); j++)
			out << matrix[i][j] << " ";
		out << "\n";
	}

	out.close();
}


static int ReadMatrix(std::vector< std::vector<double>> &matrix, std::string pathFile)
{
	std::ifstream in(pathFile);
	if (!in.is_open()) return 1;

	for (int i = 0; i < matrix.size(); i++)
	{
		for (int j = 0; j < matrix.size(); j++)
			in >> matrix[i][j];
	}

	in.close();
	return 0;
}