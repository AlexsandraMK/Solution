
#include <vector>

std::vector<double> MultMatrByVect(std::vector< std::vector<double>> matrix, std::vector<double> _vector)
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
}

double CalcScalar(std::vector<double> v1, std::vector<double> v2)
{
	int size = v1.size();

	double scalar = 0.;

	for (int i = 0; i < size; i++)
	{
		scalar += v1[i] * v2[i];
	}

	return scalar;
}

double CalcDetMatrix(std::vector< std::vector<double>> matrix)
{
	if (matrix.size() == 3)
		return matrix[0][0] * matrix[1][1] * matrix[2][2]
				+ matrix[2][0] * matrix[0][1] * matrix[1][2]
				+ matrix[1][0] * matrix[2][1] * matrix[0][2]
				- matrix[2][0] * matrix[1][1] * matrix[0][2]
				- matrix[1][0] * matrix[0][1] * matrix[2][2]
				- matrix[0][0] * matrix[1][2] * matrix[2][1];
}