#pragma once
#include "vectorInteractions.h"

class GaussSolverSLAE
{
public:
	GaussSolverSLAE(std::vector<std::vector<double>> A, std::vector<double> F);

	void Solve();

	std::vector<std::vector<double>> A;
	std::vector<double> F;
	std::vector<double> x;
};

