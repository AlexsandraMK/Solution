#include "Header.h"



int main()
{
	setlocale(LC_ALL, "Russian");
	InitialData* data = new InitialData;   // Объявили область
	InitData(data);       // Чтение исходных данных

	TimeScheme* timeScheme = new TimeScheme;
	CreateTimeScheme(*data->time_g, *timeScheme);

	SLAU* slau = new SLAU;
	CreateSLAU(data->knots.size(), *slau);

	calc_global_M(data);
	calc_global_G(data);

	// Поиск векторов весов на каждом шаге
	for (; timeScheme->time[3] <= data->time_g->end; NextScheme(*timeScheme, slau->size))
	{
		calc_global_A(data, timeScheme->time);
		calc_global_d(data, timeScheme);
		calc_first_boundary_conditions(data, timeScheme->time[3]);	// Учет первых краевых условий
		LOC(data);	// Решение СЛАУ
		write_result(*data, timeScheme->time[3]);
	}

	return 0;
}
