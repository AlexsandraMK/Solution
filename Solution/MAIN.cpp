#include "SLAU.h"
#include <locale.h>


int main()
{
	setlocale(LC_ALL, "Russian");
	InitialData* data = new InitialData();   // Объявили область

	TimeScheme* timeScheme = new TimeScheme(data);

	SLAU* slau = new SLAU(data);
	//slau->WriteResultForSolution(timeScheme->q[0], timeScheme->time[0]);
	//slau->WriteResultForSolution(timeScheme->q[1], timeScheme->time[1]);
	//slau->WriteResultForSolution(timeScheme->q[2], timeScheme->time[2]);
	slau->WriteResultForTest(timeScheme->q[0], timeScheme->time[0]);
	slau->WriteResultForTest(timeScheme->q[1], timeScheme->time[1]);
	slau->WriteResultForTest(timeScheme->q[2], timeScheme->time[2]);

	Knot* knot = new Knot();
	knot->x = 1;
	knot->y = 1;
	knot->z = 0.7;

	// Поиск векторов весов на каждом шаге
	for (; timeScheme->time[3] <= data->timeGrid->end; timeScheme->Next())
	{
		slau->SolveSLAU(data, timeScheme);
		//slau->WriteResultForSolution(slau->q, timeScheme->time[3]);
		slau->WriteResultForTest(slau->q, timeScheme->time[3]);
		timeScheme->q[3] = slau->q;

		cout << data->KEs[9]->SolveInPoint(*knot, slau->q);		//тр.призмы 1		x2 3		z2 9
	}

	return 0;
}
