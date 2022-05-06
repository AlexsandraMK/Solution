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
	//slau->WriteResultForTest(timeScheme->qx[0], timeScheme->time[0]);
	//slau->WriteResultForTest(timeScheme->qx[1], timeScheme->time[1]);
	//slau->WriteResultForTest(timeScheme->qx[2], timeScheme->time[2]);

	// Поиск векторов весов на каждом шаге
	for (; timeScheme->time[3] <= data->timeGrid->end; timeScheme->Next())
	{
		slau->SolveSLAU(data, timeScheme);
		/*slau->WriteResultForSolution(slau->q, timeScheme->time[3]);*/
		slau->WriteResultForTest(slau->q, timeScheme->time[3]);
		//slau->WriteResultForTest(slau->qy, timeScheme->time[3]);
		//slau->WriteResultForTest(slau->qz, timeScheme->time[3]);
		timeScheme->qx[3] = slau->qx;	
		timeScheme->qy[3] = slau->qy;
		timeScheme->qz[3] = slau->qz;
		/*slau->SolveInAreaForTest(data, timeScheme->time[3]);*/
		//slau->SolveInArea(data, timeScheme->time[3]);
	}

	return 0;
}
