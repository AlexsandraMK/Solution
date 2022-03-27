#include "SLAU.h"
#include <locale.h>


int main()
{
	setlocale(LC_ALL, "Russian");
	InitialData* data = new InitialData();   // Объявили область

	TimeScheme* timeScheme = new TimeScheme(data);

	SLAU* slau = new SLAU(data);
	slau->WriteResult(timeScheme->q[0], timeScheme->time[0]);
	slau->WriteResult(timeScheme->q[1], timeScheme->time[1]);
	slau->WriteResult(timeScheme->q[2], timeScheme->time[2]);
	// Поиск векторов весов на каждом шаге
	for (; timeScheme->time[3] <= data->timeGrid->end; timeScheme->Next())
	{
		slau->SolveSLAU(data, timeScheme);
		slau->WriteResult(slau->q, timeScheme->time[3]);
		timeScheme->q[3] = slau->q;
	}

	return 0;
}
