#include "SLAU.h"
#include <locale.h>


int main()
{
	setlocale(LC_ALL, "Russian");
	InitialData* data = new InitialData();   // Объявили область
	ofstream file("./res/ResultInPoint.txt");
	file.setf(ios::scientific);
	TimeScheme* timeScheme = new TimeScheme(data);

	int iPoint = 726;
	file << timeScheme->time[0] << " " << timeScheme->qx[0][iPoint] << endl;
	file << timeScheme->time[1] << " " << timeScheme->qx[1][iPoint] << endl;
	file << timeScheme->time[2] << " " << timeScheme->qx[2][iPoint] << endl;
	SLAU* slau = new SLAU(data);
	
	// Поиск векторов весов на каждом шаге
	for (; timeScheme->time[3] <= data->timeGrid->end; timeScheme->Next())
	{
		slau->SolveSLAU(data, timeScheme);
		slau->WriteResultForTest(slau->q, timeScheme->time[3]);

		timeScheme->qx[3] = slau->qx;	
		timeScheme->qy[3] = slau->qy;
		timeScheme->qz[3] = slau->qz;

		file << timeScheme->time[3] << " " << timeScheme->qx[3][iPoint] << endl;
		slau->SolveInAreaForTest(data, timeScheme->time[3]);
	}

	return 0;
}
