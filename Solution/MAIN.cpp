#include "SLAU.h"
#include <locale.h>


int main()
{
	setlocale(LC_ALL, "Russian");
	InitialData* data = new InitialData();   // Объявили область
	ofstream file("./res/ResultInPoint.txt");
	file.setf(ios::scientific);
	TimeScheme* timeScheme = new TimeScheme(data);
	file << timeScheme->time[0] << " " << timeScheme->qx[0][364] << endl;
	file << timeScheme->time[1] << " " << timeScheme->qx[1][364] << endl;
	file << timeScheme->time[2] << " " << timeScheme->qx[2][364] << endl;
	SLAU* slau = new SLAU(data);
	//slau->WriteResultForSolution(timeScheme->qx[0], timeScheme->time[0]);
	//slau->WriteResultForSolution(timeScheme->qx[1], timeScheme->time[1]);
	//slau->WriteResultForSolution(timeScheme->qx[2], timeScheme->time[2]);
	//slau->WriteResultForTest(timeScheme->q[0], timeScheme->time[0]);
	//slau->WriteResultForTest(timeScheme->q[1], timeScheme->time[1]);
	//slau->WriteResultForTest(timeScheme->q[2], timeScheme->time[2]);
	
	// Поиск векторов весов на каждом шаге
	for (; timeScheme->time[3] <= data->timeGrid->end; timeScheme->Next())
	{
		slau->SolveSLAU(data, timeScheme);
		/*slau->WriteResultForSolution(slau->q, timeScheme->time[3]);*/
		//if(timeScheme->time[3] == 3. ) 
		slau->WriteResultForTest(slau->q, timeScheme->time[3]);
		//slau->WriteResultForTest(slau->qy, timeScheme->time[3]);
		//slau->WriteResultForTest(slau->qz, timeScheme->time[3]);

		timeScheme->qx[3] = slau->qx;	
		timeScheme->qy[3] = slau->qy;
		timeScheme->qz[3] = slau->qz;

		file << timeScheme->time[3] << " " << timeScheme->qx[3][364] << endl;
		slau->SolveInAreaForTest(data, timeScheme->time[3]);
	}

	return 0;
}
