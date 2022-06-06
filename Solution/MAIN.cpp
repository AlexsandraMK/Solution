#include "SLAU.h"
#include <locale.h>

int FindPoint(vector<Knot> knots, Knot* findKnot)
{
	for (int i = 0; i < knots.size(); i++)
	{
		if (knots[i].x == findKnot->x &&
			knots[i].y == findKnot->y &&
			knots[i].z == findKnot->z) return i;
	}

	return -1;
}


int main()
{
	setlocale(LC_ALL, "Russian");
	InitialData* data = new InitialData();   // Объявили область

	int iPoint = FindPoint(data->knots, data->knotToGo);
	if (iPoint == -1)
	{
		cout << "Точка не найдена!";
		exit(1);
	}
	ofstream file("./res/ResultInPoint.txt");
	file.setf(ios::scientific);
	TimeScheme* timeScheme = new TimeScheme(data);


	file.precision(20);
	file << timeScheme->time[0] << " " << timeScheme->qx[0][iPoint] << endl;
	file << timeScheme->time[1] << " " << timeScheme->qx[1][iPoint] << endl;
	file << timeScheme->time[2] << " " << timeScheme->qx[2][iPoint] << endl;
	FEM* slau = new FEM(data);
	
	// Поиск векторов весов на каждом шаге
	for (int i = 0; timeScheme->time[3] <= data->timeGrid->end;i++, timeScheme->Next())
	{
		slau->SolveSLAU(data, timeScheme);
		slau->WriteResultForTest(slau->q, timeScheme->time[3]);

		timeScheme->qx[3] = slau->qx;	
		timeScheme->qy[3] = slau->qy;
		timeScheme->qz[3] = slau->qz;

		file << timeScheme->time[3] << " " << timeScheme->qx[3][iPoint] << " " << timeScheme->qy[3][iPoint] << " " << timeScheme->qz[3][iPoint] << endl;
		//if(i % 3200 == 0) slau->SolveInAreaForTest(data, timeScheme->time[3]);
		cout << endl << (timeScheme->time[3]- data->timeGrid->startAfterTime3) / (data->timeGrid->end-data->timeGrid->startAfterTime3) * 100 << "%" << endl;
	}

	return 0;
}
