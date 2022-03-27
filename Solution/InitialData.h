#pragma once
#include "IKE.h"
#include "Knot.h"

struct TimeGrid
{
    double start;	// Начальное время
    double end;	// Конечное время
    int nSteps;	// Количество шагов по времени
    double k;	// Множитель для подсчета следующих шагов
};


struct Bound // Граница
{
    int globalNum[4]; // Номера узла
};


class InitialData
{
public:
    vector<Knot> knots; // Массив координат узлов
    vector<IKE*> KEs;
    vector<Bound> bounds; // Массив границ
    TimeGrid* timeGrid;
    InitialData();

protected:
    void ReadCoordinates(string pathFile);
    void ReadKEs(string pathFile);
    void ReadBounds(string pathFile);
    void ReadTime(string pathFile);
};

