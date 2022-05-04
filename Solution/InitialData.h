#pragma once
#include "IKE.h"
#include "Knot.h"

struct TimeGrid
{
    double start;
    double startAfterTime3;	// Начальное время
    double end;	// Конечное время
    int nStepsAfterTime3;	// Количество шагов по времени
    double kAfterTime3;	// Множитель для подсчета следующих шагов
};


struct Bound // Граница
{
    int globalNum[4]; // Номера узла
};

struct Coeff
{
    double density; // плотность
    double Yung;    // модуль Юнга
    double Poisson; // коэффициент Пуассона
};


class InitialData
{
public:
    vector<Knot> knots; // Массив координат узлов
    vector<IKE*> KEs;
    vector<Bound> bounds; // Массив границ
    vector<Coeff> coeffs;
    TimeGrid* timeGrid;
    InitialData();

    void ReadCoefficients(string pathFile);

protected:
    void ReadCoordinates(string pathFile);
    void ReadKEs(string pathFile);
    void ReadBounds(string pathFile);
    void ReadTime(string pathFile);
};

