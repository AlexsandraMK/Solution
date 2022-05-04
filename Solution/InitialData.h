#pragma once
#include "IKE.h"
#include "Knot.h"

struct TimeGrid
{
    double start;	// ��������� �����
    double end;	// �������� �����
    int nSteps;	// ���������� ����� �� �������
    double k;	// ��������� ��� �������� ��������� �����
};


struct Bound // �������
{
    int globalNum[4]; // ������ ����
};

struct Coeff
{
    double density; // ���������
    double Yung;    // ������ ����
    double Poisson; // ����������� ��������
};


class InitialData
{
public:
    vector<Knot> knots; // ������ ��������� �����
    vector<IKE*> KEs;
    vector<Bound> bounds; // ������ ������
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

