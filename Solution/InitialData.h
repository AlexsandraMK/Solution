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


class InitialData
{
public:
    vector<Knot> knots; // ������ ��������� �����
    vector<IKE*> KEs;
    vector<Bound> bounds; // ������ ������
    TimeGrid* timeGrid;
    InitialData();

protected:
    void ReadCoordinates(string pathFile);
    void ReadKEs(string pathFile);
    void ReadBounds(string pathFile);
    void ReadTime(string pathFile);
};

