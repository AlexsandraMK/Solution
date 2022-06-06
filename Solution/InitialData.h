#pragma once
#include "IKE.h"
#include "Knot.h"

struct TimeGrid
{
    double start;
    double timeToGo;
    double startAfterTime3;	// ��������� �����
    double end;	// �������� �����
    int nStepsAfterTime3;	// ���������� ����� �� �������
    double kAfterTime3;	// ��������� ��� �������� ��������� �����
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
    Knot* knotToGo;
    Knot* uToGo;
    TimeGrid* timeGrid;
    InitialData();

    void ReadCoefficients(string pathFile);
    double ChangeCoordinates(vector<double> q);
protected:
    void ReadCoordinates(string pathFile);
    void ReadKEs(string pathFile);
    void ReadBounds(string pathFile);
    void ReadTime(string pathFile);
    void ReadPointToGo(string pathFile);
    
};

