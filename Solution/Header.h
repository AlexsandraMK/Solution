#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include "LocalSolution.h"
#include <vector>



struct Bound // �������
{
    int globalNum; // ������ ����
};

struct TimeGrid
{
    double start;	// ��������� �����
    double end;	// �������� �����
    int nSteps;	// ���������� ����� �� �������
    double k;	// ��������� ��� �������� ��������� �����
};

struct TimeScheme {
    double time[4]; // 4 ���� �� �������
    double** q;
    double k;
    double h;   
};

//TimeGrid* TIME_GRID;	// ����� �� �������

struct InitialData // �������� ������
{
    vector<Knot> knots; // ������ ��������� �����
    vector<LocalArea> locals;
    vector<Bound> bounds; // ������ ������

    int num_locals; // ����������  ��������� ���������
    int num_bounds; // ���������� ������
    TimeGrid* time_g;
};

struct SLAU {
    int size;
    double** global_M; // ���������� M
    double** global_G; // ���������� G
    double** global_A; // ���������� A
    double* global_b;
    double* global_d;   // ���������� ������ ������ �����
    double* q;  // ������� � ����� (������ �����)
    double* u;
};



void InitData(InitialData* form); // ������ �������� ������
void calc_global_matrix(InitialData* form); // ���������� ���������� �������

void calc_global_M(InitialData* form);
void calc_global_G(InitialData* form);
void calc_global_A(InitialData* form, double t_j, double t_j1, double t_j2, double t_j3);

void calc_global_F(InitialData* form, double time);	 // ���������� ����������� ������� ������ �����
void calc_global_d(InitialData* form, double t_j, double t_j1, double t_j2, double t_j3, double* q_3, double* q_2, double* q_1);

void calc_first_boundary_conditions(InitialData* form, double time); // ���� ������ ������� �������
double calc_u(Knot coord, double time); // ������� ���������
double calc_f(Knot coord, double time);
void write(InitialData form, double time); // �����
void write_result(InitialData form, double time);
void write_matrix(InitialData form); // ����� �������
void write_vector(InitialData form); // ����� �������
void mult_matr_by_vect(int size, double** matrix, double* vector, double* result); // ��������� ������� �� ������
double get_u(Knot point, InitialData form); // ��������� �������� � ������������ �����
double scalar(int size, double* v1, double* v2); // ��������� ��������� ��������
void LOC(InitialData* form); // ��������-����������� �����


//std::function<double(double, double, double, int, int, double*, double*, double*)> Mij;

double Integrate(function<double(double, double, double, int, int, double*, double*, double*)> f, int i, int j, double* x, double* y, double* z);

int mu(int index);

int v(int index);

int nu(int index);

double W(int index, double alpha);

double d_phi(int index, int what, double ksi, double etta, double tetha);

double prime_by_var(int what, double* Knot, double ksi, double etta, double tetha);

double phi(int index, double ksi, double etta, double tetha);

double det_Jacobian(double* x, double* y, double* z, double ksi, double etta, double tetha);

void CreateTimeScheme(TimeGrid& grid, TimeScheme& scheme);

void CreateSLAU(int sizeSLAU, SLAU& slau);

void NextScheme(TimeScheme& scheme, int nKnots);

