#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include "LocalSolution.h"
#include <vector>



struct Bound // Граница
{
    int globalNum; // Номера узла
};

struct TimeGrid
{
    double start;	// Начальное время
    double end;	// Конечное время
    int nSteps;	// Количество шагов по времени
    double k;	// Множитель для подсчета следующих шагов
};

struct TimeScheme {
    double time[4]; // 4 слоя по времени
    double** q;
    double k;
    double h;   
};

//TimeGrid* TIME_GRID;	// сетка по времени

struct InitialData // Исходные данные
{
    vector<Knot> knots; // Массив координат узлов
    vector<LocalArea> locals;
    vector<Bound> bounds; // Массив границ

    int num_locals; // Количество  локальных элементов
    int num_bounds; // Количество границ
    TimeGrid* time_g;
};

struct SLAU {
    int size;
    double** global_M; // Глобальная M
    double** global_G; // Глобальная G
    double** global_A; // Глобальная A
    double* global_b;
    double* global_d;   // Глобальный вектор правой части
    double* q;  // Решение в узлах (вектор весов)
    double* u;
};



void InitData(InitialData* form); // Чтение исходных данных
void calc_global_matrix(InitialData* form); // Вычисление глобальной матрицы

void calc_global_M(InitialData* form);
void calc_global_G(InitialData* form);
void calc_global_A(InitialData* form, double t_j, double t_j1, double t_j2, double t_j3);

void calc_global_F(InitialData* form, double time);	 // Вычисление глобального вектора правой части
void calc_global_d(InitialData* form, double t_j, double t_j1, double t_j2, double t_j3, double* q_3, double* q_2, double* q_1);

void calc_first_boundary_conditions(InitialData* form, double time); // Учет первых краевых условий
double calc_u(Knot coord, double time); // Решение уравнения
double calc_f(Knot coord, double time);
void write(InitialData form, double time); // Вывод
void write_result(InitialData form, double time);
void write_matrix(InitialData form); // Вывод матрицы
void write_vector(InitialData form); // Вывод вектора
void mult_matr_by_vect(int size, double** matrix, double* vector, double* result); // Умножение матрицы на вектор
double get_u(Knot point, InitialData form); // Получение значения в произвольной точке
double scalar(int size, double* v1, double* v2); // Скалярное умножение векторов
void LOC(InitialData* form); // Локально-Оптимальная Схема


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

