#include "Header.h"
#include <iostream>
#include <string>

string static str_local_area, str_cross;

double calc_f(Knot coord, double time)   // Вычисление f
{
    //return 0.; // Не зависит от пространства и времени
    //return 0.; // Не зависит от времени
    //return 1.; // Не зависит от пространства
    //return 1.; // Зависит и от пространства и от времени
    //return 2. + 2. * time; // Зависит от времени квадратично
    //return 6 * time + 3. * time * time; // Зависит от времени кубически
    return 12. * time * time + 4. * time * time * time;
}

double calc_u(Knot coord, double time)   // Вычисление решения (u)
{
    //return 1.; // Не зависит от пространства и времени
    //return coord.x_y_z[0]; // Не зависит от времени
    //return time; // Не зависит от пространства
    //return coord.x_y_z[0] + time; // Зависит и от пространства и от времени
    //return time * time; // Зависит от времени квадратично
    //return time * time * time; // Зависит от времени кубически
    return time * time * time * time;
}

void time(Knot coord, double time) {

}

// Чтение координат узлов
void ReadCross(string pathFile, vector<Knot>& knots)
{
    ifstream in(pathFile);
    int nKnots;
    in >> nKnots;
    knots.resize(nKnots);
    for (int i = 0; i < nKnots; i++)  in >> knots[i].x >> knots[i].y >> knots[i].z;
    in.close();
}

// Чтение конечных элементов
void ReadLocals(string pathFile,    // Путь к файлу
                int N_KNOTS,        // Количество узлов в КЭ
                vector<Knot> knots, // Узлы
                vector<LocalArea>& locals)  // КЭ-ты
                
{
    int firstIndLocal = locals.size();

    ifstream in(pathFile);
    int nLocalArea;
    in >> nLocalArea;     // Считываем количество конечных элементов
    locals.resize(nLocalArea + locals.size()); // Инициализируем размерность массива, который хранит локальные области

    for (int indLocalArea = firstIndLocal; indLocalArea < nLocalArea; indLocalArea++)
    {
        locals[indLocalArea].globalNum = new int[N_KNOTS] {};
        locals[indLocalArea].x = new double[N_KNOTS]{};
        locals[indLocalArea].y = new double[N_KNOTS]{};
        locals[indLocalArea].z = new double[N_KNOTS]{};
        locals[indLocalArea].size = N_KNOTS;
        for (int j = 0; j < N_KNOTS; j++)
        {
            int globalNum;
            in >> globalNum;
            locals[indLocalArea].globalNum[j] = globalNum;   // локальные узлы конечного элемента
            locals[indLocalArea].x[j] = knots[globalNum].x;
            locals[indLocalArea].y[j] = knots[globalNum].y;
            locals[indLocalArea].z[j] = knots[globalNum].z;
        }

        in >> locals[indLocalArea].lambda;    // лямбда
        in >> locals[indLocalArea].sigma;    // гамма
        in >> locals[indLocalArea].hi;    // гамма
    }

    in.close();
}

void ReadBound(string pathFile, vector<Bound>& bounds)
{
    ifstream in(pathFile);
	int nBounds;
	in >> nBounds;
	bounds.resize(nBounds);
	for (int i = 0; i < nBounds; i++) in >> bounds[i].globalNum;
	in.close();
}

// Чтение сетки по времени
void ReadTime(string pathFile, TimeGrid& time)
{
    ifstream in(pathFile);
    in >> time.start >> time.end >> time.k >> time.nSteps;
    in.close();
}

void InitData(InitialData* data) //функция считывания исходных данных(кроме краевых)
{
    ReadCross("cross.txt", data->knots);

    ReadLocals("Triangle.txt", N_TRIANGLE_PRIZM_KNOTS, data->knots, data->locals);
    ReadLocals("Shestigran.txt", N_HEXAGON_KNOTS, data->knots, data->locals);

    ReadBound("EdgeBoundary_1.txt", data->bounds);

    ReadTime("grid_time.txt", *data->time_g);
}

void mult_matr_by_vect(int size, double** matrix, double* vector, double* result)
{
    for (int i = 0; i < size; i++)
    {
        double sum = 0;
        for (int j = 0; j < size; j++)
            sum += matrix[i][j] * vector[j];
        result[i] = sum;
    }
}

void mult_matr_by_vect(int size, double matrix[8][8], double vector[8], double result[8])
{
    for (int i = 0; i < size; i++)
    {
        double sum = 0;
        for (int j = 0; j < size; j++)
            sum += matrix[i][j] * vector[j];
        result[i] = sum;
    }
}



void write(InitialData form, double time) //функция вывода в консоль
{
    cout << "Параметры узлов и первые краевые :" << endl;
    cout << " _____________________________________ " << endl;
    cout.setf(ios::left);
    cout.width(15);
    cout << "| № элемента " << "  | ";
    cout.width(5);
    cout << "x" << "| ";
    cout.width(5);
    cout << "y" << "| ";
    cout.width(5);
    cout << "z" << "|" << endl;

    cout << "|----------------|------|------|------|" << endl;

    for (int i = 0; i < form.knots.size(); i++) {     // Заполняем координаты узлов для области
        cout << "| ";
        cout.width(15);
        cout << i + 1 << "| ";
        cout.width(5);
        cout << form.knots[i].x << "| " << form.knots[i].y << "| " << form.knots[i].z;
        cout << endl;
    }

    cout << endl << "Конечные элементы:" << endl;
    cout << " _____________________________________________________________________________ " << endl;
    cout.setf(ios::left);
    cout.width(15);
    cout << "| № элемента " << "  | ";
    cout.width(25);
    cout << "Узлы" << "| ";
    cout.width(15);
    cout << "lambda" << "| ";
    cout.width(15);
    cout << "gamma" << "|" << endl;
    cout << "|----------------|--------------------------|----------------|----------------|" << endl;

    for (int i = 0; i < form.num_locals; i++)
    {
        for (int j = 0; j < 8; j++) {
            str_local_area += to_string(form.locals[i].globalNum[j]) + " ";
        }
        cout << "| ";
        cout.width(15);
        cout << i+1 << "| ";
        cout.width(24);
        cout << str_local_area << " ";  // локальные узлы конечного эле-мента
        cout << "|";
        cout.width(15);
        cout << form.locals[i].lambda << " | ";
        cout.width(15);
        cout << form.locals[i].sigma << "| ";
        cout << endl;
        str_local_area = "";

    }     
}

void write_result(SLAU& slau, double time) //функция вывода в консоль
{
    cout << endl << "ВРЕМЯ: " << time << endl;
    cout << endl << "Результат в узлах (веса):" << endl;
    cout << " ___________________________________________________________________ " << endl;

    cout.setf(ios::left);
    cout.width(15);
    cout << "| № элемента " << "  | ";
    cout.width(15);
    cout << "u*" << "| ";
    cout.width(15);
    cout << "u" << "| ";
    cout.width(15);
    cout << "|u-u*|" << "|" << endl;
    cout << "|----------------|----------------|----------------|----------------|" << endl;

    for (int i = 0; i < slau.size; i++)
    {
        cout.setf(ios::left);
        cout << "| ";
        cout.width(15);
        cout << i + 1 << "| ";
        cout.width(15);
        cout << slau.q[i] << "| ";
        cout.width(15);
        cout << slau.u[i] << "| ";
        cout.width(15);
        cout << fabs(slau.q[i] - slau.u[i]) << "| " << endl;
    }
}

double scalar(int size, double* v1, double* v2)
{
    double sum = 0;
    for (int i = 0; i < size; i++)
        sum += v1[i] * v2[i];
    return sum;
}

void LOC(SLAU* slau)
{
    int i;
    double nvzk = 0., alfa = 0., beta = 0., skp = 0., eps = 9.999999682655226e-030;
    double* z, * r, * p, * f;
    slau->q = new double[slau->size]{};
    z = new double[slau->size]{};
    r = new double[slau->size]{};
    p = new double[slau->size]{};
    f = new double[slau->size]{};
    double lastnvzk;

    mult_matr_by_vect(slau->size, slau->global_A, slau->q, f);

    for (i = 0; i < slau->size; i++)
        z[i] = r[i] = slau->global_d[i] - f[i];

    mult_matr_by_vect(slau->size, slau->global_A, z, p);
    nvzk = sqrt(scalar(slau->size, r, r)) / sqrt(scalar(slau->size, slau->global_d, slau->global_d));

    for (int k = 1; k < 100000 && nvzk > eps; k++)
    {
        lastnvzk = nvzk;
        skp = scalar(slau->size, p, p);
        alfa = scalar(slau->size, p, r) / skp;

        for (i = 0; i < slau->size; i++)
        {
            slau->q[i] += alfa * z[i];
            r[i] -= alfa * p[i];
        }

        mult_matr_by_vect(slau->size, slau->global_A, r, f);
        beta = -scalar(slau->size, p, f) / skp;

        for (i = 0; i < slau->size; i++)
        {
            z[i] = r[i] + beta * z[i];
            p[i] = f[i] + beta * p[i];
        }

        nvzk = sqrt(scalar(slau->size, r, r)) / sqrt(scalar(slau->size, slau->global_d, slau->global_d));
    }
}

double get_u(Knot point, InitialData form) // Получение значения в произвольной точ-ке
{
    double x = point.x;
    double y = point.y;
    double z = point.z;

    int i;
    int igl[8]{};
    for (i = 0; i < form.num_locals; i++) // Определяем в каком конечном элементе точка
    {
        // Определяем глобальные узлы для локального элемента
        for (int j = 0; j < 8; j++)
            igl[j] = form.locals[i].globalNum[j] - 1;

        // Если указанные точки входят в область, то мы нашли, в каком элементе лежит точка
        if (x >= form.knots[igl[0]].x && x <= form.knots[igl[1]].x &&
            y >= form.knots[igl[0]].y && y <= form.knots[igl[2]].y &&
            z >= form.knots[igl[0]].z && z <= form.knots[igl[4]].z)
            break;
    }

    if (i == form.num_locals) return -1;

    double hx = form.locals[i].h_x;
    double hy = form.locals[i].h_y;
    double hz = form.locals[i].h_z;

    double X1 = (form.knots[igl[1]].x - x) / hx;
    double X2 = (x - form.knots[igl[0]].x) / hx;
    double Y1 = (form.knots[igl[2]].y - y) / hy;
    double Y2 = (y - form.knots[igl[0]].y) / hy;
    double Z1 = (form.knots[igl[4]].z - z) / hz;
    double Z2 = (z - form.knots[igl[0]].z) / hz;

    double u =
        form.q[igl[0]] * X1 * Y1 * Z1 +
        form.q[igl[1]] * X2 * Y1 * Z1 +
        form.q[igl[2]] * X1 * Y2 * Z1 +
        form.q[igl[3]] * X2 * Y2 * Z1 +
        form.q[igl[4]] * X1 * Y1 * Z2 +
        form.q[igl[5]] * X2 * Y1 * Z2 +
        form.q[igl[6]] * X1 * Y2 * Z2 +
        form.q[igl[7]] * X2 * Y2 * Z2;

    return u;
}

void calc_global_matrix(InitialData form, SLAU& slau) // Вычисление глобальной матрицы
{
    for (int i = 0; i < form.num_locals; i++)
    {
        double M[8][8]{};
        calc_local_M(form.locals[i], M);

        double G[8][8]{};
        calc_local_G(form.locals[i], G);

        double A[8][8]{};
        for (int j = 0; j < 8; j++)
            for (int k = 0; k < 8; k++)
                slau.global_A[form.locals[i].globalNum[j]][form.locals[i].globalNum[k]] += A[j][k] = M[j][k] + G[j][k];
    }
}

void calc_global_M(InitialData& form, SLAU& slau) // Вычисление глобальной M
{
    for (int i = 0; i < form.num_locals; i++)
    {
        int size = form.locals[i].size;
        double** M = new double*[size];
        for (int j = 0; j < size; j++) M[j] = new double[size];

        calc_local_M(form.locals[i], M);

        for (int j = 0; j < size; j++)
            for (int k = 0; k < size; k++)
                slau.global_M[form.locals[i].globalNum[j]][form.locals[i].globalNum[k]] += M[j][k];
    }
}
void calc_global_G(InitialData& form, SLAU& slau) // Вычисление глобальной G
{
    for (int i = 0; i < form.num_locals; i++)
    {
        double G[8][8]{};
        calc_local_G(form.locals[i], G);

        for (int j = 0; j < 8; j++)
            for (int k = 0; k < 8; k++)
                slau.global_G[form.locals[i].globalNum[j] - 1][form.locals[i].globalNum[k] - 1] += G[j][k];
    }
}

void calc_global_A(InitialData* form, SLAU& slau, double t_j, double t_j1, double t_j2, double t_j3) // Вычисление глобальной A
{
   for (int i = 0; i < form->num_locals; i++)
   {
      for (int j = 0; j < 8; j++)
      {
         int ind_1 = form->locals[i].globalNum[j] - 1;
         for (int k = 0; k < 8; k++)
         {
            int ind_2 = form->locals[i].globalNum[k] - 1;
            slau.global_A[ind_1][ind_2] =
               form->locals[i].hi * slau.global_M[ind_1][ind_2] *
                   2 * ((t_j - t_j1) + (t_j - t_j2) + (t_j - t_j3))
                   / ((t_j - t_j3) * (t_j - t_j2) * (t_j - t_j1)) +
               form->locals[i].sigma * slau.global_M[ind_1][ind_2]
                   * ((t_j - t_j2) * (t_j - t_j1) + (t_j - t_j3) * (t_j - t_j1) + (t_j - t_j3) * (t_j - t_j2))
                   / ((t_j - t_j3) * (t_j - t_j2) * (t_j - t_j1))
               + slau.global_G[ind_1][ind_2];
         }

      }

   }


   /* for (int j = 0; j < data->nKnots; j++)
        for (int k = 0; k < data->nKnots; k++)
            data->global_matrix[j][k] = data->global_M[j][k]* 
            ((t_j- t_j2)*(t_j- t_j1)+ (t_j - t_j3)*(t_j - t_j1)+ (t_j - t_j3)*(t_j - t_j2))/
            ((t_j - t_j3)*(t_j - t_j2)*(t_j - t_j1)) 
            + data->global_G[j][k];*/
}
void matrix_mult_vector(double** matrix, double* vector, int MAX, double* result) // Вычисление матрицы
{
    for (int i = 0; i < MAX; i++)
        for (int j = 0; j < MAX; j++)
            result[i] += matrix[j][i] * vector[i];
}

void calc_global_d(InitialData* form, SLAU& slau, double t_j, double t_j1, double t_j2, double t_j3, double* q_3, double* q_2, double* q_1) // Вычисление глобальной A
{

    calc_global_F(form, time);	// Вычисление глобального вектора правой части
    double* Mq_j3, * Mq_j2, * Mq_j1;
    Mq_j3 = new double[form->knots.size()]{};
    Mq_j2 = new double[form->knots.size()]{};
    Mq_j1 = new double[form->knots.size()]{};

    matrix_mult_vector(slau.global_M, q_3, form->knots.size(), Mq_j3);
    matrix_mult_vector(slau.global_M, q_2, form->knots.size(), Mq_j2);
    matrix_mult_vector(slau.global_M, q_1, form->knots.size(), Mq_j1);

   for (int i = 0; i < form->num_locals; i++)
   {
      for (int j = 0; j < 8; j++)
      {
         int ind_1 = form->locals[i].globalNum[j] - 1;
         slau.global_d[ind_1] = slau.global_b[ind_1]
            - form->locals[i].sigma * Mq_j3[ind_1] * ((t_j - t_j2) * (t_j - t_j1)) / ((t_j3 - t_j2) * (t_j3 - t_j1) * (t_j3 - t_j))
            - form->locals[i].sigma * Mq_j2[ind_1] * ((t_j - t_j3) * (t_j - t_j1)) / ((t_j2 - t_j3) * (t_j2 - t_j1) * (t_j2 - t_j))
            - form->locals[i].sigma * Mq_j1[ind_1] * ((t_j - t_j3) * (t_j - t_j2)) / ((t_j1 - t_j3) * (t_j1 - t_j2) * (t_j1 - t_j))
            - form->locals[i].hi * Mq_j3[ind_1] * 2 * ((t_j-t_j2)+(t_j-t_j1) )/ ((t_j3 - t_j) * (t_j3 - t_j2) * (t_j3 - t_j1))
            - form->locals[i].hi * Mq_j2[ind_1] * 2 * ((t_j - t_j3) + (t_j - t_j1)) / ((t_j2 - t_j3) * (t_j2 - t_j1) * (t_j2 - t_j))
            - form->locals[i].hi * Mq_j1[ind_1] * 2 * ((t_j - t_j2) + (t_j - t_j3)) / ((t_j1 - t_j3) * (t_j1 - t_j2) * (t_j1 - t_j));
      }
   }

 /*   for (int indLocalArea = 0; indLocalArea < data->nKnots; indLocalArea++)
        data->global_vector[indLocalArea] = data->global_b[indLocalArea] -
        Mq_j3[indLocalArea] * ((t_j - t_j2) * (t_j - t_j1)) / ((t_j3 - t_j2) * (t_j3 - t_j1)* (t_j3 - t_j)) -
        Mq_j2[indLocalArea] * ((t_j - t_j3) * (t_j - t_j1)) / ((t_j2 - t_j3) * (t_j2 - t_j1) * (t_j2 - t_j)) -
        Mq_j1[indLocalArea] * ((t_j - t_j3) * (t_j - t_j2)) / ((t_j1 - t_j3) * (t_j1 - t_j2) * (t_j1 - t_j));*/

}

void calc_global_F(InitialData* form, SLAU& slau, double time) // Вычисление глобального вектора правой части
{
    for (int i = 0; i < form->num_locals; i++)
    {
        double F[8]{};
        calc_local_F(form->knots, form->locals[i], F, time);

        for (int j = 0; j < 8; j++)
        {
            slau.global_b[form->locals[i].globalNum[j] - 1] += F[j];
        }
    }
}

void calc_first_boundary_conditions(InitialData* form, SLAU& slau, double time) //учет краевых условий
{
    for (int i = 0; i < form->num_bounds; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            int global_num_coord = form->bounds[i].globalNum[j] - 1;
            for (int k = 0; k < form->knots.size(); k++)
            {
                slau.global_A[global_num_coord][k] = 0;
            }
            slau.global_A[global_num_coord][global_num_coord] = 1;
            slau.global_d[global_num_coord] = calc_u(form->knots[global_num_coord], time);
        }

    }
}










void CreateScheme(TimeScheme& scheme, InitialData& data)
{
    // Подсчет начального шага
    scheme.k = data.time_g->k;
    scheme.h = data.time_g->end - data.time_g->start;
    scheme.h = (scheme.k == 1) ? scheme.h / data.time_g->nSteps : scheme.h * (1. - data.time_g->k) / (1. - pow(data.time_g->k, data.time_g->nSteps));


    double SLAUsize = data.knots.size();

    scheme.time[0] = data.time_g->start;
    scheme.time[1] = scheme.time[0] + scheme.h;
    scheme.time[2] = scheme.time[1] + scheme.h * scheme.k;
    scheme.time[3] = scheme.time[2] + scheme.h * scheme.k * scheme.k;
    
    for (int j = 0; j < 3; j++)
    {
        double time = scheme.time[j];
        scheme.q[j] = new double[SLAUsize] {};
        for (int i = 0; i < SLAUsize; i++) scheme.q[0][i] = calc_u(data.knots[i], time);
    }

    scheme.h = scheme.h * scheme.k * scheme.k;
}

void NextScheme(TimeScheme& scheme, int nKnots)
{
    for (int i = 0; i < 3; i++)
    {
        scheme.time[i] = scheme.time[i + 1];
        scheme.q[i] = scheme.q[i + 1];
    }

    scheme.time[3] = scheme.time[2] + scheme.h;
    scheme.h *= scheme.k;
    scheme.q[3] = new double[nKnots];
}


void CreateSLAU(int sizeSLAU, SLAU& slau)
{
    slau.size = sizeSLAU;
    slau.global_A = new double* [slau.size]{};  // Инициализируем размерность глобальной матрицы
    slau.global_M = new double* [slau.size]{};  // Инициализируем размерность глобальной матрицы
    slau.global_G = new double* [slau.size]{};  // Инициализируем размерность глобальной матрицы
   
    for (int i = 0; i < slau.size; i++) {
        slau.global_A[i] = new double[slau.size]{};
        slau.global_M[i] = new double[slau.size]{};
        slau.global_G[i] = new double[slau.size]{};
    }

    slau.global_d = new double[slau.size]{};   // Инициализируем размерность глобального вектора
    slau.q = new double[slau.size]{};   // Инициализируем размерность вектора весов
    slau.global_b = new double[slau.size]{};
}

