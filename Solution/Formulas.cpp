#include "Header.h"
#include <iostream>
#include <string>

string static str_local_area, str_cross;

double calc_f(Knot coord, double time)   // ���������� f
{
    //return 0.; // �� ������� �� ������������ � �������
    //return 0.; // �� ������� �� �������
    //return 1.; // �� ������� �� ������������
    //return 1.; // ������� � �� ������������ � �� �������
    //return 2. + 2. * time; // ������� �� ������� �����������
    //return 6 * time + 3. * time * time; // ������� �� ������� ���������
    return 12. * time * time + 4. * time * time * time;
}

double calc_u(Knot coord, double time)   // ���������� ������� (u)
{
    //return 1.; // �� ������� �� ������������ � �������
    //return coord.x_y_z[0]; // �� ������� �� �������
    //return time; // �� ������� �� ������������
    //return coord.x_y_z[0] + time; // ������� � �� ������������ � �� �������
    //return time * time; // ������� �� ������� �����������
    //return time * time * time; // ������� �� ������� ���������
    return time * time * time * time;
}

void time(Knot coord, double time) {

}


function<double(double, double, double, int, int, double*, double*, double*)> Mij = [](double ksi, double etta, double tetha, int i, int j, double* x, double* y, double* z)
{
    return phi(i, ksi, etta, tetha) * phi(j, ksi, etta, tetha) * det_Jacobian(x, y, z, ksi, etta, tetha);
};

double* calc_grad(int index, double ksi, double etta, double tetha)
{
    return new double[3]{ d_phi(index,1,ksi,etta,tetha), d_phi(index,2,ksi,etta,tetha), d_phi(index,3,ksi,etta,tetha) };
}

function<double(double, double, double, int, int, double*, double*, double*)> Gij = [](double ksi, double etta, double tetha, int i, int j, double* x, double* y, double* z)
{
  

    double* J_grad_i = new double[3]{};
    double* J_grad_j = new double[3]{};

    double** reversed_Jacobian = new double* [3]
    {
        new double[3]{},
        new double[3]{},
        new double[3]{}
    };


    double** Jacobian = new double* [3]
    {
        new double[3]{prime_by_var(1,x, ksi, etta, tetha),prime_by_var(1,y, ksi, etta, tetha),prime_by_var(1,z, ksi, etta, tetha)},
        new double[3]{prime_by_var(2,x, ksi,etta,tetha) ,prime_by_var(2,y, ksi,etta,tetha),prime_by_var(2,z, ksi,etta,tetha)},
        new double[3]{prime_by_var(3,x, ksi,etta,tetha),prime_by_var(3,y, ksi,etta,tetha),prime_by_var(3,z, ksi,etta,tetha)}
    };

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            double min[4]{};
            int k = 0;
            for (int im = 0; im < 3; im++)
                for (int jm = 0; jm < 3; jm++)
                {
                    if (im != i && jm != j)
                        min[k++] = Jacobian[im][jm];
                }
            reversed_Jacobian[j][i] = pow(-1, i + j + 2) * (min[0] * min[3] - min[1] * min[2]);
        }
    }
        
      
  
    mult_matr_by_vect(3, reversed_Jacobian, calc_grad(i, ksi, etta, tetha), J_grad_i);
    mult_matr_by_vect(3, reversed_Jacobian, calc_grad(j, ksi, etta, tetha), J_grad_j);

    return scalar(3, J_grad_i, J_grad_j) / det_Jacobian(x,y,z,ksi,etta,tetha);
};



// ������ ��������� �����
void ReadCross(string pathFile, vector<Knot>& knots)
{
    ifstream in(pathFile);
    int nKnots;
    in >> nKnots;
    knots.resize(nKnots);
    for (int i = 0; i < nKnots; i++)  in >> knots[i].x >> knots[i].y >> knots[i].z;
    in.close();
}

// ������ �������� ���������
void ReadLocals(string pathFile,    // ���� � �����
                int N_KNOTS,        // ���������� ����� � ��
                vector<Knot> knots, // ����
                vector<LocalArea>& locals)  // ��-��
                
{
    int firstIndLocal = locals.size();

    ifstream in(pathFile);
    int nLocalArea;
    in >> nLocalArea;     // ��������� ���������� �������� ���������
    locals.resize(nLocalArea + locals.size()); // �������������� ����������� �������, ������� ������ ��������� �������

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
            locals[indLocalArea].globalNum[j] = globalNum;   // ��������� ���� ��������� ��������
            locals[indLocalArea].x[j] = knots[globalNum].x;
            locals[indLocalArea].y[j] = knots[globalNum].y;
            locals[indLocalArea].z[j] = knots[globalNum].z;
        }

        in >> locals[indLocalArea].lambda;    // ������
        in >> locals[indLocalArea].sigma;    // �����
        in >> locals[indLocalArea].hi;    // �����
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

// ������ ����� �� �������
void ReadTime(string pathFile, TimeGrid& time)
{
    ifstream in(pathFile);
    in >> time.start >> time.end >> time.q >> time.nSteps;
    in.close();
}

void InitData(InitialData* data) //������� ���������� �������� ������(����� �������)
{
    ReadCross("cross.txt", data->knots);

    ReadLocals("Triangle.txt", N_TRIANGLE_PRIZM_KNOTS, data->knots, data->locals);
    ReadLocals("Shestigran.txt", N_HEXAGON_KNOTS, data->knots, data->locals);

    ReadBound("EdgeBoundary_1.txt", data->bounds);

    ReadTime("grid_time.txt", *data->time_g);

    // ��������� � ������ �����

    data->global_matrix = new double* [data->knots.size()]{};  // �������������� ����������� ���������� �������
    data->global_M = new double* [data->knots.size()]{};  // �������������� ����������� ���������� �������
    data->global_G = new double* [data->knots.size()]{};  // �������������� ����������� ���������� �������
    //---data->global_A = new double* [data->nKnots]{};  // �������������� ����������� ���������� �������
    for (int i = 0; i < data->knots.size(); i++) {
        data->global_matrix[i] = new double[data->knots.size()]{};
        data->global_M[i] = new double[data->knots.size()]{};
        data->global_G[i] = new double[data->knots.size()]{};
        //---data->global_A[indLocalArea] = new double[data->nKnots]{};
    }

    data->global_vector = new double[data->knots.size()]{};   // �������������� ����������� ����������� �������
    data->q = new double[data->knots.size()]{};   // �������������� ����������� ������� �����
    data->global_b = new double[data->knots.size()]{};
 
    // ������� ���������� ����
    data->time_g->h0 = data->time_g->end - data->time_g->start;

    if (data->time_g->q != 1)
    {
        // ��� �������������
        data->time_g->h0 *= (1. - data->time_g->q) / (1. - pow(data->time_g->q, data->time_g->nSteps));
    }
    else
    {
        // ��� �����������
        data->time_g->h0 /= data->time_g->nSteps;
    }
    
}

void calc_local_M(LocalArea local, double M[8][8]) // ���������� ��������� ������� �����
{
    //double M_1[8][8] =
    //{
    //    {8, 4, 4, 2, 4, 2, 2, 1},
    //    {4, 8, 2, 4, 2, 4, 1, 2},
    //    {4, 2, 8, 4, 2, 1, 4, 2},
    //    {2, 4, 4, 8, 1, 2, 2, 4},
    //    {4, 2, 2, 1, 8, 4, 4, 2},
    //    {2, 4, 1, 2, 4, 8, 2, 4},
    //    {2, 1, 4, 2, 4, 2, 8, 4},
    //    {1, 2, 2, 4, 2, 4, 4, 8}
    //};

    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
            M[i][j] = Integrate(Mij,i,j, local.x, local.y, local.z);
}

void calc_local_G(LocalArea local, double G[8][8]) // ���������� ��������� ������� ���������
{
    //double G_x[8][8] =  // ������� ���������(x)
    //{
    //    {4, -4, 2, -2, 2, -2, 1, -1},
    //    {-4, 4, -2, 2, -2, 2, -1, 1},
    //    {2, -2, 4, -4, 1, -1, 2, -2},
    //    {-2, 2, -4, 4, -1, 1, -2, 2},
    //    {2, -2, 1, -1, 4, -4, 2, -2},
    //    {-2, 2, -1, 1, -4, 4, -2, 2},
    //    {1, -1, 2, -2, 2, -2, 4, -4},
    //    {-1, 1, -2, 2, -2, 2, -4, 4}
    //};

    //double G_y[8][8] =  // ������� ���������(y)
    //{
    //    {4, 2, -4, -2, 2, 1, -2, -1},
    //    {2, 4, -2, -4, 1, 2, -1, -2},
    //    {-4, -2, 4, 2, -2, -1, 2, 1},
    //    {-2, -4, 2, 4, -1, -2, 1, 2},
    //    {2, 1, -2, -1, 4, 2, -4, -2},
    //    {1, 2, -1, -2, 2, 4, -2, -4},
    //    {-2, -1, 2, 1, -4, -2, 4, 2},
    //    {-1, -2, 1, 2, -2, -4, 2, 4}
    //};

    //double G_z[8][8] =  // ������� ���������(z)
    //{
    //    {4, 2, 2, 1, -4, -2, -2, -1},
    //    {2, 4, 1, 2, -2, -4, -1, -2},
    //    {2, 1, 4, 2, -2, -1, -4, -2},
    //    {1, 2, 2, 4, -1, -2, -2, -4},
    //    {-4, -2, -2, -1, 4, 2, 2, 1},
    //    {-2, -4, -1, -2, 2, 4, 1, 2},
    //    {-2, -1, -4, -2, 2, 1, 4, 2},
    //    {-1, -2, -2, -4, 1, 2, 2, 4}
    //};

    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
            G[i][j] = local.lambda * Integrate(Gij,i,j,local.x,local.y,local.z);
            /*local.h_x * local.h_y * local.h_z / 36 *
            (
                1 / (local.h_x * local.h_x) * G_x[indLocalArea][j] +
                1 / (local.h_y * local.h_y) * G_y[indLocalArea][j] +
                1 / (local.h_z * local.h_z) * G_z[indLocalArea][j]
                );*/
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

void calc_local_F(vector<Knot> knots, LocalArea local, double F[8], double time) // ���������� ���������� ������� ������ �����
{
    double M_1[8][8]{};
    calc_local_M(local, M_1);

 /*   double M_1[8][8] =
    {
        {8, 4, 4, 2, 4, 2, 2, 1},
        {4, 8, 2, 4, 2, 4, 1, 2},
        {4, 2, 8, 4, 2, 1, 4, 2},
        {2, 4, 4, 8, 1, 2, 2, 4},
        {4, 2, 2, 1, 8, 4, 4, 2},
        {2, 4, 1, 2, 4, 8, 2, 4},
        {2, 1, 4, 2, 4, 2, 8, 4},
        {1, 2, 2, 4, 2, 4, 4, 8}
    };*/

    double f[8]{};

    for (int i = 0; i < 8; i++)
        f[i] = calc_f(knots[local.globalNum[i] - 1], time);

    mult_matr_by_vect(8, M_1, f, F);

    for (int i = 0; i < 8; i++)
        F[i] *= local.h_x * local.h_y * local.h_z / 216;
       
}

void write(InitialData form, double time) //������� ������ � �������
{
    cout << "��������� ����� � ������ ������� :" << endl;
    cout << " _____________________________________ " << endl;
    cout.setf(ios::left);
    cout.width(15);
    cout << "| � �������� " << "  | ";
    cout.width(5);
    cout << "x" << "| ";
    cout.width(5);
    cout << "y" << "| ";
    cout.width(5);
    cout << "z" << "|" << endl;

    cout << "|----------------|------|------|------|" << endl;

    for (int i = 0; i < form.knots.size(); i++) {     // ��������� ���������� ����� ��� �������
        cout << "| ";
        cout.width(15);
        cout << i + 1 << "| ";
        cout.width(5);
        cout << form.knots[i].x << "| " << form.knots[i].y << "| " << form.knots[i].z;
        cout << endl;
    }

    cout << endl << "�������� ��������:" << endl;
    cout << " _____________________________________________________________________________ " << endl;
    cout.setf(ios::left);
    cout.width(15);
    cout << "| � �������� " << "  | ";
    cout.width(25);
    cout << "����" << "| ";
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
        cout << str_local_area << " ";  // ��������� ���� ��������� ���-�����
        cout << "|";
        cout.width(15);
        cout << form.locals[i].lambda << " | ";
        cout.width(15);
        cout << form.locals[i].sigma << "| ";
        cout << endl;
        str_local_area = "";

    }     
}
void write_result(InitialData form, double time) //������� ������ � �������
{
    cout << endl << "�����: " << time << endl;
    cout << endl << "��������� � ����� (����):" << endl;
    cout << " ___________________________________________________________________ " << endl;

    cout.setf(ios::left);
    cout.width(15);
    cout << "| � �������� " << "  | ";
    cout.width(15);
    cout << "u*" << "| ";
    cout.width(15);
    cout << "u" << "| ";
    cout.width(15);
    cout << "|u-u*|" << "|" << endl;
    cout << "|----------------|----------------|----------------|----------------|" << endl;

    for (int i = 0; i < form.knots.size(); i++)
    {
        double u = calc_u(form.knots[i], time);

        cout.setf(ios::left);
        cout << "| ";
        cout.width(15);
        cout << i + 1 << "| ";
        cout.width(15);
        cout << form.q[i] << "| ";
        cout.width(15);
        cout << u << "| ";
        cout.width(15);
        cout << fabs(form.q[i] - u) << "| " << endl;
    }
}

double scalar(int size, double* v1, double* v2)
{
    double sum = 0;
    for (int i = 0; i < size; i++)
        sum += v1[i] * v2[i];
    return sum;
}

void LOC(InitialData* form)
{
    int i;
    double nvzk = 0., alfa = 0., beta = 0., skp = 0., eps = 9.999999682655226e-030;
    double* z, * r, * p, * f;
    form->q = new double[form->knots.size()]{};
    z = new double[form->knots.size()]{};
    r = new double[form->knots.size()]{};
    p = new double[form->knots.size()]{};
    f = new double[form->knots.size()]{};
    double lastnvzk;

    mult_matr_by_vect(form->knots.size(), form->global_matrix, form->q, f);

    for (i = 0; i < form->knots.size(); i++)
        z[i] = r[i] = form->global_vector[i] - f[i];

    mult_matr_by_vect(form->knots.size(), form->global_matrix, z, p);
    nvzk = sqrt(scalar(form->knots.size(), r, r)) / sqrt(scalar(form->knots.size(), form->global_vector, form->global_vector));

    for (int k = 1; k < 100000 && nvzk > eps; k++)
    {
        lastnvzk = nvzk;
        skp = scalar(form->knots.size(), p, p);
        alfa = scalar(form->knots.size(), p, r) / skp;

        for (i = 0; i < form->knots.size(); i++)
        {
            form->q[i] += alfa * z[i];
            r[i] -= alfa * p[i];
        }

        mult_matr_by_vect(form->knots.size(), form->global_matrix, r, f);
        beta = -scalar(form->knots.size(), p, f) / skp;

        for (i = 0; i < form->knots.size(); i++)
        {
            z[i] = r[i] + beta * z[i];
            p[i] = f[i] + beta * p[i];
        }

        nvzk = sqrt(scalar(form->knots.size(), r, r)) / sqrt(scalar(form->knots.size(), form->global_vector, form->global_vector));
    }
}

double get_u(Knot point, InitialData form) // ��������� �������� � ������������ ���-��
{
    double x = point.x;
    double y = point.y;
    double z = point.z;

    int i;
    int igl[8]{};
    for (i = 0; i < form.num_locals; i++) // ���������� � ����� �������� �������� �����
    {
        // ���������� ���������� ���� ��� ���������� ��������
        for (int j = 0; j < 8; j++)
            igl[j] = form.locals[i].globalNum[j] - 1;

        // ���� ��������� ����� ������ � �������, �� �� �����, � ����� �������� ����� �����
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

void calc_global_matrix(InitialData* form) // ���������� ���������� �������
{
    for (int i = 0; i < form->num_locals; i++)
    {
        double M[8][8]{};
        calc_local_M(form->locals[i], M);

        double G[8][8]{};
        calc_local_G(form->locals[i], G);

        double A[8][8]{};
        for (int j = 0; j < 8; j++)
            for (int k = 0; k < 8; k++)
                form->global_matrix[form->locals[i].globalNum[j]][form->locals[i].globalNum[k]] += A[j][k] = M[j][k] + G[j][k];
    }
}

void calc_global_M(InitialData* form) // ���������� ���������� M
{
    for (int i = 0; i < form->num_locals; i++)
    {
        double M[8][8]{};
        calc_local_M(form->locals[i], M);

        for (int j = 0; j < 8; j++)
            for (int k = 0; k < 8; k++)
                form->global_M[form->locals[i].globalNum[j]][form->locals[i].globalNum[k]] += M[j][k];
    }
}
void calc_global_G(InitialData* form) // ���������� ���������� G
{
    for (int i = 0; i < form->num_locals; i++)
    {
        double G[8][8]{};
        calc_local_G(form->locals[i], G);

        for (int j = 0; j < 8; j++)
            for (int k = 0; k < 8; k++)
                form->global_G[form->locals[i].globalNum[j] - 1][form->locals[i].globalNum[k] - 1] += G[j][k];
    }
}

void calc_global_A(InitialData* form, double t_j, double t_j1, double t_j2, double t_j3) // ���������� ���������� A
{
   for (int i = 0; i < form->num_locals; i++)
   {
      for (int j = 0; j < 8; j++)
      {
         int ind_1 = form->locals[i].globalNum[j] - 1;
         for (int k = 0; k < 8; k++)
         {
            int ind_2 = form->locals[i].globalNum[k] - 1;
            form->global_matrix[ind_1][ind_2] =
               form->locals[i].hi * form->global_M[ind_1][ind_2] *
                   2 * ((t_j - t_j1) + (t_j - t_j2) + (t_j - t_j3))
                   / ((t_j - t_j3) * (t_j - t_j2) * (t_j - t_j1)) +
               form->locals[i].sigma * form->global_M[ind_1][ind_2]
                   * ((t_j - t_j2) * (t_j - t_j1) + (t_j - t_j3) * (t_j - t_j1) + (t_j - t_j3) * (t_j - t_j2))
                   / ((t_j - t_j3) * (t_j - t_j2) * (t_j - t_j1))
               + form->global_G[ind_1][ind_2];
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
void matrix_mult_vector(double** matrix, double* vector, int MAX, double* result) // ���������� �������
{
    for (int i = 0; i < MAX; i++)
        for (int j = 0; j < MAX; j++)
            result[i] += matrix[j][i] * vector[i];
}

void calc_global_d(InitialData* form, double t_j, double t_j1, double t_j2, double t_j3, double* q_3, double* q_2, double* q_1) // ���������� ���������� A
{
    double* Mq_j3, * Mq_j2, * Mq_j1;
    Mq_j3 = new double[form->knots.size()]{};
    Mq_j2 = new double[form->knots.size()]{};
    Mq_j1 = new double[form->knots.size()]{};

    matrix_mult_vector(form->global_M, q_3, form->knots.size(), Mq_j3);
    matrix_mult_vector(form->global_M, q_2, form->knots.size(), Mq_j2);
    matrix_mult_vector(form->global_M, q_1, form->knots.size(), Mq_j1);

   for (int i = 0; i < form->num_locals; i++)
   {
      for (int j = 0; j < 8; j++)
      {
         int ind_1 = form->locals[i].globalNum[j] - 1;
         form->global_vector[ind_1] = form->global_b[ind_1] 
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

void calc_global_F(InitialData* form, double time) // ���������� ����������� ������� ������ �����
{
    for (int i = 0; i < form->num_locals; i++)
    {
        double F[8]{};
        calc_local_F(form->knots, form->locals[i], F, time);

        for (int j = 0; j < 8; j++)
        {
            form->global_b[form->locals[i].globalNum[j] - 1] += F[j];
        }
    }
}

void calc_first_boundary_conditions(InitialData* form, double time) //���� ������� �������
{
    for (int i = 0; i < form->num_bounds; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            int global_num_coord = form->bounds[i].globalNum[j] - 1;
            for (int k = 0; k < form->knots.size(); k++)
            {
                form->global_matrix[global_num_coord][k] = 0;
            }
            form->global_matrix[global_num_coord][global_num_coord] = 1;
            form->global_vector[global_num_coord] = calc_u(form->knots[global_num_coord], time);
        }

    }
}


int mu(int index)
{
    return ((index - 1) % 2) + 1;
}

int v(int index)
{
    return ((index - 1) / 2) % 2 + 1;
}

int nu(int index)
{
    return (index - 1) / 4 + 1;
}

double W(int index, double alpha)
{
    switch (index)
    {
    case 1: return 1. - alpha;
    case 2: return alpha;
    }
}

double d_phi(int index, int what, double ksi, double etta, double tetha)
{
    double d_phi = 0.;
    switch (what)
    {
    case 1:    // ksi
    {
        d_phi = W(v(index), etta) * W(nu(index), tetha);
        if (mu(index) == 1) d_phi *= -1;
        break;
    }
    case 2:     // etta
    {
        d_phi = W(mu(index), ksi) * W(nu(index), tetha);
        if (v(index) == 1) d_phi *= -1;
        break;
    }
    case 3:     // tetha
    {
        d_phi = W(mu(index), ksi) * W(v(index), etta);
        if (nu(index) == 1) d_phi *= -1;
        break;
    }
    }

    return d_phi;
}

double prime_by_var(int what, double* Knot, double ksi, double etta, double tetha)
{
    double x = 0.;
    for (int i = 1; i <= 8; i++)
    {
        x += Knot[i-1] * d_phi(i, what, ksi, etta, tetha);
    }
    return x;
}

double phi(int index, double ksi, double etta, double tetha)
{
    return  W(mu(index), ksi) * W(v(index), etta) * W(nu(index), tetha);
}



double det_Jacobian(double* x, double* y, double* z, double ksi, double etta, double tetha)
{
    return prime_by_var(1, x, ksi, etta, tetha) * prime_by_var(2, y, ksi, etta, tetha) * prime_by_var(3, z, ksi, etta, tetha)
        + prime_by_var(3, x, ksi, etta, tetha) * prime_by_var(1, y, ksi, etta, tetha) * prime_by_var(2, z, ksi, etta, tetha)
        + prime_by_var(2, x, ksi, etta, tetha) * prime_by_var(3, y, ksi, etta, tetha) * prime_by_var(1, z, ksi, etta, tetha)
        - prime_by_var(3, x, ksi, etta, tetha) * prime_by_var(2, y, ksi, etta, tetha) * prime_by_var(1, z, ksi, etta, tetha)
        - prime_by_var(1, x, ksi, etta, tetha) * prime_by_var(3, y, ksi, etta, tetha) * prime_by_var(2, z, ksi, etta, tetha)
        - prime_by_var(2, x, ksi, etta, tetha) * prime_by_var(1, y, ksi, etta, tetha) * prime_by_var(3, z, ksi, etta, tetha);
}




double Integrate(function<double(double, double, double, int, int, double*, double*, double*)> f, int i, int j, double* x, double* y, double* z)
{
   const int nKnot = 5; // ���������� ����� ������

   double xj[nKnot] = { -sqrt(5. + 2. * (sqrt(10. / 7.))) / 3., -sqrt(5. - 2. * (sqrt(10. / 7.))) / 3.,	// �������� ��������� �� �����
              0. , sqrt(5. - 2. * (sqrt(10. / 7.))) / 3. , sqrt(5. + 2. * (sqrt(10. / 7.))) / 3. };

   double qj[nKnot] = { (322. - 13. * sqrt(70.)) / 900., (322. + 13. * sqrt(70.)) / 900., 128. / 225.,	// ���� �����
                  (322. + 13. * sqrt(70.)) / 900., (322. - 13. * sqrt(70.)) / 900. };
   // ����
   double hX = 1. / 2.,
      hY = 1. / 2.,
      hZ = 1. / 2.,
      // ������
      cX = 1. / 2.,
      cY = 1. / 2.,
      cZ = 1. / 2.;

   double result = 0.;
   for (int ix = 0; ix < nKnot; ix++)
      for (int iy = 0; iy < nKnot; iy++)
         for (int iz = 0; iz < nKnot; iz++)
            result += qj[ix] * qj[iy] * qj[iz] * (f(cX + xj[ix] * hX, cY + xj[iy] * hY, cZ + xj[iz] * hZ, i + 1, j + 1, x, y, z));
   return 1. * 1. * 1. * result / 8.; // ���������������
}


