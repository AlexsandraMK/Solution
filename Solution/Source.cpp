#include "Header.h"



int main()
{
	setlocale(LC_ALL, "Russian");
	InitialData* form = new InitialData;   // �������� �������
	InitData(form);       // ������ �������� ������
	//calc_global_matrix(form);	// ���������� ���������� �������
	double time = form->time_g->start;

	double* q0 = new double[form->knots.size()]{};
	for (int i = 0; i < form->knots.size(); i++)
	{
		q0[i] = calc_u(form->knots[i], time);
		form->q[i] = q0[i];
	}
	write_result(*form, time);
	time += form->time_g->h0;

	double* q1 = new double[form->knots.size()]{};
	for (int i = 0; i < form->knots.size(); i++)
	{
		q1[i] = calc_u(form->knots[i], time);
		form->q[i] = q1[i];
	}
	write_result(*form, time);
	time += form->time_g->h0 * form->time_g->q;

	double* q2 = new double[form->knots.size()]{};
	for (int i = 0; i < form->knots.size(); i++)
	{
		q2[i] = calc_u(form->knots[i], time);
		form->q[i] = q2[i];
	}
	write_result(*form, time);
	calc_global_M(form);
	calc_global_G(form);

	double* q0_ = q0;
	double* q1_ = q1;
	double* q2_ = q2;

	double t0 = form->time_g->start;
	double t1 = form->time_g->start + form->time_g->h0;
	double t2 = form->time_g->h0 * form->time_g->q + t1;

	double h = form->time_g->h0 * form->time_g->q * form->time_g->q;
	
	// ����� �������� ����� �� ������ ����
	for (double time = t2 + h; time <= form->time_g->end;
		h *= form->time_g->q,
		t0 = t1, t1 = t2, t2 = time,
		q0_ = q1_,
		q1_ = q2_,
		q2_ = form->q,
		form->q = new double[form->knots.size()]{},
		time += h)
	{

		form->global_b = new double[form->knots.size()]{};
		form->global_vector = new double[form->knots.size()]{};
		calc_global_A(form, time, t2, t1, t0);
		calc_global_F(form, time);	// ���������� ����������� ������� ������ �����
		calc_global_d(form, time, t2, t1, t0, q0_, q1_, q2_);
		calc_first_boundary_conditions(form, time);	// ���� ������ ������� �������
		LOC(form);	// ������� ����
		write_result(*form, time);
	}

	return 0;
}
