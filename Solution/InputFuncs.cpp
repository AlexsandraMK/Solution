#include "InputFuncs.h"
#include <cmath>

/// <summary>
/// ��������� f
/// </summary>
/// <param name="coord">���������� �����</param>
/// <param name="time">�����</param>
double GetF(Knot coord, double time)
{
    //return 0.; // �� ������� �� ������������ � �������
    //return 0; // �� ������� �� �������
    //return 1.; // �� ������� �� ������������
    //return 1.; // ������� � �� ������������ � �� �������
    //return 2. + 2. * time; // ������� �� ������� �����������
    //return 6 * time + 3. * time * time; // ������� �� ������� ���������
    //return 12. * time * time + 4. * time * time * time; // ������ ������ �����������
    //return 2. /*+ coord.x * coord.x*/;
    //return 6. * coord.x;
    //return 12. * coord.x * coord.x;
    //return 30. * coord.x * coord.x * coord.x * coord.x;

    return 0;
}


/// <summary>
/// ��������� u
/// </summary>
/// <param name="coord">���������� �����</param>
/// <param name="time">�����</param>

double GetU(Knot coord, double time)
{
    //return 1.; // �� ������� �� ������������ � �������
    //return coord.x; // �� ������� �� �������
    //return time; // �� ������� �� ������������
    //return coord.x + time; // ������� � �� ������������ � �� �������
    //return time * time; // ������� �� ������� �����������
    //return time * time * time; // ������� �� ������� ���������
    //return time * time * time * time; // ������ ������ �����������
    //return coord.x * coord.x;
    //return coord.x * coord.x* coord.x;
    //return coord.x * coord.x * coord.x * coord.x;
    //return coord.x * coord.x * coord.x * coord.x * coord.x * coord.x;
    if (time == 4 && 
        coord.x == 0.857143 &&
        coord.y == 0. && 	
        coord.z == 2.5) return sin(time + 2 * (coord.x + coord.y));
    return 0;
}
