#include "InputFuncs.h"
#include <cmath>

/// <summary>
/// ��������� f
/// </summary>
/// <param name="coord">���������� �����</param>
/// <param name="time">�����</param>
double GetF(Knot coord, double time)
{
    // �� ������� ������!!!!
                

    // ���������������          
    //return 0.; // �� ������� �� ������������ � �������                                     //u:   //return 1.; // �� ������� �� ������������ � �������
    //return 0; // �� ������� �� �������                                                            //return coord.x; // �� ������� �� �������
    //return 1.; // �� ������� �� ������������                                                      //return time; // �� ������� �� ������������
    //return 1.; // ������� � �� ������������ � �� �������                                          //return coord.x + time; // ������� � �� ������������ � �� �������
    //return 2.; // ������� �� ������� �����������                                      //return time * time; // ������� �� ������� �����������
    //return 6 * time; // ������� �� ������� ���������                           //return time * time * time; // ������� �� ������� ���������
    //return 12. * time * time; // ������ ������ �����������              //return time * time * time * time; // ������ ������ �����������
    //return -2.;                                                                                 //return coord.x * coord.x;
    //return -6. * coord.x;                                                                         //return coord.x * coord.x* coord.x;
    //return -12. * coord.x * coord.x;                                                              //return coord.x * coord.x * coord.x * coord.x;
    
                                                                                                    
    // ������������� ������ 
    //return 1.;
    //return coord.x;
    //return -2. + coord.x * coord.x;                                                               //return coord.x * coord.x;
    //return -6. * coord.x + coord.x * coord.x * coord.x;   
    //return -12. * coord.x* coord.x + coord.x * coord.x * coord.x* coord.x;   
    //return -2. + coord.y * coord.y;                                                               //return coord.y * coord.y;
    //return -2. + coord.z * coord.z;                                                               //return coord.z * coord.z;
    //return -6. * coord.z + coord.z * coord.z * coord.z;                                             //return coord.z * coord.z * coord.z;
    //return -12. * coord.z * coord.z + coord.z * coord.z * coord.z * coord.z;                        //return coord.z * coord.z * coord.z;

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
    //return coord.x * coord.x * coord.x;
    //return coord.x * coord.x * coord.x * coord.x;
    //return coord.y * coord.y;
    //return coord.z * coord.z;
    //return coord.z * coord.z * coord.z;
    //return coord.z * coord.z * coord.z * coord.z;
    //return coord.x * coord.x* coord.x;
    //return coord.x * coord.x * coord.x * coord.x;
    //return coord.x * coord.x * coord.x * coord.x * coord.x * coord.x;
    if (time == 1.25e-06 &&
        coord.x == 0 && //0.005
        coord.y == 0 &&
        coord.z == 2.5) return 10/**std::sin(2*3.14 / 5 *time)*/;
    return 0;
}
