#include "InputFuncs.h"
#include <cmath>

double GetFx(Knot coord, double time)
{
    return 0.;  // 1, x+y+z, t, solution, x^4, x^3
    //return -1.; // x^2 + y^2 + z^2
    //return -0.5; // x^2 + y^2
    //return -1.5 * coord.y - 1.5 * coord.z; // x^3 + y^3 + z^3
    //return -3. * coord.y - 3. * coord.z; // x^4 + y^4 + z^4
    //return -3. * coord.z; // z^4
    //return 2.; // t^2
    //return 6. * time;   // t^3
    //return 12. * time * time;   //t^4

    // �������� �� �������������
    //return 1;
    /*return coord.x + coord.y;*/
    return 1. + coord.x * coord.x + coord.y * coord.y ;
    //return 0;

}

double GetFy(Knot coord, double time)
{
    return 0.;  // 1, x+y+z, t, solution
    //return -1.; // x^2 + y^2 + z^2
    //return -0.5; // x^2 
    //return -1.5 * coord.x - 1.5 * coord.z; // x^3 + y^3 + z^3
    //return -1.5 * coord.x; // x^3
    //return -3. * coord.x - 3. * coord.z; // x^4 + y^4 + z^4
    //return -3. * coord.x; // x^4
    //return -3. * coord.z; // z^4
    //return 2.; // t^2
    //return 6. * time;   // t^3
    return 12. * time * time;   //t^4


    // �������� �� �������������
    //return 1;
    /*return coord.x + coord.y;*/
    //return coord.x + coord.y + coord.z;
    return 1. + coord.x * coord.x + coord.y * coord.y;
    return 0;
    return 0;
}

double GetFz(Knot coord, double time)
{
    return 0.;  // 1, x+y+z, t, solution
    //return -1.; // x^2 + y^2 + z^2, // x^2 + y^2 
    //return -0.5; // x^2 
    //return -1.5 * coord.x - 1.5 * coord.y; // x^3 + y^3 + z^3
    //return -1.5 * coord.x; // x^3
    //return -3. * coord.x - 3. * coord.y; // x^4 + y^4 + z^4
    //return -3. * coord.x; // x^4
    //return 2.; // t^2
    //return 6. * time;   // t^3
    return 12. * time * time;   //t^4


    // �������� �� �������������
    //return 1;
    /*return coord.x + coord.y;*/
    //return coord.x + coord.y + coord.z;
    return 1. + coord.x * coord.x + coord.y * coord.y;
    return 0;
    return 0;
}





double GetUx(Knot coord, double time, double timeToGo)
{
    Knot* knot_to_go = new Knot(0.0266667, 0., 2.5);
    if (time == timeToGo && (
        coord.x == knot_to_go->x && coord.y == knot_to_go->y && coord.z == knot_to_go->z)
        )
        return 0.04 - knot_to_go->x;
    else return 0.;

    //return 1.;
    //return coord.x + coord.y;
    //return coord.x + coord.y + coord.z;
    //return coord.x * coord.x;
    //return coord.x * coord.x  + coord.y * coord.y;
    //return coord.x * coord.x + coord.y * coord.y + coord.z * coord.z;
    //return coord.x * coord.x * coord.x + coord.y * coord.y * coord.y + coord.z * coord.z * coord.z;
    //return coord.x * coord.x * coord.x;
    //return coord.x * coord.x * coord.x * coord.x + coord.y * coord.y * coord.y * coord.y + coord.z * coord.z * coord.z * coord.z;
    //return coord.x * coord.x * coord.x * coord.x;
    //return coord.z * coord.z * coord.z * coord.z;
    //return coord.x * coord.x * coord.x * coord.x * coord.x;
    //return time ;
    //return time * time;
    //return time * time * time;
    //return time * time * time * time;

    return coord.x * time;
    
};
double GetUy(Knot coord, double time, double timeToGo)
{
    return 0.;
    //return 1.;
    //return coord.x + coord.y;
    //return coord.x + coord.y + coord.z;
    //return coord.x * coord.x;
    //return coord.x * coord.x + coord.y * coord.y;
    //return coord.x * coord.x + coord.y * coord.y + coord.z * coord.z;
    //return coord.x * coord.x ;
    //return coord.x * coord.x * coord.x + coord.y * coord.y * coord.y + coord.z * coord.z * coord.z;
    //return coord.x * coord.x * coord.x;
    //return coord.x * coord.x * coord.x * coord.x + coord.y * coord.y * coord.y * coord.y + coord.z * coord.z * coord.z * coord.z;
    //return coord.x * coord.x * coord.x * coord.x;
    //return coord.z * coord.z * coord.z * coord.z;
    //return coord.y * coord.y * coord.y * coord.y;
    //return coord.y * coord.y * coord.y * coord.y * coord.y;
    //return time;
    //return time * time;
    //return time * time * time;
    return time * time * time * time;
};

double GetUz(Knot coord, double time, double timeToGo)
{
    return 0.;
    //return 1.;
    //return coord.x + coord.y;
    //return coord.x + coord.y + coord.z;
    //return coord.x * coord.x;
    //return coord.x * coord.x + coord.y * coord.y;
    //return coord.x * coord.x + coord.y * coord.y + coord.z * coord.z;
    //return coord.x * coord.x ;
    //return coord.x * coord.x * coord.x + coord.y * coord.y * coord.y + coord.z * coord.z * coord.z;
    //return coord.x * coord.x * coord.x;
    //return coord.x * coord.x * coord.x * coord.x + coord.y * coord.y * coord.y * coord.y + coord.z * coord.z * coord.z * coord.z;
    //return coord.x * coord.x * coord.x * coord.x;
    //return coord.z * coord.z * coord.z * coord.z;
    //return coord.z * coord.z * coord.z * coord.z * coord.z;
    //return time;
    //return time * time;
    //return time * time * time;
    return time * time * time * time;
};


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
    //return 2. + 2.*time; // ������� �� ������� �����������                                      //return time * time; // ������� �� ������� �����������
    //return 6 * time + 3.*time*time; // ������� �� ������� ���������                           //return time * time * time; // ������� �� ������� ���������
    return 12. * time * time + 4 * time * time * time; // ������ ������ �����������              //return time * time * time * time; // ������ ������ �����������
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
    return time * time * time * time; // ������ ������ �����������
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