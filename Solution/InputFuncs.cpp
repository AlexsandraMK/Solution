#include "InputFuncs.h"
#include <cmath>

/// <summary>
/// Вычисляет f
/// </summary>
/// <param name="coord">Координата точки</param>
/// <param name="time">Время</param>
double GetF(Knot coord, double time)
{
    // НЕ ЗАБЫВАЙ МИНУСЫ!!!!
                

    // Гиперболическая          
    //return 0.; // Не зависит от пространства и времени                                     //u:   //return 1.; // Не зависит от пространства и времени
    //return 0; // Не зависит от времени                                                            //return coord.x; // Не зависит от времени
    //return 1.; // Не зависит от пространства                                                      //return time; // Не зависит от пространства
    //return 1.; // Зависит и от пространства и от времени                                          //return coord.x + time; // Зависит и от пространства и от времени
    //return 2. + 2.*time; // Зависит от времени квадратично                                      //return time * time; // Зависит от времени квадратично
    //return 6 * time + 3.*time*time; // Зависит от времени кубически                           //return time * time * time; // Зависит от времени кубически
    return 12. * time * time + 4*time*time*time; // Должна давать погрешность              //return time * time * time * time; // Должна давать погрешность
    //return -2.;                                                                                 //return coord.x * coord.x;
    //return -6. * coord.x;                                                                         //return coord.x * coord.x* coord.x;
    //return -12. * coord.x * coord.x;                                                              //return coord.x * coord.x * coord.x * coord.x;
    
                                                                                                    
    // Эллиптическая задача 
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

double GetFx(Knot coord, double time)
{
    return 0.;

    return 1.;

    //return 0.;
    //return 2.;
    //return 6. * time;
    return 12. * time * time;
}

double GetFy(Knot coord, double time)
{
    return 0.;

    return 1.;

    //return 0.;
    //return 2.;
    //return 6. * time;
    return 12. * time * time; 
}

double GetFz(Knot coord, double time)
{
    return 0.;

    return 1.;

    //return 0.;
    //return 2.;
    //return 6. * time;
    return 12. * time * time;
}



/// <summary>
/// Вычисляет u
/// </summary>
/// <param name="coord">Координата точки</param>
/// <param name="time">Время</param>
double GetU(Knot coord, double time)
{
    //return 1.; // Не зависит от пространства и времени
    //return coord.x; // Не зависит от времени
    //return time; // Не зависит от пространства
    //return coord.x + time; // Зависит и от пространства и от времени
    //return time * time; // Зависит от времени квадратично
    //return time * time * time; // Зависит от времени кубически
    return time * time * time * time; // Должна давать погрешность
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

double GetUx(Knot coord, double time)
{
    Knot* center = new Knot(0., 0., 2.5);


    if (time == 1. && (
        coord.x == -0.28 && coord.y == 3.42901e-17 && coord.z == 2.5)
        )
        return 1.5*10e+14;
    else return 0.;

    //return 1.;
    //return coord.x + coord.y + coord.z;

    return coord.x * coord.x + coord.y * coord.y + coord.z * coord.z;

    //return time ;
    //return time * time;
    //return time * time * time;
    return time * time * time * time;
    
};
double GetUy(Knot coord, double time)
{
    return 0.;
    //return 1.;
    //return coord.x + coord.y+ coord.z;

    return coord.x*coord.x + coord.y*coord.y + coord.z*coord.z;

    //return time;
    //return time * time;
    //return time * time * time;
    return time * time * time * time;
};
double GetUz(Knot coord, double time) 
{
    return 0.;
    //return 1.;
    //return coord.x + coord.y + coord.z;

    return coord.x * coord.x + coord.y * coord.y + coord.z * coord.z;

    //return time;
    //return time * time;
    //return time * time * time;
    return time * time * time * time;
};


