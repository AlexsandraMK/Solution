#include "InputFuncs.h"
#include <cmath>

/// <summary>
/// Вычисляет f
/// </summary>
/// <param name="coord">Координата точки</param>
/// <param name="time">Время</param>
double GetF(Knot coord, double time)
{
    //return 0.; // Не зависит от пространства и времени
    //return 0; // Не зависит от времени
    //return 1.; // Не зависит от пространства
    //return 1.; // Зависит и от пространства и от времени
    //return 2. + 2. * time; // Зависит от времени квадратично
    //return 6 * time + 3. * time * time; // Зависит от времени кубически
    //return 12. * time * time + 4. * time * time * time; // Должна давать погрешность
    //return 2. /*+ coord.x * coord.x*/;
    //return 6. * coord.x;
    //return 12. * coord.x * coord.x;
    //return 30. * coord.x * coord.x * coord.x * coord.x;

    return 0;
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
    //return time * time * time * time; // Должна давать погрешность
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
