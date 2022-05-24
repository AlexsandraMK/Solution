#pragma once
#include "Knot.h"

double GetF(Knot coord, double time);
double GetFx(Knot coord, double time);
double GetFy(Knot coord, double time);
double GetFz(Knot coord, double time);


double GetU(Knot coord, double time);
double GetUx(Knot coord, double time, double timeToGo);
double GetUy(Knot coord, double time, double timeToGo);
double GetUz(Knot coord, double time, double timeToGo);

enum axis
{
    x, y, z
};
