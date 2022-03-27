#pragma once

#include "InitialData.h"

class TimeScheme
{
public:
    vector<double> time; 
    vector<vector<double>> q;
    TimeScheme(InitialData* data);
    void Next();

protected:
    double k;
    double h;
};

