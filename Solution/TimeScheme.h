#pragma once

#include "InitialData.h"

class TimeScheme
{
public:
    vector<double> time; 
    vector<vector<double>> q;
    vector<vector<double>> qx;
    vector<vector<double>> qy;
    vector<vector<double>> qz;
    TimeScheme(InitialData* data);
    void Next();

protected:
    double k;
    double h;
};

