#include "TimeScheme.h"
#include "InputFuncs.h"


TimeScheme::TimeScheme(InitialData* data)
{
    TimeGrid* timeGrid = data->timeGrid;

    // Подсчет начального шага
    k = timeGrid->kAfterTime3;
    h = timeGrid->end - timeGrid->startAfterTime3;
    h = (k == 1) ? h / timeGrid->nStepsAfterTime3 : h * (1. - k) / (1. - pow(k, timeGrid->nStepsAfterTime3));


    double SLAUsize = data->knots.size();
    time.resize(4);
    time[0] = timeGrid->start;
    time[1] = time[0] + (timeGrid->startAfterTime3 - timeGrid->start) / 2.;
    time[2] = timeGrid->startAfterTime3;
    time[3] = time[2] + h;

    q.resize(4);
    for (int j = 0; j < 3; j++)
    {
        double timeJ = time[j];
        q[j].resize(SLAUsize);
        for (int i = 0; i < SLAUsize; i++) q[j][i] = GetU(data->knots[i], timeJ);
    }
}

void TimeScheme::Next()
{
    for (int i = 0; i < 3; i++)
    {
        time[i] = time[i + 1];
        q[i] = q[i + 1];
    }
    h *= k;
    time[3] = time[2] + h;

    vector<double> qNext;
    qNext.resize(q[0].size());
    q[3] = qNext;
}