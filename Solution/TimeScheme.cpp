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

    qx.resize(4);
    for (int j = 0; j < 3; j++)
    {
        double timeJ = time[j];
        qx[j].resize(SLAUsize);
        for (int i = 0; i < SLAUsize; i++) qx[j][i] = GetUx(data->knots[i], timeJ);
    }

    qy.resize(4);
    for (int j = 0; j < 3; j++)
    {
        double timeJ = time[j];
        qy[j].resize(SLAUsize);
        for (int i = 0; i < SLAUsize; i++) qy[j][i] = GetUy(data->knots[i], timeJ);
    }

    qz.resize(4);
    for (int j = 0; j < 3; j++)
    {
        double timeJ = time[j];
        qz[j].resize(SLAUsize);
        for (int i = 0; i < SLAUsize; i++) qz[j][i] = GetUz(data->knots[i], timeJ);
    }
}

void TimeScheme::Next()
{
    for (int i = 0; i < 3; i++)
    {
        time[i] = time[i + 1];
        //q[i] = q[i + 1];
        qx[i] = qx[i + 1];
        qy[i] = qy[i + 1];
        qz[i] = qz[i + 1];
    }


    h *= k;
    time[3] = time[2] + h;

    vector<double> qNext, qxNext, qyNext, qzNext;
    //qNext.resize(q[0].size());
    //q[3] = qNext;
    qxNext.resize(qx[0].size());
    qx[3] = qxNext;
    qyNext.resize(qy[0].size());
    qy[3] = qyNext;
    qzNext.resize(qz[0].size());
    qz[3] = qzNext;
}

void TimeScheme::ChangeTimeScheme(vector<double> q)
{
    double num_knots = qx[0].size();

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < num_knots; j++)
        {
            qx[i][j] -= q[j * 3];
            qy[i][j] -= q[j * 3 + 1];
            qz[i][j] -= q[j * 3 + 2];
        }
    }

}
