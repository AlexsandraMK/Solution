#include "InitialData.h"
#include <fstream>
#include <iostream>

/// <summary>
/// Инициализация задачи (пространственной и временной сеток)
/// </summary>
/// <param name="pathFile">Имя файла</param>
/// <return name="timeGrid">Сетка по времени</return>
InitialData::InitialData()
{
    ReadCoordinates("cross.txt");
    ReadCoefficients("coeffs.txt");
    //ReadKEs("TriangularPrism.txt");
    ReadKEs("Hexagon.txt");

    ReadBounds("FirstBounds.txt");

    timeGrid = new TimeGrid();
    ReadTime("grid_time.txt");
}

/// <summary>
/// Чтение коэффициентов
/// </summary>
/// <param name="pathFile">Имя файла</param>
/// <return name="knots">Координты узлов</return>
void InitialData::ReadCoefficients(string pathFile)
{
    ifstream in(pathFile);
    int nCoeffs;
    in >> nCoeffs;
    coeffs.resize(nCoeffs);
    for (int i = 0; i < nCoeffs; i++)  in >> coeffs[i].density >> coeffs[i].Yung >> coeffs[i].Poisson;
    in.close();
}


/// <summary>
/// Чтение координат узлов
/// </summary>
/// <param name="pathFile">Имя файла</param>
/// <return name="knots">Координты узлов</return>
void InitialData::ReadCoordinates(string pathFile)
{
    ifstream in(pathFile);
    int nKnots;
    in >> nKnots;
    knots.resize(nKnots);
    for (int i = 0; i < nKnots; i++)  in >> knots[i].x >> knots[i].y >> knots[i].z;
    in.close();
}

/// <summary>
/// Чтение конечных элементов
/// </summary>
/// <param name="pathFile">Имя файла</param>
/// <param name="knots">Координаты узлов</param>
/// <return name="KEs">Конечные элементы</return>
void InitialData::ReadKEs(string pathFile) // Узлы                
{
    int indFirstKE = KEs.size();

    ifstream in(pathFile);
    int nKEsInFile;
    in >> nKEsInFile;     // Считываем количество конечных элементов
    KEs.resize(nKEsInFile + KEs.size()); // Инициализируем размерность массива, который хранит локальные области

    IKE* ke;
    pathFile.erase(pathFile.find_last_of('.'));



    for (int iKE = indFirstKE; iKE < KEs.size(); iKE++)
    {
        if (pathFile == TriangularPrism::ToString())
            ke = new TriangularPrism();
        else
            if (pathFile == Hexagon::ToString())
                ke = new Hexagon();
            else
            {
                cout << "Error";
                exit(0);
            }
        int globalNum;
        int numKnots = ke->GetCountKnots();
        for (int j = 0; j < numKnots; j++)
        {
            in >> globalNum;
            ke->SetGlobalKnotNum(globalNum, knots[globalNum]);
        }

        //in >> ke->lambda;    // лямбда
        //in >> ke->sigma;     // гамма
        //in >> ke->hi;        // хи
        in >> ke->iCoeff;


        KEs[iKE] = ke;
    }
    in.close();
}


/// <summary>
/// Чтение границ
/// </summary>
/// <param name="pathFile">Имя файла</param>
/// <return name="bounds">Границы</return>
void InitialData::ReadBounds(string pathFile)
{
    ifstream in(pathFile);
    int nBounds;
    in >> nBounds;
    bounds.resize(nBounds);
    int boundSize = 4;
    for (int i = 0; i < nBounds; i++)
    {
        for (int j = 0; j < boundSize; j++)
            in >> bounds[i].globalNum[j];
    }
    in.close();
}

/// <summary>
/// Чтение сетки по времени
/// </summary>
/// <param name="pathFile">Имя файла</param>
/// <return name="timeGrid">Сетка по времени</return>
void InitialData::ReadTime(string pathFile)
{
    ifstream in(pathFile);
    in >> timeGrid->start >> timeGrid->startAfterTime3 >> timeGrid->end >> timeGrid->kAfterTime3 >> timeGrid->nStepsAfterTime3;
    in.close();
}


