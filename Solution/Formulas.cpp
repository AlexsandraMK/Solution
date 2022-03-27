



//void write(InitialData form, double time) //функция вывода в консоль
//{
//    cout << "Параметры узлов и первые краевые :" << endl;
//    cout << " _____________________________________ " << endl;
//    cout.setf(ios::left);
//    cout.width(15);
//    cout << "| № элемента " << "  | ";
//    cout.width(5);
//    cout << "x" << "| ";
//    cout.width(5);
//    cout << "y" << "| ";
//    cout.width(5);
//    cout << "z" << "|" << endl;
//
//    cout << "|----------------|------|------|------|" << endl;
//
//    for (int i = 0; i < form.knots.size(); i++) {     // Заполняем координаты узлов для области
//        cout << "| ";
//        cout.width(15);
//        cout << i + 1 << "| ";
//        cout.width(5);
//        cout << form.knots[i].x << "| " << form.knots[i].y << "| " << form.knots[i].z;
//        cout << endl;
//    }
//
//    cout << endl << "Конечные элементы:" << endl;
//    cout << " _____________________________________________________________________________ " << endl;
//    cout.setf(ios::left);
//    cout.width(15);
//    cout << "| № элемента " << "  | ";
//    cout.width(25);
//    cout << "Узлы" << "| ";
//    cout.width(15);
//    cout << "lambda" << "| ";
//    cout.width(15);
//    cout << "gamma" << "|" << endl;
//    cout << "|----------------|--------------------------|----------------|----------------|" << endl;
//
//    for (int i = 0; i < form.num_locals; i++)
//    {
//        for (int j = 0; j < 8; j++) {
//            str_local_area += to_string(form.KEs[i].globalNum[j]) + " ";
//        }
//        cout << "| ";
//        cout.width(15);
//        cout << i+1 << "| ";
//        cout.width(24);
//        cout << str_local_area << " ";  // локальные узлы конечного эле-мента
//        cout << "|";
//        cout.width(15);
//        cout << form.KEs[i].lambda << " | ";
//        cout.width(15);
//        cout << form.KEs[i].sigma << "| ";
//        cout << endl;
//        str_local_area = "";
//
//    }     
//}




//
//double get_u(Knot point, InitialData form) // Получение значения в произвольной точ-ке
//{
//    double x = point.x;
//    double y = point.y;
//    double z = point.z;
//
//    int i;
//    int igl[8]{};
//    for (i = 0; i < form.num_locals; i++) // Определяем в каком конечном элементе точка
//    {
//        // Определяем глобальные узлы для локального элемента
//        for (int j = 0; j < 8; j++)
//            igl[j] = form.KEs[i].globalNum[j] - 1;
//
//        // Если указанные точки входят в область, то мы нашли, в каком элементе лежит точка
//        if (x >= form.knots[igl[0]].x && x <= form.knots[igl[1]].x &&
//            y >= form.knots[igl[0]].y && y <= form.knots[igl[2]].y &&
//            z >= form.knots[igl[0]].z && z <= form.knots[igl[4]].z)
//            break;
//    }
//
//    if (i == form.num_locals) return -1;
//
//    double hx = form.KEs[i].h_x;
//    double hy = form.KEs[i].h_y;
//    double hz = form.KEs[i].h_z;
//
//    double X1 = (form.knots[igl[1]].x - x) / hx;
//    double X2 = (x - form.knots[igl[0]].x) / hx;
//    double Y1 = (form.knots[igl[2]].y - y) / hy;
//    double Y2 = (y - form.knots[igl[0]].y) / hy;
//    double Z1 = (form.knots[igl[4]].z - z) / hz;
//    double Z2 = (z - form.knots[igl[0]].z) / hz;
//
//    double u =
//        form.q[igl[0]] * X1 * Y1 * Z1 +
//        form.q[igl[1]] * X2 * Y1 * Z1 +
//        form.q[igl[2]] * X1 * Y2 * Z1 +
//        form.q[igl[3]] * X2 * Y2 * Z1 +
//        form.q[igl[4]] * X1 * Y1 * Z2 +
//        form.q[igl[5]] * X2 * Y1 * Z2 +
//        form.q[igl[6]] * X1 * Y2 * Z2 +
//        form.q[igl[7]] * X2 * Y2 * Z2;
//
//    return u;
//}










