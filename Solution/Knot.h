#pragma once


struct Knot // Координаты узлов
{
	double x, y, z;

	Knot(double x_, double y_, double z_)
	{
		x = x_;
		y = y_;
		z = z_;
	}

	Knot()
	{
		x = 0;
		y = 0;
		z = 0;
	}
};