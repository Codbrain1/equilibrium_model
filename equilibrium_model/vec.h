#pragma once
#include<math.h>
class vec
{
public:
	double x, y, z;
	vec() noexcept : x(0),y(0),z(0) {}

	vec(double _x, double _y, double _z) noexcept: x(_x),y(_y),z(_z) {}
	vec(const vec& temp) = default;

	vec operator*(const double temp) const noexcept
	{
		return vec(x * temp, y * temp, z * temp);
	}

	vec operator+(const vec& temp) const
	{
		return vec(x + temp.x, y + temp.y, z + temp.z);
	}
	vec operator-(const vec& temp) const
	{
		return vec(x - temp.x, y - temp.y, z - temp.z);
	}
	vec operator/(const double temp) const
	{
		return vec(x / temp, y / temp, z / temp);
	}
	vec operator*(const vec& temp)const		//vector product
	{
		return vec(y * temp.z - z * temp.y, z * temp.x - x * temp.z, x * temp.y - y * temp.x);
	}
	inline double module()const { return sqrt(x * x + y * y + z * z); }

	/// <summary>
	/// return module square
	/// </summary>
	/// <returns></returns>
	inline double module_2() const noexcept { return x * x + y * y + z * z; }

};
//vec operator*(const double c, const vec& _vec)
//{
//	return _vec * c;
//}