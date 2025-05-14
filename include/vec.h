#pragma once
#include<math.h>
class vec
{
public:
	long double x, y, z;
	vec() noexcept : x(0.0L),y(0.0L),z(0.0L) {}

	vec(long double _x, long double _y, long double _z) noexcept: x(_x),y(_y),z(_z) {}
	
	vec(const vec& temp) = default;

	vec operator*(const long double temp) const noexcept
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
	vec operator/(const long double temp) const
	{
		return vec(x / temp, y / temp, z / temp);
	}
	vec operator*(const vec& temp)const	
	{
		return vec(y * temp.z - z * temp.y, z * temp.x - x * temp.z, x * temp.y - y * temp.x);
	}
	friend vec operator*(long double scalar, const vec& v) {
		return vec(scalar * v.x, scalar * v.y, scalar * v.z);
	}
	inline long double module()const { return sqrtl(x * x + y * y + z * z); }

	inline long double module_2() const noexcept { return x * x + y * y + z * z; }

};
