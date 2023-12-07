/*
	Copyright (c) 2021, Moscow Center for Diagnostics & Telemedicine
	All rights reserved.
	This file is licensed under BSD-3-Clause license. See LICENSE file for details.
*/
XRAD_BEGIN



//--------------------------------------------------------------
//
//	private division algorithm
//

template<class PT, class ST>
template<class PT1, class ST1>
complex_sample<PT, ST> complex_sample<PT, ST>::complex_division_algorithm(const complex_sample<PT, ST> x, const complex_sample<PT1, ST1> y) const
{
	self result;
	part_type epsilon, denominator;

	if(norma(y.im)<norma(y.re))
	{
		epsilon = y.im/y.re;
		denominator = y.re + y.im*epsilon;
		result.re = (x.re + x.im*epsilon)/denominator;
		result.im = (x.im - x.re*epsilon)/denominator;
	}
	else
	{
		epsilon = y.re/y.im;
		denominator = y.im+y.re*epsilon;
		result.re = (x.im + x.re*epsilon)/denominator;
		result.im = (-x.re + x.im*epsilon)/denominator;
	}

	return result;
}



//--------------------------------------------------------------
//
//	arithmetics (member functions)
//
//--------------------------------------------------------------



template<class PT, class ST>
template<class PT1, class ST1>
inline complex_sample<PT, ST> &complex_sample<PT, ST>::operator += (const complex_sample<PT1, ST1> &y)
{
	re += y.re;
	im += y.im;
	return (*this);
}

//--------------------------------------------------------------

template<class PT, class ST>
template<class PT1, class ST1>
inline complex_sample<PT, ST> &complex_sample<PT, ST>::operator -= (const complex_sample<PT1, ST1> &y)
{
	re -= y.re;
	im -= y.im;
	return (*this);
}

//--------------------------------------------------------------

template<class PT, class ST>
template<class PT1, class ST1>
inline complex_sample<PT, ST> &complex_sample<PT, ST>::operator *= (const complex_sample<PT1, ST1> &c)
{
	part_type new_im = re*c.im + im*c.re;
	re = re*c.re - im*c.im;
	im = new_im;
	return (*this);
}

//--------------------------------------------------------------

template<class PT, class ST>
template<class PT1, class ST1>
inline complex_sample<PT, ST> &complex_sample<PT, ST>::operator %= (const complex_sample<PT1, ST1> &y)
{
	part_type new_im = im*y.re - re*y.im;
	re = re*y.re + im*y.im;
	im = new_im;
	return (*this);
}

template<class PT, class ST>
template<class PT1, class ST1>
inline	complex_sample<PT, ST> &complex_sample<PT, ST>::operator /= (const complex_sample<PT1, ST1> &z)
{
	self buffer((*this));
	return (*this) = complex_division_algorithm(buffer, z);
}



//--------------------------------------------------------------
//
//	access
//
//--------------------------------------------------------------



template<class PT, class ST>
inline PT &real(complex_sample<PT, ST> &x) {
	return x.re;
}

template<class PT, class ST>
inline const PT &real(const complex_sample<PT, ST> &x) {
	return x.re;
}

template<class PT, class ST>
inline PT &imag(complex_sample<PT, ST> &x) {
	return x.im;
}

template<class PT, class ST>
inline const PT &imag(const complex_sample<PT, ST> &x) {
	return x.im;
}



//--------------------------------------------------------------
//
//	функции "норма" и "аргумент"
//
//	возвращаемый тип всегда с максимальной точностью.
//	иначе при возврате любого из complex_sample<Tn> возможна такая ошибка:
//	polar(1,1) преобразует результат к complex_sample<int>
//	кнс. 16 мая 2011
//
//--------------------------------------------------------------



template<class PT, class ST>
inline double cabs(const complex_sample<PT, ST> &x)
{
	return norma(x);
}

template<class PT, class ST>
inline double cabs2(const complex_sample<PT, ST> &x)
{
	return quadratic_norma(x);
}

template<class PT, class ST>
inline double complex_to_decibel(const complex_sample<PT, ST> &x)
{
	return power_to_decibel(quadratic_norma(x));
}

template<class PT, class ST>
inline double arg(const complex_sample<PT, ST> &x)
{
	if(x.re == 0 && x.im == 0)
		return 0;
	else
		return std::atan2(double(x.im), double(x.re));
}

template<class T1, class T2>
inline complex_sample<double, double> polar(T1 mag, T2 arg)
{
	return complex_sample<double, double>(mag*cos(double(arg)), mag*sin(double(arg)));
}

template<class T1, class T2>
inline complex_sample<double, double> polard(T1 mag, T2 arg)
{
	return complex_sample<double, double>(mag*cos(double(arg)/degrees_per_radian()), mag*sin(double(arg)/degrees_per_radian()));
}

template<class T>
inline complex_sample<double, double> phasor(T arg)
{
	return complex_sample<double, double>(cos(double(arg)), sin(double(arg)));
}

template<class T>
inline complex_sample<double, double> phasord(T arg)
{
	return complex_sample<double, double>(cos(double(arg)/degrees_per_radian()), sin(double(arg)/degrees_per_radian()));
}



//--------------------------------------------------------------
//
// arithmetic (not member)
//
//--------------------------------------------------------------



template<class PT, class ST>
inline complex_sample<PT, ST> operator + (double x, const complex_sample<PT, ST> &y)
{
	return y+x;
}

template<class PT, class ST>
inline complex_sample<PT, ST> operator - (double x, const complex_sample<PT, ST> &y)
{
	return complex_sample<PT, ST>(x-y.re, y.im);
}

template<class PT, class ST>
inline complex_sample<PT, ST> operator * (double x, const complex_sample<PT, ST> &y)
{
	return y*x;
}

template<class PT, class ST>
inline complex_sample<PT, ST> operator % (double x, const complex_sample<PT, ST> &y)
{
	return complex_sample<PT, ST>(x*y.re, -x*y.im);
}

template<class PT, class ST>
inline complex_sample<PT, ST> operator / (double x, const complex_sample<PT, ST> &y)
{
	return complex_sample<PT, ST>(x)/y;
}



//--------------------------------------------------------------

// действия вида (*this)=x*y с комплексными числами
template<class PT, class ST>
template<class PT1, class ST1, class PT2, class ST2>
complex_sample<PT, ST>	&complex_sample<PT, ST>::add(const complex_sample<PT1, ST1> &x, const complex_sample<PT2, ST2> &y)
{
	re=x.re+y.re;
	im=x.im+y.im;
	return *this;
}

template<class PT, class ST>
template<class PT1, class ST1, class PT2, class ST2>
complex_sample<PT, ST>	&complex_sample<PT, ST>::add_i(const complex_sample<PT1, ST1> &x, const complex_sample<PT2, ST2> &y)
{
	// result = x+iy
	re=x.re-y.im;
	im=x.im+y.re;
	return *this;
}

template<class PT, class ST>
template<class PT1, class ST1, class PT2, class ST2>
complex_sample<PT, ST>	&complex_sample<PT, ST>::subtract(const complex_sample<PT1, ST1> &x, const complex_sample<PT2, ST2> &y)
{
	re=x.re-y.re;
	im=x.im-y.im;
	return *this;
}

template<class PT, class ST>
template<class PT1, class ST1, class PT2, class ST2>
complex_sample<PT, ST>	&complex_sample<PT, ST>::subtract_i(const complex_sample<PT1, ST1> &x, const complex_sample<PT2, ST2> &y)
{
	// result = x-iy
	re=x.re+y.im;
	im=x.im-y.re;
	return *this;
}

template<class PT, class ST>
template<class PT1, class ST1, class PT2, class ST2>
complex_sample<PT, ST>	&complex_sample<PT, ST>::multiply(const complex_sample<PT1, ST1> &x, const complex_sample<PT2, ST2> &y)
{
	im = x.re*y.im + x.im*y.re;
	re = x.re*y.re - x.im*y.im;
	return (*this);
}

template<class PT, class ST>
template<class PT1, class ST1, class PT2, class ST2>
complex_sample<PT, ST>	&complex_sample<PT, ST>::multiply_conj(const complex_sample<PT1, ST1> &x, const complex_sample<PT2, ST2> &y)
{
	im = x.im*y.re - x.re*y.im;
	re = x.re*y.re + x.im*y.im;
	return (*this);
}

template<class PT, class ST>
template<class PT1, class ST1, class PT2, class ST2>
complex_sample<PT, ST>	&complex_sample<PT, ST>::divide(const complex_sample<PT1, ST1> &x, const complex_sample<PT2, ST2> &y)
{
	return *this = complex_division_algorithm(x, y);
}

template<class PT, class ST>
template<class PT1, class ST1, class PT2, class ST2>
complex_sample<PT, ST>	&complex_sample<PT, ST>::divide_conj(const complex_sample<PT1, ST1> &x, const complex_sample<PT2, ST2> &y)
{
	return *this = complex_division_algorithm(x, ~y);
}



//--------------------------------------------------------------
// действия вида (*this)+=x*y с комплексными числами
template<class PT, class ST>
template<class PT1, class ST1, class PT2, class ST2>
complex_sample<PT, ST>	&complex_sample<PT, ST>::add_multiply(const complex_sample<PT1, ST1> &x, const complex_sample<PT2, ST2> &y)
{
	im += x.re*y.im + x.im*y.re;
	re += x.re*y.re - x.im*y.im;
	return (*this);
}

template<class PT, class ST>
template<class PT1, class ST1, class PT2, class ST2>
complex_sample<PT, ST>	&complex_sample<PT, ST>::add_multiply_conj(const complex_sample<PT1, ST1> &x, const complex_sample<PT2, ST2> &y)
{
	im += x.im*y.re - x.re*y.im;
	re += x.re*y.re + x.im*y.im;
	return (*this);
}

template<class PT, class ST>
template<class PT1, class ST1, class PT2, class ST2>
complex_sample<PT, ST>	&complex_sample<PT, ST>::add_divide(const complex_sample<PT1, ST1> &x, const complex_sample<PT2, ST2> &y)
{
	return *this += complex_division_algorithm(x, y);
}

template<class PT, class ST>
template<class PT1, class ST1, class PT2, class ST2>
complex_sample<PT, ST>	&complex_sample<PT, ST>::add_divide_conj(const complex_sample<PT1, ST1> &x, const complex_sample<PT2, ST2> &y)
{
	return *this += complex_division_algorithm(x, ~y);
}

template<class PT, class ST>
template<class PT1, class ST1, class PT2, class ST2>
complex_sample<PT, ST>	&complex_sample<PT, ST>::subtract_multiply(const complex_sample<PT1, ST1> &x, const complex_sample<PT2, ST2> &y)
{
	im -= x.re*y.im + x.im*y.re;
	re -= x.re*y.re - x.im*y.im;
	return (*this);
}

template<class PT, class ST>
template<class PT1, class ST1, class PT2, class ST2>
complex_sample<PT, ST>	&complex_sample<PT, ST>::subtract_multiply_conj(const complex_sample<PT1, ST1> &x, const complex_sample<PT2, ST2> &y)
{
	im -= x.im*y.re - x.re*y.im;
	re -= x.re*y.re + x.im*y.im;
	return (*this);
}

template<class PT, class ST>
template<class PT1, class ST1, class PT2, class ST2>
complex_sample<PT, ST>	&complex_sample<PT, ST>::subtract_divide(const complex_sample<PT1, ST1> &x, const complex_sample<PT2, ST2> &y)
{
	return *this -= complex_division_algorithm(x, y);
}

template<class PT, class ST>
template<class PT1, class ST1, class PT2, class ST2>
complex_sample<PT, ST>	&complex_sample<PT, ST>::subtract_divide_conj(const complex_sample<PT1, ST1> &x, const complex_sample<PT2, ST2> &y)
{
	return *this -= complex_division_algorithm(x, ~y);
}



//--------------------------------------------------------------
// действия вида (*this)=x*y с комплексным и скаляром
template<class PT, class ST>
template<class PT1, class ST1>
complex_sample<PT, ST>	&complex_sample<PT, ST>::multiply(const complex_sample<PT1, ST1> &x, const scalar_type &a){
	re=x.re*a; im=x.im*a; return *this;
}

template<class PT, class ST>
template<class PT1, class ST1>
complex_sample<PT, ST>	&complex_sample<PT, ST>::divide(const complex_sample<PT1, ST1> &x, const scalar_type &a){
	re=x.re/a; im=x.im/a; return *this;
}



//--------------------------------------------------------------
// действия вида (*this)+=x*y с комплексным и скаляром
template<class PT, class ST>
template<class PT1, class ST1>
complex_sample<PT, ST>	&complex_sample<PT, ST>::add_multiply(const complex_sample<PT1, ST1> &x, const scalar_type &a){
	re+=x.re*a; im+=x.im*a; return *this;
}

template<class PT, class ST>
template<class PT1, class ST1>
complex_sample<PT, ST>	&complex_sample<PT, ST>::add_divide(const complex_sample<PT1, ST1> &x, const scalar_type &a){
	re+=x.re/a; im+=x.im/a; return *this;
}

template<class PT, class ST>
template<class PT1, class ST1>
complex_sample<PT, ST>	&complex_sample<PT, ST>::subtract_multiply(const complex_sample<PT1, ST1> &x, const scalar_type &a){
	re-=x.re*a; im-=x.im*a; return *this;
}

template<class PT, class ST>
template<class PT1, class ST1>
complex_sample<PT, ST>	&complex_sample<PT, ST>::subtract_divide(const complex_sample<PT1, ST1> &x, const scalar_type &a){
	re-=x.re/a; im-=x.im/a; return *this;
}



//--------------------------------------------------------------
// взвешенное сложение двух массивов, (*this) = x*a1 + y*a2
template<class PT, class ST>
template<class PT1, class ST1, class PT2, class ST2>
complex_sample<PT, ST>	&complex_sample<PT, ST>::mix(const complex_sample<PT1, ST1> &x, const complex_sample<PT2, ST2> &y, scalar_type a1, scalar_type a2){
	re=x.re*a1+y.re*a2; im=x.im*a1+y.im*a2; return *this;
}



//--------------------------------------------------------------

XRAD_END
