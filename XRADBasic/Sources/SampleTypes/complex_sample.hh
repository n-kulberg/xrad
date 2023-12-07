/*
	Copyright (c) 2021, Moscow Center for Diagnostics & Telemedicine
	All rights reserved.
	This file is licensed under BSD-3-Clause license. See LICENSE file for details.
*/

namespace xrad
{
	//--------------------------------------------------------------
	//
	//	Private division algorithm
	//
	template<class PT, class ST>
	template<class PT1, class ST1>
	complex_sample<PT, ST> complex_sample<PT, ST>::complex_division_algorithm(const complex_sample<PT, ST> x, const complex_sample<PT1, ST1> y) const
	{
		self result;
		part_type epsilon, denominator;

		if (norma(y.im) < norma(y.re))
		{
			epsilon = y.im / y.re;
			denominator = y.re + y.im * epsilon;
			result.re = (x.re + x.im * epsilon) / denominator;
			result.im = (x.im - x.re * epsilon) / denominator;
		}
		else
		{
			epsilon = y.re / y.im;
			denominator = y.im + y.re * epsilon;
			result.re = (x.im + x.re * epsilon) / denominator;
			result.im = (-x.re + x.im * epsilon) / denominator;
		}

		return result;
	}


	//--------------------------------------------------------------
	//
	//	Arithmetics (member functions)
	//

	template<class PT, class ST>
	template<class PT1, class ST1>
	inline complex_sample<PT, ST>& complex_sample<PT, ST>::operator += (const complex_sample<PT1, ST1>& y)
	{
		re += static_cast<PT> (y.re);
		im += static_cast<PT> (y.im);
	
		return (*this);
	}


	template<class PT, class ST>
	template<class PT1, class ST1>
	inline complex_sample<PT, ST>& complex_sample<PT, ST>::operator -= (const complex_sample<PT1, ST1>& y)
	{
		re -= static_cast<PT> (y.re);
		im -= static_cast<PT> (y.im);
		
		return (*this);
	}


	template<class PT, class ST>
	template<class PT1, class ST1>
	inline complex_sample<PT, ST>& complex_sample<PT, ST>::operator *= (const complex_sample<PT1, ST1>& c)
	{
		part_type new_im = re * static_cast<PT> (c.im) + im * static_cast<PT> (c.re);
		re = re * static_cast<PT> (c.re) - im * static_cast<PT> (c.im);
		im = new_im;
		
		return (*this);
	}


	template<class PT, class ST>
	template<class PT1, class ST1>
	inline complex_sample<PT, ST>& complex_sample<PT, ST>::operator %= (const complex_sample<PT1, ST1>& y)
	{
		part_type new_im = im * static_cast<PT> (y.re) - re * static_cast<PT> (y.im);
		re = re * static_cast<PT> (y.re) + im * static_cast<PT> (y.im);
		im = new_im;
		
		return (*this);
	}

	template<class PT, class ST>
	template<class PT1, class ST1>
	inline	complex_sample<PT, ST>& complex_sample<PT, ST>::operator /= (const complex_sample<PT1, ST1>& z)
	{
		self buffer((*this));
		return (*this) = complex_division_algorithm(buffer, z);
	}


	//--------------------------------------------------------------
	//
	//	Access
	//

	template<class PT, class ST>
	inline PT& real(complex_sample<PT, ST>& x) 
	{
		return x.re;
	}

	template<class PT, class ST>
	inline const PT& real(const complex_sample<PT, ST>& x) 
	{
		return x.re;
	}

	template<class PT, class ST>
	inline PT& imag(complex_sample<PT, ST>& x) 
	{
		return x.im;
	}

	template<class PT, class ST>
	inline const PT& imag(const complex_sample<PT, ST>& x) 
	{
		return x.im;
	}


	//--------------------------------------------------------------
	//
	//	"Norma" (absolute value) and "argument" (phase)
	//
	//	return type is always double (max available precision)
	//

	template<class PT, class ST>
	inline double cabs(const complex_sample<PT, ST>& x)
	{
		return norma(x);
	}

	template<class PT, class ST>
	inline double cabs2(const complex_sample<PT, ST>& x)
	{
		return quadratic_norma(x);
	}

	template<class PT, class ST>
	inline double complex_to_decibel(const complex_sample<PT, ST>& x)
	{
		return power_to_decibel(quadratic_norma(x));
	}

	template<class PT, class ST>
	inline double arg(const complex_sample<PT, ST>& x)
	{
		if (x.re == 0 && x.im == 0)
			return 0;
		else
			return std::atan2(double(x.im), double(x.re));
	}

	template<class T1, class T2>
	inline complex_sample<double, double> polar(T1 mag, T2 arg)
	{
		return complex_sample<double, double>(mag * cos(double(arg)), mag * sin(double(arg)));
	}

	template<class T1, class T2>
	inline complex_sample<double, double> polard(T1 mag, T2 arg)
	{
		return complex_sample<double, double>(mag * cos(double(arg) / degrees_per_radian()), mag * sin(double(arg) / degrees_per_radian()));
	}

	template<class T>
	inline complex_sample<double, double> phasor(T arg)
	{
		return complex_sample<double, double>(cos(double(arg)), sin(double(arg)));
	}

	template<class T>
	inline complex_sample<double, double> phasord(T arg)
	{
		return complex_sample<double, double>(cos(double(arg) / degrees_per_radian()), sin(double(arg) / degrees_per_radian()));
	}


	//--------------------------------------------------------------
	//
	// Arithmetic (not member)
	//
	template<class PT, class ST>
	inline complex_sample<PT, ST> operator + (double x, const complex_sample<PT, ST>& y)
	{
		return y + static_cast<PT> (x);
	}

	template<class PT, class ST>
	inline complex_sample<PT, ST> operator - (double x, const complex_sample<PT, ST>& y)
	{
		return complex_sample<PT, ST>(static_cast<PT> (x) - y.re, y.im);
	}

	// In all multiply/divide-by-double operations by  do not cast 'x' to 'PT', always multiply by double
	template<class PT, class ST>
	inline complex_sample<PT, ST> operator * (double x, const complex_sample<PT, ST>& y)
	{
		return y * x; 
	}

	template<class PT, class ST>
	inline complex_sample<PT, ST> operator % (double x, const complex_sample<PT, ST>& y)
	{
		return complex_sample<PT, ST>(x * y.re, -x * y.im);
	}

	template<class PT, class ST>
	inline complex_sample<PT, ST> operator / (double x, const complex_sample<PT, ST>& y)
	{
		return complex_sample<double, ST>(x) / y;
	}


	//--------------------------------------------------------------

	template<class PT, class ST>
	template<class PT1, class ST1, class PT2, class ST2>
	complex_sample<PT, ST>& complex_sample<PT, ST>::add(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y)
	{
		re = static_cast<PT> (x.re) + static_cast<PT> (y.re);
		im = static_cast<PT> (x.im) + static_cast<PT>  (y.im);
		return *this;
	}

	template<class PT, class ST>
	template<class PT1, class ST1, class PT2, class ST2>
	complex_sample<PT, ST>& complex_sample<PT, ST>::add_i(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y)
	{
		// result = x+iy
		re = static_cast<PT> (x.re) - static_cast<PT> (y.im);
		im = static_cast<PT> (x.im) + static_cast<PT> (y.re);
		return *this;
	}

	template<class PT, class ST>
	template<class PT1, class ST1, class PT2, class ST2>
	complex_sample<PT, ST>& complex_sample<PT, ST>::subtract(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y)
	{
		re = static_cast<PT> (x.re) - static_cast<PT> (y.re);
		im = static_cast<PT> (x.im) - static_cast<PT> (y.im);
		return *this;
	}

	template<class PT, class ST>
	template<class PT1, class ST1, class PT2, class ST2>
	complex_sample<PT, ST>& complex_sample<PT, ST>::subtract_i(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y)
	{
		// result = x-iy
		re = static_cast<PT> (x.re) + static_cast<PT> (y.im);
		im = static_cast<PT> (x.im) - static_cast<PT> (y.re);
		return *this;
	}

	template<class PT, class ST>
	template<class PT1, class ST1, class PT2, class ST2>
	complex_sample<PT, ST>& complex_sample<PT, ST>::multiply(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y)
	{
		im = static_cast<PT> (x.re * y.im) + static_cast<PT> (x.im * y.re);
		re = static_cast<PT> (x.re * y.re) - static_cast<PT> (x.im * y.im);
		return (*this);
	}

	template<class PT, class ST>
	template<class PT1, class ST1, class PT2, class ST2>
	complex_sample<PT, ST>& complex_sample<PT, ST>::multiply_conj(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y)
	{
		im = static_cast<PT> (x.im * y.re) - static_cast<PT> (x.re * y.im);
		re = static_cast<PT> (x.re * y.re) + static_cast<PT> (x.im * y.im);
		return (*this);
	}

	template<class PT, class ST>
	template<class PT1, class ST1, class PT2, class ST2>
	complex_sample<PT, ST>& complex_sample<PT, ST>::divide(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y)
	{
		return *this = complex_division_algorithm(x, y);
	}

	template<class PT, class ST>
	template<class PT1, class ST1, class PT2, class ST2>
	complex_sample<PT, ST>& complex_sample<PT, ST>::divide_conj(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y)
	{
		return *this = complex_division_algorithm(x, ~y);
	}


	//--------------------------------------------------------------
	// Actions like (*this)+=x*y with complex numbers
	template<class PT, class ST>
	template<class PT1, class ST1, class PT2, class ST2>
	complex_sample<PT, ST>& complex_sample<PT, ST>::add_multiply(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y)
	{
		im += static_cast<PT> (x.re * y.im) + static_cast<PT> (x.im * y.re);
		re += static_cast<PT> (x.re * y.re) - static_cast<PT> (x.im * y.im);
		return (*this);
	}

	template<class PT, class ST>
	template<class PT1, class ST1, class PT2, class ST2>
	complex_sample<PT, ST>& complex_sample<PT, ST>::add_multiply_conj(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y)
	{
		im += static_cast<PT> (x.im * y.re) - static_cast<PT> (x.re * y.im);
		re += static_cast<PT> (x.re * y.re) + static_cast<PT> (x.im * y.im);
		return (*this);
	}

	template<class PT, class ST>
	template<class PT1, class ST1, class PT2, class ST2>
	complex_sample<PT, ST>& complex_sample<PT, ST>::add_divide(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y)
	{
		return *this += complex_division_algorithm(x, y);
	}

	template<class PT, class ST>
	template<class PT1, class ST1, class PT2, class ST2>
	complex_sample<PT, ST>& complex_sample<PT, ST>::add_divide_conj(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y)
	{
		return *this += complex_division_algorithm(x, ~y);
	}

	template<class PT, class ST>
	template<class PT1, class ST1, class PT2, class ST2>
	complex_sample<PT, ST>& complex_sample<PT, ST>::subtract_multiply(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y)
	{
		im -= static_cast<PT> (x.re * y.im) + static_cast<PT> (x.im * y.re);
		re -= static_cast<PT> (x.re * y.re) - static_cast<PT> (x.im * y.im);
		return (*this);
	}

	template<class PT, class ST>
	template<class PT1, class ST1, class PT2, class ST2>
	complex_sample<PT, ST>& complex_sample<PT, ST>::subtract_multiply_conj(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y)
	{
		im -= static_cast<PT> (x.im * y.re) - static_cast<PT> (x.re * y.im);
		re -= static_cast<PT> (x.re * y.re) + static_cast<PT> (x.im * y.im);
		return (*this);
	}

	template<class PT, class ST>
	template<class PT1, class ST1, class PT2, class ST2>
	complex_sample<PT, ST>& complex_sample<PT, ST>::subtract_divide(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y)
	{
		return *this -= complex_division_algorithm(x, y);
	}

	template<class PT, class ST>
	template<class PT1, class ST1, class PT2, class ST2>
	complex_sample<PT, ST>& complex_sample<PT, ST>::subtract_divide_conj(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y)
	{
		return *this -= complex_division_algorithm(x, ~y);
	}


	//--------------------------------------------------------------
	// Actions like (*this)=x*y with complex and scalar
	template<class PT, class ST>
	template<class PT1, class ST1>
	complex_sample<PT, ST>& complex_sample<PT, ST>::multiply(const complex_sample<PT1, ST1>& x, const scalar_type& a)
	{
		re = static_cast<PT> (x.re * a);
		im = static_cast<PT> (x.im * a);

		return *this;
	}

	template<class PT, class ST>
	template<class PT1, class ST1>
	complex_sample<PT, ST>& complex_sample<PT, ST>::divide(const complex_sample<PT1, ST1>& x, const scalar_type& a) 
	{
		re = static_cast<PT> (x.re / a); 
		im = static_cast<PT> (x.im / a);

		return *this;
	}


	//--------------------------------------------------------------
	// Actions like (*this)+=x*y with complex and scalar
	template<class PT, class ST>
	template<class PT1, class ST1>
	complex_sample<PT, ST>& complex_sample<PT, ST>::add_multiply(const complex_sample<PT1, ST1>& x, const scalar_type& a)
	{
		re += static_cast<PT> (x.re * a); 
		im += static_cast<PT> (x.im * a); 
		
		return *this;
	}

	template<class PT, class ST>
	template<class PT1, class ST1>
	complex_sample<PT, ST>& complex_sample<PT, ST>::add_divide(const complex_sample<PT1, ST1>& x, const scalar_type& a) 
	{
		re += static_cast<PT> (x.re / a); 
		im += static_cast<PT> (x.im / a); 
		
		return *this;
	}

	template<class PT, class ST>
	template<class PT1, class ST1>
	complex_sample<PT, ST>& complex_sample<PT, ST>::subtract_multiply(const complex_sample<PT1, ST1>& x, const scalar_type& a) {
		re -= static_cast<PT> (x.re * a); 
		im -= static_cast<PT> (x.im * a); 
		return *this;
	}

	template<class PT, class ST>
	template<class PT1, class ST1>
	complex_sample<PT, ST>& complex_sample<PT, ST>::subtract_divide(const complex_sample<PT1, ST1>& x, const scalar_type& a) {
		re -= static_cast<PT> (x.re / a); 
		im -= static_cast<PT> (x.im / a); 
		return *this;
	}


	//--------------------------------------------------------------
	// Weighted sum (*this) = x*a1 + y*a2
	template<class PT, class ST>
	template<class PT1, class ST1, class PT2, class ST2>
	complex_sample<PT, ST>& complex_sample<PT, ST>::mix(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y, scalar_type a1, scalar_type a2)
	{
		re = static_cast<PT> (x.re * a1 + y.re * a2); 
		im = static_cast<PT> (x.im * a1 + y.im * a2);

		return *this;
	}

} // namespace xrad