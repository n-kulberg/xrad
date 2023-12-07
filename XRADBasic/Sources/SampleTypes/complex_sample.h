/*
	Copyright (c) 2021, Moscow Center for Diagnostics & Telemedicine
	All rights reserved.
	This file is licensed under BSD-3-Clause license. See LICENSE file for details.
*/
// \file complex_sample.h
//
// A part of XRAD
// Complex numbers library for scientific calculations.
//
//--------------------------------------------------------------
#ifndef XRAD__File_complex_sample_h
#define XRAD__File_complex_sample_h
//--------------------------------------------------------------

#ifndef XRAD_STANDALONE_COMPLEX
#define XRAD_STANDALONE_COMPLEX 1
#endif //XRAD_STANDALONE_COMPLEX


#if !XRAD_STANDALONE_COMPLEX
	#include <XRADBasic/Sources/Core/Config.h>
	//#include <XRADBasic/Sources/Core/BasicMacros.h>
	#include <XRADBasic/Sources/Core/NumberTraits.h>
	#include "HomomorphSamples.h"
	#include <type_traits>
#else
namespace xrad 
{
	constexpr double degrees_per_radian() { return 180. / 3.14159265358979323846; }
	constexpr double radians_per_degree() { return 3.14159265358979323846 / 180.; }
	template<class T> double norma(const T &x) { return abs(x); }
}//namespace xrad
#endif //XRAD_STANDALONE_COMPLEX

namespace xrad {


	enum class complex_number_part
	{
		real,
		imag
	};


	/*!
		\brief Complex-valued number

		PT is the type of complex number, must be nonconst scalar numeric type
		ST is the type of scalar multiplier

		'operator ~' is  complex conjugation (like unary '-' operator)
		'operator %", "%=' mean multiplication by conjugated number
	*/
	template<class PT, class ST>
	class complex_sample
	{
		// prohibit use of 'complex_sample<const PT, ST>'
		static_assert(!std::is_const<PT>::value, "Error: Using complex_sample<PT, ST> with const PT type. PT must not be const. Use const complex_sample<PT, ST> instead.");
	public:
		typedef	PT part_type;
		typedef	ST scalar_type;
		typedef complex_sample<PT, ST> self;

	public:
		PT re, im;

	public:
		//	constructors

		complex_sample() {}
		// No initialization in default constructor

		explicit complex_sample(PT r) { re = r; im = 0; }
		// 'explicit' prohibits accidental appearance of complex numbers instead of real ones

		complex_sample(PT r, PT i) { re = r; im = i; }

		template<class PT1, class ST1>
		complex_sample(const complex_sample<PT1, ST1>& c) { re = c.re; im = c.im; }

		// Assignments
		template<class PT1, class ST1>
		complex_sample& operator = (const complex_sample<PT1, ST1>& y) { re = y.re; im = y.im; return *this; }

		complex_sample<PT, ST>& operator = (PT y) { re = y; im = 0; return *this; }

		// Unary arithmetics
		complex_sample<PT, ST> operator-() const { return complex_sample<PT, ST>(-re, -im); }
		complex_sample<PT, ST> operator~() const { return complex_sample<PT, ST>(re, -im); }

		complex_sample<PT, ST>& operator++() { ++re; return *this; }
		complex_sample<PT, ST>& operator--() { --re; return *this; }

		complex_sample<PT, ST> operator++(int) { complex_sample<PT, ST> result(*this); re++; return result; }
		complex_sample<PT, ST> operator--(int) { complex_sample<PT, ST> result(*this); re--; return result; }

		// Assignment arithmetic
		template<class PT1, class ST1>
		complex_sample<PT, ST>& operator += (const complex_sample<PT1, ST1>& y);

		template<class PT1, class ST1>
		complex_sample<PT, ST>& operator -= (const complex_sample<PT1, ST1>& y);

		template<class PT1, class ST1>
		complex_sample<PT, ST>& operator *= (const complex_sample<PT1, ST1>& c);

		template<class PT1, class ST1>
		complex_sample<PT, ST>& operator %= (const complex_sample<PT1, ST1>& y);

		template<class PT1, class ST1>
		complex_sample<PT, ST>& operator /= (const complex_sample<PT1, ST1>& y);

		complex_sample<PT, ST>& operator += (PT y) { re += y; return *this; }
		complex_sample<PT, ST>& operator -= (PT y) { re -= y; return *this; }
		complex_sample<PT, ST>& operator *= (ST y) { re *= y; im *= y; return *this; }
		complex_sample<PT, ST>& operator /= (ST y) { re /= y; im /= y; return *this; }

		// Arithmetic
		template<class PT1, class ST1>
		complex_sample<PT, ST> operator + (const complex_sample<PT1, ST1>& y) const { return complex_sample<PT, ST>(*this) += y; }
		template<class PT1, class ST1>
		complex_sample<PT, ST> operator - (const complex_sample<PT1, ST1>& y) const { return complex_sample<PT, ST>(*this) -= y; }
		template<class PT1, class ST1>
		complex_sample<PT, ST> operator * (const complex_sample<PT1, ST1>& y) const { return complex_sample<PT, ST>(*this) *= y; }
		template<class PT1, class ST1>
		complex_sample<PT, ST> operator % (const complex_sample<PT1, ST1>& y) const { return complex_sample<PT, ST>(*this) %= y; }
		template<class PT1, class ST1>
		complex_sample<PT, ST> operator / (const complex_sample<PT1, ST1>& y) const { return complex_sample<PT, ST>(*this) /= y; }
		//
		complex_sample<PT, ST> operator + (PT y) const { return complex_sample<PT, ST>(*this) += y; }
		complex_sample<PT, ST> operator - (PT y) const { return complex_sample<PT, ST>(*this) -= y; }
		complex_sample<PT, ST> operator * (ST y) const { return complex_sample<PT, ST>(*this) *= y; }
		complex_sample<PT, ST> operator / (ST y) const { return complex_sample<PT, ST>(*this) /= y; }




		template<class PT1, class ST1>
		complex_sample<PT, ST>& conjugate(const complex_sample<PT1, ST1>& x) { re = x.re; im = -im.re; return *this; }

		// Ternary actions, result of binary operation is placed to *this.
		// Calling 'a.add(b,c)' avoids creating buffer variables in operations like 'a=b+c'.
		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& add(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y);
		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& add_i(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y); // x + i*y

		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& subtract(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y);
		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& subtract_i(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y); // x - i*y

		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& multiply(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y);
		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& multiply_conj(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y);

		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& divide(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y);
		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& divide_conj(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y);

		// Ternary actions like (*this)+=x*y with two complex numbers
		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& add_multiply(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y);
		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& add_multiply_conj(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y);

		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& add_divide(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y);
		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& add_divide_conj(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y);

		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& subtract_multiply(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y);
		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& subtract_multiply_conj(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y);

		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& subtract_divide(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y);
		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& subtract_divide_conj(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y);

		// Ternary actions like (*this)=x*y with complex and scalar number
		template<class PT1, class ST1>
		complex_sample<PT, ST>& multiply(const complex_sample<PT1, ST1>& x, const scalar_type& a);

		template<class PT1, class ST1>
		complex_sample<PT, ST>& divide(const complex_sample<PT1, ST1>& x, const scalar_type& a);

		// ternary actions (*this)+=x*y with complex and scalar number
		template<class PT1, class ST1>
		complex_sample<PT, ST>& add_multiply(const complex_sample<PT1, ST1>& x, const scalar_type& a);

		template<class PT1, class ST1>
		complex_sample<PT, ST>& add_divide(const complex_sample<PT1, ST1>& x, const scalar_type& a);

		template<class PT1, class ST1>
		complex_sample<PT, ST>& subtract_multiply(const complex_sample<PT1, ST1>& x, const scalar_type& a);

		template<class PT1, class ST1>
		complex_sample<PT, ST>& subtract_divide(const complex_sample<PT1, ST1>& x, const scalar_type& a);

		// weighted sum, (*this) = x*a1 + y*a2
		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& mix(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y, scalar_type a1, scalar_type a2);

		// compare
		template<class PT1, class ST1>
		bool operator == (const complex_sample<PT1, ST1>& x) const { return (re == static_cast <PT> (x.re) && im == static_cast <PT> (x.im)); }
		template<class PT1, class ST1>
		bool operator != (const complex_sample<PT1, ST1>& x) const { return (re != static_cast <PT> (x.re) || im != static_cast <PT> (x.im)); }
		template<class PT1, class ST1>
		bool operator < (const complex_sample<PT1, ST1>& x) const { return (re * re + (im * im) < static_cast <PT> (x.re* x.re + x.im * x.im)); }
		template<class PT1, class ST1>
		bool operator > (const complex_sample<PT1, ST1>& x) const { return (re * re + im * im > static_cast <PT> (x.re * x.re + x.im * x.im)); }
		template<class PT1, class ST1>
		bool operator <= (const complex_sample<PT1, ST1>& x) const { return (re * re + im * im <= static_cast <PT> (x.re * x.re + x.im * x.im)); }
		template<class PT1, class ST1>
		bool operator >= (const complex_sample<PT1, ST1>& x) const { return (re * re + im * im >= static_cast <PT> (x.re * x.re + x.im * x.im)); }
		//
		bool operator == (PT x) const { return (re == x && !im); }
		bool operator != (PT x) const { return (re != x || im); }
		bool operator < (PT x) const { return ((re * re + im * im) < static_cast <PT> (norma (x*x))); }
		bool operator > (PT x) const { return ((re * re + im * im) > static_cast <PT> (norma (x*x))); }
		bool operator <= (PT x) const { return ((re * re + im * im) <= static_cast <PT> (norma (x*x))); }
		bool operator >= (PT x) const { return ((re * re + im * im) >= static_cast <PT> (norma (x*x))); }


		//	Distance in memory between components of the complex numbers (in samples)
		//	In most cases it is equal to 1, but fields alignment inside the structure may result in greater value
		//	Used for separate extraction of real/imaginary parts in arrays.
		static ptrdiff_t parts_distance()
		{
			enum
			{
				self_size = sizeof(self),
				part_size = sizeof(part_type),
				distance = self_size / part_size
			};
			static_assert(distance * part_size == self_size, "complex_sample layout problem.");

			return distance;
		}
	private:
		template<class PT1, class ST1>
		complex_sample<PT, ST> complex_division_algorithm(const complex_sample<PT, ST> x, const complex_sample<PT1, ST1> y) const;
	};



	//--------------------------------------------------------------
	//
	// Number traits definitions. see comment in xrad/number_traits.h (ignore when using standalone)
	//
	template<class PT, class ST>
	size_t	n_components(const complex_sample<PT, ST>&) { return 2; }

	template<class PT, class ST>
	PT& component(complex_sample<PT, ST>& x, size_t n)
	{
		if (!n) return x.re;
		return x.im;
	}

	template<class PT, class ST>
	const PT& component(const complex_sample<PT, ST>& x, size_t n)
	{
		if (!n) return x.re;
		return x.im;
	}

	template<class PT, class ST>
	double	norma(const complex_sample<PT, ST>& x) { return std::hypot(x.re, x.im); }

	template<class PT, class ST>
	double	fast_norma(const complex_sample<PT, ST>& x) { return (x.re * x.re + x.im * x.im); }

	template<class PT, class ST>
	double	quadratic_norma(const complex_sample<PT, ST>& x) { return (x.re * x.re + x.im * x.im); }

	template<class PT, class ST>
	complex_sample<PT, ST> zero_value(const complex_sample<PT, ST>&) { return complex_sample<PT, ST>(0); }\

		template<class PT, class ST>
	void make_zero(complex_sample<PT, ST>& value) { value = complex_sample<PT, ST>(0); }


	//--------------------------------------------------------------
	//
	//	conjugated multiplication used inside the scalar product operator for complex-valued arrays
	//
	template<class PT, class PT1, class PT2, class ST, class ST1, class ST2>
	void	scalar_product_action(complex_sample<PT, ST>& result, const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y)
	{
		// multiplication of a complex 'x' by a conjugated complex 'y'
		result.add_multiply_conj(x, y);
	}

	template<class PT, class PT1, class ST, class ST1, class ST2>
	void	scalar_product_action(complex_sample<PT, ST>& result, const complex_sample<PT1, ST1>& x, const ST2& y)
	{
		// multiplication of a complex 'x' by a scalar 'y' (for templates compatibility)
		result.add_multiply(x, y);
	}

	template<class PT, class PT1, class ST, class ST1, class ST2>
	void	scalar_product_action(complex_sample<PT, ST>& result, const ST2& x, const complex_sample<PT1, ST1>& y)
	{
		// multiplication of a scalar 'x' by conjugated complex 'y' (for templates compatibility)
		result.add_multiply(~y, x);
	}

	template<class PT, class ST, class T1, class T2>
	void	scalar_product_action(complex_sample<PT, ST>& result, const T1& x, const T2& y)
	{
		// multiplication of two real values, result stored to a complex number (for templates compatibility)
		result += x * y;
	}


	//--------------------------------------------------------------
	//
	//	Utility
	//

	template<class T, class ST>
	inline double	amplitude_to_decibel(const complex_sample<T, ST> a)
	{
		return 10. * log10(cabs2(a));
	}

	template<class T, class ST>
	inline double	power_to_decibel(const complex_sample<T, ST> a)
	{
		return 5. * log10(cabs2(a));//5*(...cabs2) save a little on calculating the square root
	}

	//--------------------------------------------------------------
	//
	//	Predefined types
	//
	using complexF32 = complex_sample<float, double>;
	using complexF64 = complex_sample<double, double>;
	using complexI32F = complex_sample<int32_t, double>;
	using complexI16F = complex_sample<int16_t, double>;
	using complexI8F = complex_sample<int8_t, double>;
	using complexI32 = complex_sample<int32_t, int>;
	using complexI16 = complex_sample<int16_t, int>;
	using complexI8 = complex_sample<int8_t, int>;
	// No definitions for 'unsigned' types (nonsense)


#if !XRAD_STANDALONE_COMPLEX

	template<class PT, class ST>
	number_complexity_e complexity_e(const complex_sample<PT, ST>&) { return number_complexity_e::complex; }

	template<class PT, class ST>
	const number_complexity::complex* complexity_t(const complex_sample<PT, ST>&) { return nullptr; }

	//--------------------------------------------------------------

	//! \addtogroup gr_FloatingAnalog
	//! @{

	template<class T, class ST> struct FloatingAnalog32<complex_sample<T, ST>, typename enable_if<is_arithmetic_but_bool<T>::value>::type> { typedef complexF32 type; };
	template<class T, class ST> struct FloatingAnalog64<complex_sample<T, ST>, typename enable_if<is_arithmetic_but_bool<T>::value>::type> { typedef complexF64 type; };

	//! @} <!-- ^group gr_FloatingAnalog -->
	//! \addtogroup gr_ReducedWidth
	//! @{

	template <class T, class ST>
	struct	ReducedWidth<complex_sample<T, ST>>
	{
		using type = complex_sample<typename ReducedWidth<T>::type, ST>;
	};

	//! @} <!-- ^group gr_ReducedWidth -->

	//--------------------------------------------------------------

	check_if_number_traits_defined(complexF64)

#endif //!XRAD_STANDALONE_COMPLEX


} //namespace xrad

#include "complex_sample.hh"

//--------------------------------------------------------------
#endif // XRAD__File_complex_sample_h
