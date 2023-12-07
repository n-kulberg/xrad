/*
	Copyright (c) 2021, Moscow Center for Diagnostics & Telemedicine
	All rights reserved.
	This file is licensed under BSD-3-Clause license. See LICENSE file for details.
*/
// file complex_sample.h
//
// A part of XRAD
// Complex numbers library for scientific calculations.
//
// 2014 KNS typedefs complexF/complexD --> complexF32/complexF64
// 2013 KNS class complexF --> template<class PT, class ST> complex_sample{}
// 1999 Modified by ACS (Алексей Борисович Елизаров) (added constructor with SNoInit)
// 19?? class complexF created by Nicholas S. Kulberg (Кульберг Николай Сергеевич)
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
#endif //XRAD_STANDALONE_COMPLEX

namespace xrad {


	//--------------------------------------------------------------
	//
	// шаблон complex_sample<scalar>
	// "~" комплексное сопряжение
	// "%", "%=" умножение на комплексно сопряженное
	//
	//--------------------------------------------------------------


	enum class complex_number_part
	{
		real,
		imag
	};


	/*!
		\brief Комплексное число

		PT должен быть неконстантным.
		Использование complex_sample<const PT, ST> считаем неправильным.
		Там, где требуется константность, следует использовать const complex_sample<PT, ST>.
	*/
	template<class PT, class ST>
	class complex_sample
	{
		// Запрещаем использование complex_sample<const PT, ST> явным образом.
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
		// конструктор по умолчанию -- без инициализации

		explicit complex_sample(PT r) { re = r; im = 0; }
		// explicit обязательно во избежание
		// случайных появлений комплексных чисел вместо действительных

		complex_sample(PT r, PT i) { re = r; im = i; }

		template<class PT1, class ST1>
		complex_sample(const complex_sample<PT1, ST1>& c) { re = c.re; im = c.im; }

		// assignments

		template<class PT1, class ST1>
		complex_sample& operator = (const complex_sample<PT1, ST1>& y) { re = y.re; im = y.im; return *this; }

		complex_sample<PT, ST>& operator = (PT y) { re = y; im = 0; return *this; }

		// unary arithmetics

		complex_sample<PT, ST> operator-() const { return complex_sample<PT, ST>(-re, -im); }
		complex_sample<PT, ST> operator~() const { return complex_sample<PT, ST>(re, -im); }

		complex_sample<PT, ST>& operator++() { ++re; return *this; }
		complex_sample<PT, ST>& operator--() { --re; return *this; }

		complex_sample<PT, ST> operator++(int) { complex_sample<PT, ST> result(*this); re++; return result; }
		complex_sample<PT, ST> operator--(int) { complex_sample<PT, ST> result(*this); re--; return result; }

		// assignment arithmetic

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
		//
		complex_sample<PT, ST>& operator += (PT y) { re += y; return *this; }
		complex_sample<PT, ST>& operator -= (PT y) { re -= y; return *this; }
		complex_sample<PT, ST>& operator *= (ST y) { re *= y; im *= y; return *this; }
		complex_sample<PT, ST>& operator /= (ST y) { re /= y; im /= y; return *this; }
		//

		// arithmetic

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



		// тернарные действия, результат над двумя аргументами помещается в *this. позволяет избежать
		// создания буферных переменных при операциях вида a=b+c;

		template<class PT1, class ST1>
		complex_sample<PT, ST>& conjugate(const complex_sample<PT1, ST1>& x) { re = x.re; im = -im.re; return *this; }

		// действия вида (*this)=x*y с комплексными числами
		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& add(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y);
		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& add_i(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y);

		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& subtract(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y);
		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& subtract_i(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y);

		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& multiply(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y);
		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& multiply_conj(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y);

		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& divide(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y);
		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& divide_conj(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y);



		// действия вида (*this)+=x*y с комплексными числами
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



		// действия вида (*this)=x*y с комплексным и скаляром
		template<class PT1, class ST1>
		complex_sample<PT, ST>& multiply(const complex_sample<PT1, ST1>& x, const scalar_type& a);

		template<class PT1, class ST1>
		complex_sample<PT, ST>& divide(const complex_sample<PT1, ST1>& x, const scalar_type& a);



		// действия вида (*this)+=x*y с комплексным и скаляром
		template<class PT1, class ST1>
		complex_sample<PT, ST>& add_multiply(const complex_sample<PT1, ST1>& x, const scalar_type& a);

		template<class PT1, class ST1>
		complex_sample<PT, ST>& add_divide(const complex_sample<PT1, ST1>& x, const scalar_type& a);

		template<class PT1, class ST1>
		complex_sample<PT, ST>& subtract_multiply(const complex_sample<PT1, ST1>& x, const scalar_type& a);

		template<class PT1, class ST1>
		complex_sample<PT, ST>& subtract_divide(const complex_sample<PT1, ST1>& x, const scalar_type& a);



		// взвешенное сложение двух массивов, (*this) = x*a1 + y*a2
		template<class PT1, class ST1, class PT2, class ST2>
		complex_sample<PT, ST>& mix(const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y, scalar_type a1, scalar_type a2);



		// compare

		template<class PT1, class ST1>
		bool operator == (const complex_sample<PT1, ST1>& x) const { return (re == x.re && im == x.im); }
		template<class PT1, class ST1>
		bool operator != (const complex_sample<PT1, ST1>& x) const { return (re != x.re || im != x.im); }
		template<class PT1, class ST1>
		bool operator < (const complex_sample<PT1, ST1>& x) const { return (re * re + im * im < x.re* x.re + x.im * x.im); }
		template<class PT1, class ST1>
		bool operator > (const complex_sample<PT1, ST1>& x) const { return (re * re + im * im > x.re * x.re + x.im * x.im); }
		template<class PT1, class ST1>
		bool operator <= (const complex_sample<PT1, ST1>& x) const { return (re * re + im * im <= x.re * x.re + x.im * x.im); }
		template<class PT1, class ST1>
		bool operator >= (const complex_sample<PT1, ST1>& x) const { return (re * re + im * im >= x.re * x.re + x.im * x.im); }
		//
		bool operator == (PT x) const { return (re == x && !im); }
		bool operator != (PT x) const { return (re != x || im); }
		bool operator < (PT x) const { return ((re * re + im * im) < x * fabs(x)); }
		bool operator > (PT x) const { return ((re * re + im * im) > x * fabs(x)); }
		bool operator <= (PT x) const { return ((re * re + im * im) <= x * fabs(x)); }
		bool operator >= (PT x) const { return ((re * re + im * im) >= x * fabs(x)); }



		//	определение расстояния между компонентами комплексного числа (вспомогательная функция)
		//	чаще всего оно равно 1, но выравнивание полей внутри структуры может дать и большее значение
		//
		//	функция используется при вычислении шага действительной или мнимой компоненты комплексного массива MathFunctionC.
		//	писать новый шаг как 2*step() нехорошо, т.к., в общем случае не гарантируется, что sizeof(complex_sample<T>) == 2*sizeof(T).
		//	теоретически это расстояние может оказаться даже некратным размеру T (sizeof(self) != n*sizeof(T));
		//	в таком случае использовать функции real(MathFunctionC) и imag(MathFunctionC) нельзя.
		//	в этом случае при попытке использовать эти функции возникнет ошибка компилятора.
		static ptrdiff_t parts_distance()
		{
			enum
			{
				self_size = sizeof(self),
				part_size = sizeof(part_type),
				distance = self_size / part_size
			};
			static_assert(distance * part_size == self_size, "complex_sample layout problem.");
			// если компилятор здесь выдаст ошибку, см. комментарий выше

			return distance;
		}
	private:
		template<class PT1, class ST1>
		complex_sample<PT, ST> complex_division_algorithm(const complex_sample<PT, ST> x, const complex_sample<PT1, ST1> y) const;
	};



	//--------------------------------------------------------------
	//
	// number traits definitions. see comment in NumberTraits.h
	//

	template<class PT, class ST>
	number_complexity_e complexity_e(const complex_sample<PT, ST>&) { return number_complexity_e::complex; }

	template<class PT, class ST>
	const number_complexity::complex* complexity_t(const complex_sample<PT, ST>&) { return nullptr; }

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

	//
	//--------------------------------------------------------------



	//--------------------------------------------------------------
	//
	//	сопряженное умножение для вычисления
	//	скалярного произведения. при использовании комплексных
	//	векторов. занимает место универсального шаблона, объявленного
	//	в файле NumberTraits.h
	//
	template<class PT, class PT1, class PT2, class ST, class ST1, class ST2>
	void	scalar_product_action(complex_sample<PT, ST>& result, const complex_sample<PT1, ST1>& x, const complex_sample<PT2, ST2>& y)
	{
		// сопряженное перемножение компонент комплексных векторов
		result.add_multiply_conj(x, y);
	}

	template<class PT, class PT1, class ST, class ST1, class ST2>
	void	scalar_product_action(complex_sample<PT, ST>& result, const complex_sample<PT1, ST1>& x, const ST2& y)
	{
		// умножение комплексного числа на действительное
		result.add_multiply(x, y);
	}

	template<class PT, class PT1, class ST, class ST1, class ST2>
	void	scalar_product_action(complex_sample<PT, ST>& result, const ST2& x, const complex_sample<PT1, ST1>& y)
	{
		// умножение действительного числа на комплексное сопряженное
		result.add_multiply(~y, x);
	}

	template<class PT, class ST, class T1, class T2>
	void	scalar_product_action(complex_sample<PT, ST>& result, const T1& x, const T2& y)
	{
		// перемножение действительных чисел, запись в комплексный результат. такое тоже иногда бывает нужно
		result += x * y;
	}

	//
	//--------------------------------------------------------------



	template<class T, class ST>
	inline double	amplitude_to_decibel(const complex_sample<T, ST> a)
	{
		return 10. * log10(cabs2(a));//5*(...cabs2) небольшая экономия на вычислении квадратного корня
	}

	template<class T, class ST>
	inline double	power_to_decibel(const complex_sample<T, ST> a)
	{
		return 5. * log10(cabs2(a));//5*(...cabs2) небольшая экономия на вычислении квадратного корня
	}



	//--------------------------------------------------------------



	typedef	complex_sample<float, double> complexF32;
	typedef	complex_sample<double, double> complexF64;

	typedef	complex_sample<int32_t, double> complexI32F;
	typedef	complex_sample<int16_t, double> complexI16F;
	typedef	complex_sample<int8_t, double> complexI8F;

	typedef	complex_sample<int32_t, int> complexI32;
	typedef	complex_sample<int16_t, int> complexI16;
	typedef	complex_sample<int8_t, int> complexI8;
	// для unsigned определений не задаем, нонсенс


#if !XRAD_STANDALONE_COMPLEX

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
