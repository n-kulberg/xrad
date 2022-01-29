/*
	Copyright (c) 2021, Moscow Center for Diagnostics & Telemedicine
	All rights reserved.
	This file is licensed under BSD-3-Clause license. See LICENSE file for details.
*/
#include "pre.h"
#include "WindowFunction.h"

XRAD_BEGIN

//--------------------------------------------------------------
//
//	Параметрические окна
//
namespace window_function_private
{


	//--------------------------------------------------------------

	//! \brief Считается по формуле exp(-x*x/(2*w);
	//! где x по всему массиву пробегает значения [-1;1]
	class	gauss_window : public wf_generator
	{
		const double width;
		double	calculate(double x) const override { return gauss(2. * x - 1., width); }
	public:
		gauss_window(double w = sqrt(0.5)) : width(w) {}
	};


	// \brief Гибрид плоского и косинусного окна
	class	tukey_window : public wf_generator
	{
		const double flatness;
		double	calculate(double x) const override
		{
			if (!in_range(x, 0, 1)) return 0;
			// Написано вчерне, есть недочеты

			double	argument = 2. * x - 1.;
			//	double	x = fabs(2.*double(n)/double(N-1) - 1.);//четная функция


			if (argument <= flatness) return 1.;
			else
			{
				double	dx = (argument - flatness) / (1. - flatness);
				return (1. + cos(pi() * dx)) / 2.;
			}
		}

	public:
		tukey_window(double f = 0.5) : flatness(f) {}
	};


	//--------------------------------------------------------------
	//
	//	High- and moderate-resolution windows
	//
	//--------------------------------------------------------------



	struct	constant_window : public wf_generator
	{
		double	calculate(double x) const override { return in_range(x, 0, 1) ? 1 : 0; }
		using wf_generator::operator();
	};

	//--------------------------------------------------------------
	//
	//	Окна среднего и выского разрешения
	//

	/*!
		\brief Треугольное окно

		Can be seen as the convolution of two half-sized rectangular windows,
		giving it a main lobe width of twice the width of a regular rectangular window.
		The nearest lobe is -26 dB down from the main lobe.
	*/
	struct	triangular_window : public wf_generator
	{
		double	calculate(double x) const
		{
			if (!in_range(x, 0, 1)) return 0;
			return 1. - fabs(2 * x - 1.);
		}
		using wf_generator::operator();
	};

	using	bartlett_window = triangular_window;

	/*!
		\brief Класс косинусных окон

		Объединяет в один класс cos2_window, hamming_window.
		Видимо, следует подобным образом
		соединить окна Наттола-Блэкмена-Харриса, идущие ниже.
	*/
	struct	raised_cosine_window : public wf_generator
	{
	protected:
		const double	a0, a1;
		raised_cosine_window(double in_a0, double in_a1) :a0(in_a0), a1(in_a1) {};
		double	calculate(double x) const
		{
			if (!in_range(x, 0, 1)) return 0;
			return a0 - a1 * cos(2. * pi() * x);
		}
	};

	/*!
		\brief Окно cos^2 (Hann window)

		The ends of the cosine just touch zero, so the side-lobes roll off at about 18 dB per octave.
		The Hann and Hamming windows, both of which are in the family known as "raised cosine" windows,
		are respectively named after Julius von Hann and Richard Hamming.
		The term "Hanning window" is sometimes used to refer to the Hann window.
	*/
	struct	cos2_window : public raised_cosine_window
	{
		cos2_window() : raised_cosine_window(0.5, 0.5) {}
	};

	typedef	cos2_window hann_window;

	/*!
		\brief Окно Хэмминга (уточненные коэффициенты)

		The "raised cosine" with these particular coefficients
		was proposed by Richard W. Hamming. The window is optimized
		to minimize the maximum (nearest) side lobe, giving it a height
		of about one-fifth that of the Hann window,
		a raised cosine with simpler coefficients.

		Здесь используется уточненная формула. Классическая формула
		предполагает значения коэффициентов 0.54, 0.46.
	*/
	struct	hamming_window : public raised_cosine_window
	{
		hamming_window() : raised_cosine_window(0.53836, 0.46164) {}
	};



	//--------------------------------------------------------------
	//
	//	Low-resolution (high-dynamic-range) windows
	//	Окна с большим динамическим диапазоном
	//
	//--------------------------------------------------------------



	class nuttall_window : public wf_generator
	{
		double	calculate(double x) const override
		{
			if (!in_range(x, 0, 1)) return 0;
			return 0.355768 - 0.487396 * cos(2. * pi() * x) + 0.144232 * cos(4. * pi() * x) - 0.012604 * cos(6. * pi() * x);
		}
	};

	/*!
		\brief Окно Блэкмана-Харриса

		A generalization of the Hamming family, produced by adding more shifted sinc functions,
		meant to minimize side-lobe levels.
	*/
	class	blackman_harris_window : public wf_generator
	{
		double	blackman_harris_window::calculate(double x) const override
		{
			if (!in_range(x, 0, 1)) return 0;
			return 0.35875 - 0.48829 * cos(2. * pi() * x) + 0.14128 * cos(4. * pi() * x) - 0.01168 * cos(6. * pi() * x);
		}
	};

	//! \brief All side lobes at level about 100 dB (almost non-decreasing)
	class	blackman_nuttall_window : public wf_generator
	{
		double	calculate(double x) const override
		{
			if (!in_range(x, 0, 1)) return 0;
			return 0.3635819 - 0.4891775 * cos(2. * pi() * x) + 0.1365995 * cos(4. * pi() * x) - 0.0106411 * cos(6. * pi() * x);
		}
	};

	struct	flat_top_window : public wf_generator
	{
		double	flat_top_window::calculate(double x) const override
		{
			if (!in_range(x, 0, 1)) return 0;
			return 1. - 1.93 * cos(2. * pi() * x) + 1.29 * cos(4. * pi() * x) - 0.388 * cos(6. * pi() * x) + 0.032 * cos(8. * pi() * x);
		}
	};

}



string	GetWindowFunctionName(window_function_e wfe)
{
	switch (wfe)
	{
	case e_constant_window:
		return string("Constant");
	case e_triangular_window:
		return string("Triangular");
	case e_cos2_window:
		return string("1+cos(x)/2");
	case e_hamming_window:
		return string("Hamming");
	case e_nuttall_window:
		return string("Nuttall");
	case e_blackman_harris_window:
		return string("Blackman-Harris");
	case e_blackman_nuttall_window:
		return string("Blackman-Nuttall");
	case e_flat_top_window:
		return string("Flat top");
	default:
		return string("Unknown window");
	};
}

namespace window_function_private
{
	RealFunctionF64 create(const wf_generator& w_left, const wf_generator& w_right, size_t size, size_t s0, size_t s1)
	{
		RealFunctionF64	result(size, 0);
		if (!s1) s1 = result.size();
		if (s1 <= s0)
			return result;
		size_t	s = s1 - s0;
		if (s < 3) return result;

		auto it = result.begin();
		auto window_start = it + s0;
		auto window_middle = window_start + s / 2;
		auto window_end = window_start + s;
		auto it_end = result.end();


		for (; it < window_start; ++it)
			*it = 0;

		size_t	i = 0;
		for (; it < window_middle; ++it, ++i)
			*it = w_left(i, s);

		for (; it < window_end; ++it, ++i)
			*it = w_right(i, s);

		for (; it < it_end; ++it)
			*it = 0;

		return result;
	}

	unique_ptr<window_function_private::wf_generator>	GetWindowFunctionByEnum(window_function_e wfe)
	{
		switch (wfe)
		{
		case e_constant_window:
			return make_unique<constant_window>();
		case e_triangular_window:
			return make_unique<triangular_window>();
		case e_cos2_window:
			return make_unique<cos2_window>();
		case e_hamming_window:
			return make_unique<hamming_window>();
		case e_nuttall_window:
			return make_unique<nuttall_window>();
		case e_blackman_harris_window:
			return make_unique<blackman_harris_window>();
		case e_blackman_nuttall_window:
			return make_unique<blackman_nuttall_window>();
		case e_flat_top_window:
			return make_unique<flat_top_window>();
		default:
			throw invalid_argument("GetWindowFunctionByEnum, unknown argument");
		};
	}
}//namespace window_function_private

XRAD_END
