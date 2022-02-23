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

namespace window_function
{

	string	name(type wfe)
	{
		switch (wfe)
		{
		case window_function::constant:
			return string("Constant");
		case window_function::triangular:
			return string("Triangular");
		case window_function::cos2:
			return string("1+cos(x)/2");
		case window_function::hamming:
			return string("Hamming");
		case window_function::nuttall:
			return string("Nuttall");
		case window_function::blackman_harris:
			return string("Blackman-Harris");
		case window_function::blackman_nuttall:
			return string("Blackman-Nuttall");
		case window_function::flat_top:
			return string("Flat top");
		default:
			return string("Unknown window");
		};
	}

	RealFunctionF64 create(const window_function::type ew_left, const window_function::type ew_right, size_t size, size_t left_indent, size_t right_indent)
	{
		unique_ptr<internal::generator>	generator_left = internal::generator_by_enum(ew_left);
		unique_ptr<internal::generator>	generator_right = internal::generator_by_enum(ew_right);

		return internal::create_util(*generator_left, *generator_right, size, left_indent, right_indent);
	}

	RealFunctionF64 create(const window_function::type win, const size_t size, size_t left_indent, size_t right_indent)
	{
		return create(win, win, size, left_indent, right_indent);
	}

	namespace internal
	{

		//! \brief Считается по формуле exp(-x*x/(2*w);
		//! где x по всему массиву пробегает значения [-1;1]
		class	gauss_window : public generator
		{
			const double width;
			double	calculate(double x) const override { return gauss(2. * x - 1., width); }
		public:
			gauss_window(double w = sqrt(0.5)) : width(w) {}
		};


		// \brief Гибрид плоского и косинусного окна
		class	tukey_window : public generator
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


		//	High- and moderate-resolution windows

		struct	constant_window : public generator
		{
			double	calculate(double x) const override { return in_range(x, 0, 1) ? 1 : 0; }
			using generator::operator();
		};

		//	Окна среднего и выского разрешения

		/*!
			\brief Треугольное окно

			Can be seen as the convolution of two half-sized rectangular windows,
			giving it a main lobe width of twice the width of a regular rectangular window.
			The nearest lobe is -26 dB down from the main lobe.
		*/
		struct	triangular_window : public generator
		{
			double	calculate(double x) const
			{
				if (!in_range(x, 0, 1)) return 0;
				return 1. - fabs(2 * x - 1.);
			}
			using generator::operator();
		};

		using	bartlett_window = triangular_window;

		/*!
			\brief Класс косинусных окон

			Объединяет в один класс cos2_window, hamming_window.
			Видимо, следует подобным образом
			соединить окна Наттола-Блэкмена-Харриса, идущие ниже.
		*/
		struct	raised_cosine_window : public generator
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
			\brief Окно cos/2 + 1/2 (Hann window)

			The ends of the cosine just touch zero, so the side-lobes roll off at about 18 dB per octave.
			The Hann and Hamming windows, both of which are in the family known as "raised cosine" windows,
			are respectively named after Julius von Hann and Richard Hamming.
			The term "Hanning window" is sometimes used to refer to the Hann window.
		*/
		struct	cos2_window : public raised_cosine_window
		{
			cos2_window() : raised_cosine_window(0.5, 0.5) {}
		};


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



		//	Low-resolution (high-dynamic-range) windows
		//	Окна с большим динамическим диапазоном
		class nuttall_window : public generator
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
		class	blackman_harris_window : public generator
		{
			double	blackman_harris_window::calculate(double x) const override
			{
				if (!in_range(x, 0, 1)) return 0;
				return 0.35875 - 0.48829 * cos(2. * pi() * x) + 0.14128 * cos(4. * pi() * x) - 0.01168 * cos(6. * pi() * x);
			}
		};

		//! \brief All side lobes at level about 100 dB (almost non-decreasing)
		class	blackman_nuttall_window : public generator
		{
			double	calculate(double x) const override
			{
				if (!in_range(x, 0, 1)) return 0;
				return 0.3635819 - 0.4891775 * cos(2. * pi() * x) + 0.1365995 * cos(4. * pi() * x) - 0.0106411 * cos(6. * pi() * x);
			}
		};

		struct	flat_top_window : public generator
		{
			double	flat_top_window::calculate(double x) const override
			{
				if (!in_range(x, 0, 1)) return 0;
				return 1. - 1.93 * cos(2. * pi() * x) + 1.29 * cos(4. * pi() * x) - 0.388 * cos(6. * pi() * x) + 0.032 * cos(8. * pi() * x);
			}
		};

	}


	namespace internal
	{
		RealFunctionF64 create_util(const generator& generator_left, const generator& generator_right, size_t size, size_t left_indent, size_t right_indent)
		{
			RealFunctionF64	result(size, 0);
			size_t	right_boundary = result.size() - right_indent;
			if (right_boundary <= left_indent)
				return result;
			size_t	s = right_boundary - left_indent;
			if (s < 3) return result;

			auto it = result.begin();
			auto window_start = it + left_indent;
			auto window_middle = window_start + s / 2;
			auto window_end = window_start + s;
			auto it_end = result.end();


			for (; it < window_start; ++it)
				*it = 0;

			size_t	i = 0;
			for (; it < window_middle; ++it, ++i)
				*it = generator_left(i, s);

			for (; it < window_end; ++it, ++i)
				*it = generator_right(i, s);

			for (; it < it_end; ++it)
				*it = 0;

			return result;
		}

		unique_ptr<generator>	generator_by_enum(window_function::type wfe)
		{
			switch (wfe)
			{
			case window_function::constant:
				return make_unique<constant_window>();
			case window_function::triangular:
				return make_unique<triangular_window>();
			case window_function::cos2:
				return make_unique<cos2_window>();
			case window_function::hamming:
				return make_unique<hamming_window>();
			case window_function::nuttall:
				return make_unique<nuttall_window>();
			case window_function::blackman_harris:
				return make_unique<blackman_harris_window>();
			case window_function::blackman_nuttall:
				return make_unique<blackman_nuttall_window>();
			case window_function::flat_top:
				return make_unique<flat_top_window>();
			default:
				throw invalid_argument("GetWindowFunctionByEnum, unknown argument");
			};
		}
	}//namespace internal
}//namespace window_function



XRAD_END
