/*
	Copyright (c) 2021, Moscow Center for Diagnostics & Telemedicine
	All rights reserved.
	This file is licensed under BSD-3-Clause license. See LICENSE file for details.
*/
#ifndef XRAD__File_window_function_cc
#define XRAD__File_window_function_cc

XRAD_BEGIN



//--------------------------------------------------------------
//
//	Window function
//
//--------------------------------------------------------------
namespace window_function
{
namespace internal
{

class	generator
{
	//! \brief Вычисление веса окна от double на отрезке (0,1)
	virtual double	calculate(double x) const = 0;

public:
// Вычисление веса от дискретного аргумента для окна длины N.
// Возможны три варианта, о выборе подумать.
//	Вариант с https://en.wikipedia.org/wiki/Window_function : 
//	virtual double	operator()(size_t n, size_t N) const { return calculate(double(n)/double(N-1)); }
// приводит к обнулению концов отрезка на многих окнах.
// Здесь нули подразумеваются за пределами отрезка, кажется, должно быть более правильно
//	virtual double	operator()(size_t n, size_t N) const { return calculate(double(n + 1) / double(N + 1)); }
// Промежуточный вариант между 1 и 2. Для приподнятого косинуса сохраняет разбиение единицы
	virtual double	operator()(size_t n, size_t N) const { return calculate((double(n)+0.5)/double(N)); }
	virtual ~generator() {}

};




RealFunctionF64 create_util(const generator& generator_left, const generator& generator_right, size_t size, size_t left_indent = 0, size_t right_indent = 0);

template <typename ...args>
void apply_util(MathFunction<args...> &function, const generator & generator_left, const generator & generator_right, size_t left_indent, size_t right_indent)
{
	auto wf = create_util(generator_left, generator_right, function.size(), left_indent, right_indent);
	function *= wf;

}

template <typename ...args>
void apply_util(MathFunction<args...> &function, const generator & win, size_t left_indent, size_t right_indent)
{
	apply_util(function, win, win, left_indent, right_indent);
}

template <typename ...args>
void create_util(MathFunction<args...> &function, const generator & generator_left, const generator & generator_right, size_t left_indent, size_t right_indent)
{
	using value_type = typename MathFunction<args...>::value_type;
	function.fill(value_type(1));
	apply_util(function, generator_left, generator_right, left_indent, right_indent);
}

template <typename ...args>
void create_util(MathFunction<args...> &function, const generator & win, size_t left_indent, size_t right_indent)
{
	function.fill(typename MathFunction<args...>::value_type(1));
	apply_util(function, win, left_indent, right_indent);
}

unique_ptr<internal::generator> generator_by_enum(window_function::type);

}//namespace internal


template <typename ...args>
void apply(MathFunction<args...>& function, window_function::type ew_left, window_function::type ew_right, size_t left_indent, size_t right_indent)
{
	auto	generator_left = internal::generator_by_enum(ew_left);
	auto	generator_right = internal::generator_by_enum(ew_right);

	internal::apply_util(function, *generator_left, *generator_right, left_indent, right_indent);
}

template <typename ...args>
void apply(MathFunction<args...>& function, window_function::type ewin, size_t left_indent, size_t right_indent)
{
	auto	win = internal::generator_by_enum(ewin);
	internal::apply_util(function, *win, left_indent, right_indent);
}

template <typename ...args>
void create(MathFunction<args...>& function, window_function::type ew_left, window_function::type ew_right, size_t left_indent, size_t right_indent)
{
	auto	generator_left = internal::generator_by_enum(ew_left);
	auto	generator_right = internal::generator_by_enum(ew_right);

	internal::create_util(function, *generator_left, *generator_right, left_indent, right_indent);
}

template <typename ...args>
void create(MathFunction<args...>& function, window_function::type ewin, size_t left_indent, size_t right_indent)
{
	auto	win = internal::generator_by_enum(ewin);
	internal::create_util(function, *win, left_indent, right_indent);
}



}//namespace window_function


//--------------------------------------------------------------




XRAD_END

#endif //XRAD__File_window_function_cc
