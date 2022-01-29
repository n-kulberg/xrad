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
namespace window_function_private
{

class	wf_generator
{


	//! \brief Вычисление веса окна от double на отрезке (0,1)
	virtual double	calculate(double x) const = 0;

public:
// Вычисление веса от дискретного аргумента для окна длины N.
// Возможны три варианта, о выборе подумать.
//	Вариант с https://en.wikipedia.org/wiki/Window_function : { return operator()(double(n)/double(N-1)); }
// приводит к обнулению концов отрезка на многих окнах.
//! Здесь нули подразумеваются за пределами отрезка, кажется, должно быть более правильно
	virtual double	operator()(size_t n, size_t N) const { return calculate(double(n + 1) / double(N + 1)); }
// Промежуточный вариант между 1 и 2. Для приподнятого косинуса сохраняет разбиение единицы
//virtual double	operator()(size_t n, size_t N) const { return operator()((double(n)+0.5)/double(N)); }
	virtual ~wf_generator() {}

};




RealFunctionF64 create(const wf_generator& w_left, const wf_generator& w_right, size_t size, size_t s0 = 0, size_t s1 = 0);

template <typename ...args>
void apply(MathFunction<args...> &function, const wf_generator & w_left, const wf_generator & w_right, size_t s0, size_t s1)
{
	auto wf = create(w_left, w_right, function.size(), s0, s1);
	function *= wf;

}

template <typename ...args>
void apply(MathFunction<args...> &function, const wf_generator & win, size_t s0, size_t s1)
{
	apply(function, win, win, s0, s1);
}

template <typename ...args>
void create(MathFunction<args...> &function, const wf_generator & w_left, const wf_generator & w_right, size_t s0, size_t s1)
{
	using value_type = typename MathFunction<args...>::value_type;
	function.fill(value_type(1));
	apply(function, w_left, w_right, s0, s1);
}

template <typename ...args>
void create(MathFunction<args...> &function, const wf_generator & win, size_t s0, size_t s1)
{
	function.fill(typename MathFunction<args...>::value_type(1));
	apply(function, win, s0, s1);
}

unique_ptr<window_function_private::wf_generator> GetWindowFunctionByEnum(window_function_e);

}//window_function_private


//--------------------------------------------------------------



template <typename ...args>
void ApplyWindowFunction(MathFunction<args...> &function, window_function_e ew_left, window_function_e ew_right, size_t s0, size_t s1)
{
	unique_ptr<window_function_private::wf_generator>	w_left = window_function_private::GetWindowFunctionByEnum(ew_left);
	unique_ptr<window_function_private::wf_generator>	w_right = window_function_private::GetWindowFunctionByEnum(ew_right);

	window_function_private::apply(function, *w_left, *w_right, s0, s1);
}

template <typename ...args>
void ApplyWindowFunction(MathFunction<args...> &function, window_function_e ewin, size_t s0, size_t s1)
{
	unique_ptr<window_function_private::wf_generator>	win = window_function_private::GetWindowFunctionByEnum(ewin);
	window_function_private::apply(function, *win, s0, s1);
}

template <typename ...args>
void CreateWindowFunction(MathFunction<args...> &function, window_function_e ew_left, window_function_e ew_right, size_t s0, size_t s1)
{
	unique_ptr<window_function_private::wf_generator>	w_left = window_function_private::GetWindowFunctionByEnum(ew_left);
	unique_ptr<window_function_private::wf_generator>	w_right = window_function_private::GetWindowFunctionByEnum(ew_right);

	window_function_private::create(function, *w_left, *w_right, s0, s1);
}

template <typename ...args>
void CreateWindowFunction(MathFunction<args...> &function, window_function_e ewin, size_t s0, size_t s1)
{
	unique_ptr<window_function_private::wf_generator>	win = window_function_private::GetWindowFunctionByEnum(ewin);
	window_function_private::create(function, *win, s0, s1);
}



XRAD_END

#endif //XRAD__File_window_function_cc
