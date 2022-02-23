/*
	Copyright (c) 2021, Moscow Center for Diagnostics & Telemedicine
	All rights reserved.
	This file is licensed under BSD-3-Clause license. See LICENSE file for details.
*/
//	file WindowFunction.h
//--------------------------------------------------------------
#ifndef XRAD__File_window_function_h
#define XRAD__File_window_function_h
/*!
	\file
	\brief Вычисление типовых оконных функций,
	а также применение их к данным в контейнере MathFunction.
	Во всех случаях используется функтор от
	двух целых: номер индекса и размер массива
*/
//--------------------------------------------------------------

#include "XRADBasic/MathFunctionTypes.h"
#include <memory>

XRAD_BEGIN

//--------------------------------------------------------------

//--------------------------------------------------------------
//
//	Задание окна через enum
//
//--------------------------------------------------------------


namespace window_function
{
	enum	type
	{
		constant = 0,
		triangular,
		bartlett = triangular,
		cos2,
		hann = cos2,
		hamming,
		nuttall,
		blackman_harris,
		blackman_nuttall,
		flat_top,

		n_window_functions
	};

	string	name(type);


	template <typename ...args> void apply(MathFunction<args...>& original, window_function::type win_left, window_function::type win_right, size_t left_indent = 0, size_t right_indent = 0);
	template <typename ...args> void apply(MathFunction<args...>& original, window_function::type win, size_t left_indent = 0, size_t right_indent = 0);

	template <typename ...args> void create(MathFunction<args...>& original, window_function::type win_left, window_function::type win_right, size_t left_indent = 0, size_t right_indent = 0);
	template <typename ...args> void create(MathFunction<args...>& original, window_function::type win, size_t left_indent = 0, size_t right_indent = 0);

	RealFunctionF64 create(const window_function::type w_left, const window_function::type w_right, size_t size, size_t left_indent = 0, size_t right_indent = 0);
	RealFunctionF64 create(const window_function::type win, const size_t size, size_t left_indent = 0, size_t right_indent = 0);

}//namespace window_function


//	Вычисление оконной функции заданного размера
//	Возможно несимметричное окно, состоящее	из двух половинок стандартных окон.


#if 0
//Legacy definitions may be used to make work old code. 
using window_function_e = window_function::type;
#define ApplyWindowFunction window_function::apply
#define CreateWindowFunction window_function::create
#define GetWindowFunctionName window_function::name;
#endif


XRAD_END

#include "WindowFunction.hh"

//--------------------------------------------------------------
#endif //XRAD__File_window_function_h
