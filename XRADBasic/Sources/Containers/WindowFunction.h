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



enum	window_function_e
{
	e_constant_window = 0,
	e_triangular_window,
	e_bartlett_window = e_triangular_window,
	e_cos2_window,
	e_hann_window = e_cos2_window,
	e_hamming_window,
	e_nuttall_window,
	e_blackman_harris_window,
	e_blackman_nuttall_window,
	e_flat_top_window,

	n_window_functions
};

//	Вычисление оконной функции заданного размера
//	Возможно несимметричное окно, состоящее	из двух половинок стандартных окон.


template <typename ...args> void ApplyWindowFunction(MathFunction<args...>& original, window_function_e win_left, window_function_e win_right, size_t s0 = 0, size_t s1 = 0);
template <typename ...args> void ApplyWindowFunction(MathFunction<args...>& original, window_function_e win, size_t s0 = 0, size_t s1 = 0);

template <typename ...args> void CreateWindowFunction(MathFunction<args...>& original, window_function_e win_left, window_function_e win_right, size_t s0 = 0, size_t s1 = 0);
template <typename ...args> void CreateWindowFunction(MathFunction<args...>& original, window_function_e win, size_t s0 = 0, size_t s1 = 0);

string	GetWindowFunctionName(window_function_e);



XRAD_END

#include "WindowFunction.hh"

//--------------------------------------------------------------
#endif //XRAD__File_window_function_h
