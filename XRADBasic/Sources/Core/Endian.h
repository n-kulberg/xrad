﻿//	file Endian.h
//	Created by KNS on 10.11.10
//--------------------------------------------------------------
#ifndef __Endian_h
#define __Endian_h
/*!
	\addtogroup gr_CompilerSpecific
	\defgroup gr_CompilerSpecificByteOrder Разрядность и порядок байтов
	\ingroup gr_CompilerSpecific
	@{
	\file
	\brief Определение констант порядка байтов

	Внутренний файл библиотеки.

	Этот файл автоматически подключается в <CompilerSpecific.h>.
	@}

	\addtogroup gr_CompilerSpecificByteOrder
	@{
	# Разрядность

	Библиотека поддерживает размер байта 8 бит, размер int 32 и 64 бита.

	# Порядок байтов целевой платформы

	Информация о порядке байтов в различных компиляторах подается по-разному (если вообще подается).

	Здесь реализован следующий механизм. В файле <Endian.h> заданы числовые макроопределения
	XRAD_LITTLE_ENDIAN, XRAD_BIG_ENDIAN, определяющие конкретный порядок байтов. При подключении
	<CompilerSpecific.h> задается макрос XRAD_ENDIAN со значением, равным одному из заданных
	здесь значений.

	В коде, который зависит от порядка байтов, при организации ветвления должно производиться сравнение
	XRAD_ENDIAN со всеми имеющимися здесь значениями:
	- XRAD_LITTLE_ENDIAN;
	- XRAD_BIG_ENDIAN.

	При несовпадении значения ни с одним из этих должна выдаваться ошибка компиляции.

	Это требование связано с тем, что существуют архитектуры с порядком байтов, отличным
	от big endian и little endian. Поэтому в будущем при переносе на такие архитектуры
	возможно появление новых констант, задающих порядок байтов.

	Замечание. Числовые значения констант XRAD_LITTLE_ENDIAN и т.д. могут в будущем поменяться.

	# Способы автоматического определения порядка байтов на этапе компиляции

	Краткая сводка, что приходилось наблюдать:
	- на TI CCS (gcc) задано  __little_endian__ или __big_endian__ в зависимости от опций;
	- на CodeWarrior задано было __BIG_ENDIAN__;
	- на BCB ничего не было, ориентировались по #if defined(_M_IX86);
	- на MSVC то же самое;
	- gcc macosx? -- todo.

	Была попытка реализовать универсальный алгоритм определения порядка байтов
	для платформ с CHAR_BIT==8:
	задавать макроопределение XRAD_ENDIAN, которое автоматически оказывается
	равным либо XRAD_LITTLE_ENDIAN, либо XRAD_BIG_ENDIAN. Для этого константы с порядком
	байтов задавались таким образом, чтобы каждый байт двоичного представления
	содержал свой порядковый номер. А константа XRAD_ENDIAN задавалась
	особым образом. А именно:

	~~~~
	... // #if по размеру uint. Опущены варианты для разрядностей, отличных от 32.
	#elif UINT_MAX == 0xffffffff
	#define XRAD_LITTLE_ENDIAN 0x01020304
	#define XRAD_BIG_ENDIAN 0x04030201
	#define XRAD_PDP_ENDIAN 0x03040102
	#define XRAD_ENDIAN '\4\3\2\1'
	#elif
	...
	~~~~

	На GCC этот трюк работал, кажется: при задании числа строковым литералом
	порядок байтов был такой же, как и в памяти. На MSVS это не работало:
	порядок символов литерала получался такой же, как и в числе.

	@}
*/

//--------------------------------------------------------------

//! \brief Константа для сравнения с XRAD_ENDIAN: если равны, то порядок байтов целевой платформы little endian
#define XRAD_LITTLE_ENDIAN 1

//! \brief Константа для сравнения с XRAD_ENDIAN: если равны, то порядок байтов целевой платформы big endian
#define XRAD_BIG_ENDIAN 2

//--------------------------------------------------------------
#endif // __Endian_h
