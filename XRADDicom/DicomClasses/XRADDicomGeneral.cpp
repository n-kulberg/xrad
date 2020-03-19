﻿/*!
	* \file XRADDicomGeneral.cpp
	* \date 2/16/2018 2:28:02 PM
	*
	* \author kovbas
	*
	* \brief
	*
	* TODO: long description
	*
	* \note
*/
#include "pre.h"
#include "XRADDicomGeneral.h"


XRAD_BEGIN

// выбор варианта сохранения в файл
const vector<wstring> &saving_decision_options()
{
	static	vector<wstring>	options(e_saving_decision_options_exit + 1);

	options[e_save_to_new_file] = L"Save to NEW File.";
	options[e_save_to_old_file] = L"Save to OLD File.";
	options[e_saving_decision_options_exit] = L"Don't save.";

	return	options;
}

namespace Dicom
{
	map<e_dcm_hierarchy_level_t, wstring> dcm_hierarchy_level_names =
	{
		{e_dcm_hierarchy_level_t::patient, L"patient"},
		{e_dcm_hierarchy_level_t::study, L"study"},
		{e_dcm_hierarchy_level_t::series, L"series"},
		{e_dcm_hierarchy_level_t::stack, L"stack"},
		{e_dcm_hierarchy_level_t::acquisition, L"acquisition"},
		{e_dcm_hierarchy_level_t::instance, L"instance"},
		{e_dcm_hierarchy_level_t::unknown, L"unknown"},
	};

} //end namespace Dicom


XRAD_END
