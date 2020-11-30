﻿/*!
	\file
	\date 2019/10/21 10:00
	\author novik

	\brief Реализация функций для создания-записи-считывания json объектов с информацией о файлах
	в директории
*/
#include "pre.h"
#include "SingleDirectoryIndexJson.h"

#include "DicomClustering.h"
#include "DicomCatalogIndex.h"
#include "DicomFileIndexJson.h"

#include <XRADDicom/Sources/DicomClasses/Instances/ct_slice.h>
#include <XRADDicom/Sources/DicomClasses/Instances/xray_image.h>
#include <XRADDicom/Sources/DicomClasses/Instances/mr_slice.h>
#include <XRADDicom/Sources/DicomClasses/Instances/mr_slice_siemens.h>

#include <XRADDicom/XRADDicom.h>

#include <XRADBasic/ThirdParty/nlohmann/json.hpp>
#include <typeinfo>

/// директива для возможности тестирования записанных и считанных json объектов
#define TEST_JSON

XRAD_BEGIN

namespace Dicom
{



using json = nlohmann::json;

namespace
{

const map<string, string> map_header_json_type1 =
{
	{ "ID","purpose: Dicom catalog" },
	{ "version","0.0" },
	{ "type", index_file_label(index_file_type::hierarchical) }
};

const map<string, string> map_header_json_type2 =
{
	{ "ID","purpose: Dicom catalog" },
	{ "version","0.0" },
	{ "type", index_file_label(index_file_type::raw) }
};



// из списка тэгов для каждого файла сгенерировать json файл type1 (sample1.json)
void dir_info_to_json_type1(const SingleDirectoryIndex& dcmDirectoryIndex, json& json_type1)
{
	// заполнить обязательные поля о "ID", "version", "type" из map_header_json_type1
	for (auto& map_v : map_header_json_type1)
		json_type1[map_v.first] = map_v.second;

	// выделить уникальные исследования по 5-ти первым тэгам
	json json_dicom_files; // хранение информации о файлах в древовидной структуре первых 5-тэгов
	json json_not_dicom; // информация о не-DICOM файлах (идентично записи в json type 2)

	// Для определенности сохраняем элементы в порядке "возрастания" имени файла с учетом регистра
	// в utf-8.
	map<string, const DicomFileIndex*> ordered_index;
	for (auto &di: dcmDirectoryIndex)
	{
		ordered_index[convert_to_string8(di.get_file_name())] = &di;
	}

	for (auto &odi: ordered_index)
	{
		auto &file_tags = *odi.second;
		if (file_tags.is_dicom())
		{
			json json_inner_file_tag;
			if(to_json_type1_inner_block(file_tags, json_inner_file_tag))		// записать все тэги, характеризующие дайком файл во внутренний блок
			{
				json *refInnerJson = &json_dicom_files;
				for(size_t i = NFIELDS_TYPE_1; i-- > 0;)			// в json файл заполнить поля об первых 5 тэгах, характеризующих исследование
				{
					string str_tag_value = convert_to_string8(file_tags.get_dicom_tags_value(i));
					refInnerJson = &(*refInnerJson)[str_tag_value];
				}
				refInnerJson->push_back(std::move(json_inner_file_tag));
			}
		}
		else
		{
			// заполнить информацию о не dicom файлах
			json json_file_tag;
			// записать все тэги, характеризующие дайком файл, во внутренний блок
			if(to_json_type1_inner_block(file_tags, json_file_tag))
				json_not_dicom.push_back(std::move(json_file_tag));
		}
	}

	// если не пусто, записать древовидную информацию dicom файлав в тэг "dicomlist"
	if (!json_dicom_files.empty())
		json_type1["dicomlist"] = std::move(json_dicom_files);

	// записать "non Dicom" файлы, это идентично записи в json type 2
	if (!json_not_dicom.empty())
		json_type1["non Dicom"] = std::move(json_not_dicom);
}



// из списка тэгов для каждого файла сгенерировать json файл type2 (sample2.json)
void dir_info_to_json_type2(const SingleDirectoryIndex& dcmDirectoryIndex, json& json_type2)
{
	if (!dcmDirectoryIndex.size())
		return;

	// заполнить обязательные поля о "ID", "version", "type" из map_header_json_type1
	for (auto& map_v : map_header_json_type2)
	{
		json_type2[map_v.first] = map_v.second;
	}

	// Для определенности сохраняем элементы в порядке "возрастания" имени файла с учетом регистра
	// в utf-8.
	map<string, const DicomFileIndex*> ordered_index;
	for (auto &di: dcmDirectoryIndex)
	{
		ordered_index[convert_to_string8(di.get_file_name())] = &di;
	}

	for (auto &odi: ordered_index)
	{
		auto &dicom_file_Tags = *odi.second;
		json json_file_tag;
		if(to_json_type2(dicom_file_Tags, json_file_tag))	 // записать все тэги, характеризующие дайком файл
		{
			json_type2["filelist"].push_back(std::move(json_file_tag));
		}
	}
}



// загрузить json файл с деревом исследований Dicom файлов type 1
void load_json_type1_tree(SingleDirectoryIndex& result, json& json_type1_dicom_section)
{
	XRAD_ASSERT_THROW(!json_type1_dicom_section.is_null());   // проверить, есть ли информация в json объекте

	vector<string> reference_string;
	vector<vector<string>> flatten_result;
	vector<json> vec_json_dicom;
//#error
	json_parse_type1(reference_string, json_type1_dicom_section, flatten_result, vec_json_dicom);

	// для считывания этой древовидной структуры
	XRAD_ASSERT_THROW (flatten_result.size() != 0)  // если парсинг дал пустой результат

	for (size_t i = 0; i < flatten_result.size(); i++)
	{
		vector<string> fields1 = flatten_result[i];
		if (fields1.size() != NFIELDS_TYPE_1)
			break;
		DicomFileIndex fileindex_base;
		for (size_t i = NFIELDS_TYPE_1; i-- > 0;)			// в json файл заполнить поля об первых 5 тэгах, характеризующих исследование
		{
			//Dicom::tag_t tag_value = fileindex_base.m_dicom_tags[i];
			fileindex_base.set_dicom_tags_value(i, convert_to_wstring(fields1[NFIELDS_TYPE_1 - 1 - i]));
		}

		for (auto& el : vec_json_dicom[i])  // для каждого элемента array
		{
			try
			{
				DicomFileIndex fileindex = fileindex_base;
				fileindex.append(from_json_get_file_index(el, index_file_type::hierarchical));
				if(fileindex.is_dicom()) // взять остальные поля
				{
					result.add_file_index(std::move(fileindex));  // после функции move объект fileindex_from_json уже не хранит информации
				}
			}
			catch(...)
			{
			}
		}

	}
}

} // namespace



// загрузить json файл, содержащий информацию об файлах, в структуру SingleDirectoryIndex& dcmDirectoryIndex
SingleDirectoryIndex load_parse_json(const wstring& json_fname)
{
	json json_from_file = load_json(json_fname);

	for (auto& map_v : map_header_json_type1)  // проверить наличие обязательных полей = { "ID", "version", "type" };
	{
		XRAD_ASSERT_THROW(json_from_file.find(map_v.first) != json_from_file.end());
	}

	// проверить тип json файла
	index_file_type json_type = interpret_index_file_type(json_from_file["type"]);

	XRAD_ASSERT_THROW(json_type != index_file_type::unknown);

	SingleDirectoryIndex result(RemoveTrailingPathSeparator(file_path(json_fname)));
	if (json_type == index_file_type::hierarchical)
	{
		// json файл тип type 1
		// Парсинг type 1 нужно переписать полностью. Теперешний код, хотя и работает, содержит много натянутых зависимостей, которые усложняют сопровождение
		//Error("json type 1 index is not tested");

		XRAD_ASSERT_THROW(json_from_file.find("dicomlist") != json_from_file.end());		// если обязательное поле "dicom" отсутствует, то выйти
		XRAD_ASSERT_THROW(!json_from_file["dicomlist"].is_null());
		json& json_dicom_files = json_from_file.at("dicomlist");				// хранение информации о файлах в древовидной структуре первых 5-тэгов
		
		load_json_type1_tree(result, json_dicom_files);
	}

	// json файл тип type 1 или 2

	string tag_name = json_type == index_file_type::hierarchical ? "non Dicom" : "filelist";
	XRAD_ASSERT_THROW (json_from_file.find(tag_name) != json_from_file.end())		 // если обязательное поле отсутствует, то выйти

	const json& json_filelist = *json_from_file.find(tag_name);
	// из поле "filelist", считать всю инф-цию о dicom файлах
	for (size_t i = 0; i < json_filelist.size(); i++)
	{
		try
		{
			DicomFileIndex fileindex = from_json_get_file_index(json_filelist[i], index_file_type::raw);
			result.add_file_index(fileindex);
		}
		catch(...)
		{
		}
	}

	return result;
}



// записать json файл,
wstring save_to_jsons(const SingleDirectoryIndex& dcmDirectoryIndex, index_file_type json_type)
{

	// для каждого уникального кластера (директории с dicom файлами) сформировать json объект
	// json объект формируется согласно примерам образцов json файлов: sample1.json и sample2.json
	wstring wstr_json_fname;
	if (!dcmDirectoryIndex.size())
		return wstr_json_fname; // TODO: Сделать удаление файла, если нет индексируемого содержимого.

	json json_to_save;
	switch (json_type)
	{
		case index_file_type::hierarchical:
			dir_info_to_json_type1(dcmDirectoryIndex, json_to_save);
			wstr_json_fname = dcmDirectoryIndex.get_path() + wpath_separator() +
					index_filename_type1();
			break;
		case index_file_type::raw:
			dir_info_to_json_type2(dcmDirectoryIndex, json_to_save);
			wstr_json_fname = dcmDirectoryIndex.get_path() + wpath_separator() +
					index_filename_type2();
			break;
		default:
			return wstr_json_fname;
	}
	// сохранить dicom файл  в текущей уникальной директории
	//ShowText(L"Информация о записи dicom файла", ssprintf("%d число тэгов \nимя файла: \n%s", json_to_save.size(),
	//	convert_to_string(wstr_json_fname).c_str(), false)  );

	save_json(json_to_save, wstr_json_fname);
	return wstr_json_fname;
}



// проверить адекватность записи/чтения инф-ции в/из json файлов в двух форматах
bool test_write_load_json(SingleDirectoryIndex& dcmDirectoryIndex)
{
	// если некуда писать и откуда читать
	if (dcmDirectoryIndex.get_path().size() == 0)
		return true;

	try
	{
		wstring wstr_json_fname1 = save_to_jsons(dcmDirectoryIndex, index_file_type::hierarchical);
		SingleDirectoryIndex dir_index_from_json1 = load_parse_json(wstr_json_fname1); // если проблемы чтения json файла
		
		wstring wstr_json_fname2 = save_to_jsons(dcmDirectoryIndex, index_file_type::raw);
		SingleDirectoryIndex dir_index_from_json2 = load_parse_json(wstr_json_fname2); // если проблемы чтения json файла

		bool is_equal01 = dcmDirectoryIndex == dir_index_from_json1;
		bool is_equal02 = dcmDirectoryIndex == dir_index_from_json2;
		bool is_equal12 = dir_index_from_json1 == dir_index_from_json2;
		bool res = is_equal01 && is_equal02 && is_equal12;
		//if (!res)
		//	ShowText(L"1", ssprintf("writeload res = %d\n file %s", res, convert_to_string(wstr_json_fname1).c_str()));

		return res;
	}
	catch(...)
	{
		return false;
	}
}



} // namespace Dicom

XRAD_END
