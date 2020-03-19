﻿#include "pre.h"
#include <XRADDicom/DicomClasses/Indexing/DicomFileIndex.h>
#include <XRADDicom/DicomClasses/Indexing/DicomDirectoryIndex.h>
#include <XRADDicom/XRADDicom.h>
#include <typeinfo>


#include <XRADDicom/DicomClasses/instances/ct_slice.h>
#include <XRADDicom/DicomClasses/instances/xray_image.h>
#include <XRADDicom/DicomClasses/instances/mr_slice.h>
#include <XRADDicom/DicomClasses/instances/mr_slice_siemens.h>


/*!
\file
\date 2019/09/26 14:00
\author novik

\brief 	Имплементация функционала DicomDirectoryIndex - структуры для работы с файлами в одной директории

Функции работы с json файлом вынесены в файл DicomDirectoryIndexJson.cpp
*/

XRAD_BEGIN

namespace Dicom
{

bool DicomDirectoryIndex::is_equal(const DicomDirectoryIndex& a) const
{
	for (const auto& el1 : m_FilesIndex)
	{
		bool is_equal_find = false;
		for (const auto& el2 : a.m_FilesIndex)
		{
			if (el1 == el2)
			{
				is_equal_find = true;
				break;
			}
		}
		if (!is_equal_find)  // если в "a" не найдены одинаковые елементы "el1"
			return false;
	}
	return true;
}

void	DicomDirectoryIndex::add_file_index(const DicomFileIndex& dcmFileIndex)
{
	m_FilesIndex.push_back(dcmFileIndex);
}

void DicomDirectoryIndex::update()
{
	for (auto& el : m_FilesIndex)
	{
		if (!el.get_isneed_indexing())
			continue;
		DicomFileIndex current_file_tags;
		if (current_file_tags.fill_filetags_from_file(get_path(),
				convert_to_wstring(el.get_file_name())))
		{
			el = current_file_tags;			// обновили тэги
			el.set_need_indexing(false);    // сняли метку необходимости индексации
		}
	}
	// удалить индексаторы о файлах, по которым не удалось обновить информацию
	auto predicate = [](const DicomFileIndex &v) { return v.get_isneed_indexing(); };
	m_FilesIndex.erase(remove_if(m_FilesIndex.begin(), m_FilesIndex.end(), predicate), m_FilesIndex.end());
}

bool DicomDirectoryIndex::fill_from_fileinfo(const wstring &path,
		const vector<FileInfo>& file_infos)
{
	m_path = path;
	for (auto el : file_infos)
	{
		if (!may_be_dicom_filename(el.filename))
			continue;
		if (el.filename == wstring(j_name()) + L"1." + j_extension())
		{
			m_filename_json_1 = el.filename;
			continue;
		}
		if (el.filename == wstring(j_name()) + L"2." + j_extension())
		{
			m_filename_json_2 = el.filename;
			continue;
		}

		DicomFileIndex current_file_tags;
		if (current_file_tags.fill_from_fileinfo(el))
		{
			m_FilesIndex.push_back(std::move(current_file_tags));
			// после функции move объект current_file_tags уже не хранит информации
		}
	}
	return m_FilesIndex.size() || !m_filename_json_1.empty() || !!m_filename_json_2.empty();
}

} //namespace Dicom

XRAD_END
