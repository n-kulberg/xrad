﻿/*!
	* \file XRAYInstance.cpp
	* \date 4/19/2018 5:50:09 PM
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
#include "XRAYInstance.h"


XRAD_BEGIN

void XRAYInstance::set_image(const RealFunction2D_F32 &img)
{
	//m_image.realloc(img.vsize(), img.hsize());
	//m_image.CopyData(img);
	m_image.MakeCopy(img);
}

XRAD_END
