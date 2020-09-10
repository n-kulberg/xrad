﻿/*!
	\file
	\date 10/09/2020 12:19:45 PM
	\author sbp
*/
#ifndef XRayAcquisition_h__
#define XRayAcquisition_h__

#include "XRAYInstance.h"
#include "GenericImageAcquisition.h"
#include <XRADDicom/Sources/DicomClasses/Instances/LoadGenericClasses.h>

XRAD_BEGIN

class XRayAcquisition : public GenericImageAcquisition
{
public:
	XRayAcquisition(const shared_ptr<Dicom::acquisition_loader> &acquisition_loader_p);
	XRayAcquisition(const size_t elements_amount); //note этот конструктор только для создания новых сборок

	//gets
/*	virtual std::string classname() const { return "XRayAcquisition"; }

	virtual size_t n_elements() const { return m_acquisition_loader->size(); }
	void put_elements_to_instance(Dicom::instance &instance, size_t num_element) const;

	RealFunction2D_F32 get_image(size_t no) const;
	index_vector	sizes() const
	{
		Dicom::image &first_slice = dynamic_cast<Dicom::image&>(*(m_acquisition_loader->front()));
		return{ (*m_acquisition_loader).size(), first_slice.vsize(), first_slice.hsize() };
	}


	RealFunction2D_F32	slice(size_t pos) const;		// получить 2D изображение среза с номером pos
	vector<RealFunction2D_F32>	slices() const;				// получить вектор 2D изображений
	*/
};

XRAD_END

#endif // XRAYAcquisition_h__