﻿#ifndef RadonTransform_h__
#define RadonTransform_h__

/********************************************************************
	created:	2016/09/13
	created:	13:9:2016   13:35
	filename: 	q:\programs\CTDensityArtifact\sources\RadonTransform.h
	file path:	q:\programs\CTDensityArtifact\sources
	file base:	RadonTransform
	file ext:	h
	author:		kns
	
	purpose:	
*********************************************************************/


XRAD_BEGIN

void	InitRadonTransform();
void RadonTransformReverse(RealFunction2D_F32 &generated_data, const RealFunction2D_F32 &radon_data, size_t enlarge_factor, ProgressProxy pp);
void RadonTransformForward(RealFunction2D_F32 &radon_data, const RealFunction2D_F32 &original_data, size_t enlarge_factor, ProgressProxy pp);

void MakeIsotropic(RealFunction2D_F32& original);


XRAD_END

#endif // RadonTransform_h__
