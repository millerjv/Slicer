/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkImageAccumulateDiscrete.h,v $
  Date:      $Date: 2006/04/12 22:28:50 $
  Version:   $Revision: 1.21 $

=========================================================================auto=*/

#ifndef __vtkImageAccumulateDiscrete_h
#define __vtkImageAccumulateDiscrete_h

// MRML includes
#include "vtkMRML.h"

// VTK includes
#include <vtkImageToImageFilter.h>

/// \brief Histogram the first component the scalars on the image.
///
/// This filter defines a histogram on the first component of the scalars of an image.  The histogram is the count of the
/// number of voxels that falls within the intensity range of each bin in the histogram.
///
///  The bins of the histogram are defined by the dynamic range of the data. The location of the first bin is stored
///  in the origin of the output data.  This is the center of the bin. The size of the bins is stored in the spacing of the output data.
///  The number of  bins can be controlled by the user and defaults to 100.
class VTK_MRML_EXPORT vtkImageAccumulateDiscrete : public vtkImageToImageFilter
{
public:
  static vtkImageAccumulateDiscrete *New();
  vtkTypeMacro(vtkImageAccumulateDiscrete,vtkImageToImageFilter);
  void PrintSelf(ostream& os, vtkIndent indent);

  /// Set/Get the number of bins to use on the histogram.  There is rarely a need for more than 100 bins in any
  /// histogram (the default). A good rule of thumb in defined in Freedman and Diaconis (summarized in Izenman 1991),
  ///
  ///     binWidth = 2 * (IQR)*N^(-1/3)
  ///
  /// where
  ///
  ///      N is the number of samples,
  ///      IQR = invcdf(0.75) - invcdf(0.25) in the real data
  ///
  /// giving the number of bins as
  ///
  ///      NumberOfBins = DynamicRange  / binWidth
  ///
  ///  E.g. A single slice of a CT image (512x512) has a dynamic range of about 2500 HU. The IQR for the chest/abdomen
  ///  yields a binWidth around 27 HU.  Thus, the number of bins should be around 100.
  vtkSetMacro(NumberOfBins, unsigned int);
  vtkGetMacro(NumberOfBins, unsigned int);


protected:
  vtkImageAccumulateDiscrete();
  ~vtkImageAccumulateDiscrete() {};

  void ExecuteInformation(vtkImageData *input, vtkImageData *output);
  void ExecuteInformation(){this->Superclass::ExecuteInformation();};
  void ComputeInputUpdateExtent(int inExt[6], int outExt[6]);
  void ExecuteData(vtkDataObject *);

private:
  vtkImageAccumulateDiscrete(const vtkImageAccumulateDiscrete&);
  void operator=(const vtkImageAccumulateDiscrete&);

  unsigned int NumberOfBins;
};

#endif

