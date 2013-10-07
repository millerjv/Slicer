/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkImageBimodalAnalysis.h,v $
  Date:      $Date: 2006/06/14 20:44:13 $
  Version:   $Revision: 1.23 $

=========================================================================auto=*/

#ifndef __vtkImageBimodalAnalysis_h
#define __vtkImageBimodalAnalysis_h

// MRML includes
#include "vtkMRML.h"

// VTK includes
#include <vtkImageToImageFilter.h>

#define VTK_BIMODAL_MODALITY_CT 0
#define VTK_BIMODAL_MODALITY_MR 1

/// \brief Analysis bimodal histograms
///
/// This filter assumes the input comes
/// from vtkImageAccumulateDiscrete, so there.
/// \warning
/// FIXME: only works on output floating point
/// FIXME: should use vtkTemplateMacro
class VTK_MRML_EXPORT vtkImageBimodalAnalysis : public vtkImageToImageFilter
{
public:
  static vtkImageBimodalAnalysis *New();
  vtkTypeMacro(vtkImageBimodalAnalysis,vtkImageToImageFilter);
  void PrintSelf(ostream& os, vtkIndent indent);

  ///
  /// Set the type of data, if known
  vtkSetMacro(Modality, int);
  vtkGetMacro(Modality, int);
  void SetModalityToMR() {this->SetModality(VTK_BIMODAL_MODALITY_MR);};
  void SetModalityToCT() {this->SetModality(VTK_BIMODAL_MODALITY_CT);};

  ///
  /// Get stats
  vtkGetMacro(Threshold, double);
  vtkGetMacro(Window, double);
  vtkGetMacro(Level, double);
  vtkGetMacro(Min, double);
  vtkGetMacro(Max, double);

  ///
  /// Set stats
  vtkSetMacro(Threshold, double);
  vtkSetMacro(Window, double);
  vtkSetMacro(Level, double);
  vtkSetMacro(Min, double);
  vtkSetMacro(Max, double);


protected:
  vtkImageBimodalAnalysis();
  ~vtkImageBimodalAnalysis() {};

  int Modality;

  double Threshold;
  double Window;
  double Level;
  double Min;
  double Max;


  void ExecuteInformation(vtkImageData *input, vtkImageData *output);
  void ExecuteInformation(){this->Superclass::ExecuteInformation();};
  void ExecuteData(vtkDataObject *);

private:
  vtkImageBimodalAnalysis(const vtkImageBimodalAnalysis&);
  void operator=(const vtkImageBimodalAnalysis&);
};

#endif

