/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkImageBimodalAnalysis.cxx,v $
  Date:      $Date: 2006/06/14 20:44:14 $
  Version:   $Revision: 1.21 $

=========================================================================auto=*/
#include "vtkImageBimodalAnalysis.h"

#include "vtkObjectFactory.h"
#include "vtkImageData.h"

//#include <math.h>
//#include <cstdlib>

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkImageBimodalAnalysis);

//----------------------------------------------------------------------------
// Constructor sets default values
vtkImageBimodalAnalysis::vtkImageBimodalAnalysis()
{
  this->Modality  = VTK_BIMODAL_MODALITY_CT;
  this->Threshold = 0;
  this->Window    = 0;
  this->Level     = 0;
  this->Min       = 0;
  this->Max       = 0;
}

//----------------------------------------------------------------------------
void vtkImageBimodalAnalysis::ExecuteInformation(vtkImageData *, vtkImageData *output)
{
  output->SetScalarType(VTK_FLOAT);
}

//----------------------------------------------------------------------------
// This templated function executes the filter for any type of data.
template <class T>
static void vtkImageBimodalAnalysisExecute(vtkImageBimodalAnalysis *self,
                      vtkImageData *inData, T *inPtr,
                      vtkImageData *outData, float *outPtr)
{
  int x, k;
  int min0, max0, min1, max1, min2, max2;
  int noise = 1, width = 5;

  T tmp;
  vtkFloatingPointType sum, wsum;
  int ct = (self->GetModality() == VTK_BIMODAL_MODALITY_CT) ? 1 : 0;
  double centroid, noiseCentroid, trough, window, threshold, min, max;
  vtkFloatingPointType origin[3], spacing[3];

  // inData->GetExtent(min0, max0, min1, max1, min2, max2);
  // for (unsigned int ii = min0; ii <=max0; ++ii)
  //   std::cout << "[" << ii << "] = " << inPtr[ii]  << std::endl;

  // Process x dimension only
  outData->GetExtent(min0, max0, min1, max1, min2, max2);
  inData->GetOrigin(origin);
  inData->GetSpacing(spacing);


  // Zero output
  memset((void *)outPtr, 0, (max0-min0+1)*sizeof(float));

  // For CT data, ignore -32768 (min of data is lower bound of histogram, origin - spacing/2)
  //
  // Disabling for now.  The value -32768 is sometimes used to fill the corners of a CT image.  But if don't
  // really know it is CT, then this does more harm than good.
  //
  // if (ct && (fabs(origin[0] - spacing[0]/2.0 + 32768) < 1e-4))
  //   {
  //   min0 = 1;
  //   }

  // Find min (first non-zero value in histogram)
  min = x = min0;
  while (!inPtr[x] && x <= max0)
    {
    x++;
    }
  if (x <= max0)
    {
    min = x;
    }

  // Find max (last non-zero value in histogram)
  max = x = max0;
  while (!inPtr[x] && x >= min0)
    {
    x--;
    }
  if (x >= min0)
    {
    max = x;
    }

  // Smooth
  long cnt;
  for (x = min; x <= max; x++)
    {
    cnt = 0;
    for (k=0; k < width; k++)
      {
        if (x+k < max0)
        {
          outPtr[x] += (float)inPtr[x+k];
          ++cnt;
        }
      }
    if (cnt > 0)
      {
      outPtr[x] /= (double) cnt;
      }
    }

  // Find first trough of smoothed histogram
  x = min;
  trough = min-1;
  noise = 1;
  while (x < max && trough < min)
    {
    if (noise)
      {
      if (outPtr[x] > outPtr[x+1])
        {
        if (x > min)
          {
          noise = 0;
          }
        }
      }
    else
      {
      if (outPtr[x] < outPtr[x+1])
        {
        trough = x;
        }
      }
    x++;
    }

  // Compute centroid of the histogram that PRECEEDS the trough
  // (noise lobe)
  wsum = sum = 0;
  for (x=min; x <= trough; x++)
    {
    tmp = inPtr[x];
    wsum += (vtkFloatingPointType)x*tmp;
    sum  += (vtkFloatingPointType)  tmp;
    }
  if (sum)
    {
    noiseCentroid = (int)(wsum / sum);
    }
  else
    {
    noiseCentroid = trough;
    }

  // Compute centroid of the histogram that FOLLOWS the trough
  // (signal lobe, and not noise lobe)
  wsum = sum = 0;
  for (x=trough; x <= max; x++)
    {
    tmp = inPtr[x];
    wsum += (vtkFloatingPointType)x*tmp;
    sum  += (vtkFloatingPointType)  tmp;
    }
  if (sum)
    {
    centroid = (int)(wsum / sum);
    }
  else
    {
    centroid = trough;
    }

  // Threshold
  threshold = trough;

  // Compute the window as twice the width as the smaller half
  // of the signal lobe
  if (centroid - noiseCentroid < max - centroid)
    {
    window = (centroid - noiseCentroid)*2;
    }
  else
    {
    window = (max - centroid)*2;
    }

  // Record findings (in original intensities, not bins)
  self->SetThreshold(origin[0] + spacing[0]*threshold);
  self->SetMin(origin[0] + min*spacing[0] - spacing[0]/2.0);      // lower bound of first bin
  self->SetMax(origin[0] + max*spacing[0] + spacing[0]/2.0);    // upper bound of first bin
  self->SetLevel(origin[0] + centroid*spacing[0]);
  self->SetWindow(spacing[0]*window);

  //std::cout << "min bin, max bin = " << min << ", " << max << std::endl;
  //self->Print(std::cout);
}



//----------------------------------------------------------------------------
// This method is passed a input and output Data, and executes the filter
// algorithm to fill the output from the input.
// It just executes a switch statement to call the correct function for
// the Datas data types.
void vtkImageBimodalAnalysis::ExecuteData(vtkDataObject *out)
{
  vtkImageData *outData = vtkImageData::SafeDownCast(out);
  vtkImageData *inData = this->GetInput();
  void *inPtr;
  float *outPtr;

  outData->SetExtent(outData->GetWholeExtent());
  outData->AllocateScalars();

  inPtr  = inData->GetScalarPointer();
  outPtr = (float *)outData->GetScalarPointer();

  // Components turned into x, y and z
  int c = inData->GetNumberOfScalarComponents();
  if (c > 1)
    {
    vtkErrorMacro("This filter requires 1 scalar component, not " << c);
    return;
    }

  // this filter expects that output is type float.
  if (outData->GetScalarType() != VTK_FLOAT)
    {
    vtkErrorMacro(<< "ExecuteData: out ScalarType " << outData->GetScalarType()
      << " must be float\n");
    return;
    }
  switch (inData->GetScalarType())
  {
    case VTK_CHAR:
      vtkImageBimodalAnalysisExecute(this,
              inData, (char *)(inPtr), outData, outPtr);
      break;
    case VTK_UNSIGNED_CHAR:
      vtkImageBimodalAnalysisExecute(this,
              inData, (unsigned char *)(inPtr), outData, outPtr);
      break;
    case VTK_SHORT:
      vtkImageBimodalAnalysisExecute(this,
              inData, (short *)(inPtr), outData, outPtr);
      break;
    case VTK_UNSIGNED_SHORT:
      vtkImageBimodalAnalysisExecute(this,
              inData, (unsigned short *)(inPtr), outData, outPtr);
      break;
    case VTK_INT:
      vtkImageBimodalAnalysisExecute(this,
              inData, (int *)(inPtr), outData, outPtr);
      break;
    case VTK_UNSIGNED_INT:
      vtkImageBimodalAnalysisExecute(this,
              inData, (unsigned int *)(inPtr), outData, outPtr);
      break;
    case VTK_LONG:
      vtkImageBimodalAnalysisExecute(this,
              inData, (long *)(inPtr), outData, outPtr);
      break;
    case VTK_UNSIGNED_LONG:
      vtkImageBimodalAnalysisExecute(this,
              inData, (unsigned long *)(inPtr), outData, outPtr);
      break;
    case VTK_FLOAT:
      vtkImageBimodalAnalysisExecute(this,
              inData, (float *)(inPtr), outData, outPtr);
      break;
    case VTK_DOUBLE:
      vtkImageBimodalAnalysisExecute(this,
              inData, (double *)(inPtr), outData, outPtr);
      break;
    default:
      vtkErrorMacro(<< "ExecuteData: Unsupported ScalarType");
      return;
    }
}

//----------------------------------------------------------------------------
void vtkImageBimodalAnalysis::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Modality: " << this->Modality << " (" << (this->Modality == VTK_BIMODAL_MODALITY_CT ? "CT" : "MR") << ")\n";
  os << indent << "Threshold: " << this->Threshold << "\n";
  os << indent << "Window: " << this->Window << "\n";
  os << indent << "Level: " << this->Level << "\n";
  os << indent << "Min: " << this->Min << "\n";
  os << indent << "Max: " << this->Max << "\n";

}

