/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#ifndef itkChangeRegionImageFilter_h
#define itkChangeRegionImageFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{

template< typename ImageType >
class ChangeRegionImageFilter:
  public ImageToImageFilter< ImageType, ImageType >
{
public:
  /** Standard class typedefs. */
  typedef ChangeRegionImageFilter                    Self;
  typedef ImageToImageFilter< ImageType, ImageType > Superclass;
  typedef SmartPointer< Self >                       Pointer;
  typedef SmartPointer< const Self >                 ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Image related typedefs. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      ImageType::ImageDimension);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ChangeRegionImageFilter, ImageToImageFilter);

  itkGetConstReferenceMacro(Origin, typename ImageType::PointType);
  itkGetConstReferenceMacro(BufferedRegion, typename ImageType::RegionType);
  itkGetConstReferenceMacro(LargestPossibleRegion, typename ImageType::RegionType);

  void SetBufferedRegion(typename ImageType::RegionType region);
  void SetLargestPossibleRegion(typename ImageType::RegionType region);
  void SetOrigin(typename ImageType::PointType origin);

  virtual void GenerateOutputInformation() ITK_OVERRIDE;
  virtual void GenerateInputRequestedRegion() ITK_OVERRIDE;

  /** Copy the input buffer. */
  void GenerateData() ITK_OVERRIDE;

protected:
  ChangeRegionImageFilter();
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;
  virtual void VerifyInputInformation() ITK_OVERRIDE {}

private:
  ChangeRegionImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &);             // purposely not implemented

  bool m_ChangeOrigin;
  bool m_ChangeBufferedRegion;
  bool m_ChangeLargestPossibleRegion;

  typename ImageType::PointType  m_Origin;
  typename ImageType::RegionType m_BufferedRegion;
  typename ImageType::RegionType m_LargestPossibleRegion;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkChangeRegionImageFilter.hxx"
#endif

#endif
