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

#ifndef itkChangeRegionImageFilter_hxx
#define itkChangeRegionImageFilter_hxx

#include "itkChangeRegionImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkContinuousIndex.h"
#include "itkObjectFactory.h"

namespace itk
{

template< typename ImageType >
ChangeRegionImageFilter< ImageType >
::ChangeRegionImageFilter()
{
    m_ChangeOrigin = false;
    m_ChangeBufferedRegion = false;
    m_ChangeLargestPossibleRegion = false;
}

template<typename ImageType >
void
ChangeRegionImageFilter< ImageType >
::SetOrigin(typename ImageType::PointType origin)
{
    m_Origin = origin;
    m_ChangeOrigin = true;
}

template<typename ImageType >
void
ChangeRegionImageFilter< ImageType >
::SetBufferedRegion(typename ImageType::RegionType region)
{
    m_BufferedRegion = region;
    m_ChangeBufferedRegion = true;
}

template<typename ImageType >
void
ChangeRegionImageFilter< ImageType >
::SetLargestPossibleRegion(typename ImageType::RegionType region)
{
    m_LargestPossibleRegion = region;
    m_ChangeLargestPossibleRegion = true;
}

template< typename ImageType >
void
ChangeRegionImageFilter< ImageType >
::GenerateOutputInformation()
{
  Superclass::GenerateOutputInformation();

  // Get pointers to the input and output
  typename Superclass::OutputImagePointer output = this->GetOutput();
  typename Superclass::InputImagePointer input = const_cast< ImageType * >( this->GetInput() );

  // Default is to copy input's information
  output->CopyInformation(input);

  if (m_ChangeOrigin)
      output->SetOrigin(m_Origin);

  if (m_ChangeBufferedRegion) {
      output->SetBufferedRegion(m_BufferedRegion);
      output->SetRequestedRegion(m_BufferedRegion);
  }

  if (m_ChangeLargestPossibleRegion)
      output->SetLargestPossibleRegion(m_LargestPossibleRegion);
}

template< typename ImageType >
void
ChangeRegionImageFilter< ImageType >
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();

  ImageType* input = const_cast< ImageType * >( this->GetInput() );

  input->SetRequestedRegion(input->GetBufferedRegion());
}

template< typename ImageType >
void
ChangeRegionImageFilter< ImageType >
::GenerateData()
{
  ImageType* input = const_cast< ImageType * >( this->GetInput() );
  ImageType* output = this->GetOutput();

  // No need to copy the bulk data
  output->SetPixelContainer( input->GetPixelContainer() );
}

template< typename ImageType >
void
ChangeRegionImageFilter< ImageType >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Origin: " << m_Origin << std::endl;
  os << indent << "BufferedRegion: " << m_BufferedRegion << std::endl;
  os << indent << "LargestPossibleRegion: " << m_LargestPossibleRegion << std::endl;
}

} // end namespace itk

#endif
