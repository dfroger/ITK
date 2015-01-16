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
#include "itkHDF5ImageIOFactory.h"
#include "itkVectorImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkHDF5ImageIO.h"

/* Write a VectorImage of size (SZ SY,SX,SC) in a HDF5 file, then paste in
 * this file a region (allocated and filled in memory) at (IZ,IY,IX,0) of
 * size (NZ,NY,NX,SC). */

using namespace H5;

typedef float PixelType;
const int Dimension = 3;
typedef itk::VectorImage<PixelType, Dimension> ImageType;

const unsigned int SZ = 5;
const unsigned int SY = 5;
const unsigned int SX = 5;
const unsigned int SC = 2;

const unsigned int IZ = 1;
const unsigned int IY = 1;
const unsigned int IX = 1;

const unsigned int NZ = 2;
const unsigned int NY = 2;
const unsigned int NX = 2;

const std::string FILENAME = "image.hdf5";

PixelType
pixelValue(unsigned int iz, unsigned int iy, unsigned int ix, unsigned int ic)
{
  return iz*1000 + iy*100 + ix*10 + ic;
}

PixelType
pastedPixelValue(unsigned int iz, unsigned int iy, unsigned int ix, unsigned int ic)
{
  return 990000 + iz*1000 + iy*100 + ix*10 + ic;
}

bool
isPasted(unsigned int iz, unsigned int iy, unsigned int ix)
{
    return (IZ <= iz && iz < IZ+NZ &&
            IY <= iy && iy < IY+NY &&
            IX <= ix && ix < IX+NX);
}

PixelType
expectedPixelValue(unsigned int iz, unsigned int iy, unsigned int ix, unsigned int ic)
{
    if ( isPasted(iz,iy,ix) )
        return pastedPixelValue(iz, iy, ix, ic);
    else
        return pixelValue(iz, iy, ix, ic);
}

void
writeImage()
{
  // Allocate image.
  ImageType::IndexType index;
  index.Fill(0);
  ImageType::SizeType size;
  size[2] = SZ;
  size[1] = SY;
  size[0] = SX;
  ImageType::RegionType region(index,size);
  ImageType::Pointer image = ImageType::New();
  image->SetRegions(region);
  image->SetVectorLength(SC);
  image->Allocate();

  // Set image pixel values.
  typedef itk::VariableLengthVector<PixelType> VariableVectorType;
  VariableVectorType value;
  value.SetSize(SC);
  itk::ImageRegionIterator<ImageType> it(image,image->GetLargestPossibleRegion());
  for(it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
      index = it.GetIndex();
      for (unsigned int ic=0; ic < SC; ++ic)
        {
        value[ic] = pixelValue(index[2], index[1], index[0], ic);
        }
      image->SetPixel(index, value);
    }

  // Create writer.
  typedef itk::ImageFileWriter<ImageType>  WriterType;
  typedef itk::HDF5ImageIO                 ImageIOType;
  WriterType::Pointer writer = WriterType::New();
  ImageIOType::Pointer hdf5IO = ImageIOType::New();
  writer->SetImageIO(hdf5IO);
  writer->SetFileName(FILENAME);

  // Update all the pipeflow.
  writer->SetInput( image );
  writer->Update();
}

void
pasteImage()
{
   unsigned int sz,sy,sx,sc;

  // Destruct reader, so underlying HDF5 file is closed, and can by reopen by
  // writer.
  {
  // Create reader.
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::HDF5ImageIO                   ImageIOType;
  ReaderType::Pointer reader = ReaderType::New();
  ImageIOType::Pointer hdf5IO = ImageIOType::New();
  reader->SetImageIO(hdf5IO);
  reader->SetFileName(FILENAME);

  // Read image size only.
  reader->SetFileName(FILENAME);
  reader->UpdateOutputInformation();
  sz = reader->GetImageIO()->GetDimensions(2);
  sy = reader->GetImageIO()->GetDimensions(1);
  sx = reader->GetImageIO()->GetDimensions(0);
  sc = reader->GetImageIO()->GetNumberOfComponents();
  }

  // Allocate sub-region of the existing image we want to override.
  // Buffered region.
  ImageType::IndexType bindex;
  bindex[0] = IX;
  bindex[1] = IY;
  bindex[2] = IZ;
  ImageType::SizeType bsize;
  bsize[0] = NX;
  bsize[1] = NY;
  bsize[2] = NZ;
  ImageType::RegionType bregion(bindex,bsize);
  // Largest possible region.
  ImageType::IndexType lindex;
  lindex[0] = 0;
  lindex[1] = 0;
  lindex[2] = 0;
  ImageType::SizeType lsize;
  lsize[0] = sx;
  lsize[1] = sy;
  lsize[2] = sz;
  ImageType::RegionType lregion(lindex,lsize);
  // Allocate
  ImageType::Pointer image = ImageType::New();
  image->SetBufferedRegion(bregion);
  image->SetLargestPossibleRegion(lregion);
  image->SetVectorLength(sc);
  image->Allocate();

  // Set sub-region pixel
  ImageType::IndexType index;
  typedef itk::VariableLengthVector<PixelType> VariableVectorType;
  VariableVectorType value;
  value.SetSize(sc);
  for (unsigned int iz=IZ; iz < IZ+NZ; ++iz)
    {
    index[2] = iz;
    for (unsigned int iy=IY; iy < IY+NY; ++iy)
      {
      index[1] = iy;
      for (unsigned int ix=IX; ix < IX+NX; ++ix)
        {
        index[0] = ix;
        for (unsigned int ic=0; ic < sc; ++ic)
          {
          value[ic] = pastedPixelValue(iz,iy,ix,ic);
          }
        image->SetPixel(index, value);
        }
      }
    }

  // IO region
  itk::ImageIORegion ioregion(Dimension);
  ioregion.SetIndex(0,IX);
  ioregion.SetIndex(1,IY);
  ioregion.SetIndex(2,IZ);
  ioregion.SetSize(0,NX);
  ioregion.SetSize(1,NY);
  ioregion.SetSize(2,NZ);

  // Create writer.
  typedef itk::ImageFileWriter<ImageType>  WriterType;
  typedef itk::HDF5ImageIO                 ImageIOType;
  WriterType::Pointer writer = WriterType::New();
  ImageIOType::Pointer writer_hdf5IO = ImageIOType::New();
  writer->SetImageIO(writer_hdf5IO);
  writer->SetFileName(FILENAME);
  writer->SetIORegion(ioregion);

  // Update all the pipeflow
  writer->SetInput( image );
  writer->Update();
}

int
checkImage()
{
  // Read image.
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::HDF5ImageIO                ImageIOType;
  ReaderType::Pointer reader = ReaderType::New();
  ImageIOType::Pointer hdf5IO = ImageIOType::New();
  reader->SetImageIO(hdf5IO);
  reader->SetFileName(FILENAME);
  reader->Update();
  ImageType::Pointer image = reader->GetOutput();

  // Check largest possible and buffered regions.
  ImageType::RegionType expectedRegion;
  expectedRegion.SetIndex(0,0);
  expectedRegion.SetIndex(1,0);
  expectedRegion.SetIndex(2,0);
  expectedRegion.SetSize(0,5);
  expectedRegion.SetSize(1,5);
  expectedRegion.SetSize(2,5);
  if (image->GetLargestPossibleRegion() != expectedRegion)
    {
      std::cout << "Read image largest possible region: " << image->GetLargestPossibleRegion()
                << " doesn't match expectedo one: " << expectedRegion << std::endl;
    }
  if (image->GetBufferedRegion() != expectedRegion)
    {
      std::cout << "Read image buffered region: " << image->GetBufferedRegion()
                << " doesn't match expectedo one: " << expectedRegion << std::endl;
    }

  // Check image pixel values.
  itk::ImageRegionIterator<ImageType> it(image,image->GetLargestPossibleRegion());
  ImageType::IndexType index;
  typedef itk::VariableLengthVector<PixelType> VariableVectorType;
  VariableVectorType expectedValue;
  expectedValue.SetSize(SC);
  for(it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
    index = it.GetIndex();
    for (unsigned int ic=0; ic < SC; ++ic)
      {
      expectedValue[ic] = expectedPixelValue(index[2],index[1],index[0],ic);
      }
    if(image->GetPixel(index) != expectedValue)
      {
      std::cout << "At index :" << index
                << ", expected pixel value: " << expectedValue
                << ", but got: " << image->GetPixel(index)
                << std::endl;
      return EXIT_FAILURE;
      }
    }
  return EXIT_SUCCESS;
}

int
itkHDF5ImageIOPastingFromMemoryTest(int ac, char * av [])
{
  std::string prefix("");
  if(ac > 1)
    {
    prefix = *++av;
    --ac;
    itksys::SystemTools::ChangeDirectory(prefix.c_str());
    }
  itk::ObjectFactoryBase::RegisterFactory(itk::HDF5ImageIOFactory::New() );

  writeImage();
  pasteImage();
  return checkImage();
}
