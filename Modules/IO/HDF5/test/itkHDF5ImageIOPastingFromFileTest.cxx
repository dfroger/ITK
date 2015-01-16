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
#include "itkRegionOfInterestImageFilter.h"
#include "itkChangeInformationImageFilter.h"

/* Write a VectorImage of size (SZ SY,SX,SC) in a HDF5 first file.
 * Write a VectorImage of size (TZ TY,TX,SC) in a HDF5 second file.
 * Read in the second file a region at (JZ,JY,JX,0) of size (NZ,NY,NX,SC),
 * and paste it in the first file.
 *
 *     Y ^
 *       |
 *       |        TY +----------------------------+
 *       |           |                            |
 *       | SY +------+----------------+           |
 *       |    |      |                |           |
 *       | NY +   NX +   +--------+   |           |
 *       |    |      |   |        |   |           |
 *       |    |      |   |        |   |           |
 *       |    |      |   |Pasted  |   |           |
 *       | JY +   IY +   +--------+   |           |
 *       |    |      |                |           |
 *       |    |      |First           |           |
 *       |    |   OY +---+--------+---+-----------+
 *       |    |      OX  IX       NX  |           TX
 *       |    |                       |
 *       |    |Second                 |
 *       | PY +----------+--------+---+
 *       |    PX         JX       NX  SX
 *       |
 *       +----------------------------------------------------> X
 *
 * PX + JX = OX + IX
 * PY + JY = OY + IY
 * PZ + JZ = OZ + IZ
 *
 * */

typedef float PixelType;
const int Dimension = 3;
typedef itk::VectorImage<PixelType, Dimension> ImageType;

const unsigned int OZ = 3;
const unsigned int OY = 3;
const unsigned int OX = 3;

const unsigned int PZ = 2;
const unsigned int PY = 2;
const unsigned int PX = 2;

const unsigned int SZ = 6;
const unsigned int SY = 6;
const unsigned int SX = 6;
const unsigned int SC = 2;

const unsigned int TZ = 5;
const unsigned int TY = 5;
const unsigned int TX = 5;
const unsigned int TC = SC;

const unsigned int JZ = 2;
const unsigned int JY = 2;
const unsigned int JX = 2;

const unsigned int IZ = PZ + JZ - OZ;
const unsigned int IY = PY + JY - OY;
const unsigned int IX = PX + JX - OX;

const unsigned int NZ = 2;
const unsigned int NY = 2;
const unsigned int NX = 2;

const std::string FIRST_FILENAME = "first_image.hdf5";
const std::string SECOND_FILENAME = "second_image.hdf5";

static PixelType
firstImagePixelValue(unsigned int iz, unsigned int iy, unsigned int ix, unsigned int ic)
{
  return (OZ+iz)*1000 + (OY+iy)*100 + (OX+ix)*10 + ic;
}

static PixelType
secondImagePixelValue(unsigned int jz, unsigned int jy, unsigned int jx, unsigned int jc)
{
  return 990000 + (PZ+jz)*1000 + (PY+jy)*100 + (PX+jx)*10 + jc;
}

static bool
isPasted(unsigned int iz, unsigned int iy, unsigned int ix)
{
    return (IZ <= iz && iz < IZ+NZ &&
            IY <= iy && iy < IY+NY &&
            IX <= ix && ix < IX+NX);
}

static PixelType
expectedPixelValue(unsigned int iz, unsigned int iy, unsigned int ix, unsigned int ic)
{
    if ( isPasted(iz,iy,ix) )
        return secondImagePixelValue(iz, iy, ix, ic);
    else
        return firstImagePixelValue(iz, iy, ix, ic);
}

static void
writeFirstImage()
{
  // Allocate image.
  ImageType::IndexType index;
  index.Fill(0);
  ImageType::SizeType size;
  size[0] = SX;
  size[1] = SY;
  size[2] = SZ;
  ImageType::PointType origin;
  origin[0] = OX;
  origin[1] = OY;
  origin[2] = OZ;
  ImageType::RegionType region(index,size);
  ImageType::Pointer image = ImageType::New();
  image->SetRegions(region);
  image->SetOrigin( origin );
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
        value[ic] = firstImagePixelValue(index[2], index[1], index[0], ic);
        }
      image->SetPixel(index, value);
    }

  // Create writer.
  typedef itk::ImageFileWriter<ImageType>  WriterType;
  typedef itk::HDF5ImageIO                 ImageIOType;
  WriterType::Pointer writer = WriterType::New();
  ImageIOType::Pointer hdf5IO = ImageIOType::New();
  writer->SetImageIO(hdf5IO);
  writer->SetFileName(FIRST_FILENAME);

  // Update all the pipeflow.
  writer->SetInput( image );
  writer->Update();
}

static void
writeSecondImage()
{
  // Allocate image.
  ImageType::IndexType index;
  index.Fill(0);
  ImageType::SizeType size;
  size[0] = TX;
  size[1] = TY;
  size[2] = TZ;
  ImageType::PointType origin;
  origin[0] = PX;
  origin[1] = PY;
  origin[2] = PZ;
  ImageType::RegionType region(index,size);
  ImageType::Pointer image = ImageType::New();
  image->SetRegions(region);
  image->SetOrigin( origin );
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
      for (unsigned int jc=0; jc < TC; ++jc)
        {
        value[jc] = secondImagePixelValue(index[2], index[1], index[0], jc);
        }
      image->SetPixel(index, value);
    }

  // Create writer.
  typedef itk::ImageFileWriter<ImageType>  WriterType;
  typedef itk::HDF5ImageIO                 ImageIOType;
  WriterType::Pointer writer = WriterType::New();
  ImageIOType::Pointer hdf5IO = ImageIOType::New();
  writer->SetImageIO(hdf5IO);
  writer->SetFileName(SECOND_FILENAME);

  // Update all the pipeflow.
  writer->SetInput( image );
  writer->Update();
}

static void
pasteImage()
{
  /*
  // Create reader.
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(SECOND_FILENAME);

  // Adapt origin.
  typedef itk::ChangeInformationImageFilter< ImageType > ChangeInfoType;
  ChangeInfoType::Pointer changeinfo = ChangeInfoType::New();
  changeinfo->SetInput( reader->GetOutput() );
  ImageType::PointType origin;
  origin[0] = OX;
  origin[1] = OY;
  origin[2] = OZ;
  changeinfo->SetOutputOrigin( origin );
  changeinfo->ChangeOriginOn();

  //changeinfo->Update();
  //ImageType::Pointer im1 = changeinfo->GetOutput();
  //std::cout << "Origin: "  << im1->GetOrigin() << endl;
  //std::cout << "Buffered: "  << im1->GetBufferedRegion() << endl;
  //std::cout << "LargestPossible: "  << im1->GetLargestPossibleRegion() << endl;

  // Set requested region.
  ImageType::IndexType start;
  ImageType::SizeType size;
  start[0] = JX;
  start[1] = JY;
  start[2] = JZ;
  size[0] = NX;
  size[1] = NY;
  size[2] = NZ;
  typedef itk::ImageRegion<Dimension> RegionType;
  RegionType requestedRegion;
  requestedRegion.SetIndex(start);
  requestedRegion.SetSize(size);
  changeinfo->GetOutput()->SetRequestedRegion(requestedRegion);

  changeinfo->Update();
  ImageType::Pointer im0 = changeinfo->GetOutput();
  std::cout << "Origin: "  << im0->GetOrigin() << endl;
  std::cout << "Buffered: "  << im0->GetBufferedRegion() << endl;
  std::cout << "LargestPossible: "  << im0->GetLargestPossibleRegion() << endl;

  // Create desired region.
  start[0] = JX;
  start[1] = JY;
  start[2] = JZ;
  size[0] = NX;
  size[1] = NY;
  size[2] = NZ;
  ImageType::RegionType desiredRegion;
  desiredRegion.SetSize(size);
  desiredRegion.SetIndex(start);

  // Create region of interest filter.
  typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > FilterType;
  FilterType::Pointer interest = FilterType::New();
  interest->SetRegionOfInterest(desiredRegion);
  interest->SetInput(reader->GetOutput());

  interest->Update();
  ImageType::Pointer im2 = interest->GetOutput();
  std::cout << "Origin: "  << im2->GetOrigin() << std::endl;
  std::cout << "Buffered: "  << im2->GetBufferedRegion() << std::endl;
  std::cout << "LargestPossible: "  << im2->GetLargestPossibleRegion() << std::endl;

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
  writer->SetFileName(FIRST_FILENAME);
  writer->SetIORegion(ioregion);

  // Update all the pipeflow
  writer->SetInput( reader->GetOutput() );
  writer->Update();
  */
}

static int
checkImage()
{
  // Read image.
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::HDF5ImageIO                ImageIOType;
  ReaderType::Pointer reader = ReaderType::New();
  ImageIOType::Pointer hdf5IO = ImageIOType::New();
  reader->SetImageIO(hdf5IO);
  reader->SetFileName(FIRST_FILENAME);
  reader->Update();
  ImageType::Pointer image = reader->GetOutput();

  // Check largest possible and buffered regions.
  ImageType::RegionType expectedRegion;
  expectedRegion.SetIndex(0,0);
  expectedRegion.SetIndex(1,0);
  expectedRegion.SetIndex(2,0);
  expectedRegion.SetSize(0,SX);
  expectedRegion.SetSize(1,SY);
  expectedRegion.SetSize(2,SZ);
  if (image->GetLargestPossibleRegion() != expectedRegion)
    {
      std::cout << "Read image largest possible region: " << image->GetLargestPossibleRegion()
                << " doesn't match expected one: " << expectedRegion << std::endl;
    }
  if (image->GetBufferedRegion() != expectedRegion)
    {
      std::cout << "Read image buffered region: " << image->GetBufferedRegion()
                << " doesn't match expected one: " << expectedRegion << std::endl;
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
itkHDF5ImageIOPastingFromFileTest(int ac, char * av [])
{
  std::string prefix("");
  if(ac > 1)
    {
    prefix = *++av;
    --ac;
    itksys::SystemTools::ChangeDirectory(prefix.c_str());
    }
  itk::ObjectFactoryBase::RegisterFactory(itk::HDF5ImageIOFactory::New() );

  writeFirstImage();
  writeSecondImage();
  pasteImage();
  return checkImage();
}
