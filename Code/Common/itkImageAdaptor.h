/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImageAdaptor.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


  Copyright (c) 2000 National Library of Medicine
  All rights reserved.

  See COPYRIGHT.txt for copyright details.

=========================================================================*/
#ifndef __itkImageAdaptor_h
#define __itkImageAdaptor_h

#include "itkImage.h"

namespace itk
{

/**
 * \class ImageAdaptor
 * \brief Give access to partial aspects of voxels from an Image
 *
 * ImageAdaptors are templated over the ImageType and over a functor
 * that will specify what part of the pixel can be accessed
 *
 * The basic aspects of this class are the types it defines.  
 *
 * Image adaptors can be used as intermediate classes that allow
 * to send an image to a filter, specifying what part of the
 * image pixels will the filter act on.
 *
 * The TAccessor class should implement the Get and Set methods
 * as static methods. These two will specify how data can be put
 * and get from parts of each pixel. It should define the types
 * ExternalType and InternalType too.
 *
 */

template <class TImage, class TAccessor >
class ITK_EXPORT ImageAdaptor : public DataObject 
{
public:
  /**
   * Standard "Self" typedef.
   */
  typedef ImageAdaptor  Self;

  /** 
   * Smart pointer typedef support.
   */
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** 
   * Pixel typedef support. Used to declare pixel type in filters
   * or other operations.
   */
  typedef typename TAccessor::ExternalType PixelType;

  /** 
   * Pixel typedef support. Used to declare pixel type in filters
   * or other operations.
   */
  typedef typename TAccessor::InternalType InternalPixelType;


  /** 
   *  Accessor type that convert data between internal and external
   *  representations.
   */
  typedef   TAccessor   AccessorType;



  /**
   * Dimension of the image.  This enum is used by functions that are
   * templated over image type (as opposed to being templated over pixel
   * type and dimension) when they need compile time access to the dimension
   * of the image.
   */
  enum { ImageDimension = TImage::ImageDimension };
  
  /** 
   * Index typedef support. An index is used to access pixel values.
   */
  typedef typename TImage::IndexType  IndexType;

  /** 
   * Size typedef support. A size is used to define region bounds.
   */
  typedef typename TImage::SizeType  SizeType;

  /** 
   * Region typedef support. A region is used to specify a subset of an image.
   */
  typedef typename TImage::RegionType  RegionType;
  
  /** 
   * Run-time type information (and related methods).
   */
  itkTypeMacro(ImageAdaptor, ImageSource);

  /**
   * Method for creation through the object factory.
   */
  itkNewMacro(Self);  

  /** 
   * Image dimension. The dimension of an image is fixed at construction.
   */
  static unsigned int GetImageDimension() 
    { return  TImage::ImageDimension; }

  /**
   * Set the region object that defines the size and starting index
   * for the largest possible region this image could represent.  This
   * is used in determining how much memory would be needed to load an
   * entire dataset.  It is also used to determine boundary
   * conditions.
   * \sa ImageRegion, SetBufferedRegion(), SetRequestedRegion()
   */
  void SetLargestPossibleRegion(const RegionType &region);

  /**
   * Get the region object that defines the size and starting index
   * for the largest possible region this image could represent.  This
   * is used in determining how much memory would be needed to load an
   * entire dataset.  It is also used to determine boundary
   * conditions.
   * \sa ImageRegion, GetBufferedRegion(), GetRequestedRegion()
   */
  const RegionType& GetLargestPossibleRegion();

  /**
   * Set the region object that defines the size and starting index
   * of the region of the image currently load in memory. 
   * \sa ImageRegion, SetLargestPossibleRegion(), SetRequestedRegion()
   */
  void SetBufferedRegion(const RegionType &region);


  /**
   * Get the region object that defines the size and starting index
   * of the region of the image currently load in memory. 
   * \sa ImageRegion, SetLargestPossibleRegion(), SetRequestedRegion()
   */
  const RegionType& GetBufferedRegion();
  

  /**
   * Set the region object that defines the size and starting index
   * for the region of the image requested.
   * \sa ImageRegion, SetLargestPossibleRegion(), SetBufferedRegion()
   */
  void SetRequestedRegion(const RegionType &region);


  /**
   * Get the region object that defines the size and starting index
   * for the region of the image requested.
   * \sa ImageRegion, SetLargestPossibleRegion(), SetBufferedRegion()
   */
  const RegionType& GetRequestedRegion();


  /**
   * Allocate the image memory. Dimension and Size must be set a priori.
   */
  void Allocate();


  /**
   * Set a pixel.
   */
  void SetPixel(const IndexType &index, const PixelType & value)
  { TAccessor::Set( m_Image->GetPixel(index), value ); }
  
  /**
   * Get a pixel (read only version) 
   */
  const PixelType & GetPixel(const IndexType &index) const
  { return TAccessor::Get( m_Image->GetPixel(index) ); }

  /**
   * Get a pixel for editing. 
   */
  PixelType & GetPixel(const IndexType &index)
  { return TAccessor::Get( m_Image->GetPixel(index) ); }
    
  /**
   * Access a pixel. This version can be an lvalue.
   */
  PixelType & operator[](const IndexType &index)
  { return TAccessor::Get( m_Image->GetPixel(index) ); }
  
  /**
   * Access a pixel. This version can only be an rvalue.
   */
  const PixelType & operator[](const IndexType &index) const
  { return TAccessor::Get( m_Image->GetPixel(index) ); }


  /** 
   * Set the spacing (size of a pixel) of the image.
   * \sa GetSpacing()
   */
  virtual void SetSpacing( const double values[TImage::ImageDimension] );


  /** 
   * Get the spacing (size of a pixel) of the image.
   * \sa SetSpacing()
   */
  virtual const double * GetSpacing() const;


  /** 
   * Set the origin of the image.
   * \sa GetOrigin()
   */
  virtual void SetOrigin( const double values[TImage::ImageDimension] );


  /** 
   * Get the origin of the image.
   * \sa SetOrigin()
   */
  virtual const double * GetOrigin() const;


  /** 
   * PixelContainer typedef support. Used to construct a container for
   * the pixel data.
   */
  typedef typename TImage::PixelContainer        PixelContainer;
  typedef typename TImage::PixelContainerPointer PixelContainerPointer;


  /**
   * Return a pointer to the beginning of the buffer.  This is used by
   * the image iterator class.
   */
  InternalPixelType *GetBufferPointer();

  typedef InternalPixelType * InternalPixelPointerType;
  void GetBufferPointer2( InternalPixelPointerType  & );

  /**
   * Get the offset table.  
   */
  const unsigned long *GetOffsetTable() const;
  
  /**
   * Compute an offset from the beginning of the buffer for a pixel
   * at the specified index.
   */
  unsigned long ComputeOffset(const IndexType &ind) const;

  /**
   * Compute the index of the pixel at a specified offset from the
   * beginning of the buffer.
   */
  IndexType ComputeIndex(unsigned long offset) const;


  /**
   * Set Internal Image
   */
  virtual void SetImage( TImage * );

  /**
   * Delegate Modified to the Internal Image
   */
  virtual void Modified();


  /**
   * Delegate GetMTime to the Internal Image
   */
  virtual unsigned long GetMTime();


  virtual void Update();
  virtual void UpdateOutputInformation();
  virtual void SetRequestedRegionToLargestPossibleRegion();
  virtual void CopyInformation(DataObject *data);
  virtual bool RequestedRegionIsOutsideOfTheBufferedRegion();
  virtual bool VerifyRequestedRegion();
 
protected:
  ImageAdaptor();
  virtual ~ImageAdaptor();
  ImageAdaptor(const Self&);
  void operator=(const Self&);
  void PrintSelf(std::ostream& os, Indent indent);

  
private:
  
  typename TImage::Pointer   m_Image;

};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImageAdaptor.txx"
#endif
  
  

#endif

