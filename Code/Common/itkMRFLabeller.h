/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkMRFLabeller.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


  Copyright (c) 2000 National Library of Medicine
  All rights reserved.

  See COPYRIGHT.txt for copyright details.

=========================================================================*/
#ifndef _itkMRFLabeller_h
#define _itkMRFLabeller_h

//#include "itkObject.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

#include "itkFilterImageToImage.h"
#include "itkSupervisedClassifier.h"

namespace itk
{

/** \class MRFLabeller
 * \brief Implementation of a labeller object that uses Markov Random Fields
 * to classify pixels in a 3-D data set.
 * 
 * This object classifies pixels based on a 3-D Markov Random Field (MRF) 
 * model.This implementation uses the maximum a posteriori (MAP) estimates
 * for modeling the MRF. The object traverses the data set and uses the model 
 * generated by the Gaussian classifier to gets the the distance between 
 * each pixel in the data set to a set of known classes, updates the 
 * distances by evaluating the influence of its neighboring pixels (based 
 * on a 3-D MRF model) and finally, classifies each pixel to the class 
 * which has the minimum distance to that pixel (taking the neighborhood 
 * influence under consideration).
 *
 * The a classified initial labeled image is needed. It is important
 * that the number of expected classes be set before calling the 
 * classifier. In our case we have used the GaussianSupervisedClassifer to
 * generate the initial labels. This classifier requires the user to 
 * ensure that an appropriate training image set be provided. See the 
 * documentation of the classifier class for more information.
 *
 * The influence of a three-dimensional neighborhood on a given pixel's
 * classification (the MRF term) is computed by calculating a weighted
 * sum of number of class labels in a three dimensional neighborhood.
 * The basic idea of this neighborhood influence is that if a large
 * number of neighbors of a pixel are of one class, then the current
 * pixel is likely to be of the same class.
 *
 * The dimensions of the 3-D neighborhood and values of the weighting 
 * parameters are either specified by the user through the beta 
 * parameter or a default weighting table is generated during object 
 * construction. The following table shows an example of a 3x3x3 
 * neighborhood and the weighting values used. A 3 x 3 x 3 kernel
 * is used where each value is a nonnegative parameter, which encourages 
 * neighbors to be of the same class. In this example, the influence of
 * the pixels in the same slice is assigned a weight 1.7, the influence
 * of the pixels in the same location in the previous and next slice is 
 * assigned a weight 1.5, while a weight 1.3 is assigned to the influence of 
 * the north, south, east, west and diagonal pixels in the previous and next 
 * slices. 
 * \f[\begin{tabular}{ccc}
 *  \begin{tabular}{|c|c|c|}
 *   1.3 & 1.3 & 1.3 \\
 *   1.3 & 1.5 & 1.3 \\
 *   1.3 & 1.3 & 1.3 \\
 *  \end{tabular} &
 *  \begin{tabular}{|c|c|c|}
 *   1.7 & 1.7 & 1.7 \\
 *   1.7 & 0 & 1.7 \\
 *   1.7 & 1.7 & 1.7 \\
 *  \end{tabular} &
 *  \begin{tabular}{|c|c|c|}
 *   1.3 & 1.3 & 1.3 \\
 *   1.5 & 1.5 & 1.3 \\
 *   1.3 & 1.3 & 1.3 \\
 *  \end{tabular} \\
 * \end{tabular}\f]
 *
 * For minimization of the MRF labeling function the MinimizeFunctional
 * virtual method is called. For our current implementation we use the
 * the iterated conditional modes (ICM) algorithm described by Besag in the
 * paper ``On the Statistical Analysis of Dirty Pictures'' in J. Royal Stat.
 * Soc. B, Vol. 48, 1986. 
 *
 * In each iteration, the algorithm visits each pixel in turn and 
 * determines whether to update its classification by computing the influence
 * of the classification of the pixel's neighbors and of the intensity data.
 * On each iteration after the first, we reexamine the classification of a 
 * pixel only if the classification of some of its neighbors has changed
 * in the previous iteration. The pixels' classification is updated using a 
 * synchronous scheme (iteration by iteration) until the error reaches
 * less than the threshold or the number of iteration exceed the maximum set
 * number of iterations. 
 */

template <class TInputImage, class TClassifiedImage>
class ITK_EXPORT MRFLabeller : 
  public FilterImageToImage<TInputImage,TClassifiedImage>

{
public:       
  /**
   * Standard "Self" typedef.
   */
  typedef MRFLabeller   Self;

  /**
   * Standard "Superclass" typedef
   */
  typedef FilterImageToImage<TInputImage,TClassifiedImage> Superclass;

  /** 
   * Smart pointer typedef support.
   */
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** 
   * Run-time type information (and related methods).
   */
  itkTypeMacro(MRFLabeller,Object);

  /**
   * Method for creation through the object factory.
   */
  itkNewMacro(Self);

  /**
   * Type definition for the input image.
   */
  typedef typename TInputImage::Pointer              InputImageType;  

  /**
   * Type definition for the input image pixel type.
   */
  typedef typename TInputImage::PixelType            InputPixelType;

  /**
   * Type definitions for the training image.
   */
  typedef typename TClassifiedImage::Pointer         TrainingImageType;

  /**
   * Type definitions for the training image pixel type.
   */
  typedef typename TClassifiedImage::PixelType       TrainingPixelType;

  /**
   * Type definitions for the labelled image.
   * It is derived from the training image.
   */
  typedef typename TClassifiedImage::Pointer         LabelledImageType;
      
  /**
   * Type definitions for the classified image pixel type.
   * It has to be the same type as the training image.
   */
  typedef typename TClassifiedImage::PixelType       LabelledPixelType;

  /**
   * Type definitions for classifier to be used for the MRF lavbelling.
   */
  typedef Classifier<TInputImage,TClassifiedImage> ClassifierType;

  /**
   * Pointer to the classifier to be used for the MRF lavbelling.
   */
  typename ClassifierType::Pointer m_ClassifierPtr;

  typedef typename TInputImage::PixelType      InputImagePixelType;
  typedef typename TClassifiedImage::PixelType TrainingImagePixelType;
  typedef typename TClassifiedImage::PixelType LabelledImagePixelType;

  typedef
    ImageRegionIterator< InputImagePixelType, 
                         TInputImage::ImageDimension>  
      InputImageIterator;
  typedef
    ImageRegionIterator< TrainingImagePixelType, 
                         TClassifiedImage::ImageDimension> 
      LabelledImageIterator;

  typedef typename TInputImage::PixelType::VectorType 
    InputImageVectorType;

  /**
   * Set the image required for training type classifiers
   */
  void SetTrainingImage(TrainingImageType image);

  /** 
   * Set the labelled image. 
   */
  void SetLabelledImage(LabelledImageType LabelledImage);

  /** 
   * Get the labelled image. 
   */
  LabelledImageType GetLabelledImage()
  {
    return m_LabelledImage;
  }

  /**
   * Set the pointer to the classifer being used.
   */
  void SetClassifier( typename ClassifierType::Pointer ptrToClassifier );

  /**
   * Set the Number of class macro
   */
  itkSetMacro(NumClasses, unsigned int);

  /**
   * Get the Number of class macro
   */
  itkGetMacro(NumClasses, unsigned int);

  /**
   * Set the number of iteration of the Iterated Conditional Mode
   * (ICM) algorithm. A default value is set at 50 iterations.
   */
  itkSetMacro(MaxNumIter, unsigned int);

  /**
   * Set the number of iteration of the Iterated Conditional Mode
   * (ICM) algorithm.
   */
  itkGetMacro(MaxNumIter, unsigned int);

  /**
   * Set the error tollerance level which is used as a threshold
   * to quit the iterations
   */
  itkSetMacro(ErrorTollerance, double);

  /**
   * Get the error tollerance level which is used as a threshold
   * to quit the iterations
   */
  itkGetMacro(ErrorTollerance, double);


  /**
   * Set the weighting parameters (Beta Matrix). A default 3 x 3 x 3 
   * matrix is provided. However, the user is allowed to override it
   * with their choice of weights for a 3 x 3 x 3 matrix.
   */
  virtual void SetBeta( double* );

  /**
   * Get Beta matrix
   */
  double* GetBeta()
  {
    return m_Beta3x3x3;
  }

  /**
   * Set the weighting parameters (Beta Matrix). This is an overloaded
   * function allowing the users to set the Beta Matrix by providing a 
   * a 1D array of weights. Current implementation supports only a 
   * 3 x 3 x 3 kernel. The labeler needs to be extended for a different
   * kernel size.
   */
  virtual void SetBeta( double *BetaMatrix, unsigned int kernelSize );
  virtual void SetBeta( vnl_vector<double> BetaMatrix );
      
protected:
  /**
   * Constructor
   */
  MRFLabeller();

  /**
   * Destructor
   */
  ~MRFLabeller();

  /**
   * Copy constructor
   */
  MRFLabeller(const Self&) {}

  /**
   * Assignment operator
   */
  void operator=(const Self&) {}

  /**
   * Print self identity
   */      
  void PrintSelf(std::ostream& os, Indent indent);

  /**
   * Allocate memory for Labelled Images
   */
  void Allocate();

  /**
   * Minimization algorithm to be used.
   */
  virtual void MinimizeFunctional();

  virtual void GenerateData();
  virtual void GenerateInputRequestedRegion();
  virtual void EnlargeOutputRequestedRegion( DataObject * );
  virtual void GenerateOutputInformation();

private:            
  typedef typename TInputImage::SizeType InputImageSizeType;

  InputImageType         m_InputImage;
  TrainingImageType      m_TrainingImage;
  LabelledImageType      m_LabelledImage;
          
  unsigned int           m_NumClasses;
  unsigned int           m_MaxNumIter;
  unsigned int           m_KernelSize;
  unsigned int           *m_LabelStatus;
  
  double                 m_ErrorTollerance;
  double                 *m_ClassProbability; //Class liklihood
  double                 *m_Beta3x3x3;


  int                    m_ErrorCounter;
  int                    *m_Offset;
  int                    m_kWidth;
  int                    m_kHeight;
  int                    m_kDepth;
  int                    m_imgWidth;
  int                    m_imgHeight;
  int                    m_imgDepth;

  int                    *m_WidthOffset;
  int                    *m_HeightOffset;
  int                    *m_DepthOffset;

  /**
   * Apply MRF Classifier. Label the image set using Iterated Conditional 
   * Mode algorithm by J. Besag, "On statistical analysis of dirty pictures,"
   * J. Royal Stat. Soc. B, vol. 48, pp. 259-302, 1986.
   */
  void ApplyMRFLabeller();  

  //Function implementing the ICM algorithm to label the images
  void ApplyICMLabeller();

}; // class MRFLabeller


} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMRFLabeller.txx"
#endif



#endif

