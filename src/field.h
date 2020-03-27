/*! Copyright CEA, 2015-2016
 * author : Francois Lanusse < francois.lanusse@gmail.com >
 *
 * This software is a computer program whose purpose is to reconstruct mass maps
 * from weak gravitational lensing.
 *
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 *
 */

#ifndef FIELD_H
#define FIELD_H

#include <complex>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
// #include <nfft3.h>
#include <nicaea/cosmo.h>
#include <sparse2d/IM_IO.h>
#include "survey.h"


// Reference redshift used to compute the 2D convergence maps
#define Z_INF 1000


/*! Class representing a lensing field in 2D or 3D.
 *
 * Implements the transforms from and to convergence
 *
 */
class field
{

  // Field geometry
  int      npix;		      	/*!< Number of pixels */
  int      impix;             /*!< Number of square image*/
  double   pixel_size;		      	/*!< Size of the pixels, in radians */
  int      padding_size;	      	/*!< Number of pixels included in the zero padding */
  double   size;                        /*!< Effective size of the field, including padding, in radians */
  double   convert_coordinates_unit;  	/*!< Factor to convert the provided coordinates in radians */

  int      nlp;			      	/*!< Number of lens planes */
  double   zlens;			/*!< Redshift of the lens if reconstructing a single lens plane */
  std::valarray<double> zlp_low;	/*!< Lower bound of redshift bin */
  std::valarray<double> zlp_up;    	/*!< Upper bound of redshift bin */

  // FFTs
  int NpixFFT;
  double fftFactor;
  // nfft_plan** ps;
  fftwf_plan fps;
  fftwf_plan bps;
  fftwf_complex *fframe;
  fftwf_complex *bframe;

  fftwf_complex* fft_frame;

  // Survey data
  long     ngal;                /*!< Number of galaxies */
  double * shear_gamma1;        /*!< 1D Array for storing the shear for each galaxy in the survey.*/
  double * shear_gamma2;        /*!< 1D Array for storing the shear for each galaxy in the survey.*/
  double * flexion_f1;          /*!< 1D Array for storing the first flexion for each galaxy in the survey.*/
  double * flexion_f2;          /*!< 1D Array for storing the first lexion for each galaxy in the survey.*/
  double * w_f;
  double * w_e;
  bool     include_flexion;     /*!< Flag indicating whether flexion measurements are included in the reconstruction.*/

  // Auxiliary arrays
  // double * res_gamma1;          /*!< Array storing the gamma residuals for each galaxy.*/
  // double * res_gamma2;          /*!< Array storing the gamma residuals for each galaxy.*/
  double * res_f1;              /*!< Array storing the gamma residuals for each galaxy.*/
  double * res_f2;              /*!< Array storing the gamma residuals for each galaxy.*/
  double * res_conv;            /*!< Array storing the convergence for each galaxy.*/
  // double * cov;                 /*!< Array to store the covariance matrix resulting from the reduced shear.*/
  // double * lensKernel;          /*!< Array storing the conditionned lensing efficiency kernel for each galaxy.*/
  // double * lensKernelTrue;      /*!< Array storing the original lensing efficiency kernel for each galaxy.*/

  // 3D specific variables
  double r_cond;                /*!< Condition number used for the pre-conditioning matrix. */
  double * P;                   /*!< Preconditionning matrix */
  double * PP;                  /*!< Square of the preconditionning matrix */
  double * iP;                  /*!< Inverse of the preconditionning matrix */

  survey *surv;                 /*!< Reference to the survey object */

  gsl_rng *rng;                 /*!< Random number generator */

  nicaea::error **err;          /*!< NICAEA error structure.*/
  nicaea::cosmo *model;         /*!< NICAEA cosmology used for the mapping.*/

  // Auxiliary variables used in binned version

  // std::vector<double> pos1; // stores the positional x-coordinates
  // std::vector<double> pos2; // stores the positional y-coordinates
  // std::vector<int> counter; // stores the number of galaxies per pixel
  // std::vector<int> index1d;
  double * res_gamma1_map;          /*!< Array storing the gamma residuals for each pixel.*/
  double * res_gamma2_map;          /*!< Array storing the gamma residuals for each pixel.*/
  double * res_conv_map;            /*!< Array storing the convergence for each galaxy.*/
  double * covariance_map;          /*!< Array to store the covariance matrix resulting from the reduced shear.*/
  double * Sn;
  double * mask;
  /*
   * tools for the binning process
   */
  // std::vector<double> linspace(double min, double max, int n);
	// std::vector<std::complex<double>> bincount(std::vector<int> ind1d, double * value1, double * value2);
  // std::vector<double> bincount(std::vector<int> ind1d, double * value1);
	// std::vector<int> digitize(std::vector<double> values, std::vector<double> edges);

  /*
   * binned data forward backward operators
   */
  void forward_operator_map(fftwf_complex *delta);
  void adjoint_operator_map(fftwf_complex *delta);



  /*! Computes lensing kernel for each galaxy, to reconstruct a surface mass density
   *
   */
  void compute_surface_lensing_kernel();

  /*! Computes lensing kernel for each galaxy, to reconstruct a 3D density contrast.
   *
   */
  void compute_3D_lensing_kernel();

  /*! Computes the forward lensing transform from density to shear.
   *
   */
  // void forward_operator(fftwf_complex *delta);

  /*! Compute the adjoint operation.
   *
   */
  // void adjoint_operator(fftwf_complex *delta, bool preconditionning=true);

public:
  double * nshear_power0;
  double * nshear_power1;
  double * get_nshear_power0() {
        return nshear_power0;
  }
  double * get_nshear_power1() {
        return nshear_power1;
  }
  std::vector<std::complex<double>> binned_shear; // stores the binned shear as a complex field

  /*! Constructor from configuration file and survey
   *
   */
  field(boost::property_tree::ptree config, survey *surv, dblarray shearamap, dblarray shearbmap, dblarray noise_cov); // , dblarray noise_cov

  /*! Destructor */
  ~field();

  /*! Transforms gamma to kappa
  *
  */
  //void compute_kappa(std::vector<double> s1, std::vector<double> s2, bool preconditionning, double * wrec);

  /*! Return the number of pixels for a field of size npix x npix.
  *
  */
  /*
  double * get_lensKernelTrue() {
        return lensKernelTrue;
  }
  */

  /*! Return the size of the pixels, in radians
   *
   */
   double get_pixel_size() {
         return pixel_size;
   }

   // std::vector<int> get_counter(){
   //       return counter;
   // }
   // std::vector<int> get_binindex(){
   //       return index1d;
   // }

   // std::vector<std::complex<double>> get_binned_shear() {
   //       return binned_shear;
   // }

  /*! Return the number of pixels for a field of size npix x npix.
   *
   */

  int get_npix() {
        return npix;
  }


  /*! Return the number of lens planes.
   *
   */
  int get_nlp() {
        return nlp;
  }

  /*! Returns the (ra, dec) of the pixel centers in degrees.
   * \a ra and \a dec are two preallocated arrays of size NxN
   * where N is the number of pixels
   */
   // void get_pixel_coordinates(double * ra, double *dec);

  /*! Computes the gradient of the chi_2 for a given delta field.
   *
   */
  void gradient_map(fftwf_complex *delta);

  /*! Computes the gradient of the chi_2 for a given delta field and randomized measurements.
   *
   */
  void gradient_noise_map(fftwf_complex *delta);

  /*! Updates the non-linear correction factor in the covariance matrix.
   *
   */
  void update_covariance_map(fftwf_complex *delta);

  /*! Computes the spectral norm of the lensing operator
   *
   */
  double get_spectral_norm_map(int niter, double tol);

  /*! Checks that the adjoint_operator is indeed the adjoint of the forward operator.
   *
   */
  bool check_adjoint();

  /*! Combines shear and flexion components applying a minimum variance filter
   *
   */
  void combine_components(fftwf_complex * delta, fftwf_complex * delta_comb);

  /*! Copies the same combined field into shear and flexion commponents
   *
   */
  void combine_components_inverse(fftwf_complex * delta_comb, fftwf_complex * delta);

  /*! Returns the preconditioning matrix
   *
   */
  const double * get_preconditioning_matrix(){ return P; }

   /*! Computes the signal residuals
   *
   */
  // void comp_residuals(fftwf_complex *delta, std::vector<std::complex<double>> &rshear);


  /*!
   * shear binning function, implementation based on vectors
   */
  std::vector<int> compute_index_map();
  //void compute_shear_map();
  //void print_pos();
};

#endif // FIELD_H
