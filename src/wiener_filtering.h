#ifndef WIENER_FILTERING_H
#define WIENER_FILTERING_H

#include <complex>
#include <vector>

#include "survey.h"
#include "field.h"
#include "surface_reconstruction.h"
#include <sparse2d/FFTN_2D.h>
#include <sparse2d/IM_Rot.h>
//#include "fitsio.h"

#include <numeric>
#include <algorithm>
#include <sparse2d/IM_IO.h>

class wiener_filtering
{
	long ngal;
	int npix, impix;
  double pixel_size;
	double fftFactor;
	int ps1dlen;
  // copy stuff from previous classes
  surface_reconstruction *rec;
  field *f;
  survey *surv;

	// no binning will be performed here
	/*
	std::vector<double> res1; // stores the shear residuals, 1st component
  std::vector<double> res2; // stores the shear residuals, 2nd component
  std::vector<double> pos1; // stores the positional x-coordinates
  std::vector<double> pos2; // stores the positional y-coordinates
  std::vector<double> wght; // stores the weights
	*/

	std::vector<double> kwiener; // stores the wiener solution
  std::vector<std::complex<double>> shear_res; // stores the binned shear residuals
  // std::vector<double> shear; // stores the shear values for computing the covariance
	std::vector<double> shear0;
	std::vector<double> shear1;
	std::vector<int> counter; // stores the number of galaxies per pixel
	std::vector<int> index1d;

  std::vector<double> Px; //stores the power spectrum of the signal
  std::vector<double> Sn; //stores the noise covariance, diagonal matrix in space

	std::vector<std::complex<double>> wiener_residuals; // stores the wiener filtering shear residuals

	// std::vector<double> ps_residual;

  // fftw planes
  fftw_complex *frame;
	fftw_plan planf;
	fftw_plan planb;

	fftw_complex *fft_frameHf0;
	fftw_complex *fft_frameHf1;
  fftw_complex *fft_frameHb0;
	fftw_complex *fft_frameHb1;
  fftw_plan planHf0;
	fftw_plan planHf1;
  fftw_plan planHb0;
	fftw_plan planHb1;

	/*! fft shift of a vectorized image
   *
   */
	void fftshift(complex<double> *image, complex<double> *shiftedimg, int npix);

	// std::vector<double> csp(std::vector<std::complex <double>>, double rpow, int it);

	/*! compute an isotropic power spectrum
   *
   */
	//void get_isotropic_spectrum(Ifloat & Imag, fltarray & Spectrum, float Resol, float Dens);

	/*
   * Run the core wiener algorithm
   */
	std::vector<double> run_wiener_core(unsigned int Niter);

	/*
	 * Function to transform from shear to convergence
	 */
	std::vector<std::complex<double>> H_adjoint_operator(std::vector<std::complex<double>> gamma);
	std::vector<std::complex<double>> H_adjoint_operator_new(std::vector<std::complex<double>> gamma);
	/*
	 * Function to transform from convergence to shear
	 */
	std::vector<std::complex<double>> H_operator(std::vector<std::complex<double>> kappa);
	std::vector<std::complex<double>> H_operator_new(std::vector<std::complex<double>> kappa);

	/*
   * Compute the noise covariance matrix
   */
  // void compute_noise_cov(std::vector<int> bin_ind);
	// void compute_noise_cov_old(std::vector<double> data);

	/*! Compute the noise covariance matrix using radnom rotations of the shear
   *
   */
  //void compute_noise_cov_random(std::vector<double> shearmap1, std::vector<double> shearmap2, int Niter);

	/*! Get current reconstructed residuals
   *
   */
	// void get_residuals();

  /*
   * 2-dimensional binning of residuals
   */
  // std::vector<int> bin2d();

  /*
   * tools for the binning process
   */
  // std::vector<double> linspace(double min, double max, int n);
	//std::vector<std::complex<double>> bincount(std::vector<int> ind1d, std::vector<double> residual1, std::vector<double> residual);
	// std::vector<int> digitize(std::vector<double> values, std::vector<double> edges);
	std::vector<double> powerspec2d(std::vector<double> spectrum);

public:

	std::vector<std::complex<double>> compute_wiener_residuals();


  /*! Constructor
   *
   */
  wiener_filtering ( surface_reconstruction *rec, field *f, survey *surv, dblarray shearamap, dblarray shearbmap, dblarray ps_array, dblarray noise_cov);

  /*! Destructor
   *
   */
  ~wiener_filtering();

  /*
   * main iteration of wiener filtering
   */
  void main_iteration();


	// std::vector<std::complex<double>> get_wiener_residuals(){
	// 	return wiener_residuals;
	// }

};

#endif // WIENER_FILTERING_H
