#define DEBUG_FITS 1

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

#include "field.h"

// #include <random>
#include <iostream>
#include <fstream>
#include <string>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>

#ifdef DEBUG_FITS
#include <sparse2d/IM_IO.h>
#endif

#define ZMAX 10.0

#undef pi
#undef sd

#include <armadillo>

field::field(boost::property_tree::ptree config, survey *su, dblarray shearamap, dblarray shearbmap, dblarray noise_cov)
{
    surv = su;

    // Reading configuration
    double Omega_m = config.get<double>("cosmology.Omega_m", 0.25);
    double h = config.get<double>("cosmology.h", 0.70);

    // Load  the  redshift range of the reconstruction
    nlp   = config.get<double>("field.nlp",    1);
    if (nlp == 1) {
        zlens = config.get<double>("field.zlens", -1);
        if(zlens <= 0){
            std::cout << "Warning: No redshift specified for the lens, ignoring source redsfhits" << std::endl;
        }
    } else {
        zlp_low.resize(nlp);
        zlp_up.resize(nlp);
        double zlp_min = config.get<double>("field.zlp_min");
        double zlp_max = config.get<double>("field.zlp_max");
        // Creates the lensplanes, for now just regularly spaced between zmin and zmax
        for (int i = 0; i < nlp; i++) {
            zlp_low[i] = zlp_min + (zlp_max - zlp_min) * ((double) i) / ((double) nlp);
            zlp_up[i]  = zlp_min + (zlp_max - zlp_min) * ((double)(i + 1.0)) / ((double) nlp);
        }
    }

    r_cond = config.get<double>("field.r_cond", 0.1);

    npix = shearamap.nx();
    impix = npix * npix;

    std::cout << "Number of pixels : " << npix << std::endl; // " Pixel size : " << pixel_size << std::endl;

    // Check whether flexion measurements are available
    include_flexion = config.get<bool>("field.include_flexion", false);
    if (include_flexion) {
        if (! surv->get_flexion_availability()) {
            std::cout << "No flexion measurements provided, reconstructing from shear alone" << std::endl;
            include_flexion = false;
        }
    }

    // Initialize the random number generator
    const gsl_rng_type *T;
    T = gsl_rng_default;
    rng = gsl_rng_alloc(T);

    // Initialize NICAEA
    err = new nicaea::error*;
    *err = NULL;
    model = nicaea::init_parameters(Omega_m, 1.0 - Omega_m, -1.0, 0.0, NULL, 0, h, 0.044, 0.0, 0.0, 0.80, 0.96,
                                    nicaea::smith03, nicaea::eisenhu, nicaea::growth_de, nicaea::linder,
                                    nicaea::norm_s8, 0.0, err);

    res_gamma1_map = (double *) malloc(sizeof(double) * impix);
    res_gamma2_map = (double *) malloc(sizeof(double) * impix);
    res_conv_map   = (double *) malloc(sizeof(double) * impix);
    covariance_map = (double *) malloc(sizeof(double) * impix);
    Sn             = (double *) malloc(sizeof(double) * impix);
    mask           = (double *) malloc(sizeof(double) * impix);

    // initialize arrays
    int cnt = 0;
    for(int i = 0; i < impix ; i++){
      res_gamma1_map[i] = 0.;
      res_gamma2_map[i] = 0.;
      res_conv_map[i]   = 0.;
      covariance_map[i] = 1.;

      if (noise_cov(i) > 1e2){
         Sn[i] = 0.;   // pixel is masked!
         mask[i] = 0.; // pixel is masked!
         cnt++;
      }
      else{
        mask[i] = 1.;
        Sn[i] = noise_cov(i);
      }
    }
    std::cout << cnt << " zero galaxy pixels out of " << npix*npix << std::endl;

    // // check noise_cov values
    // for(int i = 0; i < impix ; i++)
    //   if(Sn[i] > 1)
    //     std::cout << "/* message */" << i << '\n';

    // get shear residuals from wiener filtering
    binned_shear.reserve(impix);
    for (int ind = 0; ind < impix; ind++){
      binned_shear.push_back({shearamap(ind), -shearbmap(ind)}); // useless here, it is reset after wiener filtering, check surface reconstruction
    }

    fframe = fftwf_alloc_complex(npix * npix);
    bframe = fftwf_alloc_complex(npix * npix);

    fps = fftwf_plan_dft_2d(npix, npix, fframe, fframe, FFTW_FORWARD, FFTW_MEASURE);
    bps = fftwf_plan_dft_2d(npix, npix, bframe, bframe, FFTW_BACKWARD, FFTW_MEASURE);

    if (bps == NULL){
      std::cout << "NULL frame allocation" << std::endl;
    }

    fft_frame = fftwf_alloc_complex(npix * npix * nlp);

    // Normalization factor for the fft
    fftFactor = 1. / ((double)(npix*npix));

    // // Initialize the lensing kernel for each galaxy
    // // TODO::this has to change length of lenskernel as the number of pixels
    // lensKernel     = (double *) malloc(sizeof(double) * impix);
    // lensKernelTrue = (double *) malloc(sizeof(double) * impix);
    // for (int i = 0; i < impix; i++) {
    //   lensKernelTrue[i] = 1.;
    //   lensKernel[i] = 1.;
    // }



#ifdef DEBUG_FITS
    // store the variance of the noisy shear
    nshear_power0 = (double *) malloc(sizeof(double) * npix*npix);
    nshear_power1 = (double *) malloc(sizeof(double) * npix*npix);
    for (int i = 0; i < npix * npix; i++) {
      nshear_power0[i] = .0;
      nshear_power1[i] = .0;
    }

    // write binned_shear to a file
    double * img11 = (double *) malloc(sizeof(double) * npix*npix);
    double * img21 = (double *) malloc(sizeof(double) * npix*npix);
    for (int i = 0; i < npix * npix; i++) {
      img11[i] = binned_shear[i].real(); // this is input to glimpse
      img21[i] = binned_shear[i].imag();
    }
    dblarray expts11;
    dblarray expts21;
    expts11.alloc(img11, npix, npix);
    expts21.alloc(img21, npix, npix);
    fits_write_dblarr("g1_map.fits", expts11);
    fits_write_dblarr("g2_map.fits", expts21); // this succeeds

    // perform a gamma to kappa transform
    fftwf_complex *delta12 = fftwf_alloc_complex(npix * npix);

    for(int i = 0; i < impix ; i++){
      res_gamma1_map[i] = real(binned_shear[i]);
      res_gamma2_map[i] = imag(binned_shear[i]);
    }

    adjoint_operator_map(delta12);

    for (long ind = 0; ind < npix * npix; ind++) {
        bframe[ind][0] = delta12[ind][0] * fftFactor;
        bframe[ind][1] = delta12[ind][1] * fftFactor;
    }

    fftwf_execute(bps);

    double * img0 = (double *) malloc(sizeof(double) * npix*npix);
    double * img1 = (double *) malloc(sizeof(double) * npix*npix);
    for (int i = 0; i < npix * npix; i++){
      img0[i] = (double)bframe[i][0];
      img1[i] = (double)bframe[i][1];
    }

    dblarray expts0;
    expts0.alloc(img0, npix, npix);
    fits_write_dblarr("gamma_to_kappa0.fits", expts0);
    dblarray expts1;
    expts1.alloc(img1, npix, npix);
    fits_write_dblarr("gamma_to_kappa1.fits", expts1);

    // perform a kappa to gamma transform

    forward_operator_map(delta12);

    // write res_gamma to a file
    double * img3 = (double *) malloc(sizeof(double) * npix*npix);
    double * img4 = (double *) malloc(sizeof(double) * npix*npix);
    for (int i = 0; i < npix * npix; i++) {
      img3[i] = res_gamma1_map[i];
      img4[i] = res_gamma2_map[i];
    }
    dblarray expts3;
    dblarray expts4;
    expts3.alloc(img3, npix, npix);
    expts4.alloc(img4, npix, npix);
    fits_write_dblarr("kg1_map.fits", expts3);
    fits_write_dblarr("kg2_map.fits", expts4); // this succeeds

    for(int i = 0; i < impix ; i++){
      res_gamma1_map[i] = .0;
      res_gamma2_map[i] = .0;
    }

#endif

}

field::~field()
{
    gsl_rng_free(rng);
    //TODO: free nicaea

    // Free data arrays
    free(shear_gamma1);
    free(shear_gamma2);
    free(w_e);
    free(res_gamma1_map);
    free(res_gamma2_map);
    free(res_conv_map);
    free(covariance_map);
    // free(lensKernel);
    // free(lensKernelTrue);

    if (include_flexion) {
        free(res_f1);
        free(res_f2);
        free(flexion_f1);
        free(flexion_f2);
        free(w_f);
    }

    // Deallocate nfft plans
    for (int i = 0; i < nlp; i++) {
        // nfft_finalize(ps[i]);
        fftwf_destroy_plan(fps);
        fftwf_destroy_plan(bps);
    }
    // free(fps);
    // free(bps);
    fftwf_free(fft_frame);
    fftwf_free(fframe);
    fftwf_free(bframe);
}


void field::combine_components(fftwf_complex *delta, fftwf_complex *delta_comb)
{

    for (int z = 0; z < nlp; z++) {
        // Computes the convergence at the position of the
        for (int y = 0; y < npix ; y++) {
            for (int x = 0; x < npix ; x++) {
                long pos = y * npix + x + z * (npix * npix);
                delta_comb[pos][0] = delta[pos][0];
                delta_comb[pos][1] = delta[pos][1];
            }
        }
        delta_comb[0][0] = .0;
        delta_comb[0][1] = .0;
    }
}


void field::combine_components_inverse(fftwf_complex *delta_comb, fftwf_complex *delta)
{

    for (int z = 0; z < nlp; z++) {
        // Computes the convergence at the position of the galaxy
        for (int y = 0; y < npix ; y++) {
            for (int x = 0; x < npix ; x++) {
                long pos = y * npix + x + z * (npix * npix);
                delta[pos][0] = delta_comb[pos][0];
                delta[pos][1] = delta_comb[pos][1];
            }
        }
    }
}


bool field::check_adjoint()
{
    fftwf_complex *delta1 = fftwf_alloc_complex( npix * npix * nlp);
    fftwf_complex *delta2 = fftwf_alloc_complex( npix * npix * nlp);

    double *test_g1 = (double *) malloc(sizeof(double) * npix * npix);
    double *test_g2 = (double *) malloc(sizeof(double) * npix * npix);

    for (long ind = 0; ind < nlp * npix * npix ; ind++) {
        delta1[ind][0] = gsl_ran_gaussian(rng, 1.0);
        delta1[ind][1] = gsl_ran_gaussian(rng, 1.0);
    }

    double result_forward   = 0;
    double result_forward2  = 0;
    double result_backward  = 0;
    double result_backward2 = 0;

    forward_operator_map(delta1);

    for (long ind = 0; ind < npix * npix ; ind++) {
        test_g1[ind] = res_gamma1_map[ind];
        test_g2[ind] = res_gamma2_map[ind];
        res_gamma1_map[ind] = gsl_ran_gaussian(rng, 1.0);
        res_gamma2_map[ind] = gsl_ran_gaussian(rng, 1.0);
    }

    adjoint_operator_map(delta2);

    for (long ind = 0; ind < npix * npix * nlp ; ind++) {
        result_forward  +=  delta1[ind][0] * delta2[ind][0] * fftFactor + delta1[ind][1] * delta2[ind][1] * fftFactor; // ()* fftFactor
        result_forward2 += -delta1[ind][1] * delta2[ind][0] * fftFactor + delta1[ind][0] * delta2[ind][1] * fftFactor; // ()* fftFactor
    }

    for (long ind = 0; ind < npix * npix; ind++) {
        result_backward  += res_gamma1_map[ind] * test_g1[ind] + res_gamma2_map[ind] * test_g2[ind];
        result_backward2 += res_gamma2_map[ind] * test_g1[ind] - res_gamma1_map[ind] * test_g2[ind];
    }
    std::cout << " Results  of check: " << result_forward << " against " << result_backward  << std::endl;
    std::cout << " Results  of check: " << result_forward2 << " against " << result_backward2  << std::endl;
    return true;
}

typedef struct {
    double w_a;
    nicaea::error **err;
    nicaea::cosmo *model;
} int_for_3d_efficiency_params;

#define EPS_GW_INT 1.0E-14
double int_for_3d_efficiency ( double aprime, void* intpar) {

    double wprime, fKwp, fKw, fKwwp, dwda;

    int_for_3d_efficiency_params *params = (int_for_3d_efficiency_params *) intpar ;
    nicaea::cosmo *self                  = params->model;
    nicaea::error **err                  = params->err;
    double w                             = params->w_a;

    double fac = 1.5/nicaea::dsqr ( R_HUBBLE ) * ( self->Omega_m+ self->Omega_nu_mass ) /aprime;

    wprime = nicaea::w ( self, aprime, 0, err );
    quitOnError ( *err, __LINE__, stderr );

    if( wprime >= w ) return 0;

    fKwp   = nicaea::f_K ( self, wprime, err );
    quitOnError ( *err, __LINE__, stderr );

    fKw   = nicaea::f_K ( self, w, err );
    quitOnError ( *err, __LINE__, stderr );

    fKwwp  = nicaea::f_K ( self, w - wprime, err );
    quitOnError ( *err, __LINE__, stderr );

    dwda   = nicaea::dwoverda ( self,aprime,err );
    quitOnError ( *err, __LINE__, stderr );

    // Prevent very small values of Gw from perturbing the integration
    if ( fKwwp < EPS_GW_INT ) fKwwp = 0.0;

    return fac * fKwp * fKwwp / fKw * dwda;
}
#undef EPS_GW_INT


typedef struct {
    redshift_distribution * redshift;
    double *x;
    double *y;
    gsl_interp *interpolator;
    gsl_interp_accel *accelerator;
} int_for_marginalisation_params;

double int_for_marginalisation(double  z, void *intpar){
    int_for_marginalisation_params *p = (int_for_marginalisation_params *) intpar;

    return p->redshift->pdf(z) * gsl_interp_eval(p->interpolator, p->x, p->y, z, p->accelerator);
}


typedef struct {
    double w_l;
    double w_inf;
    redshift_distribution *redshift;
    nicaea::error **err;
    nicaea::cosmo *model;
} int_for_sigma_params;

double int_for_sigma(double a_s, void *intpar)
{
    int_for_sigma_params *params    = (int_for_sigma_params *) intpar ;
    nicaea::cosmo *self             = params->model;
    nicaea::error **err             = params->err;
    redshift_distribution *redshift = params->redshift;
    double w_l                      = params->w_l;
    double w_inf                    = params->w_inf;

    double w_s = nicaea::w(self, a_s, 0, err);
    if (w_s - w_l <= 0) {
        return 0;
    }

    double p = 1.0 / a_s / a_s * redshift->pdf(1.0 / a_s - 1.0);

    return  p * ((w_s - w_l) * w_inf) / ((w_inf - w_l) * w_s);
}


/*
 * Function to transform from convergence to shear
 */
void field::forward_operator_map(fftwf_complex *delta)
{
	// transform from kappa to gamma in fourier space
	float freqFactor = 2.0 * M_PI / ((float) npix); // 2.0 * M_PI / ((double) npix ) ; // 1.0 / ((double) npix); // / pixel_size

  #pragma omp parallel for
  for (int z = 0; z < nlp; z++) {
    float k1, k2, k1k1, k2k2, k1k2, ksqr;
  	float c1, c2;

    for (int y = 0; y < npix ; y++) {
  		k2  = (y < npix / 2 ? y * freqFactor : (y - npix) * freqFactor); // this comes from python fft
      // k2 = (y - npix) * freqFactor;

      for (int x = 0; x < npix ; x++) {
        k1 = (x < npix / 2 ? x * freqFactor : (x - npix) * freqFactor); // this comes from python fft
        // k1 = (x - npix) * freqFactor;

        long pos = y * npix + x ;

  			k1k1 = k1 * k1;
  			k2k2 = k2 * k2;
  			ksqr = k1k1 + k2k2;
  			c1 = k1k1 - k2k2;
  			c2 = 2.0 * k1 * k2;

        if (k1 == 0 && k2 == 0) {
            bframe[pos][0] = (( delta[pos][0] * c1 + delta[pos][1] * c2) * fftFactor);
            bframe[pos][0] = ((-delta[pos][0] * c2 + delta[pos][1] * c1) * fftFactor);
        }
        else {
        bframe[pos][0] = (( delta[pos][0] * c1 + delta[pos][1] * c2)  / ksqr * fftFactor);
			  bframe[pos][1] = ((-delta[pos][0] * c2 + delta[pos][1] * c1)  / ksqr * fftFactor);
      }
  		}
  	}

  	fftwf_execute(bps);

  }

  // Apply the lensing efficiency kernel
  #pragma omp parallel for
  for (int i = 0; i < impix ; i++) {
      res_gamma1_map[i] = 0;
      res_gamma2_map[i] = 0;
      for (int z = 0; z < nlp ; z++) {
          double q = 1.; // lensKernel[i * nlp + z];
          res_gamma1_map[i] += q * bframe[i][0] ;
          res_gamma2_map[i] += q * bframe[i][1] ;
      }
  }

  // Compute the value of the field evaluated at each galaxy position
  combine_components(delta, fft_frame);
  #pragma omp parallel for
  for (int z = 0; z < nlp; z++) {
      for (int y = 0; y < npix ; y++) {
          // int ky  = (y < npix / 2 ? y + npix / 2 : y - npix / 2);
          for (int x = 0; x < npix ; x++) {
              // int kx  = (x < npix / 2 ? x + npix / 2 : x - npix / 2);
              long pos = y * npix + x + z * (npix * npix);
              bframe[pos][0] = fft_frame[pos][0] * fftFactor;
              bframe[pos][1] = fft_frame[pos][1] * fftFactor; //* fftFactor
          }
      }
      fftwf_execute(bps);
  }

  #pragma omp parallel for
  for (int i = 0; i < impix ; i++) {
    res_conv_map[i] = .0;
    for (int z = 0; z < nlp ; z++) {
        double q = 1; //lensKernelTrue[i * nlp + z];
        res_conv_map[i] += q * bframe[i][0];
    }
  }

}

void field::adjoint_operator_map(fftwf_complex *delta)
{
  // std::cout << "Adjoint operator nlp:" << nlp << " ngal:"  << ngal << " impix:" << impix <<  std::endl;
  float freqFactor = 2.0 * M_PI / ((float) npix); // 2.0 * M_PI / ((double) npix ) ; // 1.0 / ((double) npix); // / pixel_size

  // compute the fourier transform of gamma
  #pragma omp parallel for
  for (int z = 0; z < nlp; z++) {

    // apply lensing to shear values here
    for (long i = 0; i < impix ; i++) {
      fframe[i][0] = (float) res_gamma1_map[i];
      fframe[i][1] = (float) res_gamma2_map[i];
    }

    fftwf_execute(fps);

  	float k1, k2, k1k1, k2k2, k1k2, ksqr;
  	float c1, c2;

  	for (int y = 0; y < npix ; y++) {
      k2  = (y < npix / 2  ? y * freqFactor : (y - npix) * freqFactor); // this comes from python fft
      // k2 = (y - npix) * freqFactor;

  		for (int x = 0; x < npix ; x++) {
  			k1 = (x < npix / 2 ? x * freqFactor : (x - npix) * freqFactor); // this comes from python fft
        // k1 = (x - npix) * freqFactor;

        long pos = y * npix + x;

  			k1k1 = k1 * k1;
  			k2k2 = k2 * k2;
  			ksqr = k1k1 + k2k2;
  			c1 = k1k1 - k2k2;
  			c2 = 2.0 * k1 * k2;

        if (k1 == 0 && k2 == 0){
          delta[pos][0] = (fframe[pos][0] * c1 - fframe[pos][1] * c2);
          delta[pos][1] = (fframe[pos][0] * c2 + fframe[pos][1] * c1);
        }
        else {
          delta[pos][0] = (fframe[pos][0] * c1 - fframe[pos][1] * c2)  / ksqr;
          delta[pos][1] = (fframe[pos][0] * c2 + fframe[pos][1] * c1)  / ksqr;
        }
  		}
  	}
  }
}

void field::gradient_map(fftwf_complex *delta)
{
    // std::cout << "Compute the gradient" << std::endl;
    forward_operator_map(delta);

    // Compute residuals, note that we only differentiate for the second term here.
    #pragma omp parallel for
    for (int i = 0; i < impix ; i++) {
        double factor = std::max(1. - res_conv_map[i], 0.);
        res_gamma1_map[i] = mask[i] * (covariance_map[i]) * (factor * binned_shear[i].real() - res_gamma1_map[i]); //
        res_gamma2_map[i] = mask[i] * (covariance_map[i]) * (factor * binned_shear[i].imag() - res_gamma2_map[i]); // mask[i] *
    }

    adjoint_operator_map(delta);
}


void field::gradient_noise_map(fftwf_complex *delta)
{
    // // std::cout << "Inside gradient noise map" << std::endl;
    // for (long i = 0; i < impix; i++) {
    //     double theta1 = gsl_ran_flat(rng, 0, 2.0 * M_PI);
    //     double theta2 = gsl_ran_flat(rng, 0, 2.0 * M_PI);
    //
    //     res_gamma1_map[i] = sqrt(covariance_map[i]) * (binned_shear[i].real() * cos(theta1) - binned_shear[i].imag() * sin(theta1)) ; // * w_e[ind]
    //     res_gamma2_map[i] = sqrt(covariance_map[i]) * (binned_shear[i].imag() * cos(theta1) + binned_shear[i].real() * sin(theta1)) ; // * w_e[ind]
    //
    // }

    for (long ind = 0; ind < impix; ind++) {
      res_gamma1_map[ind] = sqrt(covariance_map[ind]) * sqrt(Sn[ind] / 2.) * gsl_ran_gaussian(rng,1.0);
      res_gamma2_map[ind] = sqrt(covariance_map[ind]) * sqrt(Sn[ind] / 2.) * gsl_ran_gaussian(rng,1.0);
    }

    for (long ind = 0; ind < impix; ind++) {
      nshear_power0[ind] += pow(res_gamma1_map[ind], 2.);
      nshear_power1[ind] += pow(res_gamma2_map[ind], 2.);
    }

    adjoint_operator_map(delta);


    // for (long ind = 0; ind < impix; ind++) {
    //   // if (Sn[ind] == 0.) {
    //   //   // pixel is masked
    //   //   res_gamma1_map[ind] =  gsl_ran_gaussian(rng,1.0); //(gsl_ran_gaussian(rng,1.0) > .0) ? 1. : -1. ; // 0.006819 *
    //   //   res_gamma2_map[ind] =  gsl_ran_gaussian(rng,1.0); //(gsl_ran_gaussian(rng,1.0) > .0) ? 1. : -1. ; // 0.006819 *
    //   // }
    //   // else{
    //   // }
    //
    //   double tmp = sqrt(Sn[ind] / 2.) * gsl_ran_gaussian(rng,1.0);
    //   while ( tmp > 1.  || tmp < -1. )
    //     tmp = sqrt(Sn[ind] / 2.) * gsl_ran_gaussian(rng,1.0);
    //
    //   res_gamma1_map[ind] = tmp;
    //
    //   tmp = sqrt(Sn[ind] / 2.) * gsl_ran_gaussian(rng,1.0);
    //   while ( tmp > 1  || tmp < -1 )
    //     tmp = sqrt(Sn[ind] / 2.) * gsl_ran_gaussian(rng,1.0);
    //
    //   res_gamma2_map[ind] = tmp;
    //
    // }

    // for (long ind = 0; ind < impix; ind++) {
    //   res_gamma1_map[ind] = res_gamma1_map[ind] / 2.;
    //   res_gamma2_map[ind] = res_gamma2_map[ind] / 2.;
    // }

    // // experiment to keep the average kappa coming from the roration of the noisy shear
    // for (long ind = 0; ind < npix * npix; ind++) {
    //     fframe[ind][0] = (float) Sn[ind];
    //     fframe[ind][1] = 0.;
    // }
    // fftwf_execute(fps);
    // for (long ind = 0; ind < npix * npix; ind++) {
    //   delta[ind][0] = fframe[ind][0];
    //   delta[ind][1] = fframe[ind][1];
    // }

}


double field::get_spectral_norm_map(int niter, double tol) {

    double norm=0;
    double norm_old=0;

    long ncoeff = npix*npix;

    fftwf_complex* kap = fftwf_alloc_complex(ncoeff);
    fftwf_complex* kap_tmp = fftwf_alloc_complex(ncoeff);

    norm= 0;
    for(long ind =0; ind< ncoeff; ind++) {
        kap[ind][0] = gsl_ran_gaussian(rng,1.0);
        kap[ind][1] = gsl_ran_gaussian(rng,1.0);
        norm += kap[ind][0] * kap[ind][0] + kap[ind][1] * kap[ind][1];
    }
    norm = sqrt(norm);

    for(long ind =0; ind< ncoeff; ind++) {
        kap[ind][0] /= norm;
        kap[ind][1] /= norm;
    }

    for(int k=0; k < niter; k++) {

        // Apply operator A^t A
        combine_components(kap, kap_tmp);
        forward_operator_map(kap_tmp);

        // Compute residuals, note that we only differentiate for the second term here
        for(int i=0; i < impix ; i++) {
            res_gamma1_map[i] = (covariance_map[i]) * res_gamma1_map[i]; // TODO:: check the weights in the binned case w_e[i] *
            res_gamma2_map[i] = (covariance_map[i]) * res_gamma2_map[i]; // (Sn[i] / 2.) *
        }

        // Apply adjoint operator on normalized results
        adjoint_operator_map(kap_tmp);
        combine_components_inverse(kap_tmp,kap);

        // Compute norm
        norm= 0;
        for(long ind =0; ind< ncoeff; ind++) {
            norm += kap[ind][0] * kap[ind][0] + kap[ind][1] * kap[ind][1];
        }

        norm = sqrt(norm);

        if( fabs(norm - norm_old)/norm <= tol)
            break;

        for(long ind =0; ind< ncoeff; ind++) {
            kap[ind][0] /= norm;
            kap[ind][1] /= norm;
        }

        norm_old = norm;
        //std::cout <<  "Iter : " << k << " Value " << norm <<std::endl;
        if(k == niter -1 ) std::cout << "Warning, reached maximum number of iterations field" << std::endl;
    }

    fftwf_free(kap);
    fftwf_free(kap_tmp);

    return norm*(1.0+tol);

}


void field::update_covariance_map(fftwf_complex* delta)
{

    // Compute the value of the field evaluated at each galaxy position
    combine_components(delta, fft_frame);
    #pragma omp parallel for
    for (int z = 0; z < nlp; z++) {

        //Compute reduced shear correction factor
        for (int y = 0; y < npix ; y++) {
            // int ky  = (y < npix / 2 ? y + npix / 2 : y - npix / 2);
            for (int x = 0; x < npix ; x++) {
                // int kx  = (x < npix / 2 ? x + npix / 2 : x - npix / 2);
                long pos = y * npix + x + z * (npix * npix);
                bframe[pos][0] = fft_frame[pos][0] * fftFactor;
                bframe[pos][1] = fft_frame[pos][1] * fftFactor; //* fftFactor
            }
        }
        fftwf_execute(bps);
    }

    #pragma omp parallel for
    for (int i = 0; i < impix ; i++) {
        res_conv_map[i] = 0;
        for (int z = 0; z < nlp ; z++) {
            double q = 1; // lensKernelTrue[i * nlp + z];
            // std::cout << "\t" << bframe[i][0] ;
            res_conv_map[i] += q * bframe[i][0] ;
        }
    }

    for(int i=0; i < impix ; i++) {
        double factor = std::max(1.0 - res_conv_map[i],0.3);
        // covariance_map[i] = pow(Sn[i],2) / (factor*factor);
        covariance_map[i] = 1.0 / (factor*factor);
        // covariance_map[i] = (Sn[i]== 1e-9 ? 1 : pow(Sn[i],2)) / (factor*factor);
    }

    // std::cout << "covariance_map" << '\n';
    // for(int i=0; i < 20 ; i++)
    //   std::cout << " " << std::max(1.0 - res_conv_map[i],0.3) << ' ';
    // std::cout << '\n';
    // std::cout << "Cov map ypologizetai: " << covariance_map[0] << " , " << covariance_map[0]  <<'\n';
}


// padding_size = config.get<int>("field.padding", 0);
// std::string coordinates_unit_str = config.get("field.units", "radian");
// convert_coordinates_unit = 1.0;
// if (coordinates_unit_str.find("radian") != std::string::npos) {
//     convert_coordinates_unit = 1.0;
// } else if (coordinates_unit_str.find("arcsec") != std::string::npos) {
//     convert_coordinates_unit = M_PI / 180.0 / 60.0 / 60.0;
// } else if (coordinates_unit_str.find("arcmin") != std::string::npos) {
//     convert_coordinates_unit = M_PI / 180.0 / 60.0;
// } else if (coordinates_unit_str.find("degree") != std::string::npos) {
//     convert_coordinates_unit = M_PI / 180.0;
// } else {
//     std::cout << "Unknown coordinates units, assuming radians." << std::endl;
// }
// pixel_size = config.get<double>("field.pixel_size") * convert_coordinates_unit;


// pos1.push_back(surv->get_pos1(ind));
// pos2.push_back(surv->get_pos2(ind));
// std::cout << "Binning the shear" << std::endl;
//index1d = compute_index_map(); // returns the vectorized 2d indexes
//binned_shear = bincount(index1d, shear_gamma1, shear_gamma2);

// Initialize the lensing planes, with one nfft per plane
// ps = (nfft_plan **) malloc((nlp) * sizeof(nfft_plan *));
/*
fframe = (fftw_complex **) malloc((nlp) * sizeof(fftw_complex * ));
bframe = (fftw_complex **) malloc((nlp) * sizeof(fftw_complex * ));
fps = (fftw_plan *) malloc((nlp) * sizeof(fftw_plan ));
bps = (fftw_plan *) malloc((nlp) * sizeof(fftw_plan ));
for (int i = 0; i < nlp; i++) {
  fframe[i] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * impix);
  bframe[i] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * impix);
  fps[i] = fftw_plan_dft_2d(npix, npix, fframe[i], fframe[i], FFTW_FORWARD, FFTW_MEASURE);
  bps[i] = fftw_plan_dft_2d(npix, npix, bframe[i], bframe[i], FFTW_BACKWARD, FFTW_MEASURE);
}
std::cout << "Construct field object" << std::endl;
*/

//fframe = fftw_alloc_complex(npix * npix );
//bframe = fftw_alloc_complex(npix * npix );
//fframe = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * impix);
//bframe = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * impix);
/*
for (int i = 0; i < nlp; i++) {
    // Create the nfft plan
    ps[i] = new nfft_plan;
    nfft_init_2d(ps[i], npix, npix, ngal);
    // Set up the nodes at the galaxy positions
    for (long ind = 0; ind < ngal;  ind++) {
        double ra  = surv->get_ra(ind);
        double dec = surv->get_dec(ind);
        double denom = cos(center_dec) * cos(dec) * cos(ra  - center_ra) + sin(center_dec) * sin(dec);
        double X =  cos(dec) * sin(ra  - center_ra) / denom;
        double Y = (cos(center_dec) * sin(dec) - cos(dec) * sin(center_dec) * cos(ra - center_ra)) / denom;

        double val = -0.5 + ((X) / size);
        val = val < -0.5 ? val + 1.0 : val;
        ps[i]->x[2 * ind]  = val ;
        val = -0.5 + ((Y) / size);
        val = val < -0.5 ? val + 1.0 : val;
        ps[i]->x[2 * ind + 1] = val;
    }
    // precompute psi, the entries of the matrix B
    nfft_precompute_one_psi(ps[i]);
    if (nfft_check(ps[i])) {
        std::cout << "Problem " << nfft_check(ps[i]) << std::endl;
    }
}
*/
/*
P = (double *) malloc(sizeof(double) * nlp * nlp);
PP= (double *) malloc(sizeof(double) * nlp * nlp);
iP= (double *) malloc(sizeof(double) * nlp * nlp);
if (nlp == 1) {
     // If the lens redshift wasn't provided, use unit weights
    if(zlens <= 0){
        for(long ind =0; ind < ngal*nlp; ind++){lensKernel[ind] = 1. ;}
        for(long ind =0; ind < ngal*nlp; ind++){lensKernelTrue[ind] = 1. ;}
    }else{
        // In the 2D case, the lensing kernel is just a lensing weight based on the
        // critical surface mass density.
        compute_surface_lensing_kernel();
    }

    P[0] = 1.;
    PP[0]= 1.;
    iP[0]= 1.;
} else {
    std::cout << "Starting computation of lensing kernels" <<std::endl;
    // Compute the full 3D lensing kernel, to reconstruct the 3D density contrast
    compute_3D_lensing_kernel();
    std::cout << "Done "<<std::endl;
}
*/

// Compute the ratio of shear and flexion variance if necessary
// sig_frac = 1.0;
// if (include_flexion) {
//     double shear_mean    = 0;
//     double flexion_mean  = 0;
//     double shear_sigma   = 0;
//     double flexion_sigma = 0;
//     for (int i = 0; i < ngal; i++) {
//         shear_mean   += shear_gamma1[i];
//         flexion_mean += flexion_f1[i];
//     }
//     shear_mean /= ngal;
//     flexion_mean /= ngal;
//     for (int i = 0; i < ngal; i++) {
//         shear_sigma   += pow(shear_gamma1[i] - shear_mean, 2.0);
//         flexion_sigma += pow(flexion_f1[i] - flexion_mean, 2.0);
//     }
//     shear_sigma   /= (ngal - 1.0);
//     flexion_sigma /= (ngal - 1.0);
//
//     sig_frac = flexion_sigma / shear_sigma;
// }


// void field::get_pixel_coordinates(double* ra, double* dec)
// {
//     for(int i=0; i < npix; i++){
//         double x = (i + 0.5) * pixel_size - size/2.;
//
//         for(int j=0; j < npix; j++){
//             double y = (j + 0.5) * pixel_size - size/2.;
//
//             double z = sqrt(x*x + y*y);
//             double c = atan(z);
//
//             double delta = asin(cos(c) * sin(surv->get_center_dec()) + y / z * sin(c) * cos(surv->get_center_dec()));
//
//             double denom = z * cos(surv->get_center_dec()) * cos(c) - y * sin(surv->get_center_dec()) * sin(c);
//
//             double alpha = surv->get_center_ra() + atan2(x * sin(c), denom);
//
//             ra[j * npix +  i] = alpha / ( M_PI / 180.0 );
//             dec[j * npix + i] = delta / ( M_PI / 180.0 );
//         }
//     }
// }

/*
void field::compute_3D_lensing_kernel()
{
    gsl_interp **interpolators;
    gsl_interp_accel **accelerators;
    gsl_integration_workspace **w;
    double *x;
    double **y;

    // First step, compute the lensing efficiency kernel on an interpolation table
    // to speed up the computation
    int nzsamp = ZMAX*100;
    x = (double *) malloc(sizeof(double) * nzsamp);
    y = (double **) malloc(sizeof(double*) * nlp);
    w = (gsl_integration_workspace **) malloc(sizeof(gsl_integration_workspace *) * nlp);

    for(int z=0; z < nlp; z++){
        y[z] = (double *) malloc(sizeof(double) * nzsamp);
	      w[z] = gsl_integration_workspace_alloc(2048);
    }

    // Initialize the x array of redshift sampling of the lensing efficiency kernel
    for(int i=0; i < nzsamp; i++){
        x[i] = ZMAX/((double) (nzsamp - 1) )*i;
    }

    int_for_3d_efficiency_params intpar;
    intpar.model = model;
    intpar.err   = err;

    gsl_function F;
    F.function = &int_for_3d_efficiency;

    for(int i=0; i < nzsamp; i++){
        double a   = 1.0/(1.0 + x[i]);
        intpar.w_a = nicaea::w ( model, a, 0, err ); quitOnError ( *err, __LINE__, stderr );

        F.params = (void *) &intpar;
	       #pragma omp parallel for
         for (int z=0; z < nlp; z++){

      	    double result;
      	    double abserr;
            gsl_integration_qags(&F, 1.0/(zlp_up[z] + 1), 1.0/(zlp_low[z] + 1),
                                            0, 1.0e-5, 2048, w[z], &result, &abserr);
            y[z][i] = result;
        }
    }

    // Interpolation table for each z and integrate over p(zsamp) for each galaxy
    interpolators = (gsl_interp **) malloc(sizeof(gsl_interp *)*nlp);
    accelerators  = (gsl_interp_accel **) malloc(sizeof(gsl_interp_accel *)*nlp);
    for(int z=0; z < nlp; z++){
        interpolators[z] = gsl_interp_alloc(gsl_interp_cspline, nzsamp);
        accelerators[z]  = gsl_interp_accel_alloc();
        gsl_interp_init(interpolators[z], x, y[z], nzsamp);
    }

    // Deactivate default gsl error handling
    gsl_error_handler_t * handler =  gsl_set_error_handler_off();

    // Compute lensing efficiency kernel for each galaxy by marginalising over pdf
    for (int i = 0; i < ngal; i++) {
        if(i % 100 == 0) std::cout  << "Processed " << i << "/" << ngal << " galaxies\r" << std::flush;
        redshift_distribution * redshift=surv->get_redshift(i);

	#pragma omp parallel for
  for(int z=0; z <nlp; z++){
	    double result;
	    double abserr;

	    int_for_marginalisation_params params;
	    params.x = x;
	    gsl_function G;
	    G.function = &int_for_marginalisation;
            params.y = y[z];
            params.accelerator = accelerators[z];
            params.interpolator = interpolators[z];
            params.redshift = redshift;
            G.params = (void *) &params;
            int ret_code = gsl_integration_qags(&G, std::max(0., redshift->get_zmin()),
                                     std::min(ZMAX, redshift->get_zmax()), 0, 1.0e-4, 1024, w[z], &result, &abserr);
            // If standard gsl integration fails, falls back to trapezoid method
            if(ret_code != 0){
                //std::cout << "Warning: Integration of galaxy redshift pdf did not converge. Galaxy id: " << i << "; value: " << result << std::endl;
                double a = std::max(0., redshift->get_zmin());
                double b = std::min(ZMAX, redshift->get_zmax());
                long n = 1024;

                double h = (b - a)/((double) n);
                result = 0;
                for(long it=0; it < (n - 1); it++){
                    double x1 = it*h + a;
                    double x2 = (it+1)*h + a;
                    result += int_for_marginalisation(x1, (void *) &params) + int_for_marginalisation(x2, (void *) &params);
                }
                result *= h/2.;
            }
            lensKernel[i * nlp + z] = std::max(result, 0.);

        }

    }

    // Reset default gsl error handling
    gsl_set_error_handler(NULL);
    std::cout  << "Processed " << ngal << "/" << ngal << " galaxies" << std::endl;

#ifdef DEBUG_FITS
    dblarray toto;
    toto.alloc(lensKernel, nlp, ngal, 1);
    fits_write_dblarr("lensKernel.fits",toto);
#endif

    // Free all unnecessary arrays
    for(int z=0; z < nlp; z++){
      gsl_interp_free(interpolators[z]);
      gsl_interp_accel_free(accelerators[z]);
      gsl_integration_workspace_free(w[z]);
      free(y[z]);
    }
    free(x);
    free(y);
    free(w);

    // Apply SVD regularisation to the lensing operator
    int nsmall = std::min( ngal, 20000l );
    arma::mat A ( ngal, nlp );
    arma::mat Asmall ( nsmall, nlp );
    arma::mat U;
    arma::vec s;
    arma::mat V;

    // Select random indices out of the entire survey to compute Asmall
    std::vector<int> indices;
    for (int i=0; i<ngal; ++i) indices.push_back(i);
    std::random_shuffle(indices.begin(), indices.end());

    for ( long int ind =0; ind < ngal; ind++ )
        for ( int z=0; z < nlp; z++ ) {
            if ( ind < nsmall )   Asmall (ind, z ) = lensKernel[indices[ind]*nlp + z];

            A ( ind,z ) = lensKernel[ind*nlp + z];
            lensKernelTrue[ind*nlp + z] = lensKernel[ind*nlp + z];
        }

    // Compute the pseudo-inverse of the lensing kernel so that we can have much faster convergence
    svd ( U, s, V, Asmall );

    // Compute the preconditionning matrix from s and V
    // First, regularise s
    double maxS = s ( 0 );
    s.resize ( nlp );
    for ( int i=0; i< nlp; i++ ) {
        if ( s ( i ) > r_cond*maxS ) {
            s ( i ) = 1.0/s ( i );
        }
        else {
            s ( i ) = 1.0/ ( r_cond*maxS );
        }
        if ( i >= nsmall ) {
            s ( i ) = 1.0;
        }
    }

    // Now build the preconditionning matrix
    arma::mat Pr = V * diagmat ( s ) * V.t();
    arma::mat PPr = Pr*Pr.t();

    arma::mat IP = inv ( Pr );

    // Build the conditionned Q matrix
    arma::mat AP = A * Pr;

    // Extract the preconditionning matrix and conditionned tomographic lensing operator
    for ( long int ind =0; ind < ngal; ind++ ) {
        for ( int z=0; z < nlp; z++ ) {
            lensKernel[ind*nlp + z] = AP ( ind,z );
        }
    }

    // Saves the preconditionning matrix for later
    for ( int z1=0; z1 < nlp; z1++ ) {
        for ( int z2=0; z2 < nlp; z2++ ) {
            P [z2*nlp + z1] = Pr ( z1,z2 );
            iP[z2*nlp + z1] = IP ( z1,z2 );
            PP[z2*nlp + z1] = PPr(z1,z2);
        }
    }

#ifdef DEBUG_FITS
    dblarray top;
    top.alloc(P,nlp,nlp,1);

    dblarray toip;
    toip.alloc(iP,nlp,nlp,1);

    dblarray topp;
    topp.alloc(PP,nlp,nlp,1);

    fits_write_dblarr("P.fits", top);
    fits_write_dblarr("IP.fits", toip);
    fits_write_dblarr("PP.fits", topp);
    fits_write_dblarr("QP.fits", toto);
#endif

}
*/
/*
void field::compute_surface_lensing_kernel()
{
    double result;
    double abserr;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(2048);

    gsl_function F;
    F.function = &int_for_sigma;

    double a_inf  = 1.0 / (1.0 + Z_INF);
    double a_lens = 1.0 / (1.0 + zlens);

    double w_l   = nicaea::w(model, a_lens, 0, err);
    double w_inf = nicaea::w(model, a_inf, 0, err);

    int_for_sigma_params params;
    params.model = model;
    params.err   = err;
    params.w_l   = w_l;
    params.w_inf = w_inf;

    // TODO::compute the new lensing kernel coefficients, not active for the moment
    std::cout << "Inside compute surface lensing kernel" << std::endl;
    // this is the spot where redshifts should be binned, I need the 1dindex
    double * redshifts = (double * ) malloc(sizeof(double) * ngal);
    for (long i = 0; i < surv->get_ngal(); i++) {
      redshift_distribution *redshift = surv->get_redshift(i);
      redshifts[i] = redshift->get_zmax();
    }
    // std::vector<double> binned_redshift = bincount(index1d, redshifts);

    std::cout << "Binned redshift follows:" << '\n';
    for (long i = 0; i < 200; i++) {
      std::cout << redshifts[i] << ' ';
    }
    std::cout << '\n';

    for (long i = 0; i < surv->get_ngal(); i++) {

        redshift_distribution *redshift = surv->get_redshift(i);

        // Treat the case of spectroscopic redshifts
        if (spectroscopic_redshift *specz = dynamic_cast<spectroscopic_redshift *>(redshift)) {

            double afit = 1.0 / (1.0 + specz->get_redshift());
            double w_s = nicaea::w(model, afit, 0, err);

            if (afit >= a_lens) {
                lensKernel[i] = 0;
            } else {
                lensKernel[i] = ((w_s - w_l) * w_inf) / ((w_inf - w_l) * w_s);
            }
        } else {
            params.redshift = redshift;
            F.params = (void *) &params;
            gsl_integration_qags(&F, std::max(std::max(model->a_min, a_inf), 1./(redshift->get_zmax() + 1.)),
                                     std::min(a_lens, 1./(redshift->get_zmin() + 1.)), 0, 1.0e-5, 1024, w, &result, &abserr);
            lensKernel[i] = result;
        }
        lensKernelTrue[i] = lensKernel[i];
    }

    gsl_integration_workspace_free(w);
}
*/

/*
void field::comp_residuals(fftwf_complex *delta, std::vector<std::complex<double>> &rshear)
{
    // std::cout << "Inside field comp residuals" << std::endl;

    forward_operator_map(delta);

    // Compute residuals, note that we only differentiate for the second term here.
    // std::cout << "pushing" << std::endl;
    #pragma omp parallel for
    for (int i = 0; i < impix ; i++) {
        double factor = std::max(1. - res_conv_map[i], 0.);

        // res1.push_back(sqrt(covariance_map[i]) * (factor * binned_shear[i].real() - res_gamma1_map[i])); // TODO:: * w_e[i] decide weights
        // res2.push_back(sqrt(covariance_map[i]) * (factor * binned_shear[i].imag() - res_gamma2_map[i])); // * w_e[i]
        // sqrt(covariance_map[i]) * (factor * binned_shear[i].real() - res_gamma1_map[i])

        rshear.push_back({sqrt(covariance_map[i]) * (factor * binned_shear[i].real() - res_gamma1_map[i]), sqrt(covariance_map[i]) * (factor * binned_shear[i].imag() - res_gamma2_map[i])}); //
        // std::cout << i << " ";
    }
    // std::cout << std::endl;

    // //wiener filter the initial image
    // #pragma omp parallel for
    // for (int i = 0; i < ngal ; i++) {
    //     res1.push_back(shear_gamma1[i]);
    //     res2.push_back(shear_gamma2[i]);
    // }
    //

}

*/
/*
std::vector<int> field::compute_index_map()
{

  // std::cout << "Inside compute shear map" << std::endl;
  std::vector<double>::iterator itr;
	int index;

  // find the minimum x position
	itr = std::min_element(pos1.begin(), pos1.end());
	index  = std::distance(pos1.begin(), itr);
	double xmin = pos1[index];

	// std::cout << "xmin is: " << xmin << std::endl;

	// find the maximum x position
	itr = std::max_element(pos1.begin(), pos1.end());
	index = std::distance(pos1.begin(), itr);
	double xmax = pos1[index];

  // std::cout << "xmax : " << xmax << std::endl;

	// find the minimum y position
	itr = std::min_element(pos2.begin(), pos2.end());
	index = std::distance(pos2.begin(), itr);
	double ymin = pos2[index];

	//printf("ymin element at: %d\n", index);
  //std::cout << "ymin is: " << ymin << std::endl;

	// find the maximum y position
	itr = std::max_element(pos2.begin(), pos2.end());
	index = std::distance(pos2.begin(), itr);
	double ymax = pos2[index];

	//printf("ymax element at: %d\n", index);
  //std::cout << "ymax is: " << ymax << std::endl;
	//printf("Xmin = %f, Xmax = %f\n",xmin, xmax);
	//printf("Ymin = %f, Ymax = %f\n",ymin, ymax);

	// compute the new axis
	double halfdx = (xmax - xmin) / (2 * npix - 2);
	double halfdy = (ymax - ymin) / (2 * npix - 2);
	double xlow = xmin - halfdx;
	double xhigh = xmax + halfdx;
	double ylow = ymin - halfdy;
	double yhigh = ymax + halfdy;

	//printf("xlow = %f, xhigh = %f\n",xlow, xhigh);
	// printf("ylow = %f, yhigh = %f\n",ylow, yhigh);

	// set binning edges per axis
	std::vector<double> xedges = linspace(xlow, xhigh, npix+1);
	std::vector<double> yedges = linspace(ylow, yhigh, npix+1);

  // std::cout << "xedge len " << xedges.size() << std::endl;
	// get the binning indexes
	std::vector<int> indx = digitize(pos1, xedges);
	std::vector<int> indy = digitize(pos2, yedges);

  // std::cout << "indx len " << indx.size() << std::endl;
	// create the 1d indexing vector
	std::vector<int> ind_1d;
	ind_1d.reserve(indx.size());
	std::transform(indy.begin(), indy.end(), indy.begin(), std::bind1st(std::multiplies<int>(),npix));
	std::transform(indx.begin(), indx.end(), indy.begin(), std::back_inserter(ind_1d), std::plus<int>());

	// compute the vectorized binned shear residuals
  // std::cout << "before bincount" << std::endl;
  // binned_shear = bincount(ind_1d);
  //std::cout << "after bincount" << std::endl;

	// write binned_shear to a file

  // double * img1 = (double *) malloc(sizeof(double) * npix*npix);
  // double * img2 = (double *) malloc(sizeof(double) * npix*npix);
  // for (int i = 0; i < npix * npix; i++) {
  //   img1[i] = binned_shear[i].real();
  //   img2[i] = binned_shear[i].imag();
  // }
  // dblarray expts1;
  // dblarray expts2;
  // expts1.alloc(img1, npix, npix);
  // expts2.alloc(img2, npix, npix);
  // fits_write_dblarr("../data/flagship/gamma1_map.fits", expts1);
  // fits_write_dblarr("../data/flagship/gamma2_map.fits", expts2);


  return ind_1d;
}
*/
/*
// bincount for redshifts
std::vector<double> field::bincount(std::vector<int> ind1d, double * value1)
{
  // int impix = npix * npix;
	std::vector<double> summ1(impix, 0.0);
	std::vector<double> result;
	int mycount;

	// for (unsigned i=0; i < impix; i++) {
	// 	mycount = std::count(ind1d.begin(), ind1d.end(), i);
	// 	counter[i] = mycount;
	// }

	// remove the zero values
	// for (unsigned i=0; i < impix; i++) {
	// 	if(counter[i] == 0)
	// 		counter[i] = 1;
	// }

  // std::cout << "Inside field bincount, printing counter" << std::endl;
  // for (unsigned i=0; i < 200; i++)
  //   std::cout << counter[i] << " ";
  // std::cout << std::endl;

  // summ the residual values
	for (unsigned i=0; i < ngal; i++) {
		summ1[ind1d[i]] += value1[i];
	}

	for (unsigned i=0; i < impix; i++)
		result.push_back(summ1[i]/counter[i]);

  // std::cout << "result" << result[0] << std::endl;
  // std::cout << "result" << result[1] << std::endl;
	return result;
}
*/

/*
std::vector<std::complex<double>> field::bincount(std::vector<int> ind1d, double * value1, double * value2)
{
  int impix = npix * npix;
	std::vector<double> summ1(impix, 0.0);
	std::vector<double> summ2(impix, 0.0);
	std::vector<std::complex<double>> result;

	int mycount;
	// int vl = std::distance(ind1d.begin(),ind1d.end());

  // std::cout << "inside bincount " << std::endl;

  // count the occurences
  // std::cout << "mycount: " << std::endl;
	for (unsigned i=0; i < impix; i++) {
		mycount = std::count(ind1d.begin(), ind1d.end(), i);
    // std::cout << mycount << " ";
		counter[i] = mycount;
	}
  // std::cout << std::endl;
  // std::cout << "inside bincount - mycount " << mycount << std::endl;

	// remove the zero values
	for (unsigned i=0; i < impix; i++) {
		if(counter[i] == 0)
			counter[i] = 1;
	}

  // std::cout << "Inside field bincount, printing counter" << std::endl;
  // for (unsigned i=0; i < 200; i++)
  //   std::cout << counter[i] << " ";
  // std::cout << std::endl;

  //double q = 0;
  // for (int z = 0; z < nlp ; z++) {
  //    q += lensKernel[i * nlp + z];

  // summ the residual values
	for (unsigned i=0; i < ngal; i++) {
		summ1[ind1d[i]] += value1[i];
		summ2[ind1d[i]] += value2[i];
	}
  // std::cout << "inside bincount - summ1 " << summ1[0] << std::endl;
  // std::cout << "inside bincount - summ2 " << summ2[1] << std::endl;

	for (unsigned i=0; i < impix; i++)
		result.push_back({summ1[i]/counter[i], summ2[i]/counter[i]});

  // std::cout << "result" << result[0] << std::endl;
  // std::cout << "result" << result[1] << std::endl;
	return result;
}
*/

/*
 * Return the indices of the bins to which each value in input array belongs
 */
/*
std::vector<int> field::digitize(std::vector<double> values, std::vector<double> edges)
{
	std::vector<int> result;
  int pos;
  // int vl = std::distance(values.begin(),values.end());

	for (unsigned i=0; i < ngal; i++) {
		auto lower = std::lower_bound(edges.begin(), edges.end(), values[i] );
		pos = std::distance(edges.begin(), lower);
		result.push_back(pos-1);
	}

	return result;
}
*/
/*
std::vector<double> field::linspace(double min, double max, int n)
{
	std::vector<double> result;

	int iterator = 0;

	for (int i = 0; i <= n-2; i++)
	{
		double tmp = min + i*(max-min)/(floor((double)n) - 1);
		result.insert(result.begin() + iterator, tmp);
		iterator += 1;
	}

	result.insert(result.begin() + iterator, max);
	return result;

  // void field::print_pos(){
  //   printf("Pos2 : \n");
  //   for (int i = 0; i <= 20; i++)
  //     printf("%f ", pos2[i]);
  //   printf("\n");
  //
}
*/
