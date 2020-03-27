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
#include <iostream>
#include <cmath>
#ifdef DEBUG_FITS
#include <sparse2d/IM_IO.h>
#endif

#include "surface_reconstruction.h"

using namespace std;

surface_reconstruction::surface_reconstruction(boost::property_tree::ptree config, field *fi) // , dblarray thr_array
{
    f = fi;
    // Get options from the configuration file.
    nRecIter      = config.get<int>("parameters.niter", 500);
    nRecIterDebias= config.get<int>("parameters.niter_debias", 500);
    nscales       = config.get<int>("parameters.nscales", 4);
    lambda        = config.get<double>("parameters.lambda", 4.0);
    nrandom       = config.get<int>("parameters.nrandom", 1000.0);
    nreweights    = config.get<int>("parameters.nreweights", 5);
    positivity    = config.get<bool>("parameters.positivity", false);
    double bl_reg = config.get<double>("parameters.battle_lemarie_reg", 0.1);
    double ls_reg = config.get<double>("parameters.last_scale_reg", 0.1);


    cout << "Using following reconstruction parameters :" << std::endl;
    cout << "Number of scales         : " << nscales << endl;
    cout << "Regularisation parameter : " << lambda  << endl;

    // Get characteristics of the field
    npix = f->get_npix();

    // Allocate wavelet transform
    cout << "Allocate wavelet transform, " << npix << endl;

    wav = new wavelet_transform(npix, nscales);
    cout << "Surf recon: wavelet object Allocated " << endl;

    // Effective number of wavelet frames
    nframes = wav->get_nframes();
    cout << "Surf recon: nframes " << nframes << endl;

    ncoeff = npix * npix ;   // The factor 2 is to include both shear and flexion
    nwavcoeff = npix * npix * nframes;

    // Allocating internal arrays
    kappa       = fftwf_alloc_complex(ncoeff);
    kappa_u     = fftwf_alloc_complex(ncoeff);
    kappa_rec   = fftwf_alloc_complex(ncoeff);
    kappa_old   = fftwf_alloc_complex(ncoeff);
    kappa_grad  = fftwf_alloc_complex(ncoeff);
    kappa_tmp   = fftwf_alloc_complex(ncoeff);
    kappa_trans = fftwf_alloc_complex(ncoeff);
    alpha       = (float *) malloc(sizeof(float) * nwavcoeff);
    alpha_u     = (float *) malloc(sizeof(float) * nwavcoeff);
    alpha_res   = (float *) malloc(sizeof(float) * nwavcoeff);
    alpha_tmp   = (float *) malloc(sizeof(float) * nwavcoeff);
    thresholds  = (float *) malloc(sizeof(float) * nwavcoeff);
    weights     = (float *) malloc(sizeof(float) * nwavcoeff);
    support     = (float *) malloc(sizeof(float) * nwavcoeff);

    // Initialise internal arrays
    for (long ind = 0; ind < ncoeff; ind++) {
        kappa[ind][0]     = 0;
        kappa[ind][1]     = 0;
        kappa_u[ind][0]   = 0;
        kappa_u[ind][1]   = 0;
        kappa_old[ind][0] = 0;
        kappa_old[ind][1] = 0;
        kappa_rec[ind][0] = 0;
        kappa_rec[ind][1] = 0;
        kappa_grad[ind][0] = 0;
        kappa_grad[ind][1] = 0;
        kappa_tmp[ind][0] = 0;
        kappa_tmp[ind][1] = 0;
    }
    // cout << "sdf12" << endl;
    for (long ind = 0; ind < nwavcoeff; ind++) {
        alpha[ind]     = 0;
        alpha_u[ind]   = 0;
        alpha_res[ind] = 0;
        alpha_tmp[ind] = 0;
        thresholds[ind]= 1;
        weights[ind]   = 1;
        support[ind]   = 1;
    }

    std::cout << "Surf recon: kappas initiallized " << endl;


    // Normalization factor for the fft
    fftFactor     = 1.0 / (((double)npix) * npix);
    fft_frame     = fftwf_alloc_complex(ncoeff);
    plan_forward  = fftwf_plan_dft_2d(npix, npix, fft_frame, fft_frame, FFTW_FORWARD,  FFTW_MEASURE);
    plan_backward = fftwf_plan_dft_2d(npix, npix, fft_frame, fft_frame, FFTW_BACKWARD, FFTW_MEASURE);

    // Initialize the threshold levels, with lower thresholds on larger scales
    sigma_thr = (double *) malloc(sizeof(double) * nframes);
    for (int i = 0; i < nscales - 1; i++) {
        sigma_thr[i] = lambda * sqrt(2 * log(npix / pow(2.0, i) * npix / pow(2.0, i))) / sqrt(2 * log(npix * npix));
    }
    // cout << "sdf13" << endl;
    // Special regularisation for the smooth approximation
    sigma_thr[nscales - 1] = ls_reg;

    // // Additional regularisation for additional BL frames
    // for (int i = nscales; i < nscales + 3; i++) {
    //     sigma_thr[i] = bl_reg;
    // }

    mu1 = get_spectral_norm_prox(100, 1e-7);

    // std::cout << "To mu1 ypologizetai: " << mu1 << '\n';
    // mu2 = f->get_spectral_norm(200, 1e-7);
    sig = 1.0 / mu1;
    std::cout << "sig: " << sig << std::endl;
    // tau = 0.9 / (mu2 / 2.0 + sig * mu1);


    // // read the thresholds from input array
    // for (int j = 0; j < nframes; j++) {
    //   for (long ind = 0; ind < npix * npix; ind++) {
    //       int pos = j * npix * npix + ind;
    //       thresholds[pos] = thr_array(pos) * .2;
    //   }
    // }

}

surface_reconstruction::~surface_reconstruction()
{
    std::cout << "before deleting wav object " << std::endl;
    fftwf_free(kappa);
    fftwf_free(kappa_u);
    fftwf_free(kappa_rec);
    fftwf_free(kappa_old);
    fftwf_free(kappa_grad);
    fftwf_free(kappa_tmp);
    fftwf_free(kappa_trans);

    free(sigma_thr);
    free(alpha);
    free(alpha_u);
    free(alpha_res);
    free(alpha_tmp);
    free(thresholds);
    free(weights);

    // fftw_free(fft_frame);
    fftwf_destroy_plan(plan_forward);
    fftwf_destroy_plan(plan_backward);
    std::cout << "before deleting wav object " << std::endl;
    delete wav; // can I not delete this object? corrupted size vs. prev_size error
    std::cout << "wav object deleted " << std::endl;
}

void surface_reconstruction::run_main_iteration(long int niter, bool debias)
{
    mu2 = f->get_spectral_norm_map(200, 1e-7);
    // std::cout << "Mu2 is " << mu2 << std::endl;
    // tau =  .05 / (mu2 / 2 + sig * mu1 );
    tau =  .5 / (mu2 / 2 + sig * mu1);
    std::cout << "step size main iteration : " << tau << std::endl;

#ifdef DEBUG_FITS
    dblarray rec_kappa(f->get_npix(), f->get_npix());
    char name[256];
#endif
    std::cout << "Step size : " << tau << std::endl;
    // std::cout << "mu2 : " << mu2 << std::endl;
    // std::cout << "sig : " << sig << std::endl;
    // std::cout << "mu1 : " << mu1 << std::endl;

    for (int iter = 0; iter < niter; iter++) {
        if (debias && iter % 100 == 0) {
            // tau *= .5;
            std::cout << "Iteration :" << iter << std::endl;
        }

        // Copy kappa for computing gradient step
        for (long ind = 0; ind < ncoeff; ind++) {
            kappa_grad[ind][0] = kappa[ind][0];
            kappa_grad[ind][1] = kappa[ind][1];
            kappa_old[ind][0]  = kappa[ind][0];
            kappa_old[ind][1]  = kappa[ind][1];
        }

        f->gradient_map(kappa_grad);

        // std::cout << "Kappa_grad: " << kappa_grad[0][0] << " , " << kappa_grad[0][1] << std::endl;

        // Reconstructing from wavelet coefficients
        // wav->trans_adjoint(alpha, kappa_trans);
        wav->trans_adjoint_gen2(alpha, kappa_trans);
        // std::cout << "Kappa_trans: " << kappa_trans[0][0] << " , " << kappa_trans[0][1] << std::endl;

        f->combine_components_inverse(kappa_trans, kappa_u);

        // std::cout << "Kappa_u: " << kappa_u[0][0] << " , " << kappa_u[0][1] << std::endl;

        // Updating kappa
        for (long ind = 0; ind < ncoeff; ind++) {
            kappa[ind][0] += tau * (kappa_grad[ind][0] - kappa_u[ind][0]) ;
            kappa[ind][1] += tau * (kappa_grad[ind][1] - kappa_u[ind][1]) ;
        }

        // std::cout << "Kappa: " << kappa[0][0] << " , " << kappa[0][1] << std::endl;

        // Here is the place to compute the prox of the E mode constraint
        f->combine_components(kappa, kappa_tmp);
        for (long ind = 0; ind < npix * npix; ind++) {
            fft_frame[ind][0] = kappa_tmp[ind][0] * fftFactor;
            fft_frame[ind][1] = kappa_tmp[ind][1] * fftFactor;
        }
        fftwf_execute(plan_backward);

        if (positivity) {

            for (long ind = 0; ind < npix * npix; ind++) {
                fft_frame[ind][0] = max(fft_frame[ind][0], 0.0f);
                fft_frame[ind][1] = 0;
            }
        } else {
            // std::cout << "positivity is false" << '\n';
            for (long ind = 0; ind < npix * npix; ind++) {
                fft_frame[ind][1] = 0;
            }
        }

        fftwf_execute(plan_forward);
        for (long ind = 0; ind < npix * npix; ind++) {
            kappa_tmp[ind][0] = fft_frame[ind][0];
            kappa_tmp[ind][1] = fft_frame[ind][1];
        }
        f->combine_components_inverse(kappa_tmp, kappa);
        /////////////////////////////////////////////////////////

        for (long ind = 0; ind < ncoeff; ind++) {
            kappa_tmp[ind][0] = 2 * kappa[ind][0] - kappa_old[ind][0];
            kappa_tmp[ind][1] = 2 * kappa[ind][1] - kappa_old[ind][1];
        }

        f->combine_components(kappa_tmp, kappa_trans);
        wav->transform_gen2(kappa_trans, alpha_u);
        // wav->transform(kappa_trans, alpha_u);

        for (int j = 0; j < nframes; j++) {
            for (long ind = 0; ind < npix * npix; ind++) {
                double dum = alpha[j * npix * npix + ind] + sig * alpha_u[j * npix * npix + ind];

                if (debias) {
                    if (j == (nscales - 1)) {
                        alpha[j * npix * npix + ind] = 0;
                    } else {
                        alpha[j * npix * npix + ind] = dum - dum * support[j * npix * npix + ind] ;
                    }
                } else {
                    double val = dum - copysign(max(fabs(dum) - sigma_thr[j] * thresholds[j * npix * npix + ind] * weights[j * npix * npix + ind], 0.0), dum);
                    support[j * npix * npix + ind] = fabs(val) < fabs(dum) ? 1 : 0;
                    alpha[j * npix * npix + ind] = val;
                }
            }
        }

#ifdef DEBUG_FITS
        if (debias){
          get_convergence_map(rec_kappa.buffer());
          sprintf(name, "kappa_%03d.fits", iter);
          fits_write_dblarr(name, rec_kappa);
        }
#endif
    }

#ifdef DEBUG_FITS
    // export confidence intervals using noise
    if (debias) {
      // initialise confidence intervals
      double *sglimpse = (double *) malloc(sizeof(double) * npix * npix);
      for (long ind = 0; ind < npix * npix; ind++)
        sglimpse[ind] = .0;

      // sum the noise standard deviation of the support per scale
      for (int j = 0; j < nframes; j++) {
        for (long ind = 0; ind < npix * npix; ind++) {
          sglimpse[ind] += sigma_thr[j] * thresholds[j * npix * npix + ind] * support[j * npix * npix + ind];
        }
      }
      dblarray expt1;
      expt1.alloc(sglimpse, npix, npix);
      fits_write_dblarr("sglimpse.fits", expt1);


      // export sum of wavelet coefficients energy at support set
      double *skglimpse = (double *) malloc(sizeof(double) * npix * npix);
      for (long ind = 0; ind < npix * npix; ind++)
        skglimpse[ind] = .0;
      wav->transform_gen2(kappa, alpha);
      for (int j = 0; j < nframes; j++)
        for (int y = 0; y < npix ; y++)
            for (int x = 0; x < npix ; x++)
              skglimpse[x * npix + y] += pow(alpha[j * npix * npix + x * npix + y], 2.0) * support[j * npix * npix + x * npix + y];
      dblarray expt2;
      expt2.alloc(skglimpse, npix, npix);
      fits_write_dblarr("skglimpse.fits", expt2);

      // export support
      fltarray sup;
      sup.alloc(support, npix,npix,(float)nframes);
      fits_write_fltarr("support.fits", sup);

    }
#endif

}

void surface_reconstruction::reconstruct(std::vector<std::complex<double>> rws)
{
    // std::cout << "Inside reconstruct function" << std::endl;

    // bin the shear
    // std::vector<std::complex<double>> binned_shear = f->get_binned_shear();

    //std::vector<std::complex<double>> binned_shear(npix*npix,{.0,.0});
    // or get the shear residuals from
    // int i = 0;
    // for (int y = 0; y < npix ; y++) {
    //   for (int x = 0; x < npix ; x++) {
    //     f->binned_shear[i] = {-rws[x*npix + y].real(), rws[x*npix+y].imag()};
    //     i++;
    //   }
    // }

    // for (int y = 0; y < npix ; y++) { // not working
    //     for (int x = 0; x < npix ; x++) {
    //         long pos = (npix - y - 1) * npix + (npix - x - 1);
    //         f->binned_shear[x * npix + y] = rws[pos];
    //     }
    // }


    for (int i = 0; i < npix * npix; i++)
      f->binned_shear[i] = {real(rws[i]), imag(rws[i])};

    // std::cout << "binned_shear" << rws[0] << std::endl;
    // std::cout << "Binning finished" << std::endl;

  std::cout << "Computing thresholds" << std::endl;
  compute_thresholds(nrandom);

#ifdef DEBUG_FITS
    // Saves the thresholds
    // float * thrf = (float *) malloc(sizeof(float) * nwavcoeff);
    // for (long ind = 0; ind < nwavcoeff; ind++)
    //   thrf[ind] = (float) thresholds[ind];

    fltarray thr;
    thr.alloc(thresholds, npix,npix,nframes);
    fits_write_fltarr("thresholds.fits", thr);

    // Computes a back projection
    for (long ind = 0; ind < ncoeff; ind++) {
        kappa_grad[ind][0] = 0;
        kappa_grad[ind][1] = 0;
    }
    f->gradient_map(kappa_grad);

    for (long ind = 0; ind < npix * npix; ind++) {
        fft_frame[ind][0] = kappa_grad[ind][0] * fftFactor;
        fft_frame[ind][1] = kappa_grad[ind][1] * fftFactor;
    }

    fftwf_execute(plan_backward);

    dblarray rec_kappa(f->get_npix(), f->get_npix());
    for (int y = 0; y < npix ; y++) {
        for (int x = 0; x < npix ; x++) {
            long pos = x * npix + y; // (npix - y - 1) * npix + (npix - x - 1);
            rec_kappa.buffer()[pos] = fft_frame[pos][0];
        }
    }
    fits_write_dblarr("back_proj.fits", rec_kappa);

#endif

    f->check_adjoint();
    // std::cout << "Running main iteration" << std::endl;
    run_main_iteration(nRecIter);

    // Reweighted l1 loop
    for (int i = 0; i < nreweights ; i++) {
         f->update_covariance_map(kappa);
         compute_thresholds(nrandom);

#ifdef DEBUG_FITS
   // Saves the thresholds
   char namef[256];
   float * thrf = (float *) malloc(sizeof(float) * nwavcoeff);
   for (long ind = 0; ind < nwavcoeff; ind++)
     thrf[ind] = (float) thresholds[ind];

   fltarray thr;
   thr.alloc(thrf, npix,npix,(float)nframes);
   sprintf(namef, "thresholds_%03d.fits", i);
   fits_write_fltarr(namef, thr);

   // Computes a back projection
   for (long ind = 0; ind < ncoeff; ind++) {
       kappa_grad[ind][0] = 0;
       kappa_grad[ind][1] = 0;
   }
   f->gradient_map(kappa_grad);

   for (long ind = 0; ind < npix * npix; ind++) {
       fft_frame[ind][0] = kappa_grad[ind][0] * fftFactor;
       fft_frame[ind][1] = kappa_grad[ind][1] * fftFactor;
   }

   fftwf_execute(plan_backward);

   dblarray rec_kappa(f->get_npix(), f->get_npix());
   for (int y = 0; y < npix ; y++) {
       for (int x = 0; x < npix ; x++) {
           long pos = x * npix + y; // (npix - y - 1) * npix + (npix - x - 1);
           rec_kappa.buffer()[x * npix + y] = fft_frame[pos][0];
       }
   }
   sprintf(namef, "back_proj_%03d.fits", i);
   fits_write_dblarr(namef, rec_kappa);

   // compute forward and backward operations


#endif

         compute_weights();
         run_main_iteration(nRecIter / 2);
     }

    std::cout  << "Starting debiasing " << std::endl;
    // Final debiasing step
    f->update_covariance_map(kappa);
    run_main_iteration(nRecIterDebias, true);

}

void surface_reconstruction::compute_thresholds(int niter)
{

    for (long ind = 0; ind < npix * npix * nframes; ind++) {
        thresholds[ind] = .0;
    }

    double * kvar = (double *) malloc(sizeof(double) * npix*npix);

    for (long ind = 0; ind < npix * npix ; ind++){
      kappa_tmp[ind][0] = .0;
      kappa_tmp[ind][1] = .0;
      kvar[ind] = .0;
    }

    for (int i = 0; i < niter; i++) {
        f->gradient_noise_map(kappa_rec);

        // compute the inverse transform
        for (long ind = 0; ind < npix * npix; ind++) {
            fft_frame[ind][0] = kappa_rec[ind][0] * fftFactor;
            fft_frame[ind][1] = kappa_rec[ind][1] * fftFactor;
        }
        fftwf_execute(plan_backward);
        for (int y = 0; y < npix ; y++) {
            for (int x = 0; x < npix ; x++) {
                long pos = y * npix + x ;
                kvar[pos] += pow(fft_frame[pos][0], 2.);
            }
        }

        for (long ind = 0; ind < npix * npix ; ind++){
          kappa_tmp[ind][0]  += pow(kappa_rec[ind][0], 2.) / (niter-1);
          kappa_tmp[ind][1]  += pow(kappa_rec[ind][1], 2.) / (niter-1);
        }

        f->combine_components(kappa_rec, kappa_trans);
        wav->transform_gen2(kappa_trans, alpha_tmp);

        // Compute gradient step
        for (long ind = 0; ind < npix * npix * nframes; ind++) {
            thresholds[ind] += pow(alpha_tmp[ind], 2.0) ;
        }

    }

    // printf("Number of wavelet frames is, %d\n", nframes);

    for (long n = 0; n < nframes; n++) {
        double maxThr = 0;
        for (long ind = 0; ind < npix * npix; ind++) {
            thresholds[n * npix * npix + ind] = sqrt(1.0 / ((double) niter) * thresholds[n * npix * npix + ind]);
            maxThr = thresholds[n * npix * npix + ind] > maxThr ? thresholds[n * npix * npix + ind] : maxThr;
        }
        for (long ind = 0; ind < npix * npix; ind++) {
            thresholds[n * npix * npix + ind] = max(thresholds[n * npix * npix + ind], (float) (maxThr * 0.1));
        }
    }

#ifdef DEBUG_FITS
    // // Saves the thresholds
    // fltarray thr;
    // thr.alloc(thresholds, npix,npix,(float)nframes);
    // fits_write_fltarr("thresholds.fits", thr);

    // save kappa rec
    // sprintf(name, "kappa_rec.fits", iter);

    double *krec0 = (double *) malloc(sizeof(double) * npix * npix);
    double *krec1 = (double *) malloc(sizeof(double) * npix * npix);

    for (long ind = 0; ind < npix * npix; ind++){
      krec0[ind] = kappa_tmp[ind][0];
      krec1[ind] = kappa_tmp[ind][1];
    }

    dblarray kr0;
    dblarray kr1;
    kr0.alloc(krec0, npix, npix);
    kr1.alloc(krec1, npix, npix);
    fits_write_dblarr("kappa_rec0_fourier.fits", kr0);
    fits_write_dblarr("kappa_rec1_fourier.fits", kr1);

    dblarray kr3;
    kr3.alloc(kvar, npix, npix);
    fits_write_dblarr("kappa_rec.fits", kr3);


  // save nshear_power1 in a fits file
  dblarray expts90;
  double * tmp0 = f->get_nshear_power0();
  expts90.alloc(tmp0, npix, npix);
    fits_write_dblarr("nshear_power0.fits", expts90);

  dblarray expts91;
  double * tmp1 = f->get_nshear_power1();
  expts91.alloc(tmp1, npix, npix);
  fits_write_dblarr("nshear_power1.fits", expts91); // this succeeds


#endif

}

double surface_reconstruction::get_spectral_norm_prox(int niter, double tol)
{

    double norm = 0;
    double norm_old = 0;

    // Initialise array with random numbers
    norm = 0;
    // cout << "Gradient noise map inside get_spectral norm" << endl;
    f->gradient_noise_map(kappa_tmp);

    // std::cout << "Kappa_tmp: " << kappa_tmp[0][0] << " , " << kappa_tmp[0][1] << std::endl;

    for (long ind = 0; ind < ncoeff; ind++) {
        norm += kappa_tmp[ind][0] * kappa_tmp[ind][0] + kappa_tmp[ind][1] * kappa_tmp[ind][1];
    }
    norm = sqrt(norm);
    // std::cout << "norm is " << norm << std::endl;

    // std::cout << "norm insider spectral_norm_prox: " << norm << std::endl;

    // Normalise the input
    for (long ind = 0; ind < ncoeff; ind++) {
        kappa_tmp[ind][0] /= norm;
        kappa_tmp[ind][1] /= norm;
    }

    for (int k = 0; k < niter; k++) {
        wav->transform_gen2(kappa_rec, alpha_u);
        wav->trans_adjoint_gen2(alpha_u, kappa_rec);

        // Compute norm
        for (long ind = 0; ind < ncoeff; ind++) {
            norm += kappa_tmp[ind][0] * kappa_tmp[ind][0] + kappa_tmp[ind][1] * kappa_tmp[ind][1];
        }
        norm = sqrt(norm);

        if (fabs(norm - norm_old) / norm <= tol) {
            break;
        }

        for (long ind = 0; ind < ncoeff; ind++) {
            kappa_tmp[ind][0] /= norm;
            kappa_tmp[ind][1] /= norm;
        }

        norm_old = norm;
        if(k == niter -1 ) std::cout << "Warning, reached maximum number of iterations prox" << std::endl;
    }

    return norm * (1.0 + tol);
}



void surface_reconstruction::compute_weights()
{
    for (long ind = 0; ind < npix * npix; ind++) {
        fft_frame[ind][0] = kappa[ind][0] * fftFactor;
        fft_frame[ind][1] = kappa[ind][1] * fftFactor;
    }
    fftwf_execute(plan_backward);
    for (long ind = 0; ind < npix * npix; ind++) {
        fft_frame[ind][0] = max(fft_frame[ind][0], 0.0f);
        fft_frame[ind][1] = 0;
    }

    fftwf_execute(plan_forward);
    for (long ind = 0; ind < npix * npix; ind++) {
        kappa_tmp[ind][0] = fft_frame[ind][0];
        kappa_tmp[ind][1] = fft_frame[ind][1];
    }

    wav->transform_gen2(kappa_tmp, alpha_res);

    for (long ind = 0; ind < npix * npix * nframes; ind++) {
        if (fabs(alpha_res[ind]) < lambda * thresholds[ind]) {
            weights[ind] = 1.0;
        } else {
            weights[ind] = lambda * thresholds[ind] / fabs(alpha_res[ind]);
        }
    }
}

void surface_reconstruction::get_convergence_map(double *kap)
{

    for (long ind = 0; ind < npix * npix; ind++) {
        fft_frame[ind][0] = kappa[ind][0] * fftFactor;
        fft_frame[ind][1] = kappa[ind][1] * fftFactor;
    }

    fftwf_execute(plan_backward);

    for (int y = 0; y < npix ; y++) {
        for (int x = 0; x < npix ; x++) {
            long pos = y * npix + x ; // (npix - y - 1) * npix + (npix - x - 1); // depends on flagship or mice data
            kap[y * npix + x] = fft_frame[pos][0];
        }
    }
}
