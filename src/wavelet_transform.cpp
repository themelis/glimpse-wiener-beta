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

#include <sparse2d/MR_Obj.h>
#include "wavelet_transform.h"
#include "starlet_2d.h"

wavelet_transform::wavelet_transform(int npix, int nscale, int nlp):
    npix(npix), nscale(nscale), nlp(nlp)
{

    // fftwf_complex *frame = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * npix * npix);
    // fftwf_plan plan  = fftwf_plan_dft_2d(npix, npix, frame, frame, FFTW_FORWARD,   FFTW_MEASURE);

    // We begin with starlets combined with battle_lemarie wavelets
    nframes = nscale + 3; // all starlet frames + the 3 first BL scales

    std::cout << "nframes is : " << nframes << endl;

    fframe = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * npix * npix);
    fplan  = fftwf_plan_dft_2d(npix, npix, fframe, fframe, FFTW_FORWARD,   FFTW_MEASURE);

    // frames      = new float*[nframes];
    frames_fwd  = new float*[nframes];
    frames_bkwd = new float*[nframes];

    for (int i = 0; i < nframes; i++) {
        // frames[i]      = (float *) fftwf_malloc(sizeof(float) * npix * npix);
        frames_fwd[i]  = (float *) fftwf_malloc(sizeof(float) * npix * npix);
        frames_bkwd[i] = (float *) fftwf_malloc(sizeof(float) * npix * npix);
    }

    std::cout << "create a starlet with nscales: " << nframes << endl;
    starlet_2d star(npix, npix, nframes);
    std::cout << "starlet created!" << endl;

    // compute the filter coefficients in fourier space
    dblarray image(npix, npix);
    image.init(0);
    image(npix / 2, npix / 2) = 1.0;
    dblarray alphaStar(npix, npix, nframes);

    star.transform_gen2(image.buffer(), alphaStar.buffer());
    cout << "first transform successful" << endl;

    // // save the wavelet coefficients of the delta image
    // float *alphastarout = (float *) malloc(Nx * Ny * nscales * sizeof(float));
    // for (int i = 0; i < nscales; i++) {
    //   for (long ind = 0; ind < Nx * Ny; ind++) {
    //     alphastarout[i * Nx * Ny + ind] = alphaStar.buffer()[i * Nx * Ny + ind];
    //   }
    // }
    // fltarray alphastarbuffer;
    // alphastarbuffer.alloc(alphastarout, Nx, Ny, nscales);
    // char * filename1 = (char *)"../data/einstein_alphastarbuffer.fits";
    // fits_write_fltarr(filename1, alphastarbuffer);

    // compute the fixed wavelet filter coefficients in fourier space
    for (int i = 0; i < nframes; i++) {
        for (long ind = 0; ind < npix * npix; ind++) {
            fframe[ind][0] = alphaStar.buffer()[i * npix * npix + ind];
            fframe[ind][1] = 0;
        }

        fftwf_execute(fplan);

        for (long ind = 0; ind < npix * npix; ind++) {
            frames_fwd[i][ind] = sqrt(fframe[ind][0] * fframe[ind][0] + fframe[ind][1] * fframe[ind][1])/ (npix * npix);
        }
    }

    std::cout << "frames_fwd computed !" << std::endl;

    // here we add a delta in every scale to get the coefficients of the reconstruction
    dblarray alphaStar2(npix, npix, nframes);
    dblarray alphaStar3(npix, npix, nframes);
    for (int i = 0; i < nframes; i++) {
      for (long ind = 0; ind < npix * npix; ind++) {
        alphaStar2.buffer()[i * npix * npix + ind] = image.buffer()[ind];
      }
    }

    dblarray image_out(npix, npix);
    image_out.init(0);
    star.trans_adjoint_gen2(alphaStar.buffer(), image_out.buffer(), alphaStar3.buffer());

    // char * filename5 = (char *)"../data/delta_reconstructed.fits";
    // fits_write_dblarr(filename5, image_out);

    // reset the fft buffer fframe
    for (long ind = 0; ind < npix * npix; ind++) {
        fframe[ind][0] = 0.;
        fframe[ind][1] = 0.;
    }

    // compute the fixed wavelet filter coefficients in fourier space
    for (int i = 0; i < nframes; i++) {
        for (long ind = 0; ind < npix * npix; ind++) {
            fframe[ind][0] = alphaStar3.buffer()[i * npix * npix + ind]; // image.buffer()[ind];
            fframe[ind][1] = 0.;
        }

        fftwf_execute(fplan);

        for (long ind = 0; ind < npix * npix; ind++) {
            frames_bkwd[i][ind] = sqrt(fframe[ind][0] * fframe[ind][0] + fframe[ind][1] * fframe[ind][1]) ;
        }
    }

    std::cout << "frames_bkwd computed !" << std::endl;
    // fftwf_free(frame);
    // fftwf_destroy_plan(plan);

    // fftwf_free(fframe);
    // std::cout << "fframe freed !" << std::endl;
    // fftwf_destroy_plan(fplan);
    // std::cout << "fplan destroyed !" << std::endl;

    std::cout << "arrived here, frames_bkwd is " << frames_bkwd[2][10] << "frames_fwd is " << frames_fwd[2][10] << endl;

    // Allocate batch wavelet transform either using fftw or CUDA
    fft_frame = fftwf_alloc_complex(npix * npix * nframes);

    int dimensions[2] = {npix, npix};
    int rank = 2;

#ifdef CUDA_ACC
    cufftResult ret = cufftCreate(&fft_plan);

    // Look for the number of available GPUs
    getDeviceCount(&nGPU);
    getGPUs(whichGPUs);

    std::cout << "Performing wavelet transform using " << nGPU << " GPUs" << std::endl;
    // 2 cases: Single GPU or Multiple GPUs
    if(nGPU > 1){

        ret =cufftXtSetGPUs(fft_plan , nGPU, whichGPUs);
        if(ret != 0) std::cout <<"set gpus" << ret << std::endl;

        ret =cufftMakePlanMany(fft_plan, rank, dimensions,
                        NULL, 1, npix*npix, NULL, 1, npix*npix,
                        CUFFT_C2C, nlp * nframes, worksize);
        if(ret != 0) std::cout <<"make plan " << ret << std::endl;

        ret = cufftXtMalloc(fft_plan, &d_frameXt, CUFFT_XT_FORMAT_INPLACE);
        if(ret != 0) std::cout <<"malloc " << ret << std::endl;
     }else{

        // Select GPU
        cudaSetDevice(whichGPUs[0]);

        ret = cufftMakePlanMany(fft_plan, rank, dimensions,
                        NULL, 1, npix*npix, NULL, 1, npix*npix,
                        CUFFT_C2C, nlp * nframes, worksize);
        if(ret != 0) std::cout <<"make plan " << ret << std::endl;

        cudaMalloc(&d_frame, sizeof(cufftComplex)*nlp*nframes*npix*npix);
     }
#else
    plan_forward = fftwf_plan_many_dft(rank, dimensions, nlp * nframes,
                                      fft_frame, NULL, 1, npix * npix,
                                      fft_frame, NULL, 1, npix * npix,
                                      FFTW_FORWARD, FFTW_MEASURE);
    plan_backward = fftwf_plan_many_dft(rank, dimensions, nlp * nframes,
                                       fft_frame, NULL, 1, npix * npix,
                                       fft_frame, NULL, 1, npix * npix,
                                       FFTW_BACKWARD, FFTW_MEASURE);
#endif


}

wavelet_transform::~wavelet_transform()
{

#ifdef CUDA_ACC
    if(nGPU>1){
    cufftXtFree(d_frameXt);
    }else{
     cudaFree(d_frame);
    }
    cufftDestroy(fft_plan);

#else
    fftwf_free(fft_frame);

    std::cout << "fft_frame freed !" << std::endl;
    fftwf_destroy_plan(plan_backward);
    fftwf_destroy_plan(plan_forward);

#endif

    for (int i = 0; i < nframes; i++) {
        fftwf_free(frames_fwd[i]);
        fftwf_free(frames_bkwd[i]);
    }

    delete[] frames_fwd;
    delete[] frames_bkwd;
}

// void wavelet_transform::transform(fftwf_complex *image, float *alpha)
// {
//     #pragma omp parallel
//     for (int z = 0; z < nlp; z++) {
//         for (int i = 0; i < nframes; i++) {
//             #pragma omp for
//             for (long ind = 0; ind < npix * npix; ind++) {
//                 fft_frame[ind + i * npix * npix + z * npix * npix * nframes][0] = image[ind + z * npix * npix][0] * frames[i][ind];
//                 fft_frame[ind + i * npix * npix + z * npix * npix * nframes][1] = image[ind + z * npix * npix][1] * frames[i][ind];
//             }
//         }
//     }
//
// #ifdef CUDA_ACC
//     if(nGPU>1){
//         cufftXtMemcpy(fft_plan, d_frameXt, fft_frame, CUFFT_COPY_HOST_TO_DEVICE);
//         cufftXtExecDescriptorC2C(fft_plan, d_frameXt, d_frameXt, CUFFT_INVERSE);
//         cufftXtMemcpy(fft_plan, fft_frame, d_frameXt, CUFFT_COPY_DEVICE_TO_HOST);
//     }else{
//         cudaMemcpy(d_frame, fft_frame, sizeof(cufftComplex)* npix*npix*nlp*nframes, cudaMemcpyHostToDevice);
//         cufftExecC2C(fft_plan,d_frame,d_frame, CUFFT_INVERSE);
//         cudaMemcpy(fft_frame, d_frame, sizeof(cufftComplex)* npix*npix*nlp*nframes, cudaMemcpyDeviceToHost);
//     }
// #else
//     fftwf_execute(plan_backward);
// #endif
//
//     #pragma omp parallel
//     for (int z = 0; z < nlp; z++) {
//         for (int i = 0; i < nframes; i++) {
//             #pragma omp for
//             for (long ind = 0; ind < npix * npix; ind++) {
//                 alpha[ind + i * npix * npix + z * npix * npix * nframes] = fft_frame[ind + i * npix * npix + z * npix * npix * nframes][0];
//             }
//         }
//     }
// }
//
// void wavelet_transform::trans_adjoint(float *alpha, fftwf_complex *image)
// {
//
//     #pragma omp parallel for
//     for (long ind = 0; ind < npix * npix * nlp; ind++) {
//         image[ind][0] = 0;
//         image[ind][1] = 0;
//     }
//
//     #pragma omp parallel
//     for (int z = 0; z < nlp; z++) {
//         for (int i = 0; i < nframes; i++) {
//
//             #pragma omp for
//             for (long ind = 0; ind < npix * npix; ind++) {
//                 fft_frame[ind + i * npix * npix + z * npix * npix * nframes][0] = alpha[ind + i * npix * npix + z * npix * npix * nframes];
//                 fft_frame[ind + i * npix * npix + z * npix * npix * nframes][1] = 0;
//             }
//         }
//     }
//
// #ifdef CUDA_ACC
//     if(nGPU>1){
//         cufftXtMemcpy(fft_plan, d_frameXt, fft_frame, CUFFT_COPY_HOST_TO_DEVICE);
//         cufftXtExecDescriptorC2C(fft_plan, d_frameXt, d_frameXt, CUFFT_FORWARD);
//         cufftXtMemcpy(fft_plan, fft_frame, d_frameXt, CUFFT_COPY_DEVICE_TO_HOST);
//     }else{
//         cudaMemcpy(d_frame, fft_frame, sizeof(cufftComplex)* npix*npix*nlp*nframes, cudaMemcpyHostToDevice);
//         cufftExecC2C(fft_plan,d_frame,d_frame, CUFFT_FORWARD);
//         cudaMemcpy(fft_frame, d_frame, sizeof(cufftComplex)* npix*npix*nlp*nframes, cudaMemcpyDeviceToHost);
//     }
// #else
//     fftwf_execute(plan_forward);
// #endif
//
//     #pragma omp parallel
//     for (int z = 0; z < nlp; z++) {
//         for (int i = 0; i < nframes; i++) {
//
//             #pragma omp for
//             for (long ind = 0; ind < npix * npix; ind++) {
//                 image[ind + z * npix * npix][0] += fft_frame[ind + i * npix * npix + z * npix * npix * nframes][0] * frames[i][ind];
//                 image[ind + z * npix * npix][1] += fft_frame[ind + i * npix * npix + z * npix * npix * nframes][1] * frames[i][ind];
//             }
//         }
//     }
// }



void wavelet_transform::transform_gen2(fftwf_complex *image, float *alpha)
{

  // std::cout << "Inside wavelet transform" << '\n';

  // // export image in fourier space
  // float *imf = (float *) malloc(Nx * Ny * 2 * sizeof(float));
  // for (int i = 0; i < 2; i++) {
  //   for (long ind = 0; ind < Nx * Ny; ind++) {
  //     imf[ind + i * Nx * Ny] = image[ind][0];
  //     imf[ind + i * Nx * Ny] = image[ind][1];
  //   }
  // }
  // fltarray imfs;
	// imfs.alloc(imf, Nx, Ny, 2);
	// char * filename1 = (char *)"../data/einstein_fourier.fits";
	// fits_write_fltarr(filename1, imfs);

  // // export delta coefficients in fourier space
  // float *coeff = (float *) malloc(Nx * Ny * nscales * sizeof(float));
  // for (int i = 0; i < nscales; i++) {
  //   for (long ind = 0; ind < Nx * Ny; ind++) {
  //     coeff[i * Nx * Ny + ind] = frames_fwd[i][ind];
  //   }
  // }
  // fltarray coeffs;
	// coeffs.alloc(coeff, Nx, Ny, nscales);
	// char * filename2 = (char *)"../data/frames_fwd.fits";
	// fits_write_fltarr(filename2, coeffs);

  for (int i = 0; i < nframes; i++) {
      for (long ind = 0; ind < npix * npix; ind++) {
          fft_frame[ind + i * npix * npix ][0] = image[ind][0] * frames_fwd[i][ind] ;
          fft_frame[ind + i * npix * npix ][1] = image[ind][1] * frames_fwd[i][ind] ;
      }
  }

  fftwf_execute(plan_backward);

  for (int i = 0; i < nframes; i++) {
      for (long ind = 0; ind < npix * npix; ind++) {
          alpha[ind + i * npix * npix] = fft_frame[ind + i * npix * npix][0] ;
      }
  }

}

void wavelet_transform::trans_adjoint_gen2(float *alpha, fftwf_complex *image)
{

    for (long ind = 0; ind < npix * npix; ind++) {
        image[ind][0] = 0;
        image[ind][1] = 0;
    }

    for (int i = 0; i < nframes; i++) {
        for (long ind = 0; ind < npix * npix; ind++) {
            fft_frame[ind + i * npix * npix][0] = alpha[ind + i * npix * npix];
            fft_frame[ind + i * npix * npix][1] = 0.;
        }
    }

    fftwf_execute(plan_forward);

    for (int i = 0; i < nframes; i++) {
        for (long ind = 0; ind < npix * npix; ind++) {
            image[ind][0] += fft_frame[ind + i * npix * npix][0] * frames_bkwd[i][ind];
            image[ind][1] += fft_frame[ind + i * npix * npix][1] * frames_bkwd[i][ind];
        }
    }
}
