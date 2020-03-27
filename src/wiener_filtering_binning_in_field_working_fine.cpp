
#include "wiener_filtering.h"

wiener_filtering::wiener_filtering(surface_reconstruction * reca, field *fi, survey *su, dblarray ps_array)
{
	surv = su;
	f = fi;
	rec = reca;

	ngal = su->get_ngal();
	npix = f->get_npix();
	pixel_size = f->get_pixel_size();
	impix = npix * npix;
	fftFactor = 1.0 / ((double) impix);
	ps1dlen = ps_array.nx();
	// store shear values in a vector, for computing the noise covariance
	/*shear.reserve( 2 * ngal);
	for (unsigned i=0; i < ngal; i++) {
		shear.push_back(surv->get_gamma1(i));
		shear.push_back(surv->get_gamma2(i));
	}*/
	shear0.reserve(ngal);
	shear1.reserve(ngal);
	for (int i=0; i < ngal; i++) {
		shear0.push_back(surv->get_gamma1(i));
		shear1.push_back(surv->get_gamma2(i));
	}

	// read true kappa spectrum Px from file - replaced by reading input file
	//dblarray ps_arr1;
	//fits_read_dblarr("../micev1_dec2/Ps1d_kappa_smooth.fits", ps_arr1);
	//fits_read_dblarr("../data/flagship/ps1d_pure_gaussian_simulated.fits", ps_arr1);
	//fits_read_dblarr("../xmpl/Ps1d_dec8_gaussian_component_pxtrue.fits", ps_arr1);
	//fits_read_dblarr("../xmpl/Ps1d_kappa_true.fits", ps_arr1);
	//printf("Size of dblarray is %d\n", ps_arr1.nx() );

	// transpose 1d spectrum to 2d image
	//printf("Read input spectrum:\n");
	std::vector<double> ps1d(ps1dlen,.0);
	for (int i = 0; i < ps1dlen; i++){
		ps1d[i] = ps_array(i);
		// printf("%f ", ps_array(i));
	}
	// printf("\n");

	std::vector<double> ps2d = powerspec2d(ps1d);

	// fftshift the resulting 2d power spectrum
	std::complex<double> * rimg = (std::complex<double> *) malloc(sizeof(std::complex<double>) * impix);
	std::complex<double> * shimg = (std::complex<double> *) malloc(sizeof(std::complex<double>) * impix);
	for (int i = 0; i < impix ; i++)
		rimg[i] = {ps2d[i], .0};
	fftshift(rimg, shimg, npix);

	/* read the Px from eusipco
	dblarray ps_arr1;
	fits_read_dblarr("../xmpl/Px_eusipco.fits", ps_arr1);
	Px.reserve(impix);
	for (int i=0; i < npix; i++)
		for (int j=0; j < npix; j++)
			Px.push_back(ps_arr1(i*npix+j));
	*/

	// initialize galaxy counter and noise, signal covariance matrices
	Sn.reserve(impix);
	Px.reserve(impix);
	// counter.reserve(impix);
	for (int i=0; i < impix; i++) {
		// counter.push_back(0);
		Sn.push_back(0);
		Px.push_back(shimg[i].real()); //Px.push_back(ps_arr1(i));
	}

	// get the indexes of the binning process
	index1d = f->get_binindex();
	/*
	std::cout << "Inside constructor, printing index1d" << std::endl;
	for (unsigned i=0; i < 200; i++)
		std::cout << index1d[i] << " ";
	std::cout << std::endl;
	*/

	// get the number of galaxies per pixel
	counter = f->get_counter();
	/*
	std::cout << "Inside constructor, printing counter" << std::endl;
	for (unsigned i=0; i < impix; i++)
		std::cout << counter[i] << " ";
	std::cout << std::endl;
	*/

	// allocate the shear residuals
	shear_res.reserve(impix);

	// allocate fftw plans
	frame = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * impix);
  planf = fftw_plan_dft_2d(npix, npix, frame, frame, FFTW_FORWARD, FFTW_MEASURE);
  planb = fftw_plan_dft_2d(npix, npix, frame, frame, FFTW_BACKWARD, FFTW_MEASURE);

	fft_frameHf0 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * impix);
	fft_frameHf1 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * impix);
	fft_frameHb0 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * impix);
	fft_frameHb1 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * impix);
	planHf0 = fftw_plan_dft_2d(npix, npix, fft_frameHf0, fft_frameHf0, FFTW_FORWARD, FFTW_MEASURE);
	planHf1 = fftw_plan_dft_2d(npix, npix, fft_frameHf1, fft_frameHf1, FFTW_FORWARD, FFTW_MEASURE);
	planHb0 = fftw_plan_dft_2d(npix, npix, fft_frameHb0, fft_frameHb0, FFTW_BACKWARD, FFTW_MEASURE);
	planHb1 = fftw_plan_dft_2d(npix, npix, fft_frameHb1, fft_frameHb1, FFTW_BACKWARD, FFTW_MEASURE);

	/* vectors for the shear residual*/
	/*
	res1.reserve(ngal);
	res2.reserve(ngal);
	pos1.reserve(ngal);
	pos2.reserve(ngal);
	wght.reserve(ngal);
	*/
}

wiener_filtering::~wiener_filtering()
{
	free(planf);
	free(planb);
	free(planHf0);

	free(planHb0);
	fftw_free(frame);
	fftw_free(fft_frameHf0);
	fftw_free(fft_frameHb0);
}

/*
 *  computes the diagonal noise covariance matrix of the shear
 *  can be used whether data are masked or not
 */
void wiener_filtering::compute_noise_cov(std::vector<int> bin_ind)
{

	//std::cout << "Inside compute noise cov" << endl;

	int Niter = 30; // # of random realizations

	std::vector<std::complex<double>> noise_res(impix, .0);
	for (int i=0; i < Niter; i++) {
		// permute shear values
		std::random_shuffle(shear0.begin(), shear0.end());
		std::random_shuffle(shear1.begin(), shear1.end());

		// bin the permuted shear vectors, utilizing the same binning index
		noise_res = bincount(bin_ind, shear0, shear1);
		// std::cout << "Noise res computed" << endl;
		// compute the diagonal noise covariance matrix Sn
		for (int j = 0; j < impix; j++)
			Sn[j] += pow(noise_res[j].real(), 2) + pow(noise_res[j].imag(), 2);
	}

	// std::cout << "Sn computed" << endl;
	for (int j = 0; j < impix; j++)
		if (Sn[j]==0)
			Sn[j] = 1e19;
		else
		 	Sn[j] = Sn[j] / (Niter - 1);

}

/*
 *  computes the standard deviation of the noise
 */
/*
void wiener_filtering::compute_noise_cov_old(std::vector<double> data)
{

  // calculate the mean
	double sum = std::accumulate(data.begin(), data.end(), 0.0);
	double mean = sum / (data.size()-1);

	// printf("The mean is %f\n", mean);

  // calculate the standard deviation
	std::vector<double> diff(data.size(), 0.);
	std::transform(data.begin(), data.end(), diff.begin(), [mean](double x) { return x - mean; });

	double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
	double stdev = std::sqrt(sq_sum / (data.size()-1));
	//stdev = 0.286; // stdev coming from data generation
  printf("The noise std is %f\n",stdev);
	double var = pow(stdev, 2);

	// and divide by the number of galaxies per pixel
	//printf("Noise covariance:\n");
	for (int y = 0; y < impix ; y++) {
		Sn[y] = var / counter[y];
		// (counter[y] >= 2) ? Sn[y] = stdev / (counter[y]-1) : Sn[y] = stdev;
		// Sn[y] = var / counter[y];
		// printf("%0.9f ",Sn[y]);
  }
  // printf("\n");
}
*/

/*
 * main iteration function of wiener filtering
 */
void wiener_filtering::main_iteration()
{

	// extract the shear residuals
	//get_residuals();
	f->comp_residuals(rec->get_kappa(), shear_res);
	// std::cout << "Got residuals" << std::endl;
	/*
	printf("shear res is: \n");
	for (int i = 0; i < 200; i++)
		printf("%0.9f\t", shear_res[i].imag());
	printf("\n");
	*/

	// write shear residuals to a fits file
	double * img1 = (double *) malloc(sizeof(double) * impix);
	double * img2 = (double *) malloc(sizeof(double) * impix);
	for (int i = 0; i < impix; i++) {
		img1[i] = shear_res[i].real();
		img2[i] = shear_res[i].imag();
	}
	dblarray expts1;
	dblarray expts2;
	expts1.alloc(img1, npix, npix);
	expts2.alloc(img2, npix, npix);
	fits_write_dblarr("shear_res1_map.fits", expts1);
	fits_write_dblarr("shear_res2_map.fits", expts2);


	// bin the shear residuals - no binning here, in the field instead
	// std::vector<int> ind_bin = bin2d();

	// compute the noise covariance matrix Sx of the shear
	compute_noise_cov(index1d);
	// std::cout << "Noise covariance computed" << std::endl;
	/*
	// read covariance matrix from fits file
	dblarray sn_arr1;
	fits_read_dblarr("../xmpl/Sn_diagonal_nomask_dec2.fits", sn_arr1);
	int elt1 = 0;
	for (int i = 0; i < npix; i++){
		for (int j = 0; j < npix; j++){
			Sn[j*npix+i] = sn_arr1(elt1);
			elt1++;
		}
	}
	*/

	/*
	// store covariance matrix in a fits file
	double *sncov = (double *) malloc(sizeof(double) * impix);
	for (int i = 0; i < impix; i++)
		sncov[i] = Sn[i];
	dblarray scovm;
	scovm.alloc(sncov, npix, npix);
	fits_write_dblarr("Sn.fits", scovm); // store the covariance matrix to a fits file
	*/

	// get glimpse map
	double *kglimpse = (double *) malloc(sizeof(double) * impix);
	rec->get_convergence_map(kglimpse);
	dblarray expt1;
  expt1.alloc(kglimpse, npix, npix);
  fits_write_dblarr("kglimpse.fits", expt1); // this is the same as in glimpse output

	// run wiener filtering core for the binned shear residuals
	std::vector<double> kwiener = run_wiener_core(shear_res, 300);

	// write fits file with wiener solution
	dblarray expt;
  expt.alloc(&kwiener[0], npix, npix);
  fits_write_dblarr("kwiener.fits", expt);

	// add glimpse and wiener solutions
	double *glpluswf = (double *) malloc(sizeof(double) * impix);
	for (int i = 0; i < impix; i++)
		glpluswf[i] = kwiener[i] + kglimpse[i];

	// write result to a fits file
	dblarray expt2;
  expt2.alloc(glpluswf, npix, npix);
  fits_write_dblarr("glpluswf.fits", expt2);

}

std::vector<double> wiener_filtering::powerspec2d(std::vector<double> spectrum){

	// this function will work for square images only
	std::vector<double> ps2d(impix, 0.0);

	long unsigned int nspec = spectrum.size();

	int *ni = (int *) malloc(sizeof(int) * npix);
	for (int i = 0; i < npix; i++)
		ni[i] = int(i - npix / 2.);

	int dr;
	int ind = 0;
	for (int i=0; i < npix; i++){
		for (int j=0; j < npix; j++){
			dr = int(sqrt(ni[i] * ni[i] + ni[j] * ni[j]));
			// if ((i<20) && (j<20))
			//  	printf("dr = %d and specrum[dr]=%f\n", dr, spectrum[dr]);
			if ( (dr == 0) || (dr>nspec-1)) //(dr == 0) ||
				//printf("To ix einai %d and to iy einai %d\n", i, j);
				ps2d[ind] = 1e-19;
			else
				ps2d[ind] = spectrum[dr];
			++ind;
		}
	}

	return ps2d;
}

std::vector<double> wiener_filtering::run_wiener_core(std::vector<std::complex <double>> data, int Niter)
{
	std::vector<double> result(impix, .0);
	std::vector<double> rres(impix, .0);
	std::vector<double> ires(impix, .0);

	// initiallize
	std::vector<std::complex <double>> t(impix, (std::complex<double>){{.0,.0}});
	std::vector<std::complex <double>> xg(impix, (std::complex<double>){{.0,.0}});

	// find the minimum noise variance
	std::vector<double>::iterator itr = std::min_element(Sn.begin(), Sn.end());
	int index  = std::distance(Sn.begin(), itr);
	double tau = Sn[index];

	//printf("Minimum noise variance is %0.9f\n", tau);
	double mu = .83 * tau;
	double npow = pow(mu, 2);//  * impix;
	// double tpow = pow(tau, 2) * impix;
	// printf("tau \t%0.9f, ptau\t%0.9f\n", tau, tpow);
	// printf("mu  \t%0.9f, pmu \t%0.9f\n", mu, npow);

	std::vector<double> T(impix, mu);
	std::vector<double> Pn(impix, npow); // changed from tau * impix

	// printf("Noise power estimate: %0.9f\n", npow);

	std::vector<double> Wfc(impix, .0);
	std::vector<double> TSn(impix, .0); //(std::complex<double>){.0, .0}
	std::transform(T.begin(), T.end(), Sn.begin(), TSn.begin(), std::divides<double>()); // T / Sn

	// calculate the wiener filter coefficients
	std::transform(Px.begin(), Px.end(), Pn.begin(), Wfc.begin(), std::plus<double>()); // Wfc = Px + Pn
	std::transform(Px.begin(), Px.end(), Wfc.begin(), Wfc.begin(), std::divides<double>()); // Wfc = Px / (Px + Pn)

	// Define array just for exporting data
	double * toto = (double *) malloc(sizeof(double) * impix);
	double * toti = (double *) malloc(sizeof(double) * impix);
	dblarray yop,yopi;
	yop.alloc(toto,npix,npix);
	yopi.alloc(toti,npix,npix);
	char name[256];

	double * ppo = (double *) malloc(sizeof(double) * ps1dlen);
	double * ppi = (double *) malloc(sizeof(double) * ps1dlen);
	double * ppe = (double *) malloc(sizeof(double) * ps1dlen);
	dblarray djp, jfs, gff;
	djp.alloc(ppo, 1, ps1dlen);
	jfs.alloc(ppi, 1, ps1dlen);
	gff.alloc(ppe, 1, ps1dlen);
	std::vector<double> ps1d;

	// store TSn in a fits file to plot in notebook
	double * tmp43 = (double *) malloc(sizeof(double) * impix);
	dblarray tsnexp;
	tsnexp.alloc(tmp43, 1, impix);
	for (int i=0; i < impix; i++)
		tmp43[i] = TSn[i];
	sprintf(name, "TSn.fits");
	fits_write_dblarr(name, tsnexp);

	for (int i =0; i < Niter; ++i)
	{
		//printf("-- Inside rwc, iteration number %d\n", i);
		// printf("Spectrum of xg\n");
		ps1d = csp(xg, npow, i);

		/*
		// save it to a file
		for (int i=0; i < ps1d.size(); i++)
			ppe[i] = (double ) ps1d[i];
		sprintf(name, "psx_%03d.fits", i);
		fits_write_dblarr(name, gff);
		*/
		/*
		if (i == 0){
			// debug power spectrum of t in the first iteration
			// t equation
			t = H_operator(xg); // H * xg
			std::transform(data.begin(), data.end(), t.begin(), t.begin(), std::minus<std::complex<double>>()); // y - H * xg
			std::transform(t.begin(), t.end(), Sn.begin(), t.begin(), std::divides<std::complex<double>>()); // (1 / Sn) * (y - H * xg)

			ps1d = csp(t, npow, i);
			// save it to a file
			for (int i=0; i < ps_residual.size(); i++)
				ppi[i] = (double ) ps1d[i];
			sprintf(name, "pst_%03d_1.fits", i);
			fits_write_dblarr(name, jfs);
		}
		if (i==5){
			// debug power spectrum of xg in 10th iteration
			// see how H^T H affects the spectrum
			t = H_operator(xg); // H * xg
			t = H_adjoint_operator(xg); // H^T * H * xg
			ps1d = csp(t, npow, i);
			// save it to a file
			for (int i=0; i < ps1d.size(); i++)
				ppi[i] = (double ) ps1d[i];
			sprintf(name, "pst_%03d_2.fits", i);
			fits_write_dblarr(name, jfs);
		}
		if (i==10){
			// debug power spectrum of xg in 10th iteration
			// see how H^T H affects the spectrum
			t = H_operator(xg); // H * xg
			t = H_adjoint_operator(xg); // H^T * H * xg
			ps1d = csp(t, npow, i);
			// save it to a file
			for (int i=0; i < ps1d.size(); i++)
				ppi[i] = (double ) ps1d[i];
			sprintf(name, "pst_%03d_2.fits", i);
			fits_write_dblarr(name, jfs);
		}
		*/

		// t equation
		t = H_operator(xg); // H * xg
		std::transform(data.begin(), data.end(), t.begin(), t.begin(), std::minus<std::complex<double>>()); // y - H * xg
		std::transform(TSn.begin(), TSn.end(), t.begin(), t.begin(), std::multiplies<std::complex<double>>()); // (T / Sn) * (y - H * xg)
		t = H_adjoint_operator(t); // H^T * ( (T/Sn) * (y-H*xg) )

		//printf("Spectrum of residual\n");
		ps1d = csp(t, npow, i);
		/*
		// save it to a file
		for (int i=0; i < ps1d.size(); i++)
			ppo[i] = (double) ps1d[i];
		sprintf(name, "psr_%03d.fits", i);
		fits_write_dblarr(name, djp);
		*/

		std::transform(t.begin(), t.end(), xg.begin(), t.begin(), std::plus<std::complex<double>>()); // t = xg + H^T * ( (T / Sn) * (y - H * xg) )

		for (int ind = 0; ind < impix ; ind++)
    	toti[ind] = t[ind].real(); // sqrt(frame[ind][0]*frame[ind][0] + frame[ind][1]*frame[ind][1])* fftFactor;

		// printf("Spectrum of t\n");
		ps1d = csp(t, npow, i);
		/*
		// save it to a file
		for (int i=0; i < ps1d.size(); i++)
			ppi[i] = (double ) ps1d[i];
		sprintf(name, "pst_%03d.fits", i);
		fits_write_dblarr(name, jfs);
		*/
		// compute the fft of t
		for (int ind =0; ind < impix; ++ind)
		{
			frame[ind][0] = t[ind].real();
			frame[ind][1] = t[ind].imag();
		}
		fftw_execute(planf);
		for (int ind = 0; ind < impix ; ind++)
			t[ind] = {frame[ind][0], frame[ind][1]};

		// xg equation
		std::transform(Wfc.begin(), Wfc.end(), t.begin(), xg.begin(), std::multiplies<std::complex<double>>());

		// compute the inverse fft of xg
		for (int ind =0; ind < impix; ++ind)
		{
			frame[ind][0] = xg[ind].real() * fftFactor;
			frame[ind][1] = xg[ind].imag() * fftFactor;
		}
		fftw_execute(planb);
		for (int ind = 0; ind < impix ; ind++) {
			xg[ind] = {frame[ind][0], frame[ind][1]};
    	toto[ind] = xg[ind].real();
    }

		/*
    sprintf(name, "wkappa_%03d.fits", i);
    fits_write_dblarr(name, yop);

    sprintf(name, "wt_%03d.fits", i);
    fits_write_dblarr(name, yopi);
		*/
	}

	for (int ind = 0; ind < impix ; ind++){
		result[ind] = (double) (real(xg[ind]));
		rres[ind] = (double) (real(xg[ind]));
		ires[ind] = (double) (imag(xg[ind]));
	}

	// export result to fits files
	dblarray expmg1, expmg2;
  expmg1.alloc(&rres[0], npix, npix);
  fits_write_dblarr("rres.fits", expmg1);
	expmg2.alloc(&ires[0], npix, npix);
  fits_write_dblarr("ires.fits", expmg2);

	return result;
}

// function to compute the power spectrum of the auxiliary variable t
std::vector<double>  wiener_filtering::csp(std::vector<std::complex <double>> data, double rpow, int it){
	// steps:
	// a\ estimate the 1d power density of data

	for (int ind =0; ind < impix; ++ind)
	{
		frame[ind][0] = data[ind].real();
		frame[ind][1] = data[ind].imag();
		// (it < 100) ? frame[ind][1] = data[ind].imag() : frame[ind][1] = .0; // data[ind].imag() / 3.0;
	}
	fftw_execute(planf);

	std::complex<double> * rimg = (std::complex<double> *) malloc(impix * sizeof(std::complex<double>));
	std::complex<double> * shimg = (std::complex<double> *) malloc(impix * sizeof(std::complex<double>));

	for (int ind =0; ind < impix; ++ind)
		rimg[ind] = {(pow(frame[ind][0], 2.0) + pow(frame[ind][1], 2.0)), 0.0};

	fftshift(rimg, shimg, npix);

	// compute 1d radial profile

	// create indices
	std::vector<int> r(impix, 0.0);
	int center = (int) npix / 2.0;

	//printf("The center is: %d\n", center);
	for (int y = 0; y < npix ; y++)
		for (int x = 0; x < npix ; x++)
			r[y * npix + x] = (int) sqrt(pow(x - center, 2.0) + pow(y - center, 2.0));

	// find the maximum of r
	std::vector<int>::iterator itr;
	int index;
	itr = std::max_element(r.begin(), r.end());
	index  = std::distance(r.begin(), itr);
	int rmax = r[index];
	// printf("the maximum of r is %d\n",rmax);

	// count the occurences
	std::vector<int> cnt;

	int mycount;
	for (unsigned i=0; i < rmax+1; i++) {
		mycount = std::count(r.begin(), r.end(), i);
		cnt.push_back(mycount);
	}

	for (unsigned i=0; i < rmax+1; i++) {
		if(cnt[i] == 0)
			cnt[i] = 1;
	}

	std::vector<double> spectrum(rmax + 1, 0.0);
	for (unsigned i=0; i < impix; i++) {
		spectrum[r[i]] += shimg[i].real();
	}

	for (unsigned i=0; i < spectrum.size(); i++)
		spectrum[i] = spectrum[i] / cnt[i];

	// spectrum[0] = spectrum[1] / 2.0;
	// printf("inside csp, the 1d spectrum is: \n");
	// for (int i = 0; i < spectrum.size(); i++)
	//	printf("%0.9f, ", spectrum[i]);
	// printf("\n");
	// store spectrum of xg to a vector, use Pxg = Pt -1

	/*
	// compute the mean of Pr in previous iteration
	double sum = std::accumulate(ps_residual.begin(), ps_residual.end(), 0.0);
	double mean = sum / (ps_residual.size()-1);

	std::vector<double> ps2d = powerspec2d(ps1d);

	for (int ind = 0; ind < impix ; ind++)
		rimg[ind] = {ps2d[ind], .0};
	fftshift(rimg, shimg, npix);
	for (int ind = 0; ind < impix ; ind++)
		Px[ind] = shimg[ind].real();
	*/
	return spectrum;
}

/*
void wiener_filtering::get_residuals()
{

	std::cout << "Inside get residuals" << std::endl;
	// compute the residuals
	f->comp_residuals(rec->get_kappa(), shear_res);
	std::cout << "Residuals computed" << std::endl;

	// store values in vectors to perform binning
	for (int ind = 0; ind < ngal; ind++){
		pos1.push_back(surv->get_pos1(ind));
		pos2.push_back(surv->get_pos2(ind));
		wght.push_back(surv->get_shear_weight(ind));
	}


}
*/

/*
std::vector<int> wiener_filtering::bin2d()
{
	std::vector<double>::iterator itr;
	int index;

	// find the minimum x position
	itr = std::min_element(pos1.begin(), pos1.end());
	index  = std::distance(pos1.begin(), itr);
	double xmin = pos1[index];

	//printf("itr is: %f\n", itr);
	//printf("xmin element is: %f, found at position: %d\n", xmin, index);

	// find the maximum x position
	itr = std::max_element(pos1.begin(), pos1.end());
	index = std::distance(pos1.begin(), itr);
	double xmax = pos1[index];

	//printf("xmax element at: %d\n", index);

	// find the minimum y position
	itr = std::min_element(pos2.begin(), pos2.end());
	index = std::distance(pos2.begin(), itr);
	double ymin = pos2[index];

	//printf("ymin element at: %d\n", index);

	// find the maximum y position
	itr = std::max_element(pos2.begin(), pos2.end());
	index = std::distance(pos2.begin(), itr);
	double ymax = pos2[index];

	//printf("ymax element at: %d\n", index);
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

	// get the binning indexes
	std::vector<int> indx = digitize(pos1, xedges);
	std::vector<int> indy = digitize(pos2, yedges);

	// create the 1d indexing vector
	std::vector<int> ind_1d;
	ind_1d.reserve(indx.size());
	std::transform(indy.begin(), indy.end(), indy.begin(), std::bind1st(std::multiplies<int>(),npix));
	std::transform(indx.begin(), indx.end(), indy.begin(), std::back_inserter(ind_1d), std::plus<int>());

	// compute the vectorized binned shear residuals
	shear_res = bincount(ind_1d, res1, res2);

	printf("shear res is: \n");
	for (int i = 0; i < 200; i++)
		printf("%0.9f\t", shear_res[i].imag());
	printf("\n");

	// write shear_res to a file
	double * img1 = (double *) malloc(sizeof(double) * impix);
	double * img2 = (double *) malloc(sizeof(double) * impix);
	for (int i = 0; i < impix; i++) {
		img1[i] = shear_res[i].real();
		img2[i] = shear_res[i].imag();
	}
	dblarray expts1;
	dblarray expts2;
  expts1.alloc(img1, npix, npix);
	expts2.alloc(img2, npix, npix);
  fits_write_dblarr("shear_res1_map.fits", expts1);
	fits_write_dblarr("shear_res2_map.fits", expts2);

	return ind_1d;
}
*/

std::vector<std::complex<double>> wiener_filtering::bincount(std::vector<int> ind1d, std::vector<double> residual1, std::vector<double> residual2)
{

	// cout << "Npix is " << n << endl;
	// std::cout << "Inside bincount" << endl;
	//std::vector<double> counter;
	std::vector<std::complex<double>> result;
	result.reserve(impix);
	double *summ0 = (double *) malloc(sizeof(double) * impix);
	double *summ1 = (double *) malloc(sizeof(double) * impix);

	for (int i=0; i < impix; i++) {
		summ0[i] = .0;
		summ1[i] = .0;
		result.push_back({.0,.0});
	}

	// cout << "Starting size is " << counter.size() << endl;
	// printf("Starting Gcount size is %ld\n" , counter.size());
	// int mycount;
	//std::vector<double>::iterator it = counter.begin();
	// int vl = std::distance(ind1d.begin(),ind1d.end());
	// std::cout << "Inside bincount, distance is " << ngal << std::endl;
	/*
	// count the occurences
	for (unsigned i=0; i < impix; i++) {
		mycount = std::count(ind1d.begin(), ind1d.end(), i);
		// gcount.insert(gcount.begin() + i, mycount);
		counter[i] = mycount;
	}
	// printf("After insert Gcount size is %ld\n" , counter.size());
	// remove the zero values
	for (unsigned i=0; i < impix; i++) {
		if(counter[i] == 0)
			counter[i] = 1;
	}
	*/

	/*
	std::cout << "Inside bincount, printing ind1d" << std::endl;
	for (unsigned i=0; i < 200; i++)
		std::cout << ind1d[i] << " ";
	std::cout << std::endl;
	std::cout << "Inside bincount, printing counter" << std::endl;
	for (unsigned i=0; i < 200; i++)
		std::cout << counter[i] << " ";
	std::cout << std::endl;
	std::cout << "Inside bincount, printing shear" << std::endl;
	for (unsigned i=0; i < 200; i++)
		std::cout << residual2[i] << " ";
	std::cout << std::endl;
	std::cout << "Inside bincount, printing shear" << std::endl;
	for (unsigned i=0; i < 200; i++)
		std::cout << residual1[i] << " ";
	std::cout << std::endl;
	*/
	// printf("After accessing Gcount size is %ld\n" , counter.size());

	// summ the residual values
	for (long i=0; i < ngal; i++) {
		summ0[ind1d[i]] += residual1[i];
		summ1[ind1d[i]] += residual2[i];
	}

	for (unsigned i=0; i < impix; i++)
		result[i] = {summ0[i]/counter[i], summ1[i]/counter[i]};

	/*
	printf("Inside bincount, resulting residuals: \n");
	for (int i = 0; i < 1000; i++)
		printf("%0.9f\t", result[i].imag());
	printf("\n");
	*/

	return result;
}

/*
 * Return the indices of the bins to which each value in input array belongs
 */
 /*
std::vector<int> wiener_filtering::digitize(std::vector<double> values, std::vector<double> edges)
{
	std::vector<int> result;

	int vl = std::distance(values.begin(),values.end());
	int pos;

	for (unsigned i=0; i<vl; i++) {
		auto lower = std::lower_bound(edges.begin(), edges.end(), values[i] );
		pos = std::distance(edges.begin(), lower);
		//cout << "Value is " << values[i] << "found at " << pos << endl;
		result.push_back(pos-1);
	}

	return result;
}
*/

/*
std::vector<double> wiener_filtering::linspace(double min, double max, int n)
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
}
*/
/*
 * Function to transform from convergence to shear
 */
std::vector<std::complex<double>> wiener_filtering::H_operator(std::vector<std::complex<double>> kappa)
{

	std::vector<std::complex<double>> result(impix, {.0, .0});

	// compute the fourier transform of kappa
	for (long ind = 0; ind < impix; ind++) {
		fft_frameHb0[ind][0] = kappa[ind].real() * fftFactor;
		fft_frameHb0[ind][1] = 0.0; // kappa[ind].imag() * fftFactor;
  }

	for (long ind = 0; ind < impix; ind++) {
  	fft_frameHb1[ind][0] = kappa[ind].imag() * fftFactor;
		fft_frameHb1[ind][1] = 0.0; // kappa[ind].imag() * fftFactor;
  }
	fftw_execute(planHb0);
	fftw_execute(planHb1);

	// transform from kappa to gamma in fourier space
	double freqFactor = 2.0 * M_PI / ((double) npix); // / pixel_size; // 1.0 / ((double) npix); //
	double k1, k2, k1k1, k2k2, k1k2, ksqr;
	double c1, c2;
	long ind;
	for (int y = 0; y < npix ; y++) {
		// k2 = (y - npix / 2) * freqFactor; // (y - npix / 2)
		k2 = (y < npix / 2 ? y * freqFactor : (y - npix) * freqFactor);

		// int ky  = (y < npix / 2 ? y + npix / 2 : y - npix / 2);
		for (int x = 0; x < npix ; x++) {
			// k1 = (x - npix / 2) * freqFactor; // (x - npix / 2) *
			k1 = (x < npix / 2 ? x * freqFactor : (x - npix) * freqFactor); // this comes from python fft

			// int kx  = (x < npix / 2 ? x + npix / 2 : x - npix / 2);
			// long pos = ky * npix + kx;
			ind = y * npix + x;
			if (k1 == 0 && k2 == 0) {
				fft_frameHf0[ind][0] = .0;
				fft_frameHf0[ind][1] = .0;
				fft_frameHf1[ind][0] = .0;
				fft_frameHf1[ind][1] = .0;
				continue;
			}

			k1k1 = k1 * k1;
			k2k2 = k2 * k2;
			k1k2 = k1 * k2;
			ksqr = k1k1 + k2k2;
			c1 = -k2k2 + k1k1;
			c2 = 2.0 * k1k2;

			fft_frameHf0[ind][0] = (fft_frameHb0[ind][0] * c1 + fft_frameHb1[ind][0] * c2)  / ksqr ;
			fft_frameHf0[ind][1] = (fft_frameHb0[ind][1] * c1 + fft_frameHb1[ind][1] * c2)  / ksqr ;

			fft_frameHf1[ind][0] = (-fft_frameHb0[ind][0] * c2 + fft_frameHb1[ind][0] * c1)  / ksqr ;
			fft_frameHf1[ind][1] = (-fft_frameHb0[ind][1] * c2 + fft_frameHb1[ind][1] * c1)  / ksqr ;
		}
	}

	// compute the inverse fourier transform
	fftw_execute(planHf0);
	fftw_execute(planHf1);

	// return computed shear
  for (long ind = 0; ind < npix * npix; ind++)
  	result[ind] = {fft_frameHf0[ind][0], fft_frameHf1[ind][0]};

    return result;
}

std::vector<std::complex<double>> wiener_filtering::H_adjoint_operator(std::vector<std::complex<double>> gamma)
{
	std::vector<std::complex<double>> result(impix, {.0, .0});

	// compute the fourier transform of gamma
	for (long ind = 0; ind < impix; ind++) {
		fft_frameHb0[ind][0] = gamma[ind].real() * fftFactor;
		fft_frameHb0[ind][1] = 0.0; // gamma[ind].imag() * fftFactor;
	}

	for (long ind = 0; ind < impix; ind++) {
		fft_frameHb1[ind][0] = gamma[ind].imag() * fftFactor;
		fft_frameHb1[ind][1] = 0.0; // gamma[ind].imag() * fftFactor;
	}
	fftw_execute(planHb0);
	fftw_execute(planHb1);

	// transform from shear to convergence in fourier space
	double freqFactor = 2.0 * M_PI / ((double) npix); // / pixel_size; // 1.0 / ((double) npix); //
	// printf("freqFactor is: %0.9f\n", freqFactor);

	double k1, k2, k1k1, k2k2, k1k2, ksqr;
	double c1, c2;
	long ind;
	for (int y = 0; y < npix ; y++) {
		// k2 = (y - npix / 2) * freqFactor; //
		k2  = (y < npix / 2 ? y * freqFactor : (y - npix) * freqFactor); // this comes from python fft

		// int ky  = (y < npix / 2 ? y + npix / 2 : y - npix / 2);
		for (int x = 0; x < npix ; x++) {
			// k1 = (x - npix / 2) * freqFactor; //
			k1 = (x < npix / 2 ? x * freqFactor : (x - npix) * freqFactor); // this comes from python fft

			// int kx  = (x < npix / 2 ? x + npix / 2 : x - npix / 2);
			// long pos = ky * npix + kx;
			ind = y * npix + x;
			if (k1 == 0 && k2 == 0) {
				fft_frameHf0[ind][0] = .0;
				fft_frameHf0[ind][1] = .0;
				fft_frameHf1[ind][0] = .0;
				fft_frameHf1[ind][1] = .0;
				continue;
			}
			k1k1 = k1 * k1;
			k2k2 = k2 * k2;
			k1k2 = k1 * k2;
			ksqr = k1k1 + k2k2;
			c1 = -k2k2 + k1k1;
			c2 = 2.0 * k1k2;
			fft_frameHf0[ind][0] = (fft_frameHb0[ind][0] * c1 - fft_frameHb1[ind][0] * c2)  / ksqr ;
			fft_frameHf0[ind][1] = (fft_frameHb0[ind][1] * c1 - fft_frameHb1[ind][1] * c2)  / ksqr ;

			fft_frameHf1[ind][0] = (fft_frameHb0[ind][0] * c2 + fft_frameHb1[ind][0] * c1)  / ksqr ;
			fft_frameHf1[ind][1] = (fft_frameHb0[ind][1] * c2 + fft_frameHb1[ind][1] * c1)  / ksqr ;
		}
	}

  // compute the fourier transform
  fftw_execute(planHf0);
	fftw_execute(planHf1);

	// store computed shear to output
  for (long ind = 0; ind < impix; ind++)
    result[ind] = {fft_frameHf0[ind][0], fft_frameHf1[ind][0]};

  return result;
}

/*
 * Shift the fft of an image
 */
void wiener_filtering::fftshift(complex<double> *image, complex<double> *shiftedimg, int npix)
{

  complex<double> array2d[npix][npix];
  complex<double> rdr2d[npix][npix];

  // initialize
  int ind1 = 0;
  for (int x = 0; x < npix; x++) {
    for (int y = 0; y < npix ; y++) {
        array2d[x][y] = image[ind1];
        rdr2d[x][y] = {0.0,0.0};
        ++ind1;
		}
	}

	// reorder
	for (int x = 0; x < npix /2 ; x++){
        for (int y = 0; y < npix / 2; y++)  {
            rdr2d[x + npix / 2][y + npix / 2] = array2d[x][y];
		}
	}

	for (int x = 0; x < npix /2 ; x++){
        for (int y = npix/2; y < npix; y++)  {
            rdr2d[x + npix / 2][y - npix / 2] = array2d[x][y];
		}
	}

	for (int x = npix /2; x < npix ; x++){
        for (int y = npix/2; y < npix; y++)  {
            rdr2d[x - npix / 2][y - npix / 2] = array2d[x][y];
		}
	}

	for (int x = npix /2; x < npix ; x++){
        for (int y = 0; y < npix /2; y++)  {
            rdr2d[x - npix / 2][y + npix / 2] = array2d[x][y];
		}
	}
	// store
	int ind2 = 0;
    for (int x = 0; x < npix; x++) {
        for (int y = 0; y < npix ; y++) {
            shiftedimg[ind2] = rdr2d[x][y];
            ++ind2;
		}
	}

}


/*
 *  computes the noise power
 */
/*
void wiener_filtering::pixelize_radial_spectrum()
{

k_map =  np.zeros((size, size), dtype = float)

for (i,j), val in np.ndenumerate(power_map):

    k1 = i - size/2.0
    k2 = j - size/2.0
    k_map[i, j] = (np.sqrt(k1*k1 + k2*k2))

    if k_map[i,j]==0:
        print(i,j)
        power_map[i, j] = 1e-15
    else:
        power_map[i, j] = power_function(k_map[i, j], 10.0)
}
*/

/*
 *  computes the noise power
 */
/*
void wiener_filtering::compute_noise_pow()
{
	Ifloat Image;
	Image.alloc(npix, npix, "noise_cov");

	for (int i=0; i < npix; i++)
		for (int j=0; j < npix; j++)
			Image(i,j) = (float) Sn[i];

	fltarray NSpectrum;
	get_isotropic_spectrum(Image, NSpectrum, 1.0, 1.0);
	int NpSpec = NSpectrum.n_elem() / 3;
	printf("NSpectrum: \n");
	for (int i=0; i < NpSpec; i++)
		printf(" %f ", NSpectrum(i,1));
	printf("\n");
	printf("NSpectrum nelem: %d\n", NSpectrum.n_elem());

	// write NSpectrum to a file
	double * spt = (double *) malloc(sizeof(double) * NpSpec);
	dblarray poe;
	poe.alloc(spt, NpSpec,1);
	for (int i=0; i < NpSpec; i++)
		spt[i] = (double ) NSpectrum(i,1);
	fits_write_dblarr("NSpectrum", poe);

	//float * dat = (float *) malloc (sizeof(float) * impix);
	//for (int i=0; i < impix; i++)
	//	dat[i] = (float) Sn[i];
	//Image.buffer() = dat;

}
*/



/*
 *  computes the radially averaged power spectrum
 */
 /*
void wiener_filtering::get_isotropic_spectrum(complex<double> * fftcoeft, fltarray & Spectrum, float Resol, float Dens)
{
	int i,j;
	int Nl = npix;
	int Nc = npix;
	float Np = Nl*Nc;
	Ifloat Data(Nl,Nc, "data");
	//Icomplex_f Ima1_cf (Nl, Nc, "Buffer1 conv");
	//FFTN_2D FFT;

	//FFT.fftn2d(Imag,Ima1_cf);
	for (i=0; i < Nl; i++)
	for (j=0; j < Nc; j++) Data(i,j) =  (float) (norm(Ima1_cf(i,j)) / (double) Np);
	int N2 = (Nl+1)/2*sqrt(2.);
	int NpSpec = (int) ((N2 - 1) / Resol);

	Spectrum.alloc(NpSpec,3);
	Spectrum(0,1) = Data(Nl/2,Nc/2);
	Spectrum(0,0) = 0;
	Spectrum(0,2) = 0;

	BSPLINE_DEC BD;
	BD.SamplesToCoefficients(Data);

	// Spectrum Calculation
	for (int r = 1; r < NpSpec; r++)
	{
	  float Rad = r*Resol;
	  float L = Rad / Nl;
	  int Nu=0;
	  int Np = (int) (Dens*2*PI*Rad+0.5);  // number of points on the circle of radius r
	  for (i=0; i < Np; i++)
	  {
		 float Angle = 2. *PI / (float) Np * i;
	 double x = Nc/2 + Rad * cos(Angle);
	 double y = Nl/2 + Rad * sin(Angle);
	 if ((x >= 0) && (y >= 0) && (x < Nc) && (y < Nl))
	 {
		float I = MAX(0, BD.InterpolatedValue(Data,x,y));
		Spectrum(r,1) += I;
		Nu ++;
	 }
	 // printf("r=%5.3f, A=%5.3f, x=%5.3f, y=%5.3f, P=%5.3f\n",Rad,Angle,x,y,Spectrum(r,1));
	  }
	  if (Nu > 0) Spectrum(r,1) /= (float) Nu;
	  Spectrum(r,0) = L;
	  Spectrum(r,2) = Spectrum(r,1) * L*(L+1);
	}
}
*/
/*
// Compute the noise covariance matrix using radnom rotations of the shear, not completed yet!!!
void wiener_filtering::compute_noise_cov_random(std::vector<double> shr1, std::vector<double> shr2, int Niter)
{

	// Initialize the random number generator
    const gsl_rng_type *T;
    T = gsl_rng_default;
    gsl_rng *rng = gsl_rng_alloc(T);

	double **noise_shr1 = (double **) malloc((Niter) * sizeof(double *));
	double **noise_shr2 = (double **) malloc((Niter) * sizeof(double *));

	// allocate rotated shear arrays
	for (long n = 0; n < Niter; n++) {
		noise_shr1[n] = (double *) malloc((ngal) * sizeof(double));
		noise_shr2[n] = (double *) malloc((ngal) * sizeof(double));
	}

	// compute shear rotations
	for (long n = 0; n < Niter; n++) {
		for (long ind = 0; ind < ngal; ind++) {
			double theta1 = gsl_ran_flat(rng, 0, 2.0 * M_PI);
			double theta2 = gsl_ran_flat(rng, 0, 2.0 * M_PI);

			noise_shr1[n][ind] = wght[ind] * (shr1[ind] * cos(theta1) - shr2[ind] * sin(theta1));
			noise_shr2[n][ind] = wght[ind] * (shr2[ind] * cos(theta1) + shr1[ind] * sin(theta1));

		}
	}
	// compute their mean (still in ngal space, need to bin them as well)
	// not completed yet!!!
}
*/

/*
// allocate the binned shear, supposed to be an image!
resimg = (double **) malloc(sizeof(double *) * npix);
for (int y = 0; y < npix ; y++) {
	resimg[y] = (double *) malloc(sizeof(double) * npix);
}
*/

// read the fits image of glimpse convergence
/*
Ifloat Image;
fitsstruct Header;
io_read_ima_float("../example/kappa_499.fits", Image, &Header);
printf("Read fits file: the nl is %d\n", Image.nl());
printf("Read fits file: the nc is %d\n", Image.nc());
*/

//float * dat = Image.buffer();
//for (int i=0; i < impix; i++)
//	printf(" %f ", dat[i]);
//printf("\n");
//printf("after: the nelem is %d\n", Image.elem);

/*
fltarray Spectrum;
get_isotropic_spectrum(Image, Spectrum, 1.0, 1.0);
int NpSpec = Spectrum.n_elem() / 3;
printf("Spectrum: \n");
for (int i=0; i < NpSpec; i++)
	printf(" %f ", Spectrum(i,1));
printf("\n");
printf("Spectrum nelem: %d\n",Spectrum.n_elem());

// write xng power spectrum to a file
double * spt = (double *) malloc(sizeof(double) * NpSpec);
dblarray poe;
poe.alloc(spt, NpSpec,1);
for (int i=0; i < NpSpec; i++)
	spt[i] = (double ) Spectrum(i,1);
fits_write_dblarr("Spectrum", poe);

dblarray ps_arr1;
	fits_read_dblarr("pixspectrum.fits", ps_arr1);
*/
/*
res = (double *) malloc(sizeof(double) * ngal);
pos1 = (double *) malloc(sizeof(double) * ngal);
pos2 = (double *) malloc(sizeof(double) * ngal);
wght = (double *) malloc(sizeof(double) * ngal);
*/

//#include <sparse2d/TempArray.h>
//using namespace std;

/*

*/
/*
//Compute the power spectrum of the Gaussian component
dblarray ps_arr;
fits_read_dblarr("masked_mice_2389_catalogue_small_corrected_map_power_estim.fits", ps_arr);
*/

//printf("Init: Sn size is %lu",Sn.size());
//printf("Init: Px size is %lu",Px.size());
// allocate binned shear residuals


/*
for (int y = 0; y < npix ; y++) {
	free(resimg[y]);
	}
free(resimg);
*/

/*
for (int y = 0; y < impix ; y++) {
Sn.push_back(0.5);
}
*/
// vectors to store the wiener result
//std::vector<std::complex <double>> sr1(impix, (std::complex<double>){{0,0}});
//std::vector<std::complex <double>> sr2(impix, (std::complex<double>){{0,0}});

//printf("------- Wiener filtering --------- \n");

/*
printf("After bin2d, binned shear residuals are: \n");
for (int i = 0; i < 200; i++)
	printf("%0.9f\t", shear_res[i].real());
printf("\n");
for (int i = 0; i < 200; i++)
	printf("%0.9f\t", shear_res[i].imag());
printf("\n");
*/
/*
printf("After compute_noise_cov, noise covariance is: \n");
for (int i = 0; i < 200; i++)
	printf("%0.9f\t", Sn[i]);
printf("\n");
printf("After compute_noise_cov, noise covariance size is:  %lu\n", Sn.size());
*/

/*
printf("Glimpse solution: \n");
for (int i = 0; i < 200; i++)
	printf("%f\t", kglimpse[i]);
printf("\n");
*/

// write result to a fits file
//double * imgw = (double *) malloc(impix * sizeof(double));
//for (int ind = 0; ind < impix ; ind++) {
//	imgw[ind] = kappa[ind];
// }



	//printf("Before rwc, shear_res size is : %lu\n", shear_res.size());
	/*
	// compute the powe spectrum of kglimpse
	printf("-- compute the power spectrum of the kappa \n");
	std::vector<std::complex <double>> kvector(impix,{.0,.0});
	for (int i=0; i < impix; i++) {
			kvector[i] = {kglimpse[i], .0};
	}
	cspp(kvector);
	*/

	/*
	for (int ind = 0; ind < 30 ; ind++)
		printf(" %f+i%f", shear_res[ind].real(), shear_res[ind].imag());
	*/

	// compute the noise covariance matrix


	// compute the power spectrum of the noise
	//compute_noise_pow();

	/*
	std::vector<double> resv( 2 * ngal);
	for (int y = 0; y < impix ; y++) {
		resv.push_back(shear[y].real());
		resv.push_back(shear[y].imag());
  }
  */
	/*
	double *recn = (double *) malloc(sizeof(double)* impix);
	rec->get_convergence_map(recn);
	*/

	// try to filter the glimpse result
	//dblarray recn;
  //fits_read_dblarr("../example/kappa_499.fits", recn);

	/*
	for (int ind = 0; ind < 20 ; ind++)
		printf("%f ",recn(ind));
	printf("\n");

	printf("1\n");
	//compute the Ps of glimpse map
	for (int ind =0; ind < impix; ++ind)
		{
			frame[ind][0] = recn(ind);
			frame[ind][1] = 0;
		}
		fftw_execute(planf);

	printf("2\n");

	for (int ind = 0; ind < impix ; ind++)
		Px[ind] = 1.0 / impix  * (frame[ind][0] * frame[ind][0] + frame[ind][1] * frame[ind][1] ); //abs({,frame[ind][1]}) * abs({frame[ind][0],frame[ind][1]});
	printf("3\n");
	printf("Px size is %d\n", Px.size());
	*/


	// try to wfilter the glimpse result only
	// std::vector<std::complex<double>> resk(impix, {0.0,0.0});
	// std::vector<std::complex<double>> tmp(impix, {0.0,0.0});
	// tmp = H_operator(resk); // H * xng
	// std::transform(data.begin(), data.end(), t.begin(), t.begin(), std::minus<std::complex<double>>()); // y - H * xg
	// std::vector<double> gshear(impix);

	// store glimpse solution as a vector
	// for (int ind = 0; ind < impix ; ind++) {
	//		resk[ind] = {recn(ind),0.0};
	//	gshear.push_back(recn(ind));
	//}

	/*
	// filter the residuals
	std::vector<double> kappa = run_wiener_core(shear_res, 100);

	// write result to a fits file

	double * imgw = (double *) malloc(impix * sizeof(double));
	for (int ind = 0; ind < impix ; ind++) {
		imgw[ind] = kappa[ind];
	}
	dblarray expt;
  expt.alloc(imgw, npix,npix);
  fits_write_dblarr("wiener.fits", expt);
  */

	/*
	// compute the noise covariance
	compute_noise_cov_random(vmap1, vmap2, 200); // not completed yet!!!
	*/
	/*
	printf("--- Noise covariance ---\n");
	for (auto i: Sn)
		printf("%f ",i);
	printf("\n");
	printf("Sn size is %ld\n" , Sn.size());
	*/
	//printf("moy ah ahaha for loop finished\n");

	/*
	sr1 = run_wiener_core(&vmap1, 5000);

	compute_noise_cov(&vmap2);
	sr2 = run_wiener_core(&vmap2, 5000);
	*/

	// find the minimum of sr1 sr2 just for testing
	/*
	for (int ind = 0; ind < impix ; ind++) {
		result1[ind] = sr1[ind].real();
		result2[ind] = sr2[ind].real();
	}
	*/

	/*
	std::vector<double>::iterator itr;
	int index;
	itr = std::min_element(sr1.begin(), sr1.end());
	index  = std::distance(sr1.begin(), itr);
	double tau = sr1[index];
	printf("the minimum of shear sr1 is %f\n",tau);
	itr = std::min_element(sr2.begin(), sr2.end());
	index  = std::distance(sr2.begin(), itr);
	tau = sr2[index];
	printf("the minimum of shear sr2 is %f\n",tau);

	// write wiener result to a fits file
	double * imgw = (double *) malloc(impix * sizeof(double));

	for (int ind = 0; ind < impix ; ind++) {
		imgw[ind] = sr2[ind];
	}

	fitsfile *fptr1;
	int s1 = 0;
	long * naxes = (long *) malloc(sizeof(long) * 2);
	naxes[0] = npix;
	naxes[1] = npix;
	fits_create_file(&fptr1, "../example/wiener.fits", &s1);
	fits_create_img(fptr1,  DOUBLE_IMG, 2, &naxes[0], &s1);
	fits_write_img(fptr1, TDOUBLE, 1, impix, imgw, &s1);
    fits_close_file(fptr1, &s1);


	// array for the wiener reconstruction
	double *wrec = (double *) malloc(sizeof(double) * npix );
	// get convergence map from filtered shear residuals
	f->compute_kappa(sr1, sr2, false, wrec);
	*/

		//printf("Inside powerspec2d, the spec length is %lu\n", nspec);

		/*
		printf("Inside powerspec2d, the spectrum of t is:\n");
		for (int i = 0; i < nspec; i++)
			printf(" %0.9f", spectrum[i]);
		printf("\n");
		*/
//printf("Inside powerspec2d, new index is:\n");
//for (int i = 0; i < npix; i++)
//	printf("%d\t", ni[i]);
//	printf("\n");

	/*
	printf("Inside powerspec2d, ps2d is:\n");
	for (int i = 0; i < 200; i++)
		printf("%f\t", ps2d[i]);
	printf("\n");
	*/

	/*
	Ifloat Image;
	Image.alloc(npix, npix, "t");
	int ind = 0;
	for (int i=0; i < npix; i++)
		for (int j=0; j < npix; j++) {
			Image(i,j) = (float) std::abs(data[ind]);
			++ind;
		}

	fltarray Spectrum;
	get_isotropic_spectrum(Image, Spectrum, 1.0, 1.0);
	int NpSpec = Spectrum.n_elem() / 3;

	// it would be nice to compute the mean value of Pt

	// printf("NpSpec is: %d\n", NpSpec);
	printf("Inside csp, the spectrum of t is:\n");
	for (int i = 0; i < NpSpec; i++)
		printf(" %0.9f", Spectrum(i,1));
	printf("\n");
	*/

	/*
	printf("Inside csp, the cnt is: \n");
	for (int i = 0; i < spectrum.size(); i++)
		printf(" %d", cnt[i]);
	printf("\n");
	printf("Inside csp, spectrum size is: %lu\n", spectrum.size());
	*/
	/*
	double maxpw = 0;
	for (int i = 0; i < spectrum.size(); i++)
		maxpw = spectrum[i] > maxpw ? spectrum[i] : maxpw;
		*/
		//if (Spectrum(i,1) > 1)
		//	ps1d[i] = (double) Spectrum(i,1) - 1.;
		//else
		//	ps1d[i] = (double) Spectrum(i,1) / 10.;

		// fftshift power spectra
		//std::complex<double> * rimg = (std::complex<double> *) malloc(impix * sizeof(std::complex<double>));
		//std::complex<double> * shiftedimg = (std::complex<double> *) malloc(impix * sizeof(std::complex<double>));
		//printf("In bin2d!!");
		//for (double as : res1){
		//	printf("%f ", as);
		//}
		//printf("\n");*/

		//vmap2 = bincount(ind_1d, res2);

		/*
		printf("Vmap1 size is %ld\n" , vmap1.size());
		printf("Res1 size is %ld\n" , res1.size());
		printf("Impix is %d\n" , impix);

		printf("counter!!\n");
		for (auto i: counter)
			printf("%f ",i);
		printf("\n");
		printf("Gcount size is %ld\n" , counter.size());
		*/
		/*
		printf("\nVmap1!!\n");
		for (auto i: vmap1)
			printf("%f ",i);
		printf("\n");
		printf("Vmap1 size is %ld\n" , vmap1.size());

		printf("\nVmap2!!\n");
		for (auto i: vmap2)
			printf("%f ",i);
		printf("\n");
		printf("Vmap2 size is %ld\n" , vmap2.size());
		*/
		/*
		for (int y = 0; y < npix ; y++) {
	        for (int x = 0; x < npix ; x++) {
	            long pos = (npix - y - 1) * npix + (npix - x - 1);
	            bres[x * npix + y] = 1;
	        }

	    } */
			/*
			printf("Inside get_residuals, res1 is: \n");
			for (int i = 0; i < 200; i++)
				printf("%0.9f\t", res1[i]);
			printf("\n");
			printf("Res1 size is: %lu\n", res1.size());

			printf("Inside get_residuals, res2 is: \n");
			for (int i = 0; i < 200; i++)
				printf("%0.9f\t", res2[i]);
			printf("\n");
			printf("Res2 size is: %lu\n", res2.size());
			*/

		  //print output
			/*for (long ind = 0; ind < ngal; ind++){
				printf("%f \t %f \t %f \t %f \t %f\n", res1[ind],  res2[ind], pos1[ind], pos2[ind], wght[ind]);
			}*/



			/*
			// compute the fft of t
			for (int ind =0; ind < impix; ++ind)
			{
				frame[ind][0] = t[ind].real();
				frame[ind][1] = t[ind].imag();
			}
			fftw_execute(planf);
			*/

	// printf("After csp, Px is:\n");
	// for (int i = 0; i < 200; i++)
	// 	 printf("%0.9f\t", Px[i]);
	// printf("\n");

	/*
	printf("Px before being used:\n");
	for (int y = 0; y < 200 ; y++) {
		printf("%0.9f ",Px[y]);
  }
  printf("\n");
	printf("Pt before being used:\n");
	for (int y = 0; y < 200 ; y++) {
		printf("%0.9f ",Pt[y]);
  }
  printf("\n");
	printf("Wfc before being used:\n");
	for (int y = 0; y < 200 ; y++) {
		printf("%0.9f ",Px[y]/(Px[y]+Pt[y]));
  }
  printf("\n");
	*/


	/*
	printf("real part is :\n");
	for (int y = 0; y < 200 ; y++) {
		printf("%0.9f ",xg[y].real());
	}
	printf("\n");

	printf("imag part is:\n");
	for (int y = 0; y < 200 ; y++) {
		printf("%0.9f ",xg[y].imag());
	}
	printf("\n");
	*/
	/*
	printf("Inside csp, shifted ps2d is:\n");
	for (int i = 0; i < 200; i++)
		printf("%0.9f\t", Px[i]);
	printf("\n");
	*/
	//printf("Px size is %lu\n", Px.size());
	/*
	// write result to a fits file
	dblarray expt;
	expt.alloc(&Px[0], npix, npix);
	fits_write_dblarr("Px.fits", expt);
	*/

	/*
	// initialize residual spectrum
	ps_residual.reserve(182);
	for (int i=0; i < ps_residual.size(); i++) {
		ps_residual.push_back(.0);
	}
	*/
		/*
		if (it < 50)
			for (int i = 0; i < spectrum.size(); i++)
				ps1d[i] = std::max(spectrum[i] - mean, rpow - ps_residual[i]) ;
		else
			for (int i = 0; i < spectrum.size(); i++)
				ps1d[i] = std::abs(spectrum[i] - ps_residual[i]) * pow(.99,i);
		*/
		/* good
		for (int i = 0; i < spectrum.size(); i++)
			if (i < 50)
				ps1d[i] = std::max(spectrum[i], rpow - ps_residual[i]);
			else
				ps1d[i] = .9 * ps1d[i-1] ;

		*/
		/* good
		for (int i = 0; i < spectrum.size(); i++)
			ps1d[i] = std::max(spectrum[i] - mean, mean - ps_residual[i]) * pow(.98,i);
		*/
		/* best choice so far
		for (int i = 0; i < spectrum.size(); i++)
			if (i < 50)
				ps1d[i] = std::max(spectrum[i] - rpow * pow(.9,i), 1e-6);
			else
				ps1d[i] = (std::max(spectrum[i] - rpow, rpow - ps_residual[i])) * pow(.99,i);
		*/
		/* good
		for (int i = 0; i < spectrum.size(); i++)
			ps1d[i] = (std::max(spectrum[i] - rpow, rpow - ps_residual[i])) * pow(.99,i);
		*/
		/* good choice
		for (int i = 0; i < spectrum.size(); i++)
			if (i < 50)
				ps1d[i] = std::max(spectrum[i] - rpow * pow(.9,i), 1e-6);
			else
				ps1d[i] = (std::max(spectrum[i] - rpow * pow(.9,i), rpow * pow(.9,i)- ps_residual[i])) * pow(.99,i);
		*/

		/*
		for (int i = 0; i < spectrum.size(); i++)
			ps1d[i] = spectrum[i] - rpow* pow(.5,i); // std::min(spectrum[i] , rpow * pow(.9,i));
		*/
		/*
		for (int i = 0; i < spectrum.size(); i++)
			if (i < 50)
					ps1d[i] = spectrum[i] ;
			else
				ps1d[i] = rpow * pow(.9,i); // std::max(spectrum[i] - rpow,1e-6);
		*/
		/*
		for (int i = 0; i < 50; i++)
			ps1d[i] = std::max(spectrum[i] - rpow, pxc[i]);

		for (int i = 50; i < spectrum.size(); i++)
				ps1d[i] = pxc[i];

		*/
		/*
		if (it < 50)
			for (int i = 0; i < spectrum.size(); i++)
				ps1d[i] = std::max(spectrum[i] - rpow, 1e-6);
		else
			for (int i = 0; i < spectrum.size(); i++)
				ps1d[i] = std::max(spectrum[i] - ps_residual[i], 1e-6);
		*/
		/*
		for (int i = 0; i < spectrum.size(); i++)
			if (it < 50)
				ps1d[i] = std::max(spectrum[i] - rpow, rpow - ps_residual[i]) ; // ps_arr1(i); //std::max(spectrum[i] - rpow, 1e-10); //
			else
				ps1d[i] = std::max(spectrum[i] - rpow,1e-3) ;
				*/
		/*
		for (int i = 3; i < spectrum.size(); i++)
			if (ps1d[i] == 1e-10)
				ps1d[i] = ps1d[i-1];
		*/

		/*
		if (it>100)
			for (int i = 0; i < spectrum.size(); i++)
				ps1d[i] = std::abs(spectrum[i] - ps_residual[i]); // Px almost zero since xg starts from zero
		else{
			printf("Iteration %d\n", it);
			for (int i = 0; i < spectrum.size(); i++)
				ps1d[i] = std::max(spectrum[i] - rpow, 1e-6);
		}
		*/
		/*
		if (it < 200){
			for (int i = 0; i < spectrum.size(); i++)
				ps1d[i] = std::max(spectrum[i] - rpow, 1e-6);
		}
		else{
			for (int i = 0; i < spectrum.size(); i++)
				if (spectrum[i] - rpow < 0)
					ps1d[i] = std::abs(rpow - ps_residual[i]);
		}
		*/
		/*
		if (it < 200)
			for (int i = 0; i < spectrum.size(); i++)
				ps1d[i] = std::max(spectrum[i] - rpow, 1e-6); // pxc[i]
		else
			for (int i = 0; i < spectrum.size(); i++)
				ps1d[i] = std::max(spectrum[i] - ps_residual[i], 1e-6);
		*/
		/*
		if (it<100)
			for (int i = 0; i < spectrum.size(); i++)
				ps1d[i] = std::max(spectrum[i] - rpow, 1e-10); // or use decaying power spectrum, pxc[i]
		else
			for (int i = 0; i < spectrum.size(); i++)
				if (spectrum[i] - rpow < 0)
					ps1d[i] = ps_residual[i];
		*/
		/*
			for (int i = 0; i < spectrum.size(); i++)
			// ps1d[i] = ps_arr1(i); // read spectrum from input
			// ps1d[i] = std::max(spectrum[i] - rpow, 1e-1); //  Px = [Pt - Pr]_+, Pr = mu^2 * impix
			//if (spectrum[i] < rpow)
			//	ps1d[i] = ps1d[i-1];
			//else

			// ps1d[i] = std::abs(spectrum[i] - rpow);

		else
		for (int i = 0; i < spectrum.size(); i++)
			// ps1d[i] = ps_arr1(i); // read spectrum from input
			// ps1d[i] = std::max(spectrum[i] - rpow, 1e-1); //  Px = [Pt - Pr]_+, Pr = mu^2 * impix
			if (spectrum[i] < rpow)
					printf("Iteration %d, index %\n", );
					ps1d[i] = ps1d[i-1];
			else
				ps1d[i] = std::max(spectrum[i] - rpow, 1e-10); // use decaying power spectrum, pxc[i]
			// ps1d[i] = std::abs(spectrum[i] - rpow);
			*/

		/*
		// and smooth out any discontinuity caused by thresholding
		for (int i = 3; i < spectrum.size()-3; i++)
			ps1d[i] = (ps1d[i-1] + ps1d[i+1]) / 2;
		for (int i = 3; i < spectrum.size()-3; i++)
			ps1d[i] = (ps1d[i-1] + ps1d[i+1]) / 2;
		*/
		/*
		if (it < 200)
			for (int i = 0; i < spectrum.size(); i++)
				ps1d[i] = std::max(spectrum[i] - rpow, pxc[i]);
		else
		for (int i = 0; i < spectrum.size(); i++)
			ps1d[i] = std::max(spectrum[i] - ps_residual[i], pxc[i]);
		*/
