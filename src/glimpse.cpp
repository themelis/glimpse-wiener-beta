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
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <CCfits/CCfits>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <ctime>

#include "version.h"
#include "survey.h"
#include "field.h"
#include "surface_reconstruction.h"
#include "wiener_filtering.h"
#include "density_reconstruction.h"
#include "gpu_utils.h"


namespace po = boost::program_options;
using namespace std;
using namespace CCfits;

int main(int argc, char *argv[])
{
    gsl_rng_env_setup();
    boost::property_tree::ptree pt;

    // List of GPU indices
    std::vector<int> IDlist;

    // Read command line arguments
    po::options_description generic("Options");
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "print help message")
#ifdef CUDA_ACC
    ("gpu,g", po::value< std::string >(), "comma separated list of GPUs to use (e.g: -g 0,1)")
#endif
    ;

    po::options_description positional("Arguments");
    positional.add_options()
    ("config", po::value< std::string >()->required(), "configuration file")
    // ("data",  po::value< std::string >()->required(), "survey data file")
    ("shear_map_a",  po::value< std::string >()->required(), "shear map a")
    ("shear_map_b",  po::value< std::string >()->required(), "shear map b")
    ("spectrum",  po::value< std::string >()->required(), "power spectrum data file")
    ("noise_cov",  po::value< std::string >()->required(), "noise covariance")
    // ("thresholds",  po::value< std::string >()->required(), "threshold from shear rotation")
    ("output", po::value< std::string >()->required(), "output file name");
    po::positional_options_description positionalOptions;
    positionalOptions.add("config", 1);
    // positionalOptions.add("data", 1);
    positionalOptions.add("shear_map_a", 1);
    positionalOptions.add("shear_map_b", 1);
    positionalOptions.add("spectrum", 1);
    positionalOptions.add("noise_cov", 1);
    // positionalOptions.add("thresholds", 1);
    positionalOptions.add("output", 1);

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(positional);

    po::variables_map vm;

    // Process generic options
    try {
        po::store(po::command_line_parser(argc, argv)
                  .options(generic)
                  .run(), vm);
        po::notify(vm);
    } catch (po::error &e) {
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        std::cerr << cmdline_options << std::endl;
        return 1;
    }
    if (vm.count("help")) {
        cout << cmdline_options << "\n";
        return 1;
    }

    if (vm.count("version")) {
         cout << VERSION << "\n";
         return 1;
    }

    // In case of GPU acceleration, parse the list of GPUs
#ifdef CUDA_ACC
    if (vm.count("gpu")) {

        // Check that none of the requested GPU id is larger than nGPU
        int nGpu;
        int whichGPUs[MAX_GPUS];
        cudaGetDeviceCount(&nGpu);

        std::vector<std::string> strs;
        boost::split(strs,vm["gpu"].as<std::string>(),boost::is_any_of(",;"));

        for(int i =0; i < strs.size(); i++){
            IDlist.push_back(boost::lexical_cast<int>(strs[i]));
        }

        if(IDlist.size() > MAX_GPUS){
            cout << "ERROR: Requested more GPUs than maximum number;"<< endl;
            cout << "Maximum size of GPU array " << MAX_GPUS << endl;
            return -1;
        }
        if(IDlist.size() > nGpu){
            cout << "ERROR: Requested more GPUs than available;"<< endl;
            cout << "Number of GPUs available: " << nGpu << endl;
            return -1;
        }

        for(int i=0; i < IDlist.size(); i++){

            if(IDlist[i] >= nGpu){
                cout << "ERROR: Requested GPU id not available;"<< endl;
                cout << "Maximum GPU id : " << nGpu - 1 << endl;
                return -1;
            }
            whichGPUs[i] = IDlist[i];
        }
        setWhichGPUs( IDlist.size(), whichGPUs);
    }
#endif

    // Process positional arguments
    try {
        po::store(po::command_line_parser(argc, argv)
                  .options(cmdline_options)
                  .positional(positionalOptions)
                  .run(), vm);
        po::notify(vm);
    } catch (po::error &e) {
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        std::cerr << cmdline_options << std::endl;
        return 1;
    }

    // Load the configuration file
    try {
        boost::property_tree::ini_parser::read_ini(vm["config"].as<std::string>(), pt);
    } catch (boost::property_tree::ini_parser_error e) {
        std::cout << e.what() << std::endl;
        exit(-1);
    }

    /* ________________________________________
     * read newer options from file, such as the power spectrum,
     * noise covariance, the shear map
     * ________________________________________
    */

    // read the shear map a
    const char * shaname = vm["shear_map_a"].as<std::string>().c_str();
    try {
      std::auto_ptr<FITS> psFits(new FITS(shaname, Read));
    }
    catch (FITS::CantOpen)
    {
        std::cerr << "ERROR: Cant open shear map a FITS file" << std::endl;
        return -1;
    }
    dblarray sheara_array;
    fits_read_dblarr(shaname, sheara_array);

    // read the shear map b
    const char * shbname = vm["shear_map_b"].as<std::string>().c_str();
    try {
      std::auto_ptr<FITS> psFits(new FITS(shbname, Read));
    }
    catch (FITS::CantOpen)
    {
        std::cerr << "ERROR: Cant open shear map b FITS file" << std::endl;
        return -1;
    }
    dblarray shearb_array;
    fits_read_dblarr(shbname, shearb_array);

    // read the power spectrum
    const char * psname = vm["spectrum"].as<std::string>().c_str();
    try {
      std::auto_ptr<FITS> psFits(new FITS(psname, Read));
    }
    catch (FITS::CantOpen)
    {
        std::cerr << "ERROR: Cant open power spectrum FITS file" << std::endl;
        return -1;
    }
    // read the powerspectrum, needed for wiener filtering
    dblarray ps_array;
    fits_read_dblarr(psname, ps_array);

    // read the noise covariance
    const char * nsname = vm["noise_cov"].as<std::string>().c_str();
    try {
      std::auto_ptr<FITS> psFits(new FITS(nsname, Read));
    }
    catch (FITS::CantOpen)
    {
        std::cerr << "ERROR: Cant open noise map a FITS file" << std::endl;
        return -1;
    }
    dblarray ns_array;
    fits_read_dblarr(nsname, ns_array);

    // // read the thresholds
    // const char * thrname = vm["thresholds"].as<std::string>().c_str();
    // try {
    //   std::auto_ptr<FITS> psFits(new FITS(thrname, Read));
    // }
    // catch (FITS::CantOpen)
    // {
    //     std::cerr << "ERROR: Cant open noise map a FITS file" << std::endl;
    //     return -1;
    // }
    // dblarray thr_array;
    // fits_read_dblarr(thrname, thr_array);

    // double *ncov = (double *) malloc(sizeof(double)* ns_array.nx() * ns_array.ny());
    // double *g1map = (double *) malloc(sizeof(double)* sheara_array.nx() * sheara_array.ny());
    // double *g2map = (double *) malloc(sizeof(double)* shearb_array.nx() * shearb_array.ny());

    // Create survey object and load data
    survey *surv = new survey(pt);
    // surv->load(vm["data"].as<std::string>());

    // Initialize lensing field
    field *f = new field(pt, surv, sheara_array, shearb_array, ns_array); // , ns_array

    // Open output fits file before performing the actual reconstruction
    long naxes[3];
    naxes[0] = f->get_npix();
    naxes[1] = f->get_npix();
    naxes[2] = f->get_nlp();
    long naxis = naxes[2] == 1 ? 2 : 3; // Saving surface as 2d image, not 3d cube

    std::auto_ptr<FITS> pFits(0);
    try
    {
        // overwrite existing file if the file already exists.
        std::stringstream fileName;
        fileName << "!" << vm["output"].as<std::string>();

        pFits.reset( new FITS(fileName.str() , DOUBLE_IMG , naxis , naxes));
    }
    catch (FITS::CantCreate)
    {
        std::cerr << "ERROR: Cant create output FITS file" << std::endl;
        return -1;
    }

    // Array holding the reconstruction
    // cout << "sdf0" << endl;
    double *reconstruction = (double *) malloc(sizeof(double)* f->get_npix() * f->get_npix() * f->get_nlp());
    // cout << "sdf1" << endl;
    // If CUDA is available, give the option to reconstruct in 3D
    if (f->get_nlp() > 1) {
        //density_reconstruction rec(pt, f);
        // rec.reconstruct();

        // Extracts the reconstructed array
        // rec.get_density_map(reconstruction);
    } else {
        // Initialize reconstruction object
        // cout << "sdf2" << endl;
        surface_reconstruction rec(pt, f); //, thr_array

        cout << "Glimpse: surface reconstruction object allocated " << endl;

        clock_t begin = clock();

        // rec.reconstruct();

        // clock_t end = clock();
        // double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        // printf("Elapsed time : %f\n", elapsed_secs);


        // cout << "Surface reconstruction object initialized" << endl;
        // cout << "Reconstruct returned" << endl;

        // Extracts the reconstructed array
        // rec.get_convergence_map(reconstruction);
        /*
        printf("Glimpse solution: \n");
      	for (int i = 0; i < 200; i++)
      		printf("%f\t", reconstruction[i]);
      	printf("\n");
        cout << flush;
        */

    		// run wiener filtering on the shear residuals
    		cout << "Running wiener filtering" << endl;

        // Initialize wiener filtering class
        wiener_filtering *wf = new wiener_filtering(&rec, f, surv, sheara_array, shearb_array, ps_array, ns_array);
        cout << "before Wiener completed" << endl;

        // cout << "Wiener filtering object initialized" << endl;
        wf->main_iteration();
        cout << "Wiener completed" << endl ;
        std::vector<std::complex<double>> rwshear = wf->compute_wiener_residuals();

        // cout << "print Wiener shear residuals" << endl;
      	// for (int i = 0; i < 200; i++)
      	// 	cout << rwshear[i].real() << " " << sheara_array(i) << " ";
      	// cout << endl;

        rec.reconstruct(rwshear);
        rec.get_convergence_map(reconstruction);
        cout << "Glimpse completed" << endl;
    }

    // Number of elements in the array
    // cout <<  naxes[0] << naxes[1] << naxes[2] << endl ;
    long nelements =  naxes[0] * naxes[1] * naxes[2];
    std::valarray<double> pixels(reconstruction, nelements);
    long  fpixel(1);

    // cout << "fpixel" << fpixel << endl ;
    // Write primary array
    pFits->pHDU().write(fpixel,nelements,pixels);
    // cout << "write result" << endl ;
    // std::vector<long> extAx(2, naxes[0]);
    // string raName ("RA");
    // ExtHDU* raExt = pFits->addImage(raName,DOUBLE_IMG,extAx);
    //
    // string decName ("DEC");
    // ExtHDU* decExt = pFits->addImage(decName,DOUBLE_IMG,extAx);
    //
    // double *  ra  = (double *) malloc(sizeof(double) * naxes[0] * naxes[0]);
    // double *  dec = (double *) malloc(sizeof(double) * naxes[0] * naxes[0]);
    //
    // f->get_pixel_coordinates(ra, dec);
    //
    // long nExtElements = naxes[0] * naxes[0];
    // std::valarray<double> raVals(ra, nExtElements);
    // std::valarray<double> decVals(dec, nExtElements);
    //
    // raExt->write(fpixel, nExtElements, raVals);
    // decExt->write(fpixel, nExtElements, decVals);

    free(reconstruction);

    // free(ra);
    // free(dec);
    // delete *wf;
    // delete f;
    // delete surv;

    return 0;
}
