WGlimpse

WGlimpse is a sparsity based mass-mapping algorithm. See the
[Glimpse page](http://www.cosmostat.org/software/glimpse) on the
CosmoStat website for more information. The code computes the Wiener 
estimate and then runs Glimpse-on-a-grid on the shear residuals 

## Installation

### Requirements

Glimpse depends on the following software:

* **FFTW** version 3.3 or later
* **cfitsio** and **CCFits**
* **GSL**
* **CMake** version 2.8 or later

These dependencies can easily be installed using a package manager:

* Setting up requirements on **Linux**:
  Simply use the package manager provided by your Linux distribution, for instance on Ubuntu Linux:

    $ sudo apt-get install cmake libgsl0-dev libfftw3-3  libccfits-dev libnfft3-dev

* Setting up requirements on **MacOs X**:
  The preferred installation method for the dependencies is through [MacPorts](https://www.macports.org):
    
    $ sudo port install cmake libgsl0-dev pkgconfig gsl fftw-3 nfft-3
  
  CCFits needs to be installed manually, the sources can be found [here](http://heasarc.gsfc.nasa.gov/fitsio/ccfits/).

### Compilation

Once the requirements are installed, Glimpse can be compiled by running the  following commands:

    $ cd Glimpse
    $ mkdir build
    $ cd build
    $ cmake ..
    $ make

This will create a Glimpse executable in the build directory.

## 2D Usage


All the options of the reconstruction algorithm can be specified in a *config.ini* file such as the one provided in the *example* directory.


WGlimpse can be run with the following command line:

    $ rm *.fits; ./glimpse ../data/mice/fits_files/config2d.ini ../data/mice/fits_files/mice_g1_map_0.fits ../data/mice/fits_files/mice_g2_map_0.fits ../data/mice/fits_files/mice_cosmosis_ps1d_kappa.fits ../data/mice/fits_files/mice_noisecov_map.fits ../data/mice/fits_files/kappa_mice_fast_glimpse_0.fits > ../data/mice/output.txt

Where the input arguments are

config2d.ini                   -> the config file for fast glimpse
mice_g1_map_0.fits             -> the binned g1 map
mice_g2_map_0.fits             -> the binned g2 map
mice_cosmosis_ps1d_kappa.fits  -> the 1d powerspectrum for the mice map based on the cosmosis module
mice_noisecov_map.fits         -> the noise covariance map having a DES-like mask
kappa_mice_fast_glimpse_0.fits -> the output file


