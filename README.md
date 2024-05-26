# nmrFit
is a Matlab-based tool for analysis of nuclear magnetic resonance (NMR) spectra. The program allows to load a series of spectra, typically as a function of temperature, and then perform various analysis on the entire series.

## Features
- Automatic fast fourier transform with auto phase and auto SHL option.
- Spectral moment analysis
- Fitting spectra with numerical simulations (fminsearch)
- Background correction
- Spectra multiplication, substration, shift, and normalization
- Numerical functions are defined in separate file, so one can easily add new functions.

## Recognized file formats
- Jozef Stefan Insitute (JSI) 7NMR file format
- Andraž Krajnc developed conversion tool from Bruker NMR data files to 7NMR format
- Two column ascii data file (spectre only)

## Documentation
The use of the program is very similar to [eprFit](https://github.com/apotocnik/eprFit).


## What is new?

13.7.2014
- removed bug when calculating chi2
- removed phase correction on specter controls by setting visibility to false in the fft dialog
- temporarly removed Kanizotropy+Lorentzian simulation functions
- quadrupole simulation working also on 32bit systems?
- equation for Kanizotropy corrected at the Kasy term. Now the definition is the same as in La2C3 paper.

7.7.2014
- Removed bug with Pointers
- Removed bug with mFID
- removed option show modified FID
- removed option with Error in the moment analysis
- removed old simulation functions
- under simulations static variables are by default fixed

5.5.2014
- When set range correct also spcOld, or save fOld
- added additional FFT options
- LB is working correctly with DE option
- added mFID

9.1.2014: 
- Show NSC instead of NS. Since NSC is actual number of scans and not only the setted one.
