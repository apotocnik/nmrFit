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

# Documentation
The use of the program is very similar to [eprFit](https://github.com/apotocnik/eprFit).
