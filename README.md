# ShotgunLipidomicsR
A personal utility package with a collection of functions for analysis of direct-infusion mass spectrometry data. Focus is on shotgun lipidomics data acquired on a Thermo LTQ-Orbitrap XL.

This is not a sophisticated package, but just a personal collection of functions, that can be somewhat useful at least for me. This package is intended to be used with other R packages like xcms, MALDIquant, Rdisop etc (basically, helper/wrapper functions). Typical workflow would be read mzML file into R and plot MS1 or MS2 of interest with some peak annotation. As Thermo's RAW files and mzML files from them sometimes contains "incorrect" precursor m/z values (and only to 2 decimals), due to automatic isotope pattern recognition (even if this functionality is disabled at at the acquisition stage). So this package contains a function that "correct" the precursor m/z based on target precursor isolation window listed in the files (just select the most intense peak within the window).
