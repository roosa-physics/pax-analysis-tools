# Recalibrating PAX data

When doing the initial analysis we used the built-in MASS autocalibrate tools.
This tool is great but is not built to replace a real science level calibration.



We need, additionally:
* Tools to fit the non-trivial tin x-ray escapes
* Fitting uncertainties for uncertainty quantification of the calibaration splines.


Contents
--


Presented as a "work flow"
1. Sn-Ka preparation:
    * Takes MDFGME outputs and makes a model spectrum

2. Detector_Response:
    * Convolve model spectrum with detector response

3. Fitting (to-be-added):
    * Fits the convolved spectrum to the data and extracts calibration knots

4. Apply_Recalibration (to-be-added)
    * Takes calibration knots, makes a spline
    * Reads old csv data and generates new energies
    * logs spectrum resolution and fit uncertainties
