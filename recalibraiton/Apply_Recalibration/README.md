Apply Recalibration 
--

This is pretty straight forward: 
1. I open the merged fit result file and, for each channel
  - Pull the fit results from the TGraphErorrs
  - Pull the fit pointer for the 122 peak

2. I make a calibraiton spline for each channel
  - Calulate the interpolation uncertainty and FWHM

3. I open the TTree with the old calibrated data

4. Copy that data and the recalculated energy / resolution / uncertainty to a new tree. 



Side Note 
--
Steps 1 and 2 basically set me up to do more in-depth UQ with the splines. 
I could try diffrent spline strategies, calculate drop 1's or do a sensitivity check with my knot uncertainties. 
