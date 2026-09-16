MASS - Basic 
--

Our group has a handful of simple MASS analysis codes we pass around for PAX analysis. 
The target of this notebook is to give a rough outline of the workflow we tend to use with only basic python packages + MASS. 

Checking summary data: 

![channel summary](mass_basic/run22ch9_summaries.png)


Checking Noise: 

![Run Noise](mass_basic/run22_noise.png)


And most importantly: Generating Spectra !!

![Spectrum](mass_basic/run22coadd_spectrum_nominal.png)

![Prompt Spectrum](mass_basic/run22coadd_prompt_spectrum_nominal.png)


As well as: 
- Monitoring drift correction performance
- Building the optimal filter
- Calibrating data with the mass GPR splines 


Note
--

You will need the .py files in ![Utils](Utils) to run this.

There are some bugs that were fixed there and making some of the plot functions are hidden away there aswell. 

