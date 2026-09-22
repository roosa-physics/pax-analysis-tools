Here I have some ROOT scripts I use to generate templates from MDFGME outputs that Jorge sent me.

I assume access to a CSV with columns of energy, intensity and width.

I use the maketemplates.C to dump to terminal the outputs:

$$ root maketemplates.C > template.csv

or equivalent, and clean it up after. 


Then I plot the template lines with plotlines_template.C

plotlines skips the saving process and is just nice if you want to look at the spectrum. 





* not an exact copy of the working directory