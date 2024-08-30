# MDSfit.py (Fits MDS data to estimate Kd, Rh_free and Rh_complex)

## Description 

The script fits the MDS data from Fluidity One-M to calculate Kd, Rh_free and Rh_complex.
The equation used for fitting the data, as provided by Fluidic (Fluidic Analytics, Cambridge, UK), is:

Y=Rh_free+(Rh_complex-Rh_free)*((X+n*L+Kd)-((X+n*L+Kd)^2-4*n*L*X)^0.5)/(2*n*L)

where:

- Rh_free = Hydrodynamic radius of the labeled peptide (nm)
- Rh_complex = Hydrodynamic radius of the complex between the labeled peptide and unlabeled GST-LC3B (nm)
- n = binding sites occupied on the labeled species at saturation (i.e., stoichiometry) 
- L = the labeled peptide concentration
- Kd = the equilibrium dissociation constant or affinity of the complex

Single input data is plotted as well as the global fitting of all input data provided.
The combined plot also reports Standard error of the mean for Kd, Rh_free and Rh_complex as
well as the R squared value for the goodness of the fit (the latter is also reported on single
input data plots).
Lastly, the scripts writes a summary text file with all the relevant info for each input file analyzed. 

## Requirements

- python > 3.7
- numpy
- pandas
- sys
- argparse
- scipy
- matplotlib

## Input

- one or more input csv file obtained from the Fluidity One-M instrument,

## Usage

### Activate python
`module load python/3.7/modulefile`

### Run  
python MDSfit.py [-h] -i INPUT [INPUT ...] [-o] [-c] [-u] [-l] [-n]

### Required flags
- `-i`: input file(s) (input csv file(s) containing the data as described above)

### Optional flags
- `-o`, `--output`: name of the output summary file (default: summary.txt)
- `-c`, `--color`: color of the plot as supported by matplotlib (default: black)
- `-u`, `--unlabeled`: name of the unlabeled species to use for the plots (default: unlabeled species)
- `-l`, `--labeled`: name of the labeled species to use for the plots (default: labeled species)
- `-n`, `--normalize`: flag to plot data normalized between the fitted Rh_free (set as 0%) and the
		      fitted Rh_complex (set as 100%) (default: False)
- `-s`, `--stoichiometry`: binding sites occupied on the labeled species at saturation (integer, default: 1)

## Example
Please, see the example directory and the do.sh script therein.
## Output
The scripts writes one pdf plot file for each provided input csv file and one plot for the global 
fitting of all combined data from the csv files. Moreover, it writes a summary text file
(default: summary.txt, customizable name with the `-o` option) which contains all the relevant info 
for each analyzed input file (initial parameters used for the fitting; fitted parameters and their
errors computed as one standard deviation errors on the parameters,i.e. the square root of the variance
of the parameter estimate; the goodness of fit, i.e., the R squared value) as well as for the
global fitting (same info as for the single analyses plus the standard error of the mean for all the 
single fits). 

