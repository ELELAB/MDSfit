# load python module
module load python/3.7/modulefile

# run the script to plot from three input files, setting the color to navy blue,
# labeled species name to 'OPTN wt 169-185' and unlabeled species name to 'GST-LC3B'

python ../MDSfit.py -i optn_wt_1.csv  optn_wt_2.csv  optn_wt_3.csv -c navy -u GST-LC3B -l 'OPTN wt 169-185'
