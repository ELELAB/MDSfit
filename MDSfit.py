#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import sys
import argparse
from scipy.optimize import curve_fit
from scipy.stats import sem
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.ticker as mtick

# set the pdf to be editable with external softwares
plt.rcParams['pdf.fonttype'] = 42

# Define the binding equation with the provided formula
def binding_equation(concentration, Rh_free, Rh_complex, Kd):
    
    '''
    The function takes into input concentration, Rh_free, Rh_complex and Kd 
    values and return the computed result.

    -----
    INPUT
    _____

    concentration: float, concentration of the unlabeled specie (μM)
    Rh_free: float, hydrodynamic radius of the free state (nm)
    Rh_complex: float, hydrodynamic radius of the complex (nm)
    Kd: float, equilibrium dissociation constant (μM)

    ------
    RETURN
    ------
    The function returns the result of the equation given the inputs.

    '''

    X = concentration

    return Rh_free + (Rh_complex - Rh_free) * ((X + args.stoichiometry  * float(labeled_conc) \
           + Kd) - np.sqrt((X + args.stoichiometry * float(labeled_conc) + Kd) ** 2 -  \
           4 * args.stoichiometry * float(labeled_conc) * X)) / (2 * args.stoichiometry * float(labeled_conc))

def goodness_fit(ydata,xdata):
    
    '''
    The function calulates the goodness of fit, in terms of R squared, of the
    fitting of the curve procedure.

    -----
    INPUT
    -----
    ydata: float, data coordinates on the y axis
    xdata: float, data coordinates on the x axis

    ------
    RETURN
    ------
    The function returns the R squared value.
    '''

    # first compute the residual sum of squares (ss_res)
    residuals = ydata - binding_equation(xdata, *popt)
    ss_res = np.sum(residuals ** 2)
    
    # calculate the total sum of squares (ss_tot)
    ss_tot = np.sum((ydata - np.mean(ydata)) ** 2)
    
    # finally calculate r_squared
    r_squared = 1 - (ss_res / ss_tot)
    
    return r_squared

# Argument parser
parser = argparse.ArgumentParser()

# Required arguments
requiredNamed = parser.add_argument_group('required named arguments')

requiredNamed.add_argument('-i', '--input',
                    nargs = '+', 
                    help = 'input csv files from the Fluidity One-M instrument',
                    required = True)
# Optional arguments
parser.add_argument('-o', '--output', 
                    help = 'output name of the summary file',
                    default = 'summary.txt',
                    type = str)
parser.add_argument('-c', '--color', 
                    help = 'color of the plot',
                    default = 'black',
                    type = str)
parser.add_argument('-u', '--unlabeled', 
                    help = 'name of the unlabeled species',
                    default = 'unlabeled species',
                    type = str)
parser.add_argument('-l', '--labeled', 
                    help = 'name of the labeled species',
                    default = 'labeled species',
                    type = str)
parser.add_argument('-n', '--normalize',
                    help = 'to plot normalized data',
                    default = False,
                    action = 'store_true')
parser.add_argument('-s', '--stoichiometry',
                    help = 'binding sites occupied on the labeled species at saturation',
                    default = 1,
                    type = int)
args = parser.parse_args()

# initialize combined lists for Kd, Rh_free and Rh_complex values
fit_kds = []
fit_RhFs = []
fit_RhCs = []

# Load your data from the CSV files
count = 0
for i in args.input:

    count += 1



    # remove rows not containing the radius measurement
    #data_conc = data_conc.dropna(axis = 0, subset = ['Sample Radius nm'])

    # store the whole csv file
    data = pd.read_csv(i)

    start_index = data[data['Fluidity One-M Results File'] == ('Circuit ID')].index[0]
    # table with concentrations and responses only
    data_conc = pd.read_csv(i, skiprows = start_index + 1)

    # remove rows not containing the radius measurement
    data_conc = data_conc.dropna(axis = 0, subset = ['Sample Radius nm'])

    # Access the concentrations and responses from the specified columns
    labeled_conc = (data_conc.iloc[0,2])
    concentrations = (data_conc.iloc[:, 4])  
    responses = (data_conc.iloc[:, 14])

    # Access the values caluclated by the instrument
    Rh_free_search = data[data['Fluidity One-M Results File'] == 'Binding Rh free']
    Rh_free = (float(Rh_free_search.iloc[0,1]))

    Rh_complex_search = data[data['Fluidity One-M Results File'] == 'Binding Rh complex']
    Rh_complex = (float(Rh_complex_search.iloc[0,1]))

    kd_search = data[data['Fluidity One-M Results File'] == 'Binding KD nM']
    Kd = (float(kd_search.iloc[0,1]))
    start_index = data[data['Fluidity One-M Results File'] == ('Circuit ID')].index[0]

    # append concentrations and responses to combined data 
    if count == 1:
        conc_tot = concentrations
        resp_tot = responses

    else:
        conc_tot = conc_tot.append(concentrations, 
                                   ignore_index = True)
        resp_tot = resp_tot.append(responses, 
                                   ignore_index = True)

    # Convert the concentrations and responses to NumPy arrays 
    concentrations = concentrations.values
    responses = responses.values

    # Initial guess for parameters (Rh_free, Rh_complex, Kd)
    initial_guess = (Rh_free, Rh_complex, Kd)


    # Fit the data and calculate the error
    popt, pcov = curve_fit(binding_equation, 
                           concentrations, 
                           responses, 
                           p0 = initial_guess,
                           check_finite = False)

    perr = np.sqrt(np.diag(pcov))

    # compute the R^2 goodness of fit
    r_squared = goodness_fit(responses, concentrations)

    # Access fitted parameters
    fitted_Rh_free, fitted_Rh_complex, fitted_Kd = popt

    # append fitted Kd, Rh_free and Rh_complex to combined lists
    fit_kds.append(fitted_Kd)
    fit_RhFs.append(fitted_Rh_free)
    fit_RhCs.append(fitted_Rh_complex)

    # Normalize data if requested
    if args.normalize:
        responses = (responses - fitted_Rh_free) / \
                     (fitted_Rh_complex - fitted_Rh_free)

    # Plot the results
    plt.scatter(concentrations, 
                responses, 
                c = args.color, 
                alpha = 0.3)

    x_fit = np.linspace(min(concentrations),
                        max(concentrations), 
                        10000)
    if args.normalize:
        y_fit = binding_equation(x_fit, 
                                 0, 
                                 1, 
                                 fitted_Kd)
    else:
        y_fit = binding_equation(x_fit, 
                             fitted_Rh_free, 
                             fitted_Rh_complex, 
                             fitted_Kd)

    plt.plot(x_fit, 
             y_fit, 
             color = args.color)

    plt.xlabel(args.unlabeled + ' concentration ' + u'\u03bcM')
    plt.ylabel('Sample Radius (nm)')
    if args.normalize:
        plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(xmax = 1.0))
        plt.ylabel('% Bound')
    plt.xscale('log')  # Set the x-axis to a logarithmic scale

    fit_line = mlines.Line2D([], [], 
                             color = args.color, 
                             marker = 'o',
                             markersize = 6, 
                             label = args.labeled)
    kd= mlines.Line2D([], [],
                      color = 'white', 
                      marker = '',
                      markersize = 0, 
                      label = 'Kd = ' + str('%.2f' %round(fitted_Kd, 2)) + \
                             ' ' + u'\u03bcM' + '\n' \
                             'R$^2$ = ' + str('%.2f' %round(r_squared, 2)))

    plt.legend(handles = [fit_line, kd],
               loc = 'upper left', 
               shadow = True)

    plt.title(i)
    plt.savefig('plot_' + str(count) + '.pdf',
                format = 'pdf')
    plt.close()
    with open(args.output, 'a') as o:
        o.write('*' * (12 + len(i)) + '\n' \
                'Experiment: ' + i + '\n' + '*' * (12 + len(i)) + '\n' * 2 + \
                '-' * 15 + '\n' + \
                'Initial guesses' + '\n' + \
                '-' * 15 + '\n' * 2 + \
                'Rh_free: ' + str(Rh_free) + '\n' + \
                'Rh_complex: ' + str(Rh_complex) + '\n' + \
                'Kd: ' + str(Kd) + '\n' * 2 + \
                '-' * 17 + '\n' + \
                'Fitted parameters' + '\n' + \
                '-' * 17 + '\n' * 2 + \
                'Fitted Rh_free: ' + str(fitted_Rh_free) + '\n' + \
                'Rh_free error: ' + str(perr[0]) + '\n' * 2 + \
                'Fitted Rh_complex: ' + str(fitted_Rh_complex) + '\n' + \
                'Rh_complex error: '  + str(perr[1]) + '\n' * 2 + \
                'Fitted Kd: ' + str(fitted_Kd) + '\n' + \
                'Kd error: '+ str(perr[2]) + '\n' + \
                'R$^2$ : ' + str(r_squared) + '\n' * 4)

# calculate SEM for Kd, Rh_free and Rh_complex
kd_sem = sem(fit_kds)
RhF_sem = sem(fit_RhFs)
RhC_sem = sem(fit_RhCs)

# Initial guess for parameters (Rh_free, Rh_complex, Kd)
initial_guess = (Rh_free, Rh_complex, Kd)

# Fit the data and calculate the error
popt, pcov = curve_fit(binding_equation, 
                       conc_tot, 
                       resp_tot, 
                       p0=initial_guess,
                       check_finite=False)

perr = np.sqrt(np.diag(pcov))

# compute the R^2 goodness of fit
r_squared = goodness_fit(resp_tot, conc_tot)

# Access fitted parameters
fitted_Rh_free, fitted_Rh_complex, fitted_Kd = popt

# Normalize data if requested
if args.normalize:
    resp_tot = (resp_tot - fitted_Rh_free) / \
               (fitted_Rh_complex - fitted_Rh_free) 

# Plot the results
plt.scatter(conc_tot, 
            resp_tot, 
            c = args.color, 
            alpha = 0.3)

x_fit = np.linspace(min(conc_tot), 
                    max(conc_tot), 
                    10000)
if args.normalize:
    y_fit = binding_equation(x_fit, 
                                 0, 
                                 1, 
                                 fitted_Kd)
else:
    y_fit = binding_equation(x_fit, 
                         fitted_Rh_free, 
                         fitted_Rh_complex, 
                         fitted_Kd)
plt.plot(x_fit, 
         y_fit, 
         color = args.color)

plt.xlabel(args.unlabeled + ' concentration ' + u'\u03bcM')
plt.ylabel('Sample Radius (nm)')
if args.normalize:
    plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(xmax = 1.0))
    plt.ylabel('% Bound')
plt.xscale('log')  # Set the x-axis to a logarithmic scale

fit_line = mlines.Line2D([], [], 
                         color = args.color, 
                         marker = 'o',
                         markersize = 6, 
                         label = args.labeled)
info= mlines.Line2D([], [],
                  color = 'white', 
                  marker = '',
                  markersize = 0, 
                  label = 'Kd = ' + str('%.2f' %round(fitted_Kd,2)) + \
                        u'\u00B1' + str('%.2f' %round(kd_sem,2)) + ' ' + \
                        u'\u03bcM' + '\n' + \
                        'Rh_free = ' + str('%.2f' %round(fitted_Rh_free,2)) + \
                        u'\u00B1' + str('%.2f' %round(RhF_sem,2)) + \
                        ' nm' + '\n' + \
                        'Rh_complex = ' + \
                        str('%.2f' %round(fitted_Rh_complex,3)) + \
                        u'\u00B1' + str('%.2f' %round(RhC_sem,2)) + \
                        ' nm' + '\n' + \
                        'R$^2$ = ' + str('%.2f' %round(r_squared, 2)))

plt.legend(handles = [fit_line, info],
           loc = 'upper left', 
           shadow = True)

plt.title('All replicas')
plt.savefig('all_replicas.pdf',format = 'pdf')
plt.close()
with open(args.output,'a') as o:
    o.write('*' * 21 + '\n' \
            'All replicas combined' + '\n' + \
            '*' * 21 + '\n' * 2 + \
            '-' * 15 + '\n' \
            'Initial guesses' + '\n' + \
            '-' * 15 + '\n' * 2 + \
            'Rh_free: ' + str(Rh_free) + '\n' + \
            'Rh_complex: ' + str(Rh_complex) + '\n' + \
            'Kd: ' + str(Kd) + '\n' * 2 + \
            '-' * 17 + '\n' \
            'Fitted parameters' + '\n' +\
            '-' * 17 + '\n' * 2 + \
            'Fitted Rh_free: ' + str(fitted_Rh_free) + '\n' + \
            'Rh_free SEM: ' + str(RhF_sem) + '\n' + \
            'Rh_free error: ' + str(perr[0]) + '\n' * 2 + \
            'Fitted Rh_complex: ' + str(fitted_Rh_complex) + '\n' + \
            'Rh_complex SEM: ' + str(RhC_sem) + '\n' + \
            'Rh_complex error: ' + str(perr[1]) + '\n' * 2 + \
            'Fitted Kd: ' + str(fitted_Kd) + '\n' + \
            'Kd SEM: ' + str(kd_sem) + '\n' * 2 + \
            'Kd error: ' + str(perr[2]) + '\n' + \
            'R$^2$ : ' + str(r_squared) + '\n' * 4)
