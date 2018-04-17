#! /bin/env python
import sys
from optparse import OptionParser
import copy
import matplotlib
matplotlib.use('Agg')
import pylab
import scipy.optimize
import numpy
from numpy import array
import dadi
import random

#import demographic models
import gillettii_models

#sample and set random seed
seed = random.randint(0, 2**32)
numpy.random.seed(seed)

def runModel(outFile, nuW_start, nuC_start, T_start, m12_start, m21_start):
	
        # Parse the SNP file to generate the data dictionary
        dd = dadi.Misc.make_data_dict('/home/rajivm/home_vol1/gillettii/data/gillettii_data_syn.AN24.dadi')
	
        # Extract the spectrum for ['YRI','CEU'] from that dictionary, with both
        fs = dadi.Spectrum.from_data_dict(dd, pop_ids=['WY','CO'], projections=[12,12], polarized=False)
	
	#uncomment following line to perform conventional bootstrap
	#fs = fs.sample()
	
        ns = fs.sample_sizes
        print 'sample sizes:', ns
	
        # These are the grid point settings will use for extrapolation.
        pts_l = [20,30,40]
        # suggested that the smallest grid be slightly larger than the largest sample size. But this may take a long time.              
	
	# bidirectional migration model
        func = gillettii_models.split_migration
        params = array([nuW_start, nuC_start, T_start, m12_start, m21_start])
        upper_bound = [10, 10, 10, 100, 100]
        lower_bound = [1e-5, 1e-10, 0, 0, 0]
	
        # Make the extrapolating version of the demographic model function.
        func_ex = dadi.Numerics.make_extrap_func(func)
        # Calculate the model AFS
        model = func_ex(params, ns, pts_l)        
        # Calculate likelihood of the data given the model AFS
	
        # Likelihood of the data given the model AFS.
        ll_model = dadi.Inference.ll_multinom(model, fs)
        print 'Model log-likelihood:', ll_model, "\n"
        # The optimal value of theta given the model.
        theta = dadi.Inference.optimal_sfs_scaling(model, fs)
	
        p0 = dadi.Misc.perturb_params(params, fold=1, lower_bound=lower_bound, upper_bound=upper_bound)
        print 'perturbed parameters: ', p0, "\n"
        popt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l, upper_bound=upper_bound, lower_bound=lower_bound, maxiter=None, verbose=len(params)) 
        print 'Optimized parameters:', repr(popt), "\n"
	
        #use the optimized parameters in a new model to try to get the parameters to converge
        new_model = func_ex(popt, ns, pts_l)
        ll_opt = dadi.Inference.ll_multinom(new_model, fs)
        print 'Optimized log-likelihood:', ll_opt, "\n"
	
        # Write the parameters and log-likelihood to the outFile
        s = str(seed) + '\t' + str(nuW_start) + '\t' + str(nuC_start) + '\t' + str(T_start) + '\t' + str(m12_start) + '\t' + str(m21_start) + '\t'
        for i in range(0, len(popt)):
		s += str(popt[i]) + '\t'
        s += str(ll_opt) + '\n'
	
        outFile.write(s)

#################

def mkOptionParser():
	""" Defines options and returns parser """
	
	usage = """%prog  <outFN> <nuW_start> <nuC_start> <T_start> <m12_start> <m21_start>
	%prog performs demographic inference on gillettii RNA-seq data. """
	
	parser = OptionParser(usage)
	
	return parser


def main():
	""" see usage in mkOptionParser. """
	parser = mkOptionParser()
	options, args= parser.parse_args()
	
	if len(args) != 6:
		parser.error("Incorrect number of arguments")
	
	outFN        = args[0]
	nuW_start    = float(args[1])
	nuC_start    = float(args[2])
	T_start      = float(args[3])
	m12_start    = float(args[4])
	m21_start    = float(args[5])
	
	if outFN == '-':
		outFile = sys.stdout
	else:
		outFile = open(outFN, 'a')
	
	runModel(outFile, nuW_start, nuC_start, T_start, m12_start, m21_start)


#run main
if __name__ == '__main__':
	main()
