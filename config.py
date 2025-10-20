import numpy as np

# path of radex bin file
radexBin = '/Users/shaobo/AstroTools/Radex/src/radex'
# frequency range for Radex to calculate lines
radexFreqRange = (1, 1800) # GHz
# default bounds of molecular abundance
XB = {'co': [8e-5, 2e-4], 'hcn': [1e-9, 5e-7], 'cs': [1e-9, 5e-7], 'hcop': [1e-9, 5e-7], 'ci': [1e-5, 2e-3]}


# constant & assumption
Tbg = 2.7315 #CMB
opr = 3 #Assuming ortho to para ratio equals 3
fortho = opr/(1+opr) #fraction of ortho


# path to ram disk
tmpPath = '/Volumes/RAMDisk'
# path where you put molecule data
molPath = './radex_moldata'
# path put fitting results
outPath = './out'
# path put figures
figPath = './fig'
# path of radex bin file


# number of components
nComponent = 2 # cold and warm
# suffix of output files to distinguish each session
suffix = '_6fp'


paramsConfig = {} 
paramsConfig['nH2']  = [[4.781, 3.337],	'free',	(1.5, 6.0),		r'$lg n_{H_2}$']	#density, MUST be 'free'
paramsConfig['Tkin'] = [[0.879, 2.222],	'free',	(np.log10(Tbg), 2.8),	r'$lg T_{Kin}$']	#Tkin, MUST be 'free'
paramsConfig['Xco']  = [-4.10,	'fix',	np.log10(XB['co']),	r'$lg X_{CO}$']		#abundance
paramsConfig['Xci']  = [-4.52,	'fix',	np.log10(XB['ci']),	r'$lg X_{CI}$']		#abundance
paramsConfig['Omega']= [None,	'free',	(-5, 5), 		r'$lg \Omega$']		#solid angle, shrink by one magnitude for each component
