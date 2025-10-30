# Multi-species Multi-components MCMC Radex fitting  (MMMRadex)

### IMPORTANT
# Before running, create a RAM disk
# $ diskutil erasevolume HFS+ RAMDisk $(hdiutil attach -nomount ram://1048576)
# Remember to unmount it after you finish the fitting work
# $ diskutil unmount /Volumes/RAMDisk
# Error occurs if plotting only without runMCMC, try copy all data files in molPath to tmpPath

### require Radex, available at
#https://sronpersonalpages.nl/~vdtak/radex/index.shtml#details

import os, sys, re, shutil, tempfile, logging
import multiprocessing, subprocess, emcee
# For runing the code on the clusters
import _pickle as pickle
import warnings
import numpy as np
import pandas as pd
from astropy import units as u
import astropy.constants as con
from scipy.optimize import curve_fit, minimize
from scipy.interpolate import interp1d
from scipy.stats import norm
from scipy.ndimage import gaussian_filter

# to match the cosmo in the CO paper
from astropy.cosmology import Planck18 as cosmo
#from astropy.cosmology import FlatLambdaCDM
#cosmo = FlatLambdaCDM(H0=67.8 * u.km / u.s / u.Mpc, Om0=0.308)

import matplotlib.pyplot as plt


# path of radex bin file
radexBin = '/Users/shaobo/AstroTools/Radex/src/radex'
# frequency range for Radex to calculate lines
radexFreqRange = (1, 1800) # GHz


# path to ram disk
tmpPath = '/Volumes/RAMDisk'
# path where you put molecule data
molPath = './radex_moldata'
# path put fitting results
outPath = './out'
# path put figures
figPath = './fig'
# path of radex bin file


# constant & assumption
Tbg = 2.7315 #CMB
opr = 3 #Assuming ortho to para ratio equals 3
fortho = opr/(1+opr) #fraction of ortho

ncpu = multiprocessing.cpu_count() //3*2

#pre-define
lg = np.log10
kms = u.km / u.s
Jykms = u.Jy * kms

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
logging.getLogger('matplotlib').setLevel(logging.INFO)



def fp2p(free_params, nComponent):
	#free params -> mixed params with free, share, and fixed
	global paramsConfig
	ifp = 0
	p = np.empty( (nComponent, len(paramsConfig)) )
	for c in range(nComponent):
		for i, v in enumerate(paramsConfig.values()):
			if v[1]=='free':
				p[c, i] = free_params[ifp]
				ifp += 1
			elif v[1]=='share':
				if c==0:
					p[c,i] = (free_params[ifp])
					ifp += 1
				else:
					#solid angle, shrink by two magnitude for each component```
					if i==len(paramsConfig)-1: p[c,i] = p[0,i]-2*c
					else: p[c,i] = p[0,i]
			else:
				p[c,i] = v[0][c]
	return p


def writeINP(density, Tkin, column, species, inpFile, outFile, z=0.0, dV=1.0):
	### write .inp file for radex
	with open(inpFile, 'w') as file:
		file.write(os.path.join(tmpPath, '%s.dat\n' % species)) # mol data
		file.write('%s\n' % outFile) # out file
		file.write('%s %s\n' % radexFreqRange) # freqency range
		file.write('%.3f\n' % Tkin)	# Tkin in K		
		file.write('2\n')	# collision partners: in most cases, using H2 as only partner is a good default
		file.write('o-H2\n')	# first collision partner
		file.write('%.3e\n' % (fortho*density))	# density in cm-3
		file.write('p-H2\n')	# first collision partner
		file.write('%.3e\n' % ((1-fortho)*density))	# density in cm-3
		file.write('%s\n' % (Tbg*(1+z))) # Tbg in K
		file.write('%.3e\n' % column) # column density in cm-2
		file.write('%s\n' % dV) # dV in km/s
		file.write('0\n') # run another?
		file.close()


def readOUT(outFile, dV=1.0):
	### read .out file from radex
	with open(outFile, 'r') as file:
		lines = file.readlines()
		file.close()
		### search for flux
		for i,l in enumerate(lines):
			if lines[i].split()[-1] == '(erg/cm2/s)': break
		### error ouput
		if i >= len(lines): return None
		J = []
		intFluxPerSr = []
		tau = []
		for l in lines[i+1:]:
			lsp = l.split()
			J.append( (int(lsp[0]), int(lsp[2])) )
			freq = float(lsp[4]) # GHz
			T_R = float(lsp[8]) # K
			'''
			T_R = 1*u.K
			dV = 1*kms
			freq = 1*u.GHz 
			intIntensity = np.sqrt(np.pi)/2/np.sqrt(np.log(2)) * T_R * dV
			F = 2 * con.k_B / con.c**2 * freq**2 * intIntensity
			print('f=', F.to(Jykms).value)
			'''
			result = T_R * dV * freq**2 * 32704.24161508085
			intFluxPerSr.append(result)
			tau.append(float(lsp[7]))

		return J, intFluxPerSr, tau


def runOneRadex(density, Tkin, column, species, z=0.0, dV=1.0):
	#run a radex model of one species and one component
	with tempfile.TemporaryDirectory(dir=tmpPath) as tmpdir:
		inpFile = os.path.join(tmpdir, '%s.inp' % species)	# temp input for radex
		outFile = os.path.join(tmpdir, '%s.out' % species)	# temp output for radex

		#print('oH2 %.3e pH2 %.3e T %f col %.3e %f %f %s' % (density*0.75, density*0.25, Tkin, column, dV, z, species))
		writeINP(density, Tkin, column, species, inpFile, outFile, z=z, dV=dV)
		try:
			subprocess.run('%s < %s > log' % (radexBin, inpFile), check=True, shell=True, stderr=subprocess.DEVNULL)
		except subprocess.CalledProcessError:
			print("RADEX failed with parameters:", density, Tkin, column, species)
			return None
		try:
			result = readOUT(outFile, dV=dV)
		except FileNotFoundError:
			print("RADEX failed with parameters:", density, Tkin, column, species)
			return None
		return result


'''
'p' in the form of [log Density, log Tkin, log N(species1), log N(species2), ... logOmega] * n component
'lines' in the form of {'species1': [(Jup1,Jlo1), (Jup2,Jlo2), ...], 'species2': [(Jup1,Jlo1), (Jup2,Jlo2), ...], ...}
'''
def modelMM(p, lines, nComponent=2, dvddr=1, **kws):
	# multiple molecules multiple components model
	global D_A, z, nfp, dV
	#muLensing = kws.get('muLensing', 1.0)
	intFlux = []
	tau = []
	#followingFlux = [] # the rest flux following the last detecte CO
	for c in range(nComponent):
		for i, (species, J) in enumerate(lines.items()):
			logn, logT, logX, logO = p[c, [0, 1, 2+i, -1]]
			'''
			D_A = cosmo.angular_diameter_distance(z).to(u.kpc).value
			Omega = (phySize / D_A)**2 * pi / 4 * muLensing
			N = n * phySize * X * dV
			logN = logn + logX + lgdV + 0.5 * (logOmega - log(pi/4) - logmuLensing) + logD_A
			'''
			#logN = logn + logX + lg(dV) + 0.5*(logO - lg((np.pi/4)) - lg(muLensing)) + lg(D_A) + 21.489350545222138#lg(kpc.to(cm))
			logN = logn + logX + lg(dV) - lg(dvddr) + lg(3.086e18)
			result = runOneRadex(10**logn, 10**logT, 10**logN, species, z=z, dV=dV)
			if result:
				JModel, intFluxPerSr, tauModel = result
				j_f = dict(zip(JModel, intFluxPerSr))
				intFlux += [j_f.get(j, np.nan) * 10**logO for j in J]
				j_t = dict(zip(JModel, tauModel))
				tau += [j_t.get(j, np.nan) for j in J]
			else:
				intFlux += [np.nan]*len(J)
				tau += [np.nan]*len(J)
	nObs = sum(len(js) for js in lines.values())
	intFlux = np.array(intFlux).reshape(nComponent, nObs)
	tau = np.array(tau).reshape(nComponent, nObs)
	intFlux /= 1+z

	#followingFlux = np.array(followingFlux).reshape(nComponent, len(followingFlux)//nComponent)
	return intFlux, tau#, followingFlux


def residual(fp, flux, eflux, nComponent=2, **kws):
	global lines
	p = fp2p(fp, nComponent)
	try:
		modelFlux, tau = modelMM(p, lines, nComponent=nComponent, **kws)
		### sum all component
		modelFlux = np.sum(modelFlux, axis=0)
	except ValueError:
		return np.inf
	# Some nan sometimes....
	if np.any(np.isnan(modelFlux)):
		return np.inf
	return (flux - modelFlux) / eflux


def chi2(fp, flux, eflux, nComponent=2, **kws):
	resi = residual(fp, flux, eflux, nComponent=nComponent, **kws)
	return (resi**2).sum()


def lnLike(fp, flux, eflux, limit, nComponent=2, **kws):
	""" likelihood function"""
	global lines
	p = fp2p(fp, nComponent)
	try:
		modelFlux, tau = modelMM(p, lines, nComponent=nComponent, **kws)
		### sum all component
		modelFlux = np.sum(modelFlux, axis=0)

		### Jul2: following flux should not larger than the last observed transition
		#followingFlux = np.sum(followingFlux, axis=0)
		#if np.any(followingFlux[1:] > followingFlux[0]): return -np.inf

	except ValueError:
		return -np.inf
	# have nan sometimes....
	if np.any(np.isnan(modelFlux)):
		return -np.inf

	ll = 0.0
	#get lglikelihood
	resi = flux - modelFlux
	#no error, so it is not obseved but physical limit
	mask = (limit=='u') & np.isnan(eflux)
	if (resi[mask]<0).any(): return -np.inf
	mask = (limit=='l') & np.isnan(eflux)
	if (resi[mask]>0).any(): return -np.inf
	# detected
	mask = limit=='d'
	if mask.any():
		ll += np.sum(-0.5 * (resi[mask]/eflux[mask])**2 + np.log(2*np.pi*eflux[limit=='d']**2))
	# upper limit
	mask = (limit=='u') & np.isfinite(eflux)
	if mask.any():
		ll += np.sum(norm.logcdf(resi[mask]/eflux[mask]))
	# lower limit
	mask = (limit=='l') & np.isfinite(eflux)
	if mask.any():
		ll += np.sum(norm.logsf(resi[mask]/eflux[mask]))
	return ll
	#return -0.5 * (np.sum(((flux - modelFlux)**2.0 / eflux**2.0)+np.log(eflux**2.0)))


def lnPrior(fp, bounds, nComponent=2):
	global nfp, paramsConfig, T_d

	# First Check boundaries
	if (np.any(fp > bounds[:, 1]) or np.any(fp < bounds[:, 0])):
		return -np.inf

	p = fp2p(fp, nComponent)
	# Temperature should be in ascending order
	T = p[:,1]
	if np.any(T[:-1] >= T[1:]):
		return -np.inf

	# Size shoud be in descending order
	# (colder component should be larger!)
	if paramsConfig['Omega'][1] == 'free':
		O = p[:,-1]
		if np.any(O[:-1] <= O[1:]):
			return -np.inf

	### Jul2: pressure should be in ascending order
	pressure = p[:,0] * p[:,1]
	if np.any(pressure[:-1] >= pressure[1:]):
		return -np.inf


	Tcold = 10**p[0, 1]
	return -0.5 * ((Tcold-T_d)/(1.0*T_d))**2
	# 9 < log10(N_CO/dv) - log10(n_H2) < 18
	#if np.any(parr[:,2:-1] - parr[:,0] >= 18) or np.any(parr[:,2:-1]- parr[:,0] <= 9):
	#	return -np.inf
	#if (np.any((p[2] - p[0]) >= 18.0) or np.any((p[2] - p[0]) <= 9.0) or
	#	np.any((p[6] - p[4]) >= 18.0) or np.any((p[6] - p[4]) <= 9.0)) :
	#	return -np.inf

	# Add bounds for tau, should be < 100
	# try:
	#     tau_max = model_lvg_tau(p, R)
	# except ValueError:
	#     return -np.inf
	# 
	# if (tau_max > 100):
	#     return -np.inf

	#return 0.


def lnProb(fp, flux, eflux, limit, nComponent=2, bounds=None, **kws):
	lp = lnPrior(fp, bounds, nComponent=nComponent)
	if not np.isfinite(lp):
		return -np.inf
	return lp + lnLike(fp, flux, eflux, limit, nComponent=nComponent, **kws)


def loadCatalog(filename):
	### read an ascii table to pandas frame
	df = pd.read_csv(
		filename,
		comment='#',       # Ignore lines starting with #
		skiprows=0,        # Skip the first two rows (so the third row is used as header)
		delim_whitespace=True,  # Assume the columns are separated by any whitespace
		na_values='NaN'    # Interpret 'NaN' strings as actual NaN values
	)

	### each 'value' column followed by 'err', so correct them.
	cols = []
	for i, name in enumerate(df.columns):
		if 'err' not in name:
			cols.append(name)
		else:
			cols.append('e'+df.columns[i-1])
	cols[0] = 'SOURCE'
	df.columns = cols

	### set source as row names
	df = df.set_index('SOURCE')
	return df


def decodeMolecule(text):
	match = re.match(r'^([A-Za-z0-9]+)\(\s*(\d+)-\s*(\d+)\)$', text)
	if match:
		species = match.group(1).lower()
		species = species.replace('+', 'p') # for hco+
		Jup = int(match.group(2))
		Jlo = int(match.group(3))
		return species, Jup, Jlo
	else: return None


def encodeMolecule(species, Jup, Jlo):
	return '%s(%i-%i)' % (species.upper(), Jup, Jlo)


def copyMoleculeDat(lines):
	### copy molecule dat to ram
	for s in lines.keys():
		datFile = '%s.dat' % s
		if not os.path.exists(os.path.join(tmpPath, datFile)):
			shutil.copy(os.path.join(molPath, datFile), tmpPath)


def getSource(catalog, source):
	row = catalog.T[source]
	print('Dealing with: %s\n' % source, row.T)
	z = row['z']
	T_d = row['T_dust']
	line_width = row['FWHM']

	# get lines, flux, eflux, skip if flux is nan
	lines = {}
	flux = []
	eflux = []
	limit = []
	for c in catalog.columns:
		# skip column of error
		if c.startswith('e'): continue
		mol = decodeMolecule(c)
		if mol:
			# upper/lower limit
			if type(row[c]) is str:
				if row[c].endswith('u'): limit.append('u') # upper limit
				elif row[c].endswith('l'): limit.append('l') # lower limit
				row[c] = float(row[c][:-1])
			elif np.isnan(row[c]):
				# skip nan
				continue
			else: limit.append('d') #detected

			flux.append(row[c])
			eflux.append(row['e'+c])
			# get lines
			species, Jup, Jlo = mol
			if species not in lines: lines[species] = [(Jup, Jlo),]
			else: lines[species].append((Jup, Jlo))

	flux = np.array(flux) * Jykms
	eflux = np.array(eflux) * Jykms
	limit = np.array(limit)
	return z, T_d, line_width, lines, flux, eflux, limit


def formatParams(fp, nComponent=2):
	### output a set of free params in table format
	global paramsConfig
	p = fp2p(fp, nComponent)

	columns = ['Comp']
	for k,v in paramsConfig.items():
		if v[1]=='free': columns.append('%10s(%s)' % (k, 'f'))
		elif v[1]=='share': columns.append('%10s(%s)' % (k, 's'))
		elif v[1]=='fix': columns.append('%10s(%s)' % (k, 'x'))
	print(''.join(columns))
	lines = [['%-4i' % (c+1)] for c in range(nComponent)]
	for c in range(nComponent):
		lines[c] += ['%13.3f' % v for v in p[c]]
		print(''.join(lines[c]))

	return columns, lines


def betterStep(steps, probs=None, nComponent=2, percentile=95, method='bestmean'):
	### choose a method you prefer to derive a better parameters
	global paramsConfig
	if method == 'median':
		fp = np.median(steps, axis=0)
	elif method == 'best':
		fp = steps[np.argmax(probs)]
	elif method == 'bestmean': 			
		idx = probs > np.percentile(probs, percentile)
		fp = np.mean(steps[idx], axis=0)

	print('A better step with %s is:' % method)
	formatParams(fp, nComponent)
	return fp


def loadPklFile(source, nComponent, suffix=''):
	with open(getOutPickleFilename(source, nComponent, suffix=suffix), 'rb') as pklFile:
		return pickle.load(pklFile)


def plotProb(source, nComponent=2, suffix=''):
	### plot the posterior probability to check convergence
	source, z, bounds, T_d, dV, nComponent, paramsConfig, (lines, flux, eflux, limit), (steps, lnprobs) = loadPklFile(source, nComponent, suffix=suffix)
	#with open(getOutPickleFilename(source, nComponent, suffix=suffix), 'rb') as pklFile:
	#	source, z, bounds, T_d, dV, nComponent, paramsConfig, (lines, flux, eflux), (steps, lnprobs) = pickle.load(pklFile)

	fig, ax = plt.subplots()
	ax.plot(lnprobs, '.')
	return fig, ax


def plotCorner(source, nComponent=2, bins=20, params=None, labels=None, showtop=0.5, suffix='', subset=None, **kws):
	### plot corners
	global paramsConfig, lines, dV, z

	def corners(steps, probs, bins=10, showtop=0,
		figure_kws = dict(figsize=(9,9)),
		hist_kws = dict(density=True, color = 'grey', histtype='step'),
		values = None,
		hist_value_kws = dict(pos='top', marker='v', color='tab:red'),
		hist2d_value_kws = dict(marker='+', color='tab:red', markersize=8, markeredgewidth=2),
		limits = [0.16, 0.84],
		hist_limit_kws = dict(linestyle='--', color='k', linewidth=1, alpha=0.5),
		labels = None,
		hist_label_kws = dict(fontsize=10),
		hist2d = True,
		hist2d_kws = dict(cmap='Blues'), 
		scatter = False,
		scatter_kws = dict(s=0.2, cmap='RdYlBu', alpha=0.3),
		contour = True,
		contour_smooth_sigma = 0,
		contour_kws = dict(colors='gray', levels=[0.2, 0.4, 0.6, 0.8], linewidths=1.5, alpha=0.8),
		output = 'test.out',
		):

		ndim = steps.shape[1]
		### plot corners
		def autoLimit(s):
			### extend according to s range
			ma, mi = np.nanmax(s), np.nanmin(s)
			return mi-(ma-mi)*0.1, ma+(ma-mi)*0.1
		def scientificNotationTitle(v, l, u):
			### show value, lower, and upper uncertainty in scientific notation format
			e = np.floor(np.log10(np.abs(v)))
			if e == 1: return '$%+.2f_{%+.2f}^{%+.2f}$' % (v, l, u)
			elif e == 0: return '$%+.3f_{%+.3f}^{%+.3f}$' % (v, l, u)
			else: return '$%+.3f_{%+.3f}^{%+.3f}$ x$10^{%i}$' % (v/10**e, l/10**e, u/10**e, e)

		### sort steps (or high prob steps might be put behind)
		steps = steps[np.argsort(probs)]
		probs = np.sort(probs)		
		### only plot the top% steps, clip low prob steps.
		topidx = probs > np.quantile(probs, showtop)

		### get a square figure.
		fig, ax = plt.subplots(nrows=ndim, ncols=ndim, **figure_kws)
		plt.subplots_adjust(bottom=0.09, top=0.91, left=0.09, right=0.91, wspace=0.05, hspace=0.05)

		### decide where to put marker in hist
		valuePos = hist_value_kws.pop('pos', 'bartop')
		contourLevels = np.array(contour_kws.pop('levels', [0.2, 0.4, 0.6, 0.8]))

		### output values and limits
		f = open(output, 'w')
		for i in range(ndim):
			### hist at [i,i]
			h,xe,_ = ax[i,i].hist(steps[:,i][topidx], bins=bins, **hist_kws)

			### values
			if values is not None:
				if valuePos == 'bartop': 
					upperY = h[np.searchsorted(xe, values[i], side='right')-1] * 1.1
				elif valuePos == 'top':
					upperY = np.max(h)
			ax[i,i].plot(values[i], upperY, **hist_value_kws)

			### limits
			lim = np.quantile(steps[:,i], limits)
			for v in lim:
				ax[i,i].axvline(v, **hist_limit_kws)

			### labels
			if labels is not None:
				text = '%s=%s' % (labels[i], scientificNotationTitle(values[i], min(lim)-values[i], max(lim)-values[i]))
				ax[i,i].text(0, 1.02, text, transform=ax[i,i].transAxes, ha='left', va='bottom', **hist_label_kws)

			### output
			f.write('%s %f %f %f\n' % (labels[i], values[i], min(lim)-values[i], max(lim)-values[i]))

			### hide x ticklabels except the bottom row
			if i<ndim-1: ax[i,i].xaxis.set_ticklabels([])
			### hide all y tickslables
			ax[i,i].set_yticklabels([])

			### ticks
			ax[i,i].set_autoscale_on(False)
			ax[i,i].minorticks_on()
			ax[i,i].tick_params(top=True, left=False, direction='in', labelrotation=45, length=5)
			ax[i,i].tick_params(which='minor', top=True, left=False, direction='in', length=2)
			ax[i,i].set_xlim(autoLimit(steps[:,i][topidx]))

			### plot hist2d
			for j in range(i):
				h, xe, ye = np.histogram2d(steps[:,j][topidx], steps[:,i][topidx], bins=bins)
				### hist2d
				if hist2d:
					ax[i,j].imshow(h.T, origin='lower', extent=[xe[0], xe[-1], ye[0], ye[-1]], aspect='auto', **hist2d_kws)
				### scatter
				if scatter:
					ax[i,j].scatter(steps[:,j][topidx], steps[:,i][topidx], c=probs[topidx], **scatter_kws)
				### contour
				if contour:
					if contour_smooth_sigma is not None:
						h = gaussian_filter(h, sigma=contour_smooth_sigma)
					xc = (xe[1:]+xe[:-1])/2
					yc = (ye[1:]+ye[:-1])/2
					ax[i,j].contour(xc, yc, h.T, levels=contourLevels*np.max(h), **contour_kws)

				### values
				ax[i,j].plot(values[j], values[i], **hist2d_value_kws)
				
				### hide x ticklabels except the bottom row
				if i<ndim-1: ax[i,j].xaxis.set_ticklabels([])
				### hide y ticklabels except the left column
				if j>0: ax[i,j].yaxis.set_ticklabels([])

				### ticks
				ax[i,i].set_autoscale_on(False)
				ax[i,j].minorticks_on()
				ax[i,j].set_yticks(ax[i,i].get_xticks()) ### make sure ytick is the same as the [i,i] xtick
				ax[i,j].tick_params(top=True, right=True, direction='in', labelrotation=45, length=5)
				ax[i,j].tick_params(which='minor', top=True, right=True, direction='in', length=2)
				ax[i,j].set_xlim(autoLimit(steps[:,j][topidx]))
				ax[i,j].set_ylim(autoLimit(steps[:,i][topidx]))

				### hide upper right
				ax[j,i].axes.set_axis_off()
		f.close()
		return fig, ax


	source, z, bounds, T_d, dV, nComponent, paramsConfig, (lines, flux, eflux, limit), (steps, lnprobs) = loadPklFile(source, nComponent, suffix=suffix)
	### get label
	if params is None: 
		copyMoleculeDat(lines)
		params = betterStep(steps, lnprobs, nComponent=nComponent)
		print('chi2=', chi2(params, flux.value, eflux.value, nComponent=nComponent))

	if labels is None: 
		labels = []
		for c in range(nComponent):
			for v in paramsConfig.values():
				if v[1] == 'free':
					labels.append('%s_%i' % (v[-1],c+1))
				if v[1] == 'share':
					if c==0: labels.append('%s_%i' % (v[-1],c+1))

	if subset is not None:
		steps = steps[:, subset]
		params = params[subset]
		labels = np.array(labels)[subset]

	### filter
	fig, ax = corners(steps, lnprobs, bins=bins, values=params, labels=labels, showtop=showtop, \
		output=getOutPickleFilename(source, nComponent=nComponent, suffix=suffix)+'.out')

	fig.savefig(os.path.join(figPath, '%s_corner%s.png' % (source, suffix)), format='png', bbox_inches='tight', dpi=400)
	return fig, ax


def plotSLED(source, nComponent=2, plotLines={'co':[(j,j-1) for j in range(1,11)], 'ci':[(j,j-1) for j in range(1,3)]}, 
	suffix='', params=None, speciesSep = 1.5, obsColor='k', modelColor = '#DD8822', componentColors = ['#6699EE', '#EE8888'], 
	bestN=10, bestNColor='#FFAA44', otherModel={}, otherColor=['green'], separateComponents=True, **kws):

	### plot SLED of several species
	global z, dV, paramsConfig

	source, z, bounds, T_d, dV, nComponent, paramsConfig, (lines, flux, eflux, limit), (steps, lnprobs) = loadPklFile(source, nComponent, suffix=suffix)
	#with open(getOutPickleFilename(source, nComponent, suffix=suffix), 'rb') as pklFile:
	#	source, z, bounds, T_d, dV, nComponent, paramsConfig, (lines, flux, eflux), (steps, lnprobs) = pickle.load(pklFile)

	### copy mol dat to RAM
	copyMoleculeDat(lines)
	### modify the form of lines, fill flux/eflux/limit to the dict
	def fillFluxToLines(lines, *args):
		### modify lines to [Jup, flux, eflux, ...]
		i = 0
		filledLines = lines.copy()
		for species, J in filledLines.items():
			filledLines[species] = [J,]
			for a in args:
				filledLines[species].append(a[...,i:i+len(J)])
			i += len(J)
		return filledLines
	### smooth sled to curve
	def interp(x, y):
		if len(x)>2:
			f = interp1d(x, y, kind='cubic')
			xa = np.linspace(x[0], x[-1], int((x[-1]-x[0])*50))
			return xa, f(xa)
		else: return x, y

	# fill obs
	obsLines = fillFluxToLines(lines, flux, eflux, limit)

	# fill model
	if params is None: params = betterStep(steps, lnprobs, nComponent=nComponent)
	modelFlux, modelTau = modelMM(fp2p(params, nComponent), plotLines, nComponent = nComponent)
	modelLines = fillFluxToLines(plotLines, modelFlux)

	def readModelOut(file):
		try:
			filledLines = {}
			f = open(file, 'r')
			for l in f:
				linename, *flux, sumflux = [s for s in l.split(' ') if len(s)>0]
				species, Jup, Jlo = decodeMolecule(linename)
				flux = [float(v) for v in flux]
				sumflux = float(sumflux)
				if species in filledLines:
					filledLines[species][0].append((Jup, Jlo))
					filledLines[species][1].append(flux)
				else: filledLines[species] = [[(Jup, Jlo),],[flux,]]
			f.close()
			for s in filledLines:
				filledLines[s][1] = np.array(filledLines[s][1]).T
			return filledLines
		except:
			return None
	### another sled
	otherModelLines = {label:readModelOut(file) for label, file in otherModel.items()}

	if bestN>0:
		### plot the best n results
		idx = np.argsort(lnprobs)[lnprobs.size//2::400]
		bestFlux = []
		bestLines = plotLines.copy()
		for p in steps[idx]:
			modelFlux, modelTau = modelMM(fp2p(p, nComponent), plotLines, nComponent = nComponent)
			#followingFlux = np.sum(followingFlux, axis=0)
			#if np.any(followingFlux[1:]>followingFlux[0]*0.95): print('Raised tail after last CO:\n', followingFlux)
			bestFlux.append(modelFlux)
		bestLines = fillFluxToLines(bestLines, *bestFlux)


	### plot
	fig, ax = plt.subplots(figsize=(8,6))
	plt.subplots_adjust(bottom=0.15, left=0.1, right=0.95)

	xShift = 0
	xticks = []
	xticklabels = []
	modelOut = []
	for i, (species, J) in enumerate(plotLines.items()):
		xPos = np.array(J)[:,0] + xShift - 1

		### source donot have this species detection
		if species in obsLines:
			### observe
			_x, _y, _e, _l = obsLines[species]
			_x = np.array(_x)[:,0] + xShift - 1 #Jup

			mask = _l == 'd'
			if mask.any():
				ax.errorbar(_x[mask], _y[mask].value, _e[mask].value, fmt='o', color=obsColor, capsize=0, label='$\mathrm{data}$' if i==0 else None)
			mask = (_l == 'u') & np.isfinite(_e)
			if mask.any():
				ax.errorbar(_x[mask], _y[mask].value, _e[mask].value, fmt='_', uplims=True, color=obsColor, capsize=0)
			mask = (_l == 'l') & np.isfinite(_e)
			if mask.any():
				ax.errorbar(_x[mask], _y[mask].value, _e[mask].value, fmt='_', lolims=True, color=obsColor, capsize=0)

			### model
			_x, _y = modelLines[species]
			_x = np.array(_x)[:,0] + xShift - 1

			modelOut.append(_y.copy())
			print('       '+' '.join(['  %5i  ' % v for v in _x]))

			### each component
			if (nComponent>1) and separateComponents:
				for c in range(nComponent):
					if c ==0: label = 'low-excitation comp.'
					else: labels='high-excitation comp.'
					ax.plot(*interp(_x, _y[c]), color=componentColors[c], lw=2, label=label if i==0 else None)
					print('Comp %1s:' % (c+1) + ' '.join(['%9.2e' % v for v in _y[c]]),)
			### sum all
			_ys = np.sum(_y, axis=0)
			label = '%1i-component model' % nComponent
			ax.plot(*interp(_x, _ys), color=modelColor, lw=2, label=label if i==0 else None)
			print('All   :' + ' '.join(['%9.2e' % v for v in _ys]))

			if bestN>0:
				for b in bestLines[species][1:]:
					#ax.plot(*interp(_x, b[0]), alpha=0.1, lw=1, color=componentColors[0])
					#ax.plot(*interp(_x, b[1]), alpha=0.1, lw=1, color=componentColors[1])
					ax.plot(*interp(_x, np.sum(b, axis=0)), alpha=0.15, lw=2, color=bestNColor, zorder=0)

			for o, (label, anotherModelLines) in enumerate(otherModelLines.items()):
				_x, _y = anotherModelLines[species]
				_x = np.array(_x)[:,0] + xShift - 1
				### sum all
				_ys = np.sum(_y, axis=0)
				ax.plot(*interp(_x, _ys), color=otherColor[o], lw=2, label=label if i==0 else None)

		### gather ticks and labels
		xticks += list(xPos)
		xticklabels += ['%s(%i-%i)' % (species.upper(), *j) for j in J]

		### move to the next species
		xShift = xPos[-1] + speciesSep

	ax.legend(loc='upper right', fontsize=10)

	# set tick
	ax.set_xticks(xticks)
	ax.set_xticklabels(xticklabels, rotation=45, fontsize=12, fontweight='bold')
	ax.set_yticklabels(ax.get_yticklabels(), fontsize=12, fontweight='bold')
	ax.set_ylabel('$\mathbf{S_J}$ (Jy km s$\mathbf{^{-1}}$)', fontsize=15, fontweight='bold')

	fig.savefig(os.path.join(figPath, '%s_SLED%s.png' % (source, suffix)), format='png', bbox_inches='tight', dpi=400)

	### output model flux to dat
	outfile = open(os.path.join(outPath, '%s_%1icomp_model_flux.dat' % (source, nComponent)), 'w')
	modelOut = np.hstack(modelOut)
	for i in range(len(xticklabels)):
		outfile.write('%-12s' % xticklabels[i] + ' '.join(['%10.3e' % v for v in modelOut[:,i]]) + '%10.3e' % (modelOut[:,i].sum()) + '\n')
	outfile.close()
	return fig, ax



def	initConfig(nComponent, phySize=6, muLensing=1.0):
	### modify paramsConfig according to the properties of source
	global z, D_A, paramsConfig, T_d

	### modify Tkin bound
	TkinBound = paramsConfig['Tkin'][2]
	paramsConfig['Tkin'][2] = (max(lg(Tbg*(1+z)), lg(T_d*0.5)), TkinBound[1])	### Tkin > max(Tbg*(1+z), Tdust)

	### repeat X for components
	for s in list(paramsConfig.keys())[2:-1]:
		if paramsConfig[s][1] == 'fix':
			paramsConfig[s][0] = [paramsConfig[s][0]]*nComponent

	### modify Omega
	logOmega = lg((phySize / D_A)**2 * np.pi / 4 * muLensing)
	print('Physical size of %s kpc corresponding to solid angle of %f sr' % (phySize, logOmega))
	if paramsConfig['Omega'][1] in ['free', 'share']:
		### reset bound
		omegaBound = paramsConfig['Omega'][2]
		paramsConfig['Omega'][2] = (logOmega+omegaBound[0], logOmega+omegaBound[1])
	else:
		### shrink warmer component by one magnitude
		paramsConfig['Omega'][0] = [logOmega-c*2 for c in range(nComponent)]


def initBounds(nComponent):
	### bounds of free params draw from ([log Density, log Tkin, log X(species1), log X(species2), ... logSize] * n component)
	global paramsConfig

	bounds = []
	for c in range(nComponent):
		for v in paramsConfig.values():
			if v[1]=='free': bounds.append(v[2])
			elif v[1]=='share':
				if c==0: bounds.append(v[2])

	return np.array(bounds)


def initParams(nComponent, bounds, nwalkers, flux, eflux, limit):
	### initialize params randomly
	global line, nfp, paramsConfig

	p0 = []
	while len(p0) < nwalkers:
		#uniform normalized p
		normFp = np.random.uniform(0, 1, nfp)
		fp = normFp * (bounds[:,1] - bounds[:,0]) + bounds[:, 0]

		if np.isfinite(lnProb(fp, flux, eflux, limit, nComponent=nComponent, bounds=bounds)):
			p0.append(fp)
		'''
		### sort T
		if nComponent>=2:
			normP[:, 1] = np.sort(normP[:, 1], axis=0)
		### sort S
		if paramsConfig['Omega'][1] == 'free':
			normP[:, -1] = np.sort(normP[:, -1], axis=0)[::-1]

		normP = normP.flatten()

		### test p if it is valid for radex
		try:
			intFlux, tau = modelMM(scaleP, lines, nComponent=nComponent)
			if np.isfinite(intFlux).all():
				print(p0)
				p0.append(scaleP)
			### do a minimize to find nearby local minimum, but this is very slow
			#res = minimize(chi2, scaleP, bounds=bounds, method='Powell')
			#p0.append(res.x)
		except:
			continue
		'''
	return p0


def runMCMC(catname, source,
	nComponent = 2,
	nwalkers = 32,
	nburn = 200,
	niter = 1000,
	suffix = '',
	phySize = 7):
	### run mcmcm
	global z, lines, D_A, T_d, nfp, dV, paramsConfig

	cat = loadCatalog(catname)
	z, T_d, line_width, lines, flux, eflux, limit = getSource(cat, source)
	D_A = cosmo.angular_diameter_distance(z).to(u.kpc).value

	dV = line_width	# do not care actually, fully degenerate with the column density

	# Line : Jy km/s
	# 1 Jansky=1 Jy = 10e-26 W m-2 Hz-1 = 1e-23 erg s-1 cm-2 Hz-1
	# brightness erg s-1 cm-2 Hz-1 sr-1 ~ W m-2 Hz-1 sr-1 ~ Jy sr-1

	# Radex output : erg s-1 cm-2  = 1e3 W m-2
	# pyradex output : erg s-1 cm-2 Hz-1 sr-1

	# K*km/s is an integrated surface brightness, while the value
	# in erg/s/cm2 is an isotopic flux emitted in all directions.

	# Specicic Intensity  I(nu)  W m-2 Hz-1 sr-1
	# Flux density W m-2 Hz-1 ~ 10e26 Jy : S = I * Surf./d**2

	# -> Specific (velocity integrated ) Intensity
	# I = S_nu * d**2 / surf

	# Size:
	# R_source = D_A * sqrt(size/pi)
	# This is uncorrected for lensing magification

	# H2 density : Number density of collision partners : Unit: cm-3.  ; Allowed range: 10(-3) - 10(13)
	# Column density : Unit: cm-2 : Allowed range: 10(5) - 10(25)
	# Kinetic temperature : Unit: K. Allowed range: 0.1 - 10,000
	# linewidth : Thermal + turbulent width (FWHM) of the lines.; Unit: km s-1. ; Allowed range: 10(-3) - 10(3)
	# Size : Unit sr

	### number of free parameters
	nfp = 0
	for v in paramsConfig.values():
		if v[1]=='free': nfp += nComponent
		elif v[1]=='share': nfp += 1

	### copy mol dat to RAM
	copyMoleculeDat(lines)

	### modify paramsConfig, get bounds and p0
	initConfig(nComponent, phySize=phySize, muLensing=1.0)
	bounds = initBounds(nComponent)
	p0 = initParams(nComponent, bounds, nwalkers, flux.value, eflux.value, limit)

	### check chi2 for p0
	for fp in p0:
		formatParams(fp, nComponent)
		print('chi2=', chi2(fp, flux.value, eflux.value, nComponent=nComponent))

	ndim = nfp

	# Multithread
	with multiprocessing.Pool() as pool:
		sampler = emcee.EnsembleSampler(nwalkers, ndim, lnProb,
										args=(flux.value, eflux.value, limit),
										kwargs={'bounds': bounds, 'nComponent': nComponent},
										threads=ncpu,
										moves=emcee.moves.StretchMove(0.5))

		### run first burn-in
		logger.info("1st Burn-in")
		pos, prob, state = sampler.run_mcmc(p0, nburn, progress=True)
		#sampler.reset()

		### discard stuck walkers, but there is still risk that the best walker falls into some local minimum
		probs = sampler.flatlnprobability
		steps = sampler.flatchain
		maxProbArg = np.argmax(probs)
		p0 = steps[maxProbArg]
		logger.info("Burn-in again from %s" % p0)

		### re-initialize
		p0 = p0 + 1e-4 * np.random.randn(nwalkers, ndim)

		### burn-in again
		sampler.reset()
		logger.info("2nd Burn-in")
		pos, prob, state = sampler.run_mcmc(p0, nburn, progress=True)
		sampler.reset()

		### begin running
		logger.info("Running production...")
		pos, prob, state = sampler.run_mcmc(pos, niter, progress=True)

		#pemcee = np.percentile(sampler.flatchain, [50], axis=0)[0]

	#steps, lnprobs = sampler.chain, sampler.lnprobability
	steps = sampler.flatchain
	lnprobs = sampler.flatlnprobability

	### output
	if not os.path.exists(outPath):
		os.makedirs(outPath)
	with open(getOutPickleFilename(source, nComponent, suffix=suffix), 'wb') as pklFile:
		pickle.dump((source, z, bounds, T_d, dV, nComponent, paramsConfig, (lines, flux, eflux, limit), (steps, lnprobs)), pklFile)

	logger.info('%s done' % source)
	return



def getOutPickleFilename(source, nComponent=2, suffix=''):
	### file of fitting result
	return os.path.join(outPath, '%s_%1icomp%s.pickle' % (source, nComponent, suffix))




if __name__ == '__main__':
	'''
	### Catalog Format
	SOURCE   z       D_L     T_dust  FWHM   CO(1-0)  err      CO(2-1)  err     CO(3-2)  err    CO(4-3)  err     CO(5-4)  err     CO(6-5)  err     CO(7-6)  err     CO(8-7)  err     CI(1-0)  err     CI(2-1)  err
	1338     1.234   8801    30      89.8   NaN      NaN      0.278    0.024   NaN      NaN    0.521    0.038   0.537    0.081   NaN      NaN     0.352    0.067   NaN      NaN     0.267    0.038   0.328    0.067
	2299     1.395   10238   30      143.4  0.165    0.029    0.636    0.039   NaN      NaN    2.710    0.028   2.783    0.085   NaN      NaN     1.860    0.034   1.2u     0.1     1.186    0.031   2.006    0.042
	......

	Value ends with:
	'u' is upper limit
		if its 'err' is NaN: value is not obseved, but model should below it (more like a physical limit), and will not show in SLED plot.
		if not: value is obseved but below a threshold, and will be plotted as small bar.
	'l' is lower limit
		if its 'err' is NaN: value is not obseved, but model should above it (more like a physical limit), and will not show in SLED plot.
		if not: value is obseved but above a threshold, and will be plotted as small bar.

	### Parameter Configure
	The keys of paramsConfig should be in the order:
		'nH2', 'Tkin', 'Xmol1', 'Xmol2', ..., 'Omega', 
	though the values are base-10 logarithm.

	Their values are in order of ( initial, type, bounds, label )
		intial: only valid if it is a fixed parameter
			intial value of Omega is derived by assuming a physical size of source at its D_A
		type:
			'free': free parameter for all components;
			'share': free but share the same parameter among components;
			'fix': fixed parameter for all components.
		bounds: only valid if it is a free/share parameter
			actual bounds of Omega is its initial + bounds
		label: title of corners

	### Usage
	# calculate one source like this:
		catname = '../data/flux_for2p.dat'
		source = '2299'
		runMCMC(catname, source, nComponent=nComponent, suffix=suffix, nburn=200, niter=800, phySize=3)
		plotProb(source, nComponent = nComponent, suffix=suffix)
		plotCorner(source, nComponent = nComponent, suffix=suffix)
		plotSLED(source, nComponent = nComponent, suffix=suffix)

	# or iterate the whole table like this:
		catname = '../data/flux_for2p.dat'
		cat = loadCatalog(catname)
		for source,_ in cat.iterrows():
			runMCMC(catname, source, nComponent=nComponent, suffix=suffix, nburn=200, niter=800, phySize=3)
			plotProb(source, nComponent = nComponent, suffix=suffix)
			plotCorner(source, nComponent = nComponent, suffix=suffix)
			plotSLED(source, nComponent = nComponent, suffix=suffix)

	### cf. flux_for2p.dat for table format
	'''

	# name of catalog
	catname = './flux_for2p.dat'
	# source name in catalog
	source = '7395'
	# number of components
	nComponent = 2 # cold and warm
	# suffix of output files to distinguish each session
	#suffix = '_5fp_Tdust_zcorr_limit'
	suffix = '_6fp_fixCI_Tdust_zcorr_limit'

	# default bounds of molecular abundance
	XB = {'co': [8e-5, 2e-4], 'hcn': [1e-9, 5e-7], 'cs': [1e-9, 5e-7], 'hcop': [1e-9, 5e-7], 'ci': [1e-5, 2e-3]}
	# parameter config
	paramsConfig = dict(
		nH2  = [[4.781, 3.337],	'free',	(1.5, 6.0),	r'$lg n_{H_2}$'],	#density, MUST be 'free'
		Tkin = [[0.879, 2.222],	'free',	(lg(Tbg), 2.8),	r'$lg T_{Kin}$'],	#Tkin, MUST be 'free'
		Xco  = [-4.10,	'fix',	lg(XB['co']),	r'$lg X_{CO}$'],		#abundance
		Xci  = [-4.52,	'share',lg(XB['ci']),	r'$lg X_{CI}$'],		#abundance
		Omega= [None,	'free',	(-5, 5), 			r'$lg \Omega$']		#solid angle, if share, shrink by two magnitude for each component
	)


	runMCMC(catname, source, nComponent=nComponent, suffix=suffix, nburn=50, niter=100, phySize=3)
	plotProb(source, nComponent = nComponent, suffix=suffix)
	plotCorner(source, nComponent = nComponent, suffix=suffix, showtop=0.5)
	# draw sled of a model
	plotSLED(source, nComponent = nComponent, suffix=suffix, bestN=10)
	# draw several model together
	plotSLED(source, nComponent = nComponent, suffix=suffix, bestN=0, separateComponents=True, otherModel={'1-component model':'out/%s_1comp_model_flux.dat' % source})
	plt.show()


# See :
# - https://arxiv.org/pdf/1602.01095.pdf
# - https://arxiv.org/pdf/1401.2998.pdf
# - https://arxiv.org/pdf/1212.5955.pdf
# - https://arxiv.org/abs/0809.2337
