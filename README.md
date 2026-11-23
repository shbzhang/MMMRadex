# MMMRadex
Multi-species Multi-components MCMC Radex fitting

cf. [paper](link)

# System requirements:
- [radex](https://sronpersonalpages.nl/~vdtak/radex/index.shtml#details) 
- The code has been tested and verified to run on Python 3.8, with key dependencies:
  - **numpy** 1.23.5
  - **emcee** 3.1.2
  - **matplotlib** 3.7.5

# Usage:

### Runing with RAM
create a RAM disk in shell
```
diskutil erasevolume HFS+ RAMDisk $(hdiutil attach -nomount ram://1048576)
```
Error occurs if only plotting without runMCMC fitting, try to copy all data files in molPath to tmpPath.

Remember to unmount it after you finish the fitting work
```
diskutil unmount /Volumes/RAMDisk
```

### Configure radex and path
modify parameter for radex running  in `MMMRadex.py`
```
# path of radex bin file
radexBin = '/Users/Radex/src/radex'
# frequency range for Radex to calculate lines
radexFreqRange = (1, 1800) # GHz
# default bounds of molecular abundance
XB = {'co': [8e-5, 2e-4], 'hcn': [1e-9, 5e-7], 'cs': [1e-9, 5e-7], 'hcop': [1e-9, 5e-7], 'ci': [1e-5, 2e-3]}
```
also modify the path to input and ouput files
```
# path to ram disk
tmpPath = '/Volumes/RAMDisk'
# path to put molecule data
molPath = './radex_moldata'
# path to put fitting results
outPath = './out'
# path to put output figures
figPath = './fig'
```

### Catalog
flux of sources are stored in a catalog, each flux is followed by an error column:
```
SOURCE/LINES(Jy_km/s)     z       D_L     T_dust  FWHM   CO(1-0)  err      CO(2-1)  err     CO(3-2)  err    CO(4-3)  err     CO(5-4)  err     CO(6-5)  err     CO(7-6)  err     CO(8-7)  err     CI(1-0)  err     CI(2-1)  err
2299                      1.395   9973    30      151  0.467    0.130    1.101    0.064   NaN      NaN    2.718    0.030   2.786    0.094   NaN      NaN     1.932    0.041   NaN      NaN     1.228    0.035   1.837    0.041
42925                     1.603   11839   30      501  0.237u   0.079    0.659    0.109   NaN      NaN    1.229    0.045   1.254   0.148    NaN      NaN     1.467    0.072   NaN      NaN     0.576    0.055   1.036    0.071
......
```

- Flux value ends with `u` is upper limit, while value ends with `l` is lower limit.
- a valid error should be provided after the upper/lower flux limit to calculate likelihood.
- if the following error is NaN, the flux is not an observed result, but a constraint to the model (more like a physical limit). This kind of flux will NOT show in the SLED plot.

### Custom fitting
modify parameter  in `MMMRadex.py` to custom your fitting process
```
# number of components
nComponent = 2 # cold and warm
# suffix of output files to distinguish each session
suffix = '_6fp'
```

The dict `paramsConfig` control the number of free parameters. Its keys should be in the order:
```
'nH2', 'Tkin', 'Xmol1', 'Xmol2', ..., 'Omega'
```
though the values are base-10 logarithm.

The value of each key should be a tuple in order of `( initial, type, bounds, label )`
- **initial**: the value of a parameter, only valid if type is `fixed`.
	the initial guess of `free` or `share` parameters are randomly got within bounds.
- **type**: whether a parameter is free.
	- `free`: free parameter for all components;
	- `share`: free but share the same parameter among components;
	- `fix`: fixed parameter for all components.
- **bounds**: only valid if it is a `free` or `share` parameter
- **label**: title of corners

NOTE: the Initial value of Omega is derived by assuming a physical size of source at its D_A. So the final bounds of Omega is its `initial + bounds`. If its bound is `share`, the value is for the first component, and shrink by two magnitude for each later component

An example of `paramsConfig`:
```
paramsConfig = dict(
	# density, MUST be 'free'
	nH2  = [[4.781, 3.337],	'free',	(1.5, 6.0),	r'$lg n_{H_2}$'],
	# Tkin, MUST be 'free'	
	Tkin = [[0.879, 2.222],	'free',	(lg(Tbg), 2.8),	r'$lg T_{Kin}$'],
	# CO abundance
	Xco  = [-4.10,	'fix',	lg(XB['co']),	r'$lg X_{CO}$'],
	# CI abundance
	Xci  = [-4.52,	'share',lg(XB['ci']),	r'$lg X_{CI}$'],
	# solid angle, if share, shrink by two magnitude for each component
	Omega= [None,	'free',	(-5, 5), 		r'$lg \Omega$']
	)
```

### Start a fitting session
calculate one source like this:
```
# specify a catalog
catname = 'example.dat'
# specify a source
source = '2299'
# run
runMCMC(catname, source, nComponent=nComponent, suffix=suffix, nburn=200, niter=800, phySize=3)
```

or iterate the whole table like this:
```
# specify a catalog
catname = 'example.dat'
# load catalog and run
cat = loadCatalog(catname)
for source,_ in cat.iterrows():
	runMCMC(catname, source, nComponent=nComponent, suffix=suffix, nburn=200, niter=800, phySize=3)
```

### Output figures
Three types of figures will be output to `figPath`:
- probability of steps to check convergence
```
plotProb(source, nComponent = nComponent, suffix=suffix)
```
- corner map
```
plotCorner(source, nComponent = nComponent, suffix=suffix)
```
- SLED
```
plotSLED(source, nComponent = nComponent, suffix=suffix)
```

# Installation:
```
git clone https://github.com/shbzhang/MMMRadex
```
