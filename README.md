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

#### Runing with RAM
create a RAM disk in shell
```
diskutil erasevolume HFS+ RAMDisk $(hdiutil attach -nomount ram://1048576)
```
Error occurs if only plotting without runMCMC fitting, try to copy all data files in molPath to tmpPath.

Remember to unmount it after you finish the fitting work
```
diskutil unmount /Volumes/RAMDisk
```

#### Configure radex and path
modify parameter for radex running  in `configure.py`
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

#### Catalog
flux of sources are stored in a catalog, each flux is followed by an error column:
```
SOURCE   z       D_L     T_dust  FWHM   CO(1-0)  err      CO(2-1)  err     CO(3-2)  err    CO(4-3)  err     CO(5-4)  err     CO(6-5)  err     CO(7-6)  err     CO(8-7)  err     CI(1-0)  err     CI(2-1)  err
1338     1.234   8801    30      89.8   NaN      NaN      0.278    0.024   NaN      NaN    0.521    0.038   0.537    0.081   NaN      NaN     0.352    0.067   NaN      NaN     0.267    0.038   0.328    0.067
2299     1.395   10238   30      143.4  0.165    0.029    0.636    0.039   NaN      NaN    2.710    0.028   2.783    0.085   NaN      NaN     1.860    0.034   1.2u     0.1     1.186    0.031   2.006    0.042
......
```

- Flux value ends with `u` is upper limit, while value ends with `l` is lower limit.
- a valid error should be provided after the upper/lower flux limit to calculate likelihood.
- if the following error is NaN, the flux is not an observed result, but a constraint to the model (more like a physical limit). This kind of flux will NOT show in the SLED plot.

#### Custom fitting
modify parameter  in `configure.py` to custom your fitting process
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
paramsConfig = {} 
# volume density, MUST be 'free'
paramsConfig['nH2']  = [[4.781, 3.337],	'free',	(1.5, 6.0), r'$lg n_{H_2}$']
# kinematic temperature, MUST be 'free'
paramsConfig['Tkin'] = [[0.879, 2.222],	'free',	(np.log10(Tbg), 2.8), r'$lg T_{Kin}$']
# CO abundance
paramsConfig['Xco']  = [-4.10, 'fix', np.log10(XB['co']), r'$lg X_{CO}$']
# CI abundance
paramsConfig['Xci']  = [-4.52, 'share', np.log10(XB['ci']),	r'$lg X_{CI}$']
# solid angle
paramsConfig['Omega']= [None, 'free', (-5, 5), r'$lg \Omega$']
```

#### Start a fitting session
calculate one source like this:
```
# specify a catalog
catname = '../data/flux_for2p.dat'
# specify a source
source = '2299'
# run
runMCMC(catname, source, nComponent=nComponent, suffix=suffix, nburn=200, niter=800, phySize=3)
```

or iterate the whole table like this:
```
# specify a catalog
catname = '../data/flux_for2p.dat'
# load catalog and run
cat = loadCatalog(catname)
for source,_ in cat.iterrows():
	runMCMC(catname, source, nComponent=nComponent, suffix=suffix, nburn=200, niter=800, phySize=3)
```

#### Output figures
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
