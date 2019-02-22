# FieldGenerator
--------------
Module to create a FIELD and CONFIG file for [DL_POLY Classic](https://ccpforge.cse.rl.ac.uk/gf/project/dl_poly_classic/).
The Class-Structure lets you run optimizations faster since single parts might be updated without
recalculating everything else. It also protects you from having to copy-paste hundreds of values for your force-field.

by Patrick Melix
2017/04

Limitations so far (that I know of):
 - no support for repeat-counter, frozen atom or igrp (see [DL_POLY Manual](https://ccpforge.cse.rl.ac.uk/gf/project/dl_poly_classic/docman/))
 - no support for crossterms, pmf, constraints, shell, teth
 - no support for tbp, fbp, external field, tersoff, metal potential
 - only lj and buck potentials for VDW are supported

The code contains a part for the creation of a LFMM-Parameter section (see [doi:10.1021/ic501519ah](http://dx.doi.org/10.1021/ic501519a) for explenation). The required modifications of the DL_POLY classic code are not available publicly but can be requested from the authors of the publication. The code works fine without any LFMM-Parameters and is therefore also compatible with the official DL_POLY Classic.

## Installation:
 - Run `pip3 install .` inside the main directory (use `pip3 install --user .`) if you do not have root permissions.
 - If you don't have `pip3`, add the *fieldgeneartor* folder to your `$PYTHONPATH`
 - If you want a developer installation (symlinked files) use `pip3 install --user -e .`

## Test:
- Go to the *./test/* folder and execute `./test_generator.py`
- You should see a FIELD and a CONFIG file beeing created.

## Usage:
   - Load the module and create an Instance of *FieldCreator()*, then call `.create()`
   - Change the Variables defined in `__init()__` to affect the procedures.
   - If you don't want to use `.create()` make sure you call the public functions in the right order, check `.create()`
     to see the way this module is intended to be used.
   - If you want to change some values use `.updateFF()` or `.updateLFSE()` to get an updated FIELD fast.
   - Or just execute this file like a script to use the defaults.
### Example:
This example uses the test-files provided in *./test/*
```python
from fieldgenerator import *
#create Instane
generator = FieldGenerator()
#read input and write FIELD and CONFIG
generator.create()

#update a section with new values
gen.updateFF(list(range(0,36)),'ANGLES')
```

## Explenation of files and settings involved:
While creating an instance you can change most of the defaults:
```python
FieldGenerator(systemName='systemName', workDir='', systemXYZFileName='xyz.reference',
                  xyzFileNames=['orig.xyz'], ffFileNames=['force_field.dat'],
                  nMols=[1], molNames=['MoleculeName'], dataFileName='data.pickle',
                  unit='kcal', levcfg=0, imcon=2, engcfg=0.0,
                  floatFormat='{:20.6f}', fieldFileName='FIELD',
                  configFileName='CONFIG', disablePrinting=False)
```

- **systemXYZFileName:**   Contains the entire simulation-box in near xyz format:

```
   a1 a2 a3  --first three lines box vectors (not if imcon=0, meaning no periodicity)
   b1 b2 b3
   c1 c2 c4
   fancyC x y z (vx vy vz fx fy fz)  --then for every atom its name and position
   fancyO x y z (vx vy vz fx fy fz)  --(velocities and forces are optional)
   ...
```

- **xyzFileNames:**
   List of xyz-files of the different molecules involved. Order of appearance has to match the order in the systemXYZFile. If you have a periodic structure, make sure that you place the molecule in such a way, that all bonds can be found by using the three lattice vectors from the systemXYZFile. The origin for periodicity is always (0,0,0).

- **FFFileNames:**
   List of files containing the force-field data for the molecules in xyzFileNames, same ordering as the xyzFiles.

- **nMols:** List of integers giving the amount of molecules defined by xyzFileNames[i] and FFFileNames[i] in the entire
   system (systemXYZFile).
   
   
Number of lines in systemXYZFile must therfore be: `3 (box) + sum(nMols[:]*len(xyzFileNames[:]))`

- **ffFile**: Check the *force_field.dat* in the *./test/* folder for explenation. Supports simple mathematical expressions using `simpleeval` to make switching between different functional forms easier.
- `unit='kcal', levcfg=0, imcon=2, engcfg=0.0`: Check the DLPOLY Manual, Section FIELD-file.

## Advanced usage:

To add stuff that is not yet implemented add a function or use the self.fieldOut, which is a list of lines
of the future FIELD file. If you write a function please consider sending or commiting it to me for others to use too.

Using the save() and load() functions you can save an instance to a file and restore it later.
Don't forget to update the workDir after doing that! This feature is using `pickle`.

If you want to skip the input-file reading in the beginning consider setting the appropiate values and calling the rest of the functions manually.

## Contribute:

Feel free to forck and fix/update/add stuff. Please commit your changes and make a pull request so that other people can enjoy them too.
You can also just message me: chemistry@melix.me

## Acknowledgement:

Radii for bond-guessing taken from [PLAMS](https://github.com/SCM-NV/PLAMS) (LGPL-3.0)
If your molecule is not represented well with this method try playing with the factor radiiFactor
or overwrite the function.
