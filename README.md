# solveME: fast and reliable solution of nonlinear ME models

solveME toolbox for solving nonlinear ME models.
Includes interfaces (as external modules) to the quad MINOS-based 
Fortran 90 code by Ding Ma and Michael A. Saunders at Stanford University

If you use solveME in a scientific publication, please cite:

[Yang, L., Ma, D., Ebrahim, A., Lloyd, C. J., Saunders, M. A., & Palsson, B. O. (2016). solveME: fast and reliable solution of nonlinear ME models. BMC Bioinformatics, 17(1), 391. doi:10.1186/s12859-016-1240-1](http://doi.org/10.1186/s12859-016-1240-1)

and the following for the Quad MINOS solver:

[Ma, D., Yang, L., Fleming, R.M.T., Thiele, I., Palsson, B.O., Saunders, M.A. (2017). Reliable and efficient solution of genome-scale models of Metabolism and macromolecular Expression. Scientific Reports, 7, 40863. doi:10.1038/srep40863](http://doi.org/10.1038/srep40863)

Author: Laurence Yang

Systems Biology Research Group, UCSD

## Requirements
1. Python
	- tested on Python 2.7 and 3.6
1. cobrapy
	- tested on 0.5.11
	- not yet fully compatible with cobrapy 0.6 release
1. cobrame: follow Installation instructions here: https://github.com/SBRG/cobrame
1. gfortran (>=4.6)
	- (Ubuntu) sudo apt-get install gfortran
	- (Arch) sudo pacman -S gcc-fortran
1. quadMINOS (available for academic use from Prof. Michael A. Saunders at Stanford University)

## Installation
1. Compile quadMINOS
	- cd qminos_root
	- cd minos56; make clean; make
	- cd ..
	- cd qminos56; make clean; make
1. Copy shared libraries to solveme root directory
	- cd [solveme_root]
	- cp [qminos_root]/qminos56/lib/libquadminos.a ./
	- cp [qminos_root]/minos56/lib/libminos.a ./
1. Run: python setup.py **develop**
	- or python setup.py develop --user
	- or sudo python setup.py develop to install for all users
		- (if getting an error, [undefined reference to main], try sudo python setup.py develop)

1. *Special NERSC installation steps*
	- (shell) module swap PrgEnv-intel PrgEnv-gnu
	- python setup.py **config_fc --f90exec=ftn** develop --user
		- (need to use the **ftn** gfortran Cray wrapper, or else get symbol not found error during import)
1. ***Troubleshooting common errors***
	1. ImportError: .../solveME/solvemepy/qminospy/qwarmLP.so: symbol _gfortran_transfer_real128_write, version GFORTRAN_1.4 not defined in file libgfortran.so.3 with link time reference
		- Problem: setup.py picked up the wrong libgfortran during installation. Typically, if anaconda is installed, an older version is installed, which seems to be used by default.
		- Solution (workaround): pre-load the desired libgfortran version, bypassing the anaconda one:
			- `export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libgfortran.so.3`
			- (then, proceed with setup) python setup.py develop --user
		- [Link to original solution](https://github.com/ContinuumIO/anaconda-issues/issues/686)
            

### Use qminos to solve ME models in python
### For (reduced) ME models prior prototype 44
1. Import qminospy to access the solver methods:
~~~~{.python}
from qminospy.me2 import ME_NLP
me_nlp = ME_NLP(me)
# Solve directly using LCL (linearly constrained Lagrangian) NLP method
x,stat,hs = me_nlp.solvenlp()
# Access the solution that is saved in the original minime object
sol = me.solution
sol.f
sol.x_dict
~~~~

### 24 Feb 2016: for ME models after prototype 44
1. Import 

py to access the solver methods:
~~~~{.python}
from qminospy.me1 import ME_NLP1
# The object containing solveME methods--composite that uses a ME model object 
# Provide growth_key = 'mu' for minime models,
me_nlp = ME_NLP1(me, growth_key='mu')
# Use bisection for now (until the NLP formulation is worked out for the new prototype 44
muopt, hs, xopt, cache = me_nlp.bisectmu(precision=1e-6)    
# Access the solution that is saved in the original minime object
sol = me.solution
sol.f
sol.x_dict
~~~~

### If your ME model is based on ME 1.0 code (iOL1650, iJL1678):
1. Same as above but use growth_key='growth_rate_in_per_hour'
~~~~{.python}
from qminospy.me1 import ME_NLP1
# The object containing solveME methods--composite that uses a ME model object 
me_nlp = ME_NLP1(me, growth_key='growth_rate_in_per_hour')
# Bisection
muopt, hs, xopt, cache = me_nlp.bisectmu(precision=1e-3)    
# Access the solution that is saved in the original minime object
sol = me.solution
sol.f
sol.x_dict
~~~~
