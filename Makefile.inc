## ============================================================================
## General makefile for edf, isochrones and chem
## ============================================================================
## Dependencies:
## 	1. Tact
TACTPATH=/home/ktfm2/Project_Code/tact/

##  2. GSL
GSLPATH =/home/ktfm2/Project_Code/gsl-2.6/

##  3. GalPot -- implementation in TorusMapper code
TORUSPATH = /home/jls/work/code/Torus_github/Torus/

##  4. nlohmann/json
JSONPATH=/home/ktfm2/Project_Code/json/include/nlohmann/

##  5. cfitsio, healpix
HEALPIXPATH = /opt/ioa/software/healpix/Healpix_3.30/

##  6. logging
LOGGPATH = /home/ktfm2/Project_Code/easyloggingpp-9.96.7/

##  7. google test (if want to test)
GTESTPATH = /data/jls/gtest-1.7.0/
####################################################################
## Python module
####################################################################
## These paths are only important if you want to make python modules
## Relies on boost

BOOSTINCPATH =/opt/ioa/software/boost/1.55.0/include
BOOSTLIBPATH =-L/opt/ioa/software/boost/1.55.0/lib
PYTHONINCPATH=/usr/include/python2.7/
ifeq ($(USER),jls)
	PYTHONINCPATH=/opt/ioa/software/python/2.7.8/include/python2.7
endif
PYTHONLIBPATH=-L/opt/ioa/software/python/2.7.8/lib/

## ============================================================================
## ============================================================================
## ============================================================================
CCOMPILER=g++
CCOMPILER_MPI=mpic++
OPTS = --std=c++11 -O3 -fopenmp -fPIC #-ffast-math
# #------------------------------------------------------
# # Set paths to other code
## ============================================================================
## Tact
GENERALPATH = $(TACTPATH)general/
CUBAPATH=/home/ktfm2/Project_Code/tact/general/cuba/
CUBAINCPATH=-I$(CUBAPATH)
CUBALIBPATH=-L$(CUBAPATH)
CTPATH = $(GENERALPATH)coordtransforms/
POTPATH = $(TACTPATH)pot/
AAPATH = $(TACTPATH)aa/
GTESTINC=$(GTESTPATH)/include/
GTESTLIB=$(GTESTPATH)/lib/.libs/
## ============================================================================

JSONINC=$(JSONPATH)
LOGGINC=$(LOGGPATH)src/

TORUSINCPATH = -I$(TORUSPATH)src/ -I$(TORUSPATH)src/pot/ -I$(TORUSPATH)src/utils/ -I$(TORUSPATH)WDlib/inc/
TORUSLIBPATH = -L$(TORUSPATH)obj/ -L$(TORUSPATH)WDlib/obj/ \
               -L/data/jls/libebf_c_cpp-0.0.3/lib/
HEALPIX_INC = -I$(HEALPIXPATH)include/ \
			  -I$(HEALPIXPATH)src/cxx/optimized_gcc/include/ \
			  -I$(HEALPIXPATH)src/cxx/cxxsupport/
HEALPIX_LIB = -L$(HEALPIXPATH)lib/ \
			  -L$(HEALPIXPATH)src/cxx/optimized_gcc/lib/ \
			  -Wl,-rpath=$(HEALPIXPATH)src/cxx/optimized_gcc/lib/ \
			  -lhealpix_cxx -lchealpix -lcxxsupport -lcfitsio

HDF5PATH=/usr/
#HDF5PATH=/usr/
HDF5INCPATH=$(HDF5PATH)include/hdf5/serial/
HDF5LIBPATH=-L$(HDF5PATH)lib/x86_64-linux-gnu/ -Wl,-rpath=$(HDF5PATH)lib/x86_64-linux-gnu/

INCPY = $(INC) -I$(PYTHONINCPATH) -I$(BOOSTINCPATH)
BOOSTLIB = -lboost_python -lpython2.7
LIBPY = $(PYTHONLIBPATH) $(BOOSTLIBPATH) $(BOOSTLIB)
## ============================================================================
