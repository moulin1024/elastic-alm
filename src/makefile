################################################################################
# FC
################################################################################

FC = gfortran

################################################################################
# FFLAGS
################################################################################

FFLAGS = -O3 -frecursive#-Mpreprocess# -mcmodel=medium
# CUDA = -Mcuda -lcufft -ta=nvidia
# debug
# -check bounds -debug all -traceback -fpe0 fp-stack-check -debug inline-debug-info -parallel-source-info=2 -warn all

# prof
# -pg

# NVTX uncomment next line to add nvtx markers. Make sure that CUDA_HOME is available in your env.
# USE_NVTX = 1



################################################################################
# LIBS
################################################################################

# LIBS = -L/home/moulin/anaconda3/lib -lfftw3
# INC = #-I/opt/intel/compilers_and_libraries_2019.4.243/linux/mpi/intel64/include/
################################################################################
# OBJ, PROG
################################################################################
PROG = main

DEP0 = parameters.o structural.o
# DEP1  =   dimen.o cufft_m.o
# DEP2  =   update.o convection.o derivative.o filter.o io.o pressure.o rhs.o time.o control.o
# DEP3  =   divergence.o turbine_alm_nrel.o 
# DEP4  =   sgs.o sgs_mgm.o wall.o


# ifdef USE_NVTX
# 	FFLAGS += -DUSE_NVTX -lnvToolsExt -L$(CUDA_HOME)/lib64/
# 	DEP0   += nvtx.o
# endif

DEPS  = $(DEP0)# $(DEP1) $(DEP2) $(DEP3) $(DEP4) 
OBJ   = $(DEPS) main.o


# --- link ---
$(PROG):$(OBJ) main.o
	$(FC) $(FFLAGS) $(CUDA) $(OBJ) -o $(PROG) $(INC) $(LIBS)

# --- compile ---
main.o: main.f90 $(DEPS)
	$(FC) $(FFLAGS) $(CUDA) $(INC) $(PRECISION) -c $<
$(DEP0): %.o : %.f90
	$(FC) $(FFLAGS) $(CUDA) $(INC) $(PRECISION) -c $<
# $(DEP1): %.o: %.cuf precision.o
# 	$(FC) $(FFLAGS) $(CUDA) $(INC) $(PRECISION) -c $<
# $(DEP2): %.o: %.cuf precision.o dimen.o
# 	$(FC) $(FFLAGS) $(CUDA) $(INC) $(PRECISION) -c $<
# $(DEP3): %.o: %.cuf precision.o dimen.o derivative.o control.o
# 	$(FC) $(FFLAGS) $(CUDA) $(INC) $(PRECISION) -c $<
# $(DEP4): %.o: %.cuf precision.o dimen.o cufft_m.o filter.o
# 	$(FC) $(FFLAGS) $(CUDA) $(INC) $(PRECISION) -c $<

# --- clean ---
clean:
	rm -f *.o *.bck *.bak *.mod $(PROG)
