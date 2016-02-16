########## User-definable stuff ##########
#
###Compiler and compilation options
COMP_SER = gcc
COMP_MPI = mpicc
OPTIONS = -Wall -O3
#
### Behavioural flags
#Use double precision integer (enable in general)
DEFINEFLAGS += -D_LONGIDS
#Use custom-made frequency table
DEFINEFLAGS += -D_IRREGULAR_NUTABLE
#Generate debug help. Only useful for development
DEFINEFLAGS += -D_DEBUG
#Use double precision floating point? Set to "yes" or "no"
USE_SINGLE_PRECISION = yes
#Use OMP parallelization? Set to "yes" or "no"
USE_OMP = yes
#Use MPI parallelization? Set to "yes" or "no"
USE_MPI = yes
#
###Path to libraries and headers
###If two or more of the dependencies reside in the same paths, only
###one instance is necessary.
#GSL
GSL_INC = -I/home/dmonge/include
GSL_LIB = -L/home/dmonge/lib
#FFTW
FFTW_INC = 
FFTW_LIB = 
#HEALPix
HPIX_INC = 
HPIX_LIB = 
#LIBSHARP
LSHT_INC = -I/home/dmonge/Software/libsharp/auto/include
LSHT_LIB = -L/home/dmonge/Software/libsharp/auto/lib
#cfitsio
FITS_INC = 
FITS_LIB = 
#
########## End of user-definable ##########

ifeq ($(strip $(USE_OMP)),yes)
OPTIONS += -fopenmp
DEFINEFLAGS += -D_HAVE_OMP
endif #OMP

ifeq ($(strip $(USE_MPI)),yes)
DEFINEFLAGS += -D_HAVE_MPI
COMP_PAR = $(COMP_MPI)
else #MPI
COMP_PAR = $(COMP_SER)
endif #MPI

ifeq ($(strip $(USE_SINGLE_PRECISION)),yes)
DEFINEFLAGS += -D_SPREC

ifeq ($(strip $(USE_OMP)),yes)
LIB_FFTW += -lfftw3f_omp
endif #OMP
ifeq ($(strip $(USE_MPI)),yes)
LIB_FFTW += -lfftw3f_mpi
endif #MPI
LIB_FFTW += -lfftw3f

else #SINGLE_PRECISION

ifeq ($(strip $(USE_OMP)),yes)
LIB_FFTW += -lfftw3_omp
endif #OMP
ifeq ($(strip $(USE_MPI)),yes)
LIB_FFTW += -lfftw3_mpi
endif #MPI
LIB_FFTW += -lfftw3

endif #SINGLE_PRECISION


OPTIONS += $(DEFINEFLAGS)

INC_ALL = -I./src $(GSL_INC) $(FFTW_INC) $(HPIX_INC) $(FITS_INC) $(LSHT_INC)
LIB_ALL = $(GSL_LIB) $(FFTW_LIB) $(HPIX_LIB) $(FITS_LIB) $(LSHT_LIB) -lgsl -lgslcblas $(LIB_FFTW) -lchealpix -lsharp -lfftpack -lc_utils -lcfitsio -lm

COMMONO = src/common.o src/healpix_extra.o

COMMONO_MPI = src/common.o src/healpix_extra.o
COMMONGHO = src/common_gh.o
COSMOMADO = src/cosmo_mad.o
COSMOO = src/cosmo.o
FOURIERO = src/fourier.o
GRIDO = src/grid_tools.o
PSOURCES = src/psources.o
PIXELO = src/pixelize.o
IOGHO = src/io_gh.o
USERO = src/user_defined.o
MAINGH = src/main_gh.c
OFILESGH = $(COMMONO_MPI) $(COMMONGHO) $(USERO) $(COSMOMADO) $(COSMOO) $(FOURIERO) $(GRIDO) $(PSOURCES) $(PIXELO) $(IOGHO)

COMMONFGO = src/common_fg.o
IOFGO = src/io_fg.o
SCKO = src/sck_maps.o
GSIO = src/gsync_i.o
GSQUO = src/gsync_qu.o
MAINFG = src/main_fg.c
OFILESFG = $(COMMONO) $(COMMONFGO) $(IOFGO) $(SCKO) $(GSIO) $(GSQUO)

COMMONJTO = src/common_jt.o
IOJTO = src/io_jt.o
MAINJT = src/main_jt.c
OFILESJT = $(COMMONO) $(COMMONJTO) $(IOJTO)

EXEGH = GetHI
EXEFG = ForGet
EXEJT = JoinT

default : $(EXEGH) $(EXEFG) $(EXEJT)

%.o : %.c
	$(COMP_CC) $(OPTIONS) $(INC_ALL) -c $< -o $@

$(OFILESGH) : COMP_CC := $(COMP_PAR)

$(OFILESFG) $(OFILESJT) : COMP_CC := $(COMP_SER)

$(EXEGH) : $(OFILESGH)
	$(COMP_PAR) $(OPTIONS) $(INC_ALL) $(OFILESGH) $(MAINGH) -o $(EXEGH) $(LIB_ALL)

$(EXEFG) : $(OFILESFG)
	$(COMP_SER) $(OPTIONS) $(INC_ALL) $(OFILESFG) $(MAINFG) -o $(EXEFG) $(LIB_ALL)

$(EXEJT) : $(OFILESJT)
	$(COMP_SER) $(OPTIONS) $(INC_ALL) $(OFILESJT) $(MAINJT) -o $(EXEJT) $(LIB_ALL)

clean :
	rm -f src/*.o

cleaner : 
	rm -f *~ src/*.o src/*~  $(EXEGH) $(EXEFG) $(EXEJT)
