########## User-definable stuff ##########
#
###Compiler and compilation options
COMP = gcc
COMP_MPI = mpicc
OPTIONS = -Wall -O3
#
### Behavioural flags
#Use double precision integer (enable in general)
DEFINEFLAGS += -D_LONGIDS
#Output extra information
DEFINEFLAGS += -D_VERBOSE
#Generate debug help. Only useful for development
#DEFINEFLAGS += -D_DEBUG
#Use double precision floating point? Set to "yes" or "no"
USE_SINGLE_PRECISION = yes
#Use MPI parallelization? Set to "yes" or "no"
USE_MPI = yes
#Use OpenMP parallelization? Set to "yes" or "no". This is compatible with USE_MPI = yes
USE_OMP = yes
#
###Path to libraries and headers
###If two or more of the dependencies reside in the same paths, only
###one instance is necessary.
#GSL
GSL_INC = -I/home/damonge/include
GSL_LIB = -L/home/damonge/lib
#FFTW
FFTW_INC = 
FFTW_LIB = 
#HEALPix
HPIX_INC = 
HPIX_LIB = 
#LIBSHARP
LSHT_INC =
LSHT_LIB =
#cfitsio
FITS_INC = 
FITS_LIB = 
#
########## End of user-definable ##########

ifeq ($(strip $(USE_OMP)),yes)
DEFINEFLAGS += -D_HAVE_OMP
OPTIONS += -fopenmp
endif

ifeq ($(strip $(USE_MPI)),yes)
DEFINEFLAGS += -D_HAVE_MPI
COMP_CC = $(COMP_MPI)
else
COMP_CC = $(COMP)
endif

ifeq ($(strip $(USE_SINGLE_PRECISION)),yes)
DEFINEFLAGS += -D_SPREC


LIB_FFTW = 
ifeq ($(strip $(USE_OMP)),yes)
LIB_FFTW += -lfftw3f_omp
endif #OMP
ifeq ($(strip $(USE_MPI)),yes)
LIB_FFTW += -lfftw3f_mpi
endif #MPI
LIB_FFTW += -lfftw3f

else #SINGLE_PRECISION

LIB_FFTW = 
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

COMMONO = src/common.o
EHPIXO = src/healpix_extra.o
RNGO = src/rng.o

#COMMONGHO = src/common_gh.o
#COSMOMADO = src/cosmo_mad.o
COSMOO = src/cosmo.o
FOURIERO = src/fourier.o
GRIDO = src/grid2maps.o
#PSOURCES = src/psources.o
#PIXELO = src/pixelize.o
IOGHO = src/io_gh.o
USERO = src/user_defined.o
MAINGH = src/main_gh.c
OFILESGH = $(COMMONO) $(EHPIXO) $(RNGO) $(USERO) $(COMMONGHO) $(COSMOO) $(IOGHO) $(FOURIERO) $(GRIDO)
#$(USERO) $(COSMOMADO) $(FOURIERO) $(GRIDO) $(PSOURCES) $(PIXELO)

#COMMONFGO = src/common_fg.o
#IOFGO = src/io_fg.o
#SCKO = src/sck_maps.o
#GSIO = src/gsync_i.o
#GSQUO = src/gsync_qu.o
#MAINFG = src/main_fg.c
#OFILESFG = $(COMMONO) $(EHPIXO) $(COMMONFGO) $(IOFGO) $(SCKO) $(GSIO) $(GSQUO)

#COMMONJTO = src/common_jt.o
#IOJTO = src/io_jt.o
#MAINJT = src/main_jt.c
#OFILESJT = $(COMMONO) $(EHPIXO) $(COMMONJTO) $(IOJTO)

EXEGH = GetHI
#EXEFG = ForGet
#EXEJT = JoinT

default : $(EXEGH) $(EXEFG) $(EXEJT)

%.o : %.c
	$(COMP_CC) $(OPTIONS) $(INC_ALL) -c $< -o $@

$(EXEGH) : $(OFILESGH)
	$(COMP_CC) $(OPTIONS) $(INC_ALL) $(OFILESGH) $(MAINGH) -o $(EXEGH) $(LIB_ALL)

$(EXEFG) : $(OFILESFG)
	$(COMP_CC) $(OPTIONS) $(INC_ALL) $(OFILESFG) $(MAINFG) -o $(EXEFG) $(LIB_ALL)

$(EXEJT) : $(OFILESJT)
	$(COMP_CC) $(OPTIONS) $(INC_ALL) $(OFILESJT) $(MAINJT) -o $(EXEJT) $(LIB_ALL)

clean :
	rm -f src/*.o

cleaner : 
	rm -f *~ src/*.o src/*~  $(EXEGH) $(EXEFG) $(EXEJT)
