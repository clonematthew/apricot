##############################################################
###                                                        ###
###           Make File for the setpartGAD Code            ###
###                     pcc - 27/07/2006                   ###
###                                                        ###
##############################################################

.KEEP_STATE:

.SUFFIXES: .f .f90 .o

###FC = ifort 
###FFLAGS = -O3 -ipo -fpe0 -fno-alias -check all -traceback -static-libcxa -u -warn all
FC = ftn 
FFLAGS = -O3 -fomit-frame-pointer

FOBS = kdtree2.o cell_data.o code_options.o image_arepo.o \
        write_image.o make_image_arepo.o make_sgchem_images_arepo.o \
         uclchem.o ppv.o nhnh.o

EXEC = ../imagearepo

.f90.o :
	$(FC) $(FLAGS) -c $*.f90

$(EXEC): $(FOBS)
	$(FC) $(FOBS) $(FFLAGS) -o $(EXEC)

clean:
	\rm -f *.o *.mod $(EXEC)
