##############################################################
###                                                        ###
###           Make File for the setpartGAD Code            ###
###                     pcc - 27/07/2006                   ###
###                                                        ###
##############################################################

.KEEP_STATE:

.SUFFIXES: .f .f90 .o

FC = gfortran -lhdf5_fortran
FFLAGS = -O3

FOBS = kdtree2.o cell_data.o code_options.o image_arepo.o \
        write_image.o make_image_arepo.o make_sgchem_images_arepo.o \
         uclchem.o ppv.o nhnh.o

EXEC = ../imagearepo

.f90.o :
	$(FC) $(FLAGS) -c $*.f90

$(EXEC): $(FOBS)
	$(FC) $(FOBS) $(FFLAGS) -o $(EXEC)

clean:
	\rm -f *.o kdtree2.mod kdtree2_module.mod kdtree2_precision_module.mod kdtree2_priority_queue_module.mod celldata.mod codeoptions.mod $(EXEC)
