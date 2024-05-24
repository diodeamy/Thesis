#
# Use this makefile to create libraries for some of the programs that have a complex compilation
#

# paths to libraries
INC_PATHS = -I/net/plato/data/users/cautun/Programs/stow/include
LIB_PATHS = -L/net/plato/data/users/cautun/Programs/stow/lib
CC        = g++

# paths where to create the executables and libraries
OBJ_DIRS = ./o/
BIN_DIRS = ./bin/
LIB_DIRS = ./lib/		#this is not used



#########################################################
#                                                       #
#          Do not modify below this                     #
#                                                       #
#########################################################



# compilation flags and libraries; some additional variables
CCADDFLAGS = -O3 -fopenmp -DBOOST_DISABLE_ASSERTS
INC_PATH = $(INC_PATHS) -I./src/defines/ -I./src/classes/ -I./src/miscellaneous/ -I./src/gadgetFile/ -I./src/densityFile/ -I./src/MMF_file/ -I./src/MMF_computations/ -I./src/MMF_contraction -I./src/halo_file -I./src/output/
LIB_PATH = $(LIB_PATHS) -lfftw3_threads -lfftw3 -lfftw3f_threads -lfftw3f -lm -lpthread -lboost_filesystem -lboost_program_options -lgsl -lgslcblas

OBJ_DIR = $(OBJ_DIRS)/
BIN_DIR = $(BIN_DIRS)/
LIB_DIR = $(LIB_DIRS)/


# the final executables
nexus: MMF_response MMF_threshold MMF_clean MMF_spineDirection MMF_spineProperties
	cp $(addprefix python/, infoHeader.py MMF_combine.py density_toVisit.py density_raw.py density_fromRaw.py) $(BIN_DIR)

#nexus_extra: MMF_spineDirection MMF_spineProperties

#halo: halo_countSubstructure halo_directionVsScale halo_findSubstructure halo_findDistance




CLASSES_HEADERS = $(addprefix src/classes/, array.h arrayProperties.h binarySearch.h box.h FFTW.h FFTW_float.h grid_index.h k_space.h linear_fit.h matrix.h periodic.h vector.h )
CLASSES_FILES2  = $(addprefix classes/, FFTW_float.cc)
CLASSES_FILES   = $(addprefix src/, $(CLASSES_FILES2) )

DEFINES_HEADERS = $(addprefix src/defines/, defines.h)

MISC_HEADERS    = $(addprefix src/miscellaneous/, intervalCheck.h max.h memoryAllocationCheck.h message.h min.h misc.h miscellaneous.h optionsCheck.h quicksort.h range.h rescale.h throwError.h throwWarning.h Timer.h zero.h )
MISC_FILES2     = $(addprefix miscellaneous/, miscellaneousFile.cc miscellaneousMath.cc miscellaneousInput.cc miscellaneousUserInput.cc )
MISC_FILES      = $(addprefix src/, $(MISC_FILES2) )

GADGET_HEADERS  = $(addprefix src/gadgetFile/, gadget_header.h gadgetFile.h gadgetParticles.h parameters.ih)
GADGET_FILES2   = $(addprefix gadgetFile/, gadgetHeader.cc readGadgetFile.cc writeGadgetFile.cc)
GADGET_FILES    = $(addprefix src/, $(GADGET_FILES2) )

DENSITY_HEADERS = $(addprefix src/densityFile/, density_header.h densityFile.h )
DENSITY_FILES2  = $(addprefix densityFile/, densityHeader.cc )
DENSITY_FILES   = $(addprefix src/, $(DENSITY_FILES2) )

MMF_HEADERS     = $(addprefix src/MMF_file/, MMF_file.h MMF_header.h )
MMF_FILES2      = $(addprefix MMF_file/, MMF_header.cc )
MMF_FILES       = $(addprefix src/, $(MMF_FILES2) )

MMF_COMP_HEADERS= $(addprefix src/MMF_computations/, additionalComputations.h hessian.h MMF_computations.h objects.h response.h threshold.h )
MMF_COMP_FILES2 = $(addprefix MMF_computations/, additionalComputations.cc discard.cc eigenvectors.cc filters.cc hessian.cc objects.cc response.cc threshold.cc )
MMF_COMP_FILES  = $(addprefix src/, $(MMF_COMP_FILES2) )

MMF_CONT_HEADERS= $(addprefix src/MMF_contraction/, MMF_contraction.h )
MMF_CONT_FILES2 = $(addprefix MMF_contraction/, MMF_contraction.cc )
MMF_CONT_FILES  = $(addprefix src/, $(MMF_CONT_FILES2) )

HALO_HEADERS    = $(addprefix src/halo_file/, halo_computations.h halo_file.h halo_header.h )
HALO_FILES2     = $(addprefix halo_file/, haloHeader.cc )
HALO_FILES      = $(addprefix src/, $(HALO_FILES2) )

OUTPUT_HEADERS  = $(addprefix src/output/, compareDensityOutput.h compareMMFOutput.h MMF_thresholdOutput.h output.h raw.h )
OUTPUT_FILES2   = $(addprefix output/, ascii.cc densitySlices.cc visivo.cc )
OUTPUT_FILES    = $(addprefix src/, $(OUTPUT_FILES2) )



MMF_response: set_obj_dir $(CLASSES_HEADERS) $(MISC_HEADERS) $(DEFINES_HEADERS) $(DENSITY_HEADERS) $(MMF_HEADERS) $(MMF_COMP_HEADERS) $(OUTPUT_HEADERS) $(addsuffix .o,$(addprefix $(OBJ_DIR), MMF_response.cc $(CLASSES_FILES2) $(MISC_FILES2) $(DENSITY_FILES2) $(MMF_FILES2) $(MMF_COMP_FILES2))) Makefile
	$(CC) $(CCADDFLAGS)  $(addsuffix .o,$(addprefix $(OBJ_DIR), MMF_response.cc $(CLASSES_FILES2) $(MISC_FILES2) $(DENSITY_FILES2) $(MMF_FILES2) $(MMF_COMP_FILES2))) $(LIB_PATH) -o "$(BIN_DIR)/MMF_response"

$(OBJ_DIR)MMF_response.cc.o: set_obj_dir src/MMF_response.cc $(CLASSES_HEADERS) $(MISC_HEADERS) $(DEFINES_HEADERS) $(DENSITY_HEADERS) $(MMF_HEADERS) $(MMF_COMP_HEADERS) $(OUTPUT_HEADERS) Makefile
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)MMF_response.cc.o" "src/MMF_response.cc"


MMF_threshold: set_obj_dir $(CLASSES_HEADERS) $(MISC_HEADERS) $(DEFINES_HEADERS) $(DENSITY_HEADERS) $(MMF_HEADERS) $(MMF_COMP_HEADERS) $(OUTPUT_HEADERS) $(addsuffix .o,$(addprefix $(OBJ_DIR), MMF_threshold.cc $(CLASSES_FILES2) $(MISC_FILES2) $(DENSITY_FILES2) $(MMF_FILES2) $(MMF_COMP_FILES2))) Makefile
	$(CC) $(CCADDFLAGS)  $(addsuffix .o,$(addprefix $(OBJ_DIR), MMF_threshold.cc $(CLASSES_FILES2) $(MISC_FILES2) $(DENSITY_FILES2) $(MMF_FILES2) $(MMF_COMP_FILES2))) $(LIB_PATH) -o "$(BIN_DIR)/MMF_threshold"

$(OBJ_DIR)MMF_threshold.cc.o: set_obj_dir src/MMF_threshold.cc $(CLASSES_HEADERS) $(MISC_HEADERS) $(DEFINES_HEADERS) $(DENSITY_HEADERS) $(MMF_HEADERS) $(MMF_COMP_HEADERS) $(OUTPUT_HEADERS) Makefile
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)MMF_threshold.cc.o" "src/MMF_threshold.cc"


MMF_clean: set_obj_dir $(CLASSES_HEADERS) $(MISC_HEADERS) $(DEFINES_HEADERS) $(DENSITY_HEADERS) $(MMF_HEADERS) $(MMF_COMP_HEADERS) $(OUTPUT_HEADERS) $(addsuffix .o,$(addprefix $(OBJ_DIR), MMF_clean.cc $(CLASSES_FILES2) $(MISC_FILES2) $(DENSITY_FILES2) $(MMF_FILES2) $(MMF_COMP_FILES2))) Makefile
	$(CC) $(CCADDFLAGS)  $(addsuffix .o,$(addprefix $(OBJ_DIR), MMF_clean.cc $(CLASSES_FILES2) $(MISC_FILES2) $(DENSITY_FILES2) $(MMF_FILES2) $(MMF_COMP_FILES2))) $(LIB_PATH) -o "$(BIN_DIR)/MMF_clean"

$(OBJ_DIR)MMF_clean.cc.o: set_obj_dir src/MMF_clean.cc $(CLASSES_HEADERS) $(MISC_HEADERS) $(DEFINES_HEADERS) $(DENSITY_HEADERS) $(MMF_HEADERS) $(MMF_COMP_HEADERS) $(OUTPUT_HEADERS) Makefile
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)MMF_clean.cc.o" "src/MMF_clean.cc"



MMF_spineDirection: set_obj_dir $(CLASSES_HEADERS) $(MISC_HEADERS) $(DEFINES_HEADERS) $(DENSITY_HEADERS) $(MMF_HEADERS) $(MMF_COMP_HEADERS) $(MMF_CONT_HEADERS) $(OUTPUT_HEADERS) $(addsuffix .o,$(addprefix $(OBJ_DIR), MMF_spineDirection.cc $(CLASSES_FILES2) $(MISC_FILES2) $(DENSITY_FILES2) $(MMF_FILES2) $(MMF_COMP_FILES2) $(MMF_CONT_FILES2))) Makefile
	$(CC) $(CCADDFLAGS)  $(addsuffix .o,$(addprefix $(OBJ_DIR), MMF_spineDirection.cc $(CLASSES_FILES2) $(MISC_FILES2) $(DENSITY_FILES2) $(MMF_FILES2) $(MMF_COMP_FILES2) $(MMF_CONT_FILES2))) $(LIB_PATH) -o "$(BIN_DIR)/MMF_spineDirection"

$(OBJ_DIR)MMF_spineDirection.cc.o: set_obj_dir src/MMF_spineDirection.cc $(CLASSES_HEADERS) $(MISC_HEADERS) $(DEFINES_HEADERS) $(MMF_HEADERS) $(MMF_COMP_HEADERS) $(MMF_CONT_HEADERS) $(OUTPUT_HEADERS) Makefile
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)MMF_spineDirection.cc.o" "src/MMF_spineDirection.cc"


MMF_spineProperties: set_obj_dir $(CLASSES_HEADERS) $(MISC_HEADERS) $(DEFINES_HEADERS) $(DENSITY_HEADERS) $(MMF_HEADERS) $(MMF_COMP_HEADERS) $(MMF_CONT_HEADERS) $(OUTPUT_HEADERS) $(addsuffix .o,$(addprefix $(OBJ_DIR), MMF_spineProperties.cc $(CLASSES_FILES2) $(MISC_FILES2) $(DENSITY_FILES2) $(MMF_FILES2) $(MMF_COMP_FILES2) $(MMF_CONT_FILES2))) Makefile
	$(CC) $(CCADDFLAGS)  $(addsuffix .o,$(addprefix $(OBJ_DIR), MMF_spineProperties.cc $(CLASSES_FILES2) $(MISC_FILES2) $(DENSITY_FILES2) $(MMF_FILES2) $(MMF_COMP_FILES2) $(MMF_CONT_FILES2))) $(LIB_PATH) -o "$(BIN_DIR)/MMF_spineProperties"

$(OBJ_DIR)MMF_spineProperties.cc.o: set_obj_dir src/MMF_spineProperties.cc $(CLASSES_HEADERS) $(MISC_HEADERS) $(DEFINES_HEADERS) $(MMF_HEADERS) $(MMF_COMP_HEADERS) $(MMF_CONT_HEADERS) $(OUTPUT_HEADERS) Makefile
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)MMF_spineProperties.cc.o" "src/MMF_spineProperties.cc"



halo_findSubstructure: set_obj_dir $(CLASSES_HEADERS) $(MISC_HEADERS) $(DEFINES_HEADERS) $(HALO_HEADERS) $(OUTPUT_HEADERS) $(addsuffix .o,$(addprefix $(OBJ_DIR), halo_findSubstructure.cc $(CLASSES_FILES2) $(MISC_FILES2) $(HALO_FILES2))) Makefile
	$(CC) $(CCADDFLAGS)  $(addsuffix .o,$(addprefix $(OBJ_DIR), halo_findSubstructure.cc $(CLASSES_FILES2) $(MISC_FILES2) $(HALO_FILES2))) $(LIB_PATH) -o "$(BIN_DIR)/halo_findSubstructure"

$(OBJ_DIR)halo_findSubstructure.cc.o: set_obj_dir src/halo_findSubstructure.cc $(CLASSES_HEADERS) $(MISC_HEADERS) $(DEFINES_HEADERS) $(HALO_HEADERS) $(OUTPUT_HEADERS) Makefile
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)halo_findSubstructure.cc.o" "src/halo_findSubstructure.cc"


halo_countSubstructure: set_obj_dir $(CLASSES_HEADERS) $(MISC_HEADERS) $(DEFINES_HEADERS) $(HALO_HEADERS) $(OUTPUT_HEADERS) $(addsuffix .o,$(addprefix $(OBJ_DIR), halo_countSubstructure.cc $(CLASSES_FILES2) $(MISC_FILES2) $(HALO_FILES2))) Makefile
	$(CC) $(CCADDFLAGS)  $(addsuffix .o,$(addprefix $(OBJ_DIR), halo_countSubstructure.cc $(CLASSES_FILES2) $(MISC_FILES2) $(HALO_FILES2))) $(LIB_PATH) -o "$(BIN_DIR)/halo_countSubstructure"

$(OBJ_DIR)halo_countSubstructure.cc.o: set_obj_dir src/halo_countSubstructure.cc $(CLASSES_HEADERS) $(MISC_HEADERS) $(DEFINES_HEADERS) $(HALO_HEADERS) $(OUTPUT_HEADERS) Makefile
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)halo_countSubstructure.cc.o" "src/halo_countSubstructure.cc"


halo_findDistance: set_obj_dir $(CLASSES_HEADERS) $(MISC_HEADERS) $(DEFINES_HEADERS) $(HALO_HEADERS) $(OUTPUT_HEADERS) $(addsuffix .o,$(addprefix $(OBJ_DIR), halo_findDistance.cc $(CLASSES_FILES2) $(MISC_FILES2) $(HALO_FILES2))) Makefile
	$(CC) $(CCADDFLAGS)  $(addsuffix .o,$(addprefix $(OBJ_DIR), halo_findDistance.cc $(CLASSES_FILES2) $(MISC_FILES2) $(HALO_FILES2))) $(LIB_PATH) -o "$(BIN_DIR)/halo_findDistance"

$(OBJ_DIR)halo_findDistance.cc.o: set_obj_dir src/halo_findDistance.cc $(CLASSES_HEADERS) $(MISC_HEADERS) $(DEFINES_HEADERS) $(HALO_HEADERS) $(OUTPUT_HEADERS) Makefile
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)halo_findDistance.cc.o" "src/halo_findDistance.cc"


halo_directionVsScale: set_obj_dir $(CLASSES_HEADERS) $(MISC_HEADERS) $(DEFINES_HEADERS) $(DENSITY_HEADERS) $(HALO_HEADERS) $(OUTPUT_HEADERS) $(addsuffix .o,$(addprefix $(OBJ_DIR), halo_directionVsScale.cc $(CLASSES_FILES2) $(MISC_FILES2) $(DENSITY_FILES2) $(HALO_FILES2))) Makefile
	$(CC) $(CCADDFLAGS)  $(addsuffix .o,$(addprefix $(OBJ_DIR), halo_directionVsScale.cc $(CLASSES_FILES2) $(MISC_FILES2) $(DENSITY_FILES2) $(HALO_FILES2))) $(LIB_PATH) -o "$(BIN_DIR)/halo_directionVsScale"

$(OBJ_DIR)halo_directionVsScale.cc.o: set_obj_dir src/halo_directionVsScale.cc $(CLASSES_HEADERS) $(MISC_HEADERS) $(DEFINES_HEADERS) $(DENSITY_HEADERS) $(HALO_HEADERS) $(OUTPUT_HEADERS) Makefile
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)halo_directionVsScale.cc.o" "src/halo_directionVsScale.cc"





$(OBJ_DIR)classes/FFTW_float.cc.o: $(CLASSES_HEADERS) $(CLASSES_FILES) $(DEFINES_HEADERS) $(MISC_HEADERS) Makefile
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)classes/FFTW_float.cc.o" "src/classes/FFTW_float.cc"



$(OBJ_DIR)miscellaneous/miscellaneousFile.cc.o: $(MISC_HEADERS) $(MISC_FILES) $(DEFINES_HEADERS) $(MISC_HEADERS) Makefile
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)miscellaneous/miscellaneousFile.cc.o" "src/miscellaneous/miscellaneousFile.cc" 

$(OBJ_DIR)miscellaneous/miscellaneousInput.cc.o: $(MISC_HEADERS) $(MISC_FILES) $(DEFINES_HEADERS) $(MISC_HEADERS) Makefile
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)miscellaneous/miscellaneousInput.cc.o" "src/miscellaneous/miscellaneousInput.cc" 

$(OBJ_DIR)miscellaneous/miscellaneousMath.cc.o: $(MISC_HEADERS) $(MISC_FILES) $(DEFINES_HEADERS) $(MISC_HEADERS) Makefile
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)miscellaneous/miscellaneousMath.cc.o" "src/miscellaneous/miscellaneousMath.cc" 

$(OBJ_DIR)miscellaneous/miscellaneousUserInput.cc.o: $(MISC_HEADERS) $(MISC_FILES) $(DEFINES_HEADERS) $(MISC_HEADERS) Makefile
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)miscellaneous/miscellaneousUserInput.cc.o" "src/miscellaneous/miscellaneousUserInput.cc" 


$(OBJ_DIR)densityFile/densityHeader.cc.o: $(DENSITY_HEADERS) $(MISC_FILES) $(DEFINES_HEADERS) $(MISC_HEADERS) $(GADGET_HEADERS) Makefile
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)densityFile/densityHeader.cc.o" "src/densityFile/densityHeader.cc"


$(OBJ_DIR)MMF_file/MMF_header.cc.o: $(MMF_HEADERS) $(MMF_FILES) $(DEFINES_HEADERS) $(MISC_HEADERS) $(CLASSES_HEADERS) Makefile $(DENSITY_HEADERS)
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)MMF_file/MMF_header.cc.o" "src/MMF_file/MMF_header.cc" 


$(OBJ_DIR)MMF_computations/additionalComputations.cc.o: $(MMF_COMP_HEADERS) $(MMF_COMP_FILES) $(DEFINES_HEADERS) Makefile $(MISC_HEADERS) $(CLASSES_HEADERS)
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)MMF_computations/additionalComputations.cc.o" "src/MMF_computations/additionalComputations.cc" 

$(OBJ_DIR)MMF_computations/discard.cc.o: $(MMF_COMP_HEADERS) $(MMF_COMP_FILES) $(DEFINES_HEADERS) $(MISC_HEADERS) Makefile $(CLASSES_HEADERS)
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)MMF_computations/discard.cc.o" "src/MMF_computations/discard.cc"

$(OBJ_DIR)MMF_computations/eigenvectors.cc.o: $(MMF_COMP_HEADERS) $(MMF_COMP_FILES) $(DEFINES_HEADERS) $(MISC_HEADERS) Makefile $(CLASSES_HEADERS)
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)MMF_computations/eigenvectors.cc.o" "src/MMF_computations/eigenvectors.cc" 

$(OBJ_DIR)MMF_computations/filters.cc.o: $(MMF_COMP_HEADERS) $(DEFINES_HEADERS) $(MISC_HEADERS) $(CLASSES_HEADERS) Makefile
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)MMF_computations/filters.cc.o" "src/MMF_computations/filters.cc" 

$(OBJ_DIR)MMF_computations/hessian.cc.o: $(MMF_COMP_HEADERS) $(MMF_COMP_FILES) $(DEFINES_HEADERS) $(MISC_HEADERS) Makefile $(CLASSES_HEADERS)
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)MMF_computations/hessian.cc.o" src/"MMF_computations/hessian.cc" 

$(OBJ_DIR)MMF_computations/objects.cc.o: $(MMF_COMP_HEADERS) $(MMF_COMP_FILES) $(DEFINES_HEADERS) $(MISC_HEADERS) Makefile $(CLASSES_HEADERS)
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)MMF_computations/objects.cc.o" "src/MMF_computations/objects.cc"

$(OBJ_DIR)MMF_computations/response.cc.o: $(MMF_COMP_HEADERS) $(MMF_COMP_FILES) $(DEFINES_HEADERS) $(MISC_HEADERS) Makefile $(CLASSES_HEADERS)
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)MMF_computations/response.cc.o" "src/MMF_computations/response.cc" 
    
$(OBJ_DIR)MMF_computations/threshold.cc.o: $(MMF_COMP_HEADERS) $(MMF_COMP_FILES) $(DEFINES_HEADERS) $(MISC_HEADERS) Makefile $(CLASSES_HEADERS)
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)MMF_computations/threshold.cc.o" "src/MMF_computations/threshold.cc"  


$(OBJ_DIR)MMF_contraction/MMF_contraction.cc.o: $(MMF_HEADERS) $(MMF_CONT_FILES) $(DEFINES_HEADERS) $(MISC_HEADERS) Makefile $(CLASSES_HEADERS)
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)MMF_contraction/MMF_contraction.cc.o" "src/MMF_contraction/MMF_contraction.cc"  


$(OBJ_DIR)halo_file/haloHeader.cc.o: $(HALO_HEADERS) $(HALO_FILES) $(DEFINES_HEADERS) $(MISC_HEADERS) $(CLASSES_HEADERS) Makefile
	$(CC) $(CCADDFLAGS)  -c $(INC_PATH) -o "$(OBJ_DIR)halo_file/haloHeader.cc.o" "src/halo_file/haloHeader.cc" 




set_obj_dir:
	@ mkdir -p $(OBJ_DIR)
	@ mkdir -p $(OBJ_DIR) $(OBJ_DIR)/classes $(OBJ_DIR)/miscellaneous $(OBJ_DIR)/halo_file $(OBJ_DIR)/densityFile $(OBJ_DIR)/defines $(OBJ_DIR)/MMF_file $(OBJ_DIR)/MMF_computations $(OBJ_DIR)/MMF_contraction $(OBJ_DIR)/output
	@ mkdir -p $(BIN_DIR)

clean:
	rm -f $(BIN_DIR)/*
	rm -f -r $(OBJ_DIR)/*
	

set_directories:
	@ if !(test -d $(OBJ_DIR)); \
	then mkdir $(OBJ_DIR); \
	fi
	@ if !( test -d $(BIN_DIR) ); \
	then mkdir $(BIN_DIR); \
	fi
	@ if !( test -d $(LIB_DIR) ); \
	then mkdir $(LIB_DIR); \
	fi

EXEC_FILES = `find -type f -perm -o+rx`
copy: set_directories
	mv $(EXEC_FILES) $(BIN_DIR)
