
# paths
HIGHS_PATH = $(HOME)/Documents/HiGHS
CHOLMOD_PATH = $(HOME)/Documents/SuiteSparse
METIS_PATH = $(HOME)/Documents/METIS
LOCAL_PATH = $(HOME)/local

# source files
cpp_sources = \
		mainIPM.cpp \
		Ipm.cpp \
		IpmModel.cpp \
		VectorOperations.cpp \
		Ipm_aux.cpp \
		FactorHiGHSSolver.cpp \
		CurtisReidScaling.cpp \
		../FactorHiGHS/Analyse.cpp \
		../FactorHiGHS/Auxiliary.cpp \
		../FactorHiGHS/Factorise.cpp \
		../FactorHiGHS/Numeric.cpp \
		../FactorHiGHS/Symbolic.cpp \
		../FactorHiGHS/FormatHandler.cpp \
		../FactorHiGHS/FullFormatHandler.cpp \
		../FactorHiGHS/HybridPackedFormatHandler.cpp \
		../FactorHiGHS/HybridHybridFormatHandler.cpp \
		../FactorHiGHS/DataCollector.cpp
c_sources = ../FactorHiGHS/DenseFact.c \
			../FactorHiGHS/CallAndTimeBlas.c \
			../FactorHiGHS/timing.c

# binary file name
binary_name = ipm

# object files directory
OBJDIR = obj

# compilers
#CC=clang
#CPP=clang++
 CPP = /opt/homebrew/Cellar/llvm/17.0.6_1/bin/clang++
 CC = /opt/homebrew/Cellar/llvm/17.0.6_1/bin/clang

# compiler flags
CPPFLAGS = -std=c++11 -O3 -g3 -Wno-deprecated #-fsanitize=address #ASAN_OPTIONS=detect_leaks=1
CFLAGS = -O3 -g3 #-fsanitize=address

# mess to link openmp on mac
#OPENMP_FLAGS = -Xclang -fopenmp -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib -lomp

# rpaths for shared libraries
rpaths = -rpath $(CHOLMOD_PATH)/lib/ -rpath $(LOCAL_PATH)/lib/

# includes and libraries
includes = -I$(HIGHS_PATH)/build -I$(HIGHS_PATH)/src/ -I$(METIS_PATH)/include -I$(LOCAL_PATH)/include
libs_path = -L$(HIGHS_PATH)/build/lib -L$(METIS_PATH)/build/libmetis -L$(LOCAL_PATH)/lib
libs = -lhighs -lmetis -lGKlib -llapack -lblas

# name of objects
cpp_objects = $(cpp_sources:%.cpp=$(OBJDIR)/%.o)
c_objects = $(c_sources:%.c=$(OBJDIR)/%.o)

# dependency files
dep = $(cpp_sources:%.cpp=$(OBJDIR)/%.d)




# link ipm
$(binary_name): $(cpp_objects) $(c_objects)
	@echo Linking objects into $@
	@$(CPP) $(CPPFLAGS) $(libs_path) $(libs) $(rpaths) $^ -o $@

# manage dependencies
-include $(dep)

# compile cpp
$(cpp_objects): $(OBJDIR)/%.o: %.cpp
	@echo Compiling $<
	@mkdir -p $(@D)
	@$(CPP) -MMD -c $(CPPFLAGS) $(includes) $< -o $@

# compile c
$(c_objects): $(OBJDIR)/%.o: %.c
	@echo Compiling $<
	@mkdir -p $(@D)
	@$(CC) -MMD -c $(CFLAGS) $(includes) $< -o $@


.PHONY : clean
clean: 
	rm $(OBJDIR)/*.o
	rm FactorHiGHS/*.o
	rm $(binary_name)
	rm $(OBJDIR)/*.d
	rm FactorHiGHS/*.d

