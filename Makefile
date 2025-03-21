
# paths
HIGHS_PATH = $(HOME)/Documents/HiGHS
METIS_PATH = $(HOME)/Documents/METIS
LOCAL_PATH = $(HOME)/local
BLAS_PATH  = /Library/Developer/CommandLineTools/SDKs/MacOSX14.4.sdk/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers
OPENBLAS_PATH = /opt/homebrew/Cellar/openblas/0.3.26

# source files
cpp_sources = \
		Ipm.cpp \
		IpmModel.cpp \
		VectorOperations.cpp \
		FactorHiGHSSolver.cpp \
		CurtisReidScaling.cpp \
		IpmIterate.cpp \
		../FactorHiGHS/Analyse.cpp \
		../FactorHiGHS/Auxiliary.cpp \
		../FactorHiGHS/Factorise.cpp \
		../FactorHiGHS/Numeric.cpp \
		../FactorHiGHS/Symbolic.cpp \
		../FactorHiGHS/FormatHandler.cpp \
		../FactorHiGHS/FullFormatHandler.cpp \
		../FactorHiGHS/HybridPackedFormatHandler.cpp \
		../FactorHiGHS/HybridHybridFormatHandler.cpp \
		../FactorHiGHS/PackedPackedFormatHandler.cpp \
		../FactorHiGHS/SolveHandler.cpp \
		../FactorHiGHS/FullSolveHandler.cpp \
		../FactorHiGHS/PackedSolveHandler.cpp \
		../FactorHiGHS/HybridSolveHandler.cpp \
		../FactorHiGHS/DataCollector.cpp \
		../FactorHiGHS/DenseFactKernel.cpp \
		../FactorHiGHS/DenseFactFull.cpp \
		../FactorHiGHS/DenseFactHybrid.cpp \
		../FactorHiGHS/CallAndTimeBlas.cpp \
		../FactorHiGHS/KrylovMethods.cpp \
		../FactorHiGHS/SymScaling.cpp \
		../FactorHiGHS/DgemmParallel.cpp

# binary file name
binary_name = ipm

test_name = masterNetlib
#test_name = ../FactorHiGHS/test

# object files directory
OBJDIR = obj

# compilers
#CPP=clang++
 CPP = /opt/homebrew/Cellar/llvm/17.0.6_1/bin/clang++

# compiler flags
CPPFLAGS = -std=c++11 -O3 -g3 -Wno-deprecated #-fsanitize=address #ASAN_OPTIONS=detect_leaks=1

# mess to link openmp on mac
#OPENMP_FLAGS = -Xclang -fopenmp -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib -lomp

# includes and libraries
includes = -I$(HIGHS_PATH)/build -I$(HIGHS_PATH)/src/ -I$(METIS_PATH)/include -I$(LOCAL_PATH)/include -I$(BLAS_PATH)
libs_path = -L$(HIGHS_PATH)/build/lib -L$(METIS_PATH)/build/libmetis -L$(LOCAL_PATH)/lib #-L$(OPENBLAS_PATH)/lib
libs = -lhighs -lmetis -lGKlib -lblas

# name of objects
cpp_objects = $(cpp_sources:%.cpp=$(OBJDIR)/%.o)

# dependency files
dep = $(cpp_sources:%.cpp=$(OBJDIR)/%.d)




# link ipm
$(binary_name): $(cpp_objects) $(OBJDIR)/mainIPM.o
	@echo Linking objects into $@
	@$(CPP) $(CPPFLAGS) $(libs_path) $(libs) $^ -o $@

# manage dependencies
-include $(dep)

# compile cpp
# (include Makefile as a prerequisite, to force recompilation if Makefile changes, e.g. if flags have changed)
$(cpp_objects): $(OBJDIR)/%.o: %.cpp Makefile
	@echo Compiling $<
	@mkdir -p $(@D)
	@$(CPP) -MMD -c $(CPPFLAGS) $(includes) $< -o $@

$(OBJDIR)/mainIPM.o: mainIPM.cpp Makefile
	@echo Compiling $<
	@$(CPP) -MMD -c $(CPPFLAGS) $(includes) $< -o $@

# test script
$(OBJDIR)/$(test_name).o: $(test_name).cpp Makefile
	@echo Compiling $<
	@$(CPP) -MMD -c $(CPPFLAGS) $(includes) $< -o $@

test: $(cpp_objects) $(OBJDIR)/$(test_name).o
	@echo Linking objects into $@
	@$(CPP) $(CPPFLAGS) $(libs_path) $(libs) $^ -o $@

.PHONY : clean
clean: 
	rm $(OBJDIR)/*.d
	rm FactorHiGHS/*.d
	rm $(OBJDIR)/*.o
	rm FactorHiGHS/*.o
	rm $(binary_name)
	rm test

