
# paths
HIGHS_PATH = $(HOME)/Documents/HiGHS
METIS_PATH = $(HOME)/Documents/METIS
LOCAL_PATH = $(HOME)/local

# source files
cpp_sources = \
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
		../FactorHiGHS/PackedPackedFormatHandler.cpp \
		../FactorHiGHS/DataCollector.cpp \
		../FactorHiGHS/DenseFact.cpp \
		../FactorHiGHS/CallAndTimeBlas.cpp

# binary file name
binary_name = ipm

# object files directory
OBJDIR = obj

# compilers
#CPP=clang++
 CPP = /opt/homebrew/Cellar/llvm/17.0.6_1/bin/clang++

# compiler flags
CPPFLAGS = -std=c++11 -O3 -g3 -Wno-deprecated #-fsanitize=address #ASAN_OPTIONS=detect_leaks=1

# mess to link openmp on mac
#OPENMP_FLAGS = -Xclang -fopenmp -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib -lomp

# rpaths for shared libraries
# rpaths = -rpath $(CHOLMOD_PATH)/lib/ -rpath $(LOCAL_PATH)/lib/

# includes and libraries
includes = -I$(HIGHS_PATH)/build -I$(HIGHS_PATH)/src/ -I$(METIS_PATH)/include -I$(LOCAL_PATH)/include
libs_path = -L$(HIGHS_PATH)/build/lib -L$(METIS_PATH)/build/libmetis -L$(LOCAL_PATH)/lib
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
-include $(OBJDIR)/mainIPM.d
-include $(OBJDIR)/testNetlib.d

# compile cpp
# (include Makefile as a prerequisite, to force recompilation if Makefile changes, e.g. if flags have changed)
$(cpp_objects): $(OBJDIR)/%.o: %.cpp Makefile
	@echo Compiling $<
	@mkdir -p $(@D)
	@$(CPP) -MMD -c $(CPPFLAGS) $(includes) $< -o $@

$(OBJDIR)/mainIPM.o: mainIPM.cpp Makefile
	@echo Compiling $<
	@$(CPP) -MMD -c $(CPPFLAGS) $(includes) $< -o $@

$(OBJDIR)/masterNetlib.o: masterNetlib.cpp Makefile
	@echo Compiling $<
	@$(CPP) -MMD -c $(CPPFLAGS) $(includes) $< -o $@

test: $(cpp_objects) $(OBJDIR)/masterNetlib.o
	@echo Linking objects into $@
	@$(CPP) $(CPPFLAGS) $(libs_path) $(libs) $^ -o $@

.PHONY : clean
clean: 
	rm $(OBJDIR)/*.o
	rm FactorHiGHS/*.o
	rm $(binary_name)
	rm test
	rm $(OBJDIR)/*.d
	rm FactorHiGHS/*.d

