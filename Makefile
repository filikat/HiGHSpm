
# paths
HIGHS_PATH = $(HOME)/Documents/HiGHS
CHOLMOD_PATH = $(HOME)/Documents/SuiteSparse
QDLDL_PATH = $(HOME)/Documents/qdldl
METIS_PATH = $(HOME)/Documents/METIS
LOCAL_PATH = $(HOME)/local

# source files
cpp_sources = \
		mainIPM.cpp \
		IPM_caller.cpp \
		IPM_model.cpp \
		NormalEquations.cpp \
		ConjugateGradient.cpp \
		VectorOperations.cpp \
		IPM_aux.cpp \
		Direct.cpp \
		ExperimentData.cpp \
		Metis_caller.cpp \
		VertexCover.cpp \
		Lapack_wrapper.cpp \
		SchurContribution.cpp
c_sources = hsl_ma86_wrapper.c

# binary file name
binary_name = ipm

# object files directory
OBJDIR = obj

# compilers
CC=clang
CPP=clang++

# compiler flags
CPPFLAGS = -std=c++11 -O3 -DHAVE_CHOLMOD -DHAVE_MA86 -DHAVE_MC68
CFLAGS = -O3

# mess to link openmp on mac
OPENMP_FLAGS = -Xclang -fopenmp -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib -lomp

# rpaths for shared libraries
rpaths = -rpath $(CHOLMOD_PATH)/lib/ -rpath $(QDLDL_PATH)/build/out/ -rpath $(LOCAL_PATH)/lib/

# includes and libraries
includes = -I$(HIGHS_PATH)/build -I$(HIGHS_PATH)/src/ -I$(CHOLMOD_PATH)/include -I$(QDLDL_PATH)/include -I$(METIS_PATH)/include -I$(LOCAL_PATH)/include
libs_path = -L$(HIGHS_PATH)/build/lib -L$(CHOLMOD_PATH)/lib -L$(QDLDL_PATH)/build/out -L$(METIS_PATH)/build/libmetis -L$(LOCAL_PATH)/lib
libs = -lhighs -lcholmod -lqdldl -lmetis -lhsl_ma86 -lhsl_mc68 -lfakemetis -lGKlib -llapack

# name of objects
cpp_objects = $(cpp_sources:%.cpp=$(OBJDIR)/%.o)
c_objects = $(c_sources:%.c=$(OBJDIR)/%.o)

# dependency files
dep = $(cpp_sources:%.cpp=$(OBJDIR)/%.d)




# link ipm
$(binary_name): $(cpp_objects) $(c_objects)
	@echo Linking objects into $@
	@$(CPP) $(CPPFLAGS) $(OPENMP_FLAGS) $(libs_path) $(libs) $(rpaths) $^ -o $@

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
	rm $(binary_name)
	rm $(OBJDIR)/*.d

