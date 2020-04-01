CPP= g++ --std=c++11
CPPFLAGS= -g -O3 -march=native -fopenmp
NVCC= nvcc
NVCCFLAGS= -G -std=c++11 -gencode arch=compute_50,code=sm_50
LCUDAFLAGS= -I/usr/local/cuda/include -L/usr/local/cuda/lib64  -lcudart -lcuda

all: main

mo-vrp: crossover_mutation.o parent_selection.o ga_lib.o helper_lib.o init_population.o mo-vrp.o
	$(CPP) $(CPPFLAGS) $? -o $@ ${LCUDAFLAGS}

crossover_mutation.o: crossover_mutation.cu
	$(NVCC) $(NVCCFLAGS) -c $? -o $@

ga_lib.o: ga_lib.cu
	$(NVCC) $(NVCCFLAGS) -c $? -o $@

helper_lib.o: helper_lib.cu
	$(NVCC) $(NVCCFLAGS) -c $? -o $@

init_population.o: init_population.cu
	$(NVCC) $(NVCCFLAGS) -c $? -o $@

mo-vrp.o: mo-vrp.cu
	$(NVCC) $(NVCCFLAGS) -c $? -o $@

parent_selection.o: parent_selection.cu
	$(NVCC) $(NVCCFLAGS) -c $? -o $@




clean:
	rm -rf *.o movrp
