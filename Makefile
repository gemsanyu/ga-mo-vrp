CPP= g++ --std=c++11
CPPFLAGS= -g -O3 -march=native -fopenmp
NVCC= nvcc
NVCCFLAGS= -g -G -std=c++11
LCUDAFLAGS= -I/usr/local/cuda/include -L/usr/local/cuda/lib64  -lcudart -lcuda

# all: mo-vrp

mo-vrp: helper_lib.o ga_lib.o nsga2.o mo-vrp.o
	$(CPP) $(CPPFLAGS) $? -o $@

mo-vrp.o: mo-vrp.cpp
	$(CPP) $(CPPFLAGS) -c $? -o $@

# crossover_mutation.o: crossover_mutation.cu
	# $(NVCC) $(NVCCFLAGS) -c $? -o $@

nsga2.o: nsga2.cpp
	$(CPP) $(CPPFLAGS) -c $? -o $@

ga_lib.o: ga_lib.cpp
	$(CPP) $(CPPFLAGS) -c $? -o $@

helper_lib.o: helper_lib.cpp
	$(CPP) $(CPPFLAGS) -c $? -o $@

clean:
	rm -rf *.o movrp
