CPP= g++ --std=c++11
CPPFLAGS= -g -O3 -march=native -fopenmp
NVCC= nvcc
NVCCFLAGS= -G -std=c++11
LCUDAFLAGS= -I/usr/local/cuda/include -L/usr/local/cuda/lib64  -lcudart -lcuda

all: main

mo-vrp: helper_lib.o mo-vrp.o
	$(CPP) $(CPPFLAGS) $? -o $@ ${LCUDAFLAGS}

mo-vrp.o: mo-vrp.cu
	$(NVCC) $(NVCCFLAGS) -c $? -o $@

helper_lib.o: helper_lib.cu
	$(NVCC) $(NVCCFLAGS) -c $? -o $@

clean:
	rm -rf *.o movrp
