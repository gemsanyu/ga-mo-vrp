One Way to run this

g++ --std=c++11 -O3 -march=native ga_lib.cpp helper_lib.cpp mo-vrp.cpp -o mo-vrp; ./mo-vrp config/small-distributed-1


Config Structure

small-distributed -> dataset file name
50 -> n of customer
200 -> maximum distance
30 -> maximum capacity
5000 -> maximum iteration
0.000000001 -> fitness difference threshold
50 -> Population Size
10 -> Number of parent chosen (roulette wheel count)
0.8 -> PC
0.03 -> PM
