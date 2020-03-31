One Way to run this

make clean; make -j8 mo-vrp;

./mo-vrp config/parameter-1 data/small-distributed



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
