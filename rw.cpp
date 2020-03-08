int spinRouletteWheel(double *probs, int probSize){ 
	float 	select; 	// Variable Random Number 
	int 	result ; 	// Variable Result
	int c, r, z;	 	//variable batas
	srand(time(0));
	//Random Number
	select = (double) rand()/RAND_MAX;

			// NORMALISASI NILAI FITNESS	
// diketahui fitness 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
//	for (c = 0; c< Probsize; c++){ 
//		probs [c] = (double) (c+1)/10;
//	}
	r = 0;
	c = 0;
	z = 0;
	//code selected colom
	while(c<10){
		if (select - probs[c] < 0){
			result = c;
			z++;
		}
		c++;
	}
	return result;
	
}
