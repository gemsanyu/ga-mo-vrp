#include <iostream>
#include <conio.h>
#include <cstdlib>
#include <ctime>

using namespace std;
/* run this program using the console pauser or add your own getch, system("pause") or input loop */

int main(int argc, char** argv) {
	
	float select [5], prob [10], mat [5];
	int c, r, z;
	cout.precision (2);
	cout << "Random Number:" << endl;	
	srand(time(0));
	
	for (c =0; c<5; c++){
		select [c] = (double) rand()/RAND_MAX;
		cout << select [c] << " ";
	}
	cout << endl << "\nIncreasing Probabilities / Fitness Function:" << endl;
	
	for (c = 0; c< 10; c++){
		prob [c] = (double) (c+1)/10;
		cout << prob [c] << " ";
	}
	cout << endl;
	r = 0;
	c = 0;
	z = 0;
	
	while (z<5){
		c=0;
		while(c<10){
			if (select[z]-prob[c]<0){
				mat [z]=c;
				z++;
			}
		c++;
		}
	}
	
	cout << "\nSelected Column:" << endl;
	for (c =0; c<5; c++){
		cout << mat [c] << endl;
	}
		
	return 0;
}
