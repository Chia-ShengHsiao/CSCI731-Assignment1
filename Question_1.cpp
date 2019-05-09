// ID : 23399066
// Name : Chia-Sheng Hsiao
// Quetion 1 in Final Exam 
#include<iostream>
#include<cmath> 

using namespace std;

long double r(double ri, double rj, double t, double ti, double tj) {
  return (0.5 * ri * ti + 0.5 * rj * tj)/ t;
}

double compute(double y, double t, double c) {
	double temp = (1.0 + y/4);
	double temp2 = pow(temp, t*4.0);
	return (c/4.0) / temp2;
}

double compute_2(double F, double y) {
	double temp = (1.0 + y/4);
	double temp2 = pow(temp, 8);
	return F/temp2;
}

double sum(double y) {
	//my ID is 23399066
	return compute(y, 0.25,2.0) + compute(y, 0.5,3.0) + compute(y, 0.75,3.0) + compute(y, 1.0,9.0) + 
	compute(y, 1.25,9.0) + compute(y, 1.50,0.0) + compute(y, 1.75,6.0) + compute(y, 2.0,6.0) + 
	compute_2(100.0, y);
}

int main(int argc, char** argv) {
	double F = 100.0;
	
	cout << "*****coompute d & r*****" << endl;
	double d_050 = exp(-0.5*0.0099751);
	cout << "d_050 = " << d_050 << endl;
	double d_100 = exp(-1.0*0.0161980);
	cout << "d_100 = " << d_100 << endl;
	double d_150 = exp(-1.5*0.0198461);
	cout << "d_150 = " << d_150 << endl;
	double d_200 = exp(-2.0*0.0224432);
	cout << "d_200 = " << d_200 << endl << endl;
	
	double r_025 = r(0.0, 0.99751, 0.25, 0.0, 0.5)/100;
	cout << "r_025 = " << r_025 << endl;
	double r_075 = r(0.99751, 1.61980, 0.75, 0.5, 1.0)/100.0;
	cout << "r_075 = " << r_075 << endl;
	double r_125 = r(1.61980, 1.98461, 1.25, 1.0, 1.5)/100.0;
	cout << "r_125 = " << r_125 << endl;
	double r_175 = r(1.98461, 2.24432, 1.75, 1.5, 2.0)/100.0;
	cout << "r_175 = " << r_175 << endl << endl;
	
	double d_025 = exp( -0.25* r_025 );
	cout << "d_025 = " << d_025 << endl;
	double d_075 = exp( -0.75* r_075 );
	cout << "d_075 = " << d_075 << endl;
	double d_125 = exp( -1.25* r_125 );
	cout << "d_125 = " << d_125 << endl;
	double d_175 = exp( -1.75* r_175 );
	cout << "d_175 = " << d_175 << endl;
	
	//in formula (1.1)
	double BFV = (2.0/4.0)*d_025 + (3.0/4.0)*d_050 + (3.0/4.0)*d_075 + (9.0/4.0)*d_100 + (9.0/4.0)*d_125 + (0.0/4.0)*d_150 + (6.0/4.0)*d_175 + (F+(6.0/4.0))*d_200;
	cout << endl << "In formula (1.1), BFV = " << BFV << endl;
	
	//in formula (1.2)
	double y=0.01;
	while( sum(y)-BFV>0.0000001) {
		y = y+0.00000001;
	}
	cout << "in foluma (1.2), y = " << y << endl;
	
	//in formula (1.3)
	double temp = d_025/4.0 + d_050/4.0 + d_075/4.0 + d_100/4.0 + d_125/4.0 + d_150/4.0 + d_175/4.0 + d_200/4.0;
	double c = (BFV-F*d_200)/ temp;
	cout <<  "in foluma (1.3), c = " << c << endl;
		
	return 0;
}
