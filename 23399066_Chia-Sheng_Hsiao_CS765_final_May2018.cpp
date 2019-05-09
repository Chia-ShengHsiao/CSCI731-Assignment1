// Name : Chia-Sheng Hsiao
// ID : 23399066

#include<iostream>
#include<algorithm>
#include<fstream>
#include<vector>
#include <iomanip>
#include<cmath> 

using namespace std;

//Database class
class Database {
	public:
		//Data
		double r = 0.0;
		double q = 0.0;
		
		//Method
		Database() {
			//Constructor
		}
		
		~Database() {}
};

//TreeNode class
class TreeNode {
	public:
		//Data
		double S = 0.0;
		double V = 0.0;
		double t = 0.0;
		
		//Method
		TreeNode() {
			//Constructor
		}
		
		~TreeNode() {}
};

//Derivative class
class Derivative {
	public:
		//Data
		double T = 0.0;
		
		//Method
		virtual double TerminalPayoff(double S) const {return 0;}
		
		virtual int ValuationTests(TreeNode &node) const {return 0;}
		
		virtual ~Derivative() {}
	
	protected:
		
		Derivative() {
		}
};

//BinomialTree class
class BinomialTree {
	public:
		
		BinomialTree(int n);
		
		~BinomialTree();
		
		int FairValue(int n, const Derivative * p_derivative, const Database * p_db, double S, double sigma, double t0, double &FV );
		
		int ImpliedVolatility(int n, const Derivative * p_derivative, const Database * p_db, double S, double t0, double target, double & implied_vol, int & num_iter);
		
	private:
		//Data
		int n_tree;
		TreeNode **tree_nodes;
		
		//methods
		void Clear();
		
		int Allocate(int n);
};

BinomialTree :: BinomialTree(int n) {
	n_tree = 0;
	tree_nodes = NULL;
	Allocate(n);
}

BinomialTree :: ~BinomialTree() {
	Clear();
}

void BinomialTree :: Clear() {
	if(tree_nodes == NULL) return;
	
	for(int i=0; i<n_tree; i++) {
		delete[] tree_nodes[i];
	}
	
	delete[] tree_nodes;
	tree_nodes = NULL;
}

int BinomialTree :: Allocate(int n) {
	if(n < n_tree) return 0; 
	
	//deallocate old tree
	Clear();
	
	//allocate memory
	n_tree = n;
	tree_nodes = new TreeNode*[n + 1];
	
	for(int i=0; i<=n_tree; i++) {
		tree_nodes[i] = new TreeNode[n + 1];
		for(int j=0; j<=n_tree; j++) {
			TreeNode temp_node;
			tree_nodes[i][j] = temp_node;
		}
	}
	
	return 0;
}


int BinomialTree :: FairValue(int n, const Derivative * p_derivative, const Database * p_db, double S, double sigma, double t0, double &FV) {
	
	FV = 0;
	
	//Validate the input data
	if(n<1 || S<=0 || p_derivative==NULL || p_db==NULL || p_derivative->T<=t0 || sigma<=0) return 1;
	
	double dt = (p_derivative->T - t0)/ double(n);
	double df = exp(-1 * p_db->r * dt);
	double growth = exp((p_db->r - p_db->q) * dt);
	double u = exp(sigma * sqrt(dt));
	double d = 1.0/u;
	
	double p_prob = (growth-d)/(u-d);
	double q_prob = 1.0-p_prob;
	
	if(p_prob<0.0 || p_prob>1.0) return 1;
	
	// Allocate(n) to allocate memory for the binomial tree
	Allocate(n);
	
	//set up the stock price
	TreeNode * node_tmp = tree_nodes[0];
	node_tmp[0].S = S;
	node_tmp[0].t = t0;
	
	//set Stock price
	for(int i=1; i<=n; i++) {
		double t = t0 + i*dt;
		TreeNode* prev = tree_nodes[i-1];
		node_tmp = tree_nodes[i];
        node_tmp[0].S = prev[0].S * d;
        node_tmp[0].t = t;
		
		for(int j=1; j<=n; j++) {
			node_tmp[j].S = node_tmp[j-1].S * u * u;
			node_tmp[j].t = t;
		}
	}
	
	//set terminal payoff
	int i = n;
	node_tmp = tree_nodes[i];
	for(int j=0; j<=n; j++) {
		node_tmp[j].V = p_derivative->TerminalPayoff(node_tmp[j].S);
	}
	
	//valuation loop
	for(int i=n-1; i>=0; i--) {
		node_tmp = tree_nodes[i];
		TreeNode* node_next = tree_nodes[i+1];
		for(int j=0; j<=i; j++) {
			node_tmp[j].V = df*(p_prob*node_next[j+1].V + q_prob*node_next[j].V );
			p_derivative->ValuationTests(node_tmp[j]);
		}
	}
	
	//set Fair Value
	node_tmp = tree_nodes[0];
	FV = node_tmp[0].V;
	
	return 0;
}

int BinomialTree :: ImpliedVolatility(int n, const Derivative * p_derivative, const Database * p_db, double S, double t0, double target, double & implied_vol, int & num_iter) {
	
	//constants
	const double tol = 1.0e-6;
	const int max_iter =  100;
	
	//initialize Referenced variables
	implied_vol = 0;
	num_iter = 0;
	
	//local variables
	double sigma_low = 0.01;
	double sigma_high = 2.0;
	double sigma = 0.0;
	double FV_low = 0.0;
	double FV_high = 0.0;
	double FV = 0.0;
	
	//Validation for sigma_low
	sigma = sigma_low;
	FairValue(n, p_derivative, p_db, S, sigma_low, t0, FV_low);
	double diff_FV_low = FV_low - target;
	if(abs(diff_FV_low) <= tol) {
		implied_vol = sigma;
		return 0;
	}
	
	sigma = sigma_high;
	FairValue(n, p_derivative, p_db, S, sigma_high, t0, FV_high);
	double diff_FV_high = FV_high - target;
	if(abs(diff_FV_high) <= tol) {
		implied_vol = sigma;
		return 0;
	}
	
	//Tests if the target value lies between FV_low and FV_high
	if(diff_FV_low*diff_FV_high>0) {
		implied_vol = 0;
		return 1;
	}
	
	//Iteration bisection loop
	for(num_iter=1; num_iter<max_iter; num_iter++) {
		sigma = 0.5*(sigma_low+sigma_high);
		FairValue(n, p_derivative, p_db, S, sigma, t0, FV);
		double diff_FV = FV - target;
		
		if(abs(diff_FV) < tol) {
			implied_vol = sigma;
			return 0;
		}
		else {
			if((diff_FV_low*diff_FV) > 0) {
				sigma_low = sigma;
			}
			else {
				sigma_high = sigma;
			}
		}
		
		if(abs(sigma_high-sigma_low) < tol) {
			implied_vol = sigma;
			return 0;
		}
	}
	
	implied_vol = 0;
	num_iter = max_iter;
	
	return 1;	  
}

//******************************Question 4******************************
//Option Class
class Option : public Derivative {
	public:
		//Data
		double K = 0.0;
		bool isCall = false;
		bool isAmerican = false;
		
		//Method
		Option() {
			//Constructor
		}
		
		virtual ~Option() {}
		
		virtual double TerminalPayoff(double S) const;
		
		virtual int ValuationTests(TreeNode &node) const;	
};

double Option :: TerminalPayoff(double S) const {
	if(isCall && S>K) return S-K;
	if(!isCall && S<K) return K-S;
	return 0;
}

int Option :: ValuationTests(TreeNode & node) const {
	if(isAmerican) {
		if(isCall) {
			node.V = fmax(node.V, fmax(node.S-K, 0));
		}
		else {
			node.V = fmax(node.V, fmax(K-node.S, 0));
		}
	}
	
	return 0;
}

//******************************Question 5******************************
//Straddle Class
class Straddle : public Derivative {
	public:
		//Data
		double K = 0.0;
		bool isCall;
		bool isAmerican = false;
		
		//Method
		Straddle() {
			//Constructor
		}
		
		virtual ~Straddle() {}
		
		virtual double TerminalPayoff(double S) const;
		
		virtual int ValuationTests(TreeNode &node) const;
};

double Straddle :: TerminalPayoff(double S) const {
	return abs(S-K);
}

int Straddle :: ValuationTests(TreeNode &node) const {
	if(isAmerican) {
		double payoff = abs(node.S-K);
		if(node.V<payoff) node.V = payoff;
	}
	
	return 0;
	
}

//******************************Question 6******************************
//BinaryOption class   
class BinaryOption : public Derivative {
	public:
		//Data
		double K = 0.0;
		bool isCall = false;
		bool isAmerican = false;
		
		//Method
		BinaryOption() {
			//Constructor
		}
		
		virtual ~BinaryOption() {}
		
		virtual double TerminalPayoff(double S) const;
		
		virtual int ValuationTests(TreeNode &node) const;
};

double BinaryOption :: TerminalPayoff(double S) const {
	if(isCall && S>=K) return 1;
	if(!isCall && S<K) return 1;
	return 0;
}

int BinaryOption :: ValuationTests(TreeNode &node) const {
	//only conside the European binary call
	return 0;
}

double cum_norm(double x) {
	const double root = sqrt(0.5);
	return 0.5*(1.0 + erf(x*root));
}

double d1_func(double S, double K, double r, double q, double sigma, double T, double t0) {
	double d1 = ((log(S/K) + (r-q)*(T-t0))/(sigma*sqrt(T-t0))) + 0.5*sigma*sqrt(T-t0);
	return d1;
}

double d2_func(double S, double K, double r, double q, double sigma, double T, double t0) {
	return  d1_func(S, K, r, q, sigma, T, t0 ) - sigma*sqrt(T-t0);
}

double c_BSM_binary_func(double r, double T, double t0, double d2) {
	double delta_c = exp(-1*r*(T-t0)) * cum_norm(d2);
	return delta_c;		
}

double p_BSM_binary_func(double r, double T, double t0, double d2) {
		double delta_p = exp(-1*r*(T-t0)) * cum_norm(-1*d2);	
		return delta_p;		
}

//******************************Question 7******************************
//CovertibleBond class
class ConvertibleBond : public Derivative {
	public:
		//Data
		double K = 0.0;
		double B = 0.0; //Threshold
		bool isCall = false;
		bool isAmerican = false;
		
		//Method
		ConvertibleBond() {
			//Constructor
		}
		
		virtual ~ConvertibleBond() {}
		virtual double TerminalPayoff(double S) const;
		virtual int ValuationTests(TreeNode &node) const;
};

double ConvertibleBond :: TerminalPayoff(double S) const {
	if(isCall && S>=K) return S;
	if(isCall && S<K) return K;
	return 0;
}

int ConvertibleBond :: ValuationTests(TreeNode &node) const {
	if(node.S>=B) {
		node.V = node.S;
		return 0;
	}
	
	if(isAmerican) {
		double intrinsic = 0;
		if(isCall) {
			if(node.S >= K ) intrinsic = node.S;
		}
		else intrinsic =  0;
		
		if(intrinsic > node.V) node.V = intrinsic;
	}
	
	return 0;
}

//******************************The test function for Question 4******************************

int Q4_Binomial() {
	int rc = 0;
	
	//setting the initial value in question 4
	double S = 100;
	double K = 100;
	double r = 0.05;
	double q = 0.02;
	double sigma = 0.5;
	double T = 1.0;
	double t0 = 0.0;
	
	Database db;
	db.r = r;
	db.q = q;
	
	//Europe put
	Option Eur_put;
	Eur_put.K = K;
	Eur_put.T = T;
	Eur_put.isCall = false;
	Eur_put.isAmerican = false;
	
	//American put
	Option Am_put;
	Am_put.K = K;
	Am_put.T = T;
	Am_put.isCall = false;
	Am_put.isAmerican = true;
	
	//Europe call
	Option Eur_call;
	Eur_call.K = K;
	Eur_call.T = T;
	Eur_call.isCall = true;
	Eur_call.isAmerican = false;
	
	//American call
	Option Am_call;
	Am_call.K = K;
	Am_call.T = T;
	Am_call.isCall = true;
	Am_call.isAmerican = true;
	
	double FV_Am_put = 0.0;
	double FV_Eur_put = 0.0;
	double FV_Am_call = 0.0;
	double FV_Eur_call = 0.0;
	
	S = 100;
	int n = 100;
	BinomialTree binom(n);
	rc = binom.FairValue(n, &Am_put, &db, S, sigma, t0, FV_Am_put);
	rc = binom.FairValue(n, &Eur_put, &db, S, sigma, t0, FV_Eur_put);
	rc = binom.FairValue(n, &Am_call, &db, S, sigma, t0, FV_Am_call);
	rc = binom.FairValue(n, &Eur_call, &db, S, sigma, t0, FV_Eur_call);
	
	double EP = FV_Eur_put;
	double AP = FV_Am_put;
	double EC = FV_Eur_call;
	double AC = FV_Am_call;
	
	cout << "****************************** Question 4 ******************************" << endl << endl;
	cout << "European Put = " << EP << endl << endl;
	cout << "American put = " << AP << endl << endl;
	cout << "European Call = " << EC << endl << endl;
	cout << "American Call = " << AC << endl << endl;
}

//******************************The test function for Question 5******************************

int Q5_Straddle() {
	int rc = 0;
	double r = 0.05;
	double q = 0.02;
	double T = 1.0;
	double t0 = 0.0;
	
	Database db;
	db.r = r;
	db.q = q;

	double S = 90.0;
	double K = 100.0;
	double sigma = 0.1;
	
	int n = 100;
	Straddle Amt;
	Amt.K = K;
	Amt.T = T;
	Amt.isAmerican = true;
	
	Straddle Eur;
	Eur.K = K;
	Eur.T = T;
	Eur.isAmerican = false;
	
	double FV_Amt = 0.0;
	double FV_Eur = 0.0;
	
	BinomialTree binom(n);
	rc = binom.FairValue(n, &Amt, &db, S, sigma, t0, FV_Amt);
	rc = binom.FairValue(n, &Eur, &db, S, sigma, t0, FV_Eur);
	cout << endl << "****************************** Question 5 ******************************" << endl << endl;
	cout << "The fair value of American Straddle :" << FV_Amt << endl << endl;
	cout << "The fair value of European Straddle : " << FV_Eur << endl << endl;
	
	return 0;
}

//******************************The following table for Question 6******************************
int Q6_Following_Table() {
	//initial value
	int rc = 0;
	double r = 0.0223;
	//Depend on Question 1 risk free rate r = 0.0223
	double q = 0.02;
	double T = 1.0;
	double t0 = 0.0;
	
	Database db;
	db.r = r;
	db.q = q;

	double S = 90.0;
	double K = 100.0;
	
	int n = 1000;
	double c_binomial_binary = 0.0;
	double p_binomial_binary = 0.0;
	double c_BSM_binary = 0.0;
	double p_BSM_binary = 0.0;
	
	cout << endl << "****************************** Question 6 Follow Table ******************************" << endl << endl;
	for(int i=1; i<=10; i++) {
		double sigma = i * 0.1;
		
		BinaryOption euroBinomialBinaryCallOption;
	    euroBinomialBinaryCallOption.K= K;
	    euroBinomialBinaryCallOption.T = T;
	    euroBinomialBinaryCallOption.isCall = true;
	    euroBinomialBinaryCallOption.isAmerican = false;
	    
	    
		BinaryOption euroBinomialBinaryPutOption;
	    euroBinomialBinaryPutOption.K = K;
	    euroBinomialBinaryPutOption.T = T;
	    euroBinomialBinaryPutOption.isCall = false;
	    euroBinomialBinaryPutOption.isAmerican = false;	
		
		BinomialTree binom(n);
	    rc = binom.FairValue(n, &euroBinomialBinaryCallOption, &db, S, sigma, t0, c_binomial_binary);	
	    rc = binom.FairValue(n, &euroBinomialBinaryPutOption, &db, S, sigma, t0, p_binomial_binary);	
	    
	    double d1 = d1_func(S, K, r, q, sigma, T, t0);
	    double d2 = d2_func(S, K, r, q, sigma, T, t0);
	    double c_BSM_binary = c_BSM_binary_func(r, T, t0, d2);
	    double p_BSM_binary = p_BSM_binary_func(r, T, t0, d2);
	    
	    cout << "Sigma = " << sigma;
		cout << "; C_Binomial_Binary = " << setprecision(3) << c_binomial_binary;
		cout << "; P_Binomial_Binary = " << setprecision(3) << p_binomial_binary;
	    cout << "; C_BSM_Binary = " << setprecision(3) << c_BSM_binary;
		cout << "; P_BSM_Binary = " << setprecision(3) << p_BSM_binary << endl << endl; 
	}
	
	return 0;
}

//******************************The Plot Graph Data for Question 6******************************
int Q6_Graph_Data() {
	
	int rc = 0;
	
	//create 4 output file & starting to writing
	ofstream output_call_binomial("Q6 Call Binomianl Model Data.txt");
	ofstream output_call_BSM("Q6 Call BSM Formula Data.txt");
	ofstream output_put_binomial("Q6 Put Binomial Model Data.txt");
	ofstream output_put_BSM("Q6 Put BSM Formula Data.txt");
	
	//initial value
	double r = 0.0223;
	//Depend on Question 1 risk free rate r = 0.0223
	double q = 0.02;
	double T = 1.0;
	double t0 = 0.0;
	
	Database db;
	db.r = r;
	db.q = q;

	double S = 90.0;
	double K = 100.0;
	
	int n = 1000;
	double c_binomial_binary = 0.0;
	double p_binomial_binary = 0.0;
	double c_BSM_binary = 0.0;
	double p_BSM_binary = 0.0;
	for(int i=1; i<=100; i++){
		double sigma = i*0.01;
		
		BinaryOption euroBinomialBinaryCallOption;
	    euroBinomialBinaryCallOption.K= K;
	    euroBinomialBinaryCallOption.T = T;
	    euroBinomialBinaryCallOption.isCall = true;
	    euroBinomialBinaryCallOption.isAmerican = false;
	    
	    
		BinaryOption euroBinomialBinaryPutOption;
	    euroBinomialBinaryPutOption.K = K;
	    euroBinomialBinaryPutOption.T = T;
	    euroBinomialBinaryPutOption.isCall = false;
	    euroBinomialBinaryPutOption.isAmerican = false;	
		
		BinomialTree binom(n);
	    rc = binom.FairValue(n, &euroBinomialBinaryCallOption, &db, S, sigma, t0, c_binomial_binary);	
	    rc = binom.FairValue(n, &euroBinomialBinaryPutOption, &db, S, sigma, t0, p_binomial_binary);	
	    
	    double d1 = d1_func(S, K, r, q, sigma, T, t0);
	    double d2 = d2_func(S, K, r, q, sigma, T, t0);
	    double c_BSM_binary = c_BSM_binary_func(r, T, t0, d2);
	    double p_BSM_binary = p_BSM_binary_func(r, T, t0, d2);
	    
	    output_call_binomial << std::setw(6) << sigma << " ";
	    output_call_binomial << std::setw(16) << c_binomial_binary << " " << endl;
	    output_call_BSM << std::setw(6) << sigma << " ";
	    output_call_BSM << std::setw(16) << c_BSM_binary << " " << endl;
	    output_put_binomial << std::setw(6) << sigma << " ";
	    output_put_binomial << std::setw(16) << p_binomial_binary << " " << endl;
	    output_put_BSM << std::setw(6) << sigma <<" ";
	    output_put_BSM << std::setw(16) << p_BSM_binary << " " << endl;
	}
	
	//close all file
	output_call_binomial.close();
	output_call_BSM.close();
	output_put_binomial.close();
	output_put_BSM.close();
	return 0;
}

//******************************The Plot Graph Data for Question 7******************************

int Q7_Graph_Data() {
	
	//create 1 output file & starting to writing
	ofstream output_5("Q7 Conver Bond Plot Graph Data.txt");
	
	//initial value
	int rc = 0;
	//fair value of the convertible 
	double U = 0.0;
	//market price of the con
	double M = 0.0;
	
	double r = 0.05;
	double q = 0.0;
	double T = 5.0;
	double t0 = 0.0;
	
	Database db;
	db.r = r;
	db.q = q;
	
	double S = 0.0;
	double K = 100.0;
	double sigma = 0.5;
	
	double B = 130.0;
	int n = 1000;

	
	ConvertibleBond callableAmericanConvertibleBond;
	callableAmericanConvertibleBond.K = K;
	callableAmericanConvertibleBond.T = T;
	callableAmericanConvertibleBond.B = B;
	callableAmericanConvertibleBond.isCall = true;
	callableAmericanConvertibleBond.isAmerican = true;
    
    double FV_convertible_bond = 0.0;
	BinomialTree binom(n);
	
	cout << endl << "****************************** Question 7 Plot Graph Data ******************************" << endl << endl;
	cout<< "sigma = 0.5" << endl << endl ;
	
	for(int i=1 ; i<=150; i++) {
		S = i*1.00;
		rc = binom.FairValue(n, &callableAmericanConvertibleBond, &db, S, sigma, t0, FV_convertible_bond);
		cout << std::setw(6) << S << " ";
		cout << std::setw(16) << setprecision(6) << FV_convertible_bond << " " << endl;
		output_5 << std::setw(6) << S << " ";
		output_5 << std::setw(16) << FV_convertible_bond << " " << endl;
	}
	
	//close file
	output_5.close();
	
	return 0;
}

//******************************The Gamma Trading for Question 7******************************

int Q7_Gamma_Trading() {
	
	cout << endl << "****************************** Question 7 Gamma Trading ******************************" << endl << endl;
	
	int rc = 0;
	//fair value of the convertible 
	double U = 0.0;
	//market price of the con
	double M = 0.;
	
	//initial value
	double r = 0.05;
	double q = 0.0;
	double T = 5.0;
	double t0 = 0.0;
	
	Database db;
	db.r = r;
	db.q = q;
	
	double S = 0.0;
	double K = 100.;
		
	double B = 130.0;
	int n = 1000;
	
	// The Gamma Trading (page12)
	//My ID is 23399066
    double digit_S1 = 2339.0/10000.0;
    double digit_S2 = -9066.0/10000.0;
	
	// 1. Day 0 (page13)
	t0 = 0.0;
	double S0 = 60.0;
	//(c) The barrier option market price is M0 = 90
	double M0 = 90.0; 
	double target = M0;
	ConvertibleBond convertibleBondDay0;
	convertibleBondDay0.K = K;
	convertibleBondDay0.T = T;
	convertibleBondDay0.B = B;
	convertibleBondDay0.isCall = true;
	convertibleBondDay0.isAmerican = true;
	
	double implied_volatility = 0.0;
	int num_iter = 0;
	BinomialTree binomDay0(n);
	
	//(d) Calculate the implied volatility of the convertible bond (4 decimal places).
	rc = binomDay0.ImpliedVolatility(n, &convertibleBondDay0, &db, S0, t0, target,implied_volatility, num_iter);
	cout << "implied volatility (Sigma_0) = " << setprecision(4) << implied_volatility << endl;
	
	//(f) Calculate the Delta of the convertible bond (4 decimal places)
	double U_S0Plus1 = 0.0;
	double U_S0Minus1 = 0.0;
	binomDay0.FairValue(n, &convertibleBondDay0, &db, S0 + 1, implied_volatility, t0, U_S0Plus1);
	binomDay0.FairValue(n, &convertibleBondDay0, &db, S0 - 1, implied_volatility, t0, U_S0Minus1);
	double delta0_con_bond = (U_S0Plus1 - U_S0Minus1) / 2.0;
	cout << "delta0 con bond =  " << delta0_con_bond << endl;
	
	//(j) Calculate the value of Money0 to 2 decimal places.
	double money0 = delta0_con_bond * S0 - M0;
	cout<< "Money0 on day 0 is "<< setprecision(4) << money0 << endl << endl;
	
	// 2. Day1 (page14)
	t0 = 0.01;
	double M1 = 90.2;
	double S1 = S0 + digit_S1;
	target = M1;
	ConvertibleBond convertibleBondDay1;
	convertibleBondDay1.K = K;
	convertibleBondDay1.T = T;
	convertibleBondDay1.B = B;
	convertibleBondDay1.isCall = true;
	convertibleBondDay1.isAmerican = true;
	
	implied_volatility = 0;
	num_iter = 0;
	BinomialTree binomDay1(n);
	
	//(d) Calculate the implied volatility of the convertible bond (4 decimal places).
	rc = binomDay1.ImpliedVolatility(n, &convertibleBondDay1, &db, S1, t0, target, implied_volatility, num_iter);
	cout << "implied volatility (Sigma_1) = " << setprecision(4) << implied_volatility << endl;
	
	//(f) Calculate the Delta of the convertible bond (4 decimal places)
	double U_S1Plus1 = 0 ;
	double U_S1Minus1 = 0;
	binomDay1.FairValue(n, &convertibleBondDay1, &db, S1 + 1, implied_volatility, t0, U_S1Plus1);
	binomDay1.FairValue(n, &convertibleBondDay1, &db, S1 - 1, implied_volatility, t0, U_S1Minus1);
	double delta1_con_bond = (U_S1Plus1 - U_S1Minus1) /2.0;
	cout << "delta1 con bond =  " << delta1_con_bond << endl;
	
	//(j) Calculate the value of Money1 to 2 decimal places.
	double money1 = money0 + (delta1_con_bond - delta0_con_bond) * S1;
	cout << "Money1 on day 1 is " << setprecision(4) << money1 << endl << endl;;
	
	// 3. Day2 (page15)
	t0 = 0.02;
	double M2 = 90.15;
	double S2 = S1 + digit_S2;
	target = M2;
	ConvertibleBond convertibleBondDay2;
	convertibleBondDay2.K = K;
	convertibleBondDay2.T = T;
	convertibleBondDay2.B = B;
	convertibleBondDay2.isCall = true;
	convertibleBondDay2.isAmerican = true;
	
	implied_volatility = 0;
	num_iter = 0;
	BinomialTree binomDay2(n);
		
	//(e) Calculate the implied volatility of the convertible bond (4 decimal places).
	rc = binomDay2.ImpliedVolatility(n, &convertibleBondDay2, &db, S2, t0, target, implied_volatility, num_iter);
	cout << "implied volatility (Sigma_2) = " << setprecision(4) << implied_volatility << endl;
	
	//(g) Calculate the Delta of the convertible bond (4 decimal places)
	double U_S2Plus1 = 0.0 ;
	double U_S2Minus1 = 0.0;	
	binomDay2.FairValue(n, &convertibleBondDay1, &db, S2 + 1, implied_volatility, t0, U_S2Plus1);
	binomDay2.FairValue(n, &convertibleBondDay1, &db, S2 - 1, implied_volatility, t0, U_S2Minus1);
	double delta2_con_bond = (U_S2Plus1 - U_S2Minus1) /2.0;
	cout << "delta2 con bond =  " << delta2_con_bond << endl;
	
	//(l) Calculate the value of Money2to 2 decimal places.
	double money2 = money1 + (delta2_con_bond - delta1_con_bond)* S2;
	cout << "Money2 on day 2 is " << setprecision(4) << money2 << endl;	
	
	// close out our gamma trading portfolio at the end of day 2. (page16)
	// 5. Calculate the value of the profit to 2 decimal places.
	double Profit = money2 + M2 - delta2_con_bond * S2;

	cout << endl << endl << "The total profit is " << setprecision(2) << Profit << endl << endl;
	
	return 0;
}

int main() {
	
	//test for Question 4
	Q4_Binomial();
	
	//test for Question 5
	Q5_Straddle();
	
	//Following Table for Question 6
	Q6_Following_Table();
	
	//Plot Graph Data for Question 6
	Q6_Graph_Data();
	
	//Plot Graph Data for Question 7
	Q7_Graph_Data();
	
	//Gamma Trading for Question 7
	Q7_Gamma_Trading();
	
	return 0;
	
}


