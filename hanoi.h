#include <iostream>
#include <cmath>
#include <vector>
#include <random>

using namespace std;

/* Initialize RNG. */
random_device rd;
mt19937_64 gen(rd());
uniform_int_distribution<int> intDist;
uniform_real_distribution<double> doubleDist(0.0,1.0);

class Hanoi {
private:
	int * lattice;
	vector<int> stack;
	int length;
public: 
	Hanoi( size_t );
	~Hanoi();
	void setLattice( int );
	int getSpin( int );
	void pushIndexToStack( int );
	int HN3longRange( int );
	void displayLattice();
	void displayStack();
	int getRandomIndex();
};

Hanoi::Hanoi( size_t len ){
	lattice = new int[len];
	length = len;
	for ( int i=0; i<length; i++ )
		lattice[i]=1;
}

Hanoi::~Hanoi() {
	delete [] lattice;
}

void Hanoi::setLattice(int x){
	for ( int i=0; i<length; i++ )
		lattice[i]=x;
}

int Hanoi::getSpin( int index ){
	return lattice[index];
}

void Hanoi::pushIndexToStack( int index ){
	stack.push_back( index );
}

int Hanoi::HN3longRange( int index) {
	int temp  = index;
	int hierarchy = 1;
	
	// calculates the small-world hierarchy level 
	if ( index != 0 ) {
		while( temp % 2 == 0 ) {
			hierarchy++;
			temp /= 2;
		}
	
	// determines forward or backward connectivity to 
	// small-world neighbor using formula: n = [2(j)+1]
	// if j%2==0 sw bond is foward, else backwards
	if( ((temp-1)/2) % 2 == 0 )
		return (index + pow(2,hierarchy));
	else
		return (index - pow(2,hierarchy));
	}
	else return 0;
}

void Hanoi::displayLattice() {
	for ( int i=0; i < length; i++ )
		cout << lattice[i] << '\t';
	cout << endl;
}

void Hanoi::displayStack() {
	for ( vector<int>::iterator it = stack.begin(); it != stack.end(); ++it )
		cout << *it << '\t';
	cout << endl;
}

int Hanoi::getRandomIndex() {
	return intDist(gen) % (length+1);
}