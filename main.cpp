#include <iostream>
#include "hanoi.h"

using namespace std;

int main()
{
	Hanoi x( 33 );
	for ( int i = 0; i < 33; ++i ) {
		cout << i << '\t' << x.HN3longRange( i ) << endl;	
	}
	cout << endl;
	
	x.pushIndexToStack(4);
	x.pushIndexToStack(10);
	x.pushIndexToStack(17);
	x.displayStack();

	for ( int i = 0; i < 100; ++i )
		cout << x.getRandomIndex() << '\t';
	cout << endl;

	return 0;
}