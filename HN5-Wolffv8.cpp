#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <limits>

using namespace std;
const int SIZE = 4097;  // 2^k + 1 spin
const int LIMIT = (SIZE-1) / 2;
const double LOW_T = 4.35;
const double HIGH_T = 4.30;
const double T_INC = 0.6;
const int DATA = 2;
const double PROB = 1.0;
const int DEBUG = 0;

gsl_rng *gBaseRand;       /* global rand number generator */

int main() 
{
  /******************************/
  //     FUNCTION PROTOTYPES    //
  /******************************/
  void initialize( int[ SIZE ] );                                    // creates a lattice with randomized spins
  void uniform( int[ SIZE ] );                                       // creates a lattice with all spins +1
  int ipow( int, int );                                             // same as C++ pow() function
  void print( int[ SIZE ] );                                         // prints lattice to screen
  int print2D( double[ DATA ][ 8 ], int );                           // prints measurements to file
  int printLargestCluster( int[ DATA ][ SIZE ], int );               // prints the largest cluster
  void sweep( int[ SIZE ], double, double * );             // Metropolis algorithm
  double HK( int[ SIZE ], int[ DATA ][ SIZE ], int, double );        // Hoshen-Kopelman algorithm
  double mag( int[ SIZE ] );                                         // calculates the magnetization of the lattice, M
  int min( int, int );                                               // function that returns the small of two numbers
  int bondCountHN5( int );
  int spinEnergy( int *, int );

  /******************************/
  //         VARIABLES          //
  /******************************/
  int lattice[ SIZE ] = { 0 };                  // creates lattice array
  double m;                                     // magnetization
  double e = 0;                                 // internal energy
  double *mptr;                                 // magnetization pointer
  //  double *eptr;                                 // energy pointer
  double T = LOW_T;                             // sets initial temperature at the lowest temperature
  int infinite_loop = 1;                        // used to create infinite loop
  double avg[ DATA ][ 8 ] = { 0 };              // creates average array which contains temperature, < M >, < M^2 >
  double sum[ DATA ][ 8 ] = { 0 };              // creates the sum array which contains the sum of all measurements
  double *sum_ptr;                              // initializes sum_ptr
  int count = 0;                                // sets counter to zero
  int measurements = 0;                         // set measurements to zero
  int clusterCount = 0;                         
  int largestCluster[ DATA ][ SIZE ] = { 0 };
  int *lptr;
  int i;

  //***************RNG initialization and seeding***************//
  unsigned long randSeed;  
  gBaseRand = gsl_rng_alloc(gsl_rng_taus);
  srand(time(NULL));                    /* initialization for rand() */
  randSeed = rand();                    /* returns a non-negative integer */
  gsl_rng_set (gBaseRand, randSeed);    /* seed the PRNG */

  //initialize( lattice );   // creates random lattice
  uniform( lattice );
  m = mag( lattice );      // calculates the magnetization of lattice
  mptr = &m;               // points to magnetization, m
  //  eptr = &e;               // points to energy, e
  lptr = &lattice[ 0 ];    // a pointer to the first entry is passed to spinEnergy function

  cout << "Bond sum is: " << bondCountHN5( SIZE ) << endl;

  // calculates internal energy
  e = 0;
  for ( i = 0; i < SIZE; i++ ) e += spinEnergy( lptr, i );
  e /= 2; // corrects for double-counting

  sum_ptr = &sum[ 0 ][ 0 ];   // pointer points to sum array

  for ( i = 0; i < 100000; i++ ) sweep( lattice, T, mptr );    // equilibrates lattice


  int myCount = 0;

  while( infinite_loop == 1 ) {    // enter infinite loop

    if ( DEBUG ) cout << "T: " << T << endl;

    for ( i = 0; i < 1; i++ ) sweep( lattice, T, mptr );  
    
    ofstream acfFile( "HN5-N4096-acf-.txt", ios::app );
    acfFile << myCount << '\t' << abs(*mptr) << endl;
    myCount++;
    cout << myCount << endl;
    
    // Measurements are made
    e = 0;
    for ( i = 0; i < SIZE; i++ ) e += spinEnergy( lptr, i );
    e /= 2; // corrects for double-counting

    //   cout << "T:" << T << '\n';
    //HK( lattice, largestCluster, clusterCount, T );
    //double bigCluster = 0;
    *sum_ptr += T;  sum_ptr++;                   // assigns the temperature and advances the pointer  
    *sum_ptr += abs(*mptr);  sum_ptr++;          // assigns the magnetization and advances the pointer
    *sum_ptr += pow( *mptr, 2 );  sum_ptr++;     // assigns the square of the magnetization and advances the pointer
    *sum_ptr += pow( *mptr, 4 ); sum_ptr++;      // assigns the fourth power of the magnetization (for Binder parameter)
    double bigCluster = HK( lattice, largestCluster, clusterCount, T );
    *sum_ptr += bigCluster; sum_ptr++; 
    *sum_ptr += ipow( bigCluster, 2 ); sum_ptr++;
    *sum_ptr += e; sum_ptr++;
    *sum_ptr += ipow( e, 2 ); sum_ptr++;
    
    T += T_INC;                                  // increases the temperature
    clusterCount++;

    if ( T > HIGH_T ) {
      // start from the beginning (reset)
      count++;
      T = LOW_T;
      clusterCount = 0;
      //      uniform( lattice );
      // initialize( lattice ); 
      m = mag( lattice );
      
      // this for loop iter doesn't have to be that high
      // for uniform( lattice ) initialization
      //      for( i = 0; i < 0; i++ ) sweep( lattice, T, mptr );       // equilibrates lattice
      sum_ptr = &sum[ 0 ][ 0 ];                                   // pointer points to sum array
      
      // averages are taken every 100 measurements
      if ( count % 100 == 0 ) {
	cout << "count: " << count << endl;
	
	for ( int m = 0; m < DATA; m++ ) {
	  for ( int n = 0; n < 8; n++ ) {
	    avg[ m ][ n ] = sum[ m ][ n ] / count;
	  }
	}
	measurements += 1;              // used for averaging and calculating error
	//	print2D( avg, measurements );   // write to file
	//	printLargestCluster( largestCluster, measurements );
      }
    }
  }
  return 0;
}

void initialize( int lattice[ SIZE ] )  // creates a lattice with randomly assigned spins
{
  double random_num;
  int spin;

  for ( int i = 0; i < SIZE; i++ ) {
      random_num =  gsl_rng_uniform_pos ( gBaseRand );
      spin = (random_num > 0.5 ? 1 : -1);
      lattice[ i ] = spin ;
    }
}

void uniform( int lattice[ SIZE ] )   // creates a lattice with all spins = +1
{ 
  for ( int i = 0; i < SIZE; i++ )
    lattice[ i ] = 1;
}

int ipow(int base, int exp)
{
  int result = 1;
    while (exp)
      {
      if (exp & 1)
	  result *= base;
        exp >>= 1;
        base *= base;
    }
    
    return result;
}

void print( int lattice[ SIZE ] )    // prints the value of the lattice to cout
{
  for ( int i = 0; i < SIZE; i++ ) 
    cout << lattice[ i ] << "\t";

  cout << endl;
}

int print2D( double array[ DATA ][ 8 ], int measurements )      // writes measurements to file
{
  ofstream outFile( "xxx.txt" );
  if (!outFile) {
    cout << "Cannot open file.\n";
    return 1;
  }
  outFile << "Temperature" << "\tMagnetization" << "\tMagSquared" << "\tMag^4" << "\tCluster" << "\tClusterSquared" 
	  << "\tEnergy" << "\tEnergySquared" << "\tMagError" << "\tClusterError" << "\tEnergyError\n";
  for ( int i = 0; i < DATA; i++ ){
    for ( int j = 0; j < 8; j++ ) {
      //outFile << array[ i ][ j ] << '\t';
      outFile.precision(numeric_limits<double>::digits10 + 2);
      outFile << fixed << array[ i ][ j ] << '\t';
    }
    outFile << sqrt( (1/(static_cast<double>(measurements) - 1 ) )*(array[ i ][ 2 ] - ipow(array[ i ][ 1 ],2))) << '\t';
    outFile << sqrt( (1/(static_cast<double>(measurements) - 1 ) )*(array[ i ][ 5 ] - ipow(array[ i ][ 4 ],2))) << '\t';
    outFile << sqrt( (1/(static_cast<double>(measurements) - 1 ) )*(array[ i ][ 7 ] - ipow(array[ i ][ 6 ],2))) << '\n';
  }
  return 0;
}

int printLargestCluster( int largestCluster[ DATA ][ SIZE ], int measurements )      // writes measurements to file
{
  ofstream outFile( "HN5-64-cluster-acf.txt" );
  if (!outFile) {
    cout << "Cannot open file.\n";
    return 1;
  }
  
  for ( int i = 0; i < DATA; i++ ) {
    for ( int j = 0; j < SIZE; j++ ) 
      outFile << static_cast<double>(largestCluster[ i ][ j ]) / static_cast<double>(measurements*100) << "\t";
    outFile << endl;
  }
  return 0;
}

double mag( int lattice[ SIZE ] )         // calculates magnetization, M (*note this is not magnetization per spin, m)
{
  int spin_sum = 0;

  for( int i = 0; i < SIZE; i++ )
    spin_sum += lattice[ i ];

  return spin_sum;
}

void sweep( int lattice[ SIZE ], double temperature, double *mag_ptr )
{

  double padd = 1 - exp( -2 * (1/temperature) );
  double y_padd = padd * PROB;
  int i, sp, oldspin, newspin, current, nn;
  int stack[ SIZE ] = { 0 };
  int temp = 0;
  int count = 0;
  int hier = 1;

  // choose the seed spin for the cluster,
  // put it on the stack, and flip it

  i =  gsl_rng_uniform_int ( gBaseRand, SIZE-1 );

  stack[ 0 ] = i;
  sp = 1;
  oldspin = lattice[ i ];
  newspin = -lattice[ i ];
  lattice[ i ] = newspin;  // flip the spin
  *mag_ptr += (static_cast<double>(newspin*2));

  while (sp) {
    
    count = 0;
    hier = 1;

    // pop a site off the stack 
    current = stack[--sp];

    // check the backbone neighbor forward
    if ( current < SIZE-1 ) {
      nn = current + 1;
      if ( DEBUG ) cout << "current: " << current << "\tneighbor: " << nn << endl;
      if ( lattice[ nn ] == oldspin ) {
	if ( gsl_rng_uniform_pos( gBaseRand ) < padd ) {
	  stack[sp++] = nn;
	  lattice[ nn ] = newspin;
	  *mag_ptr += (static_cast<double>(newspin*2));
	}
      }
    }

    // check the backbone neighbor backward
    if ( current > 0 ) {
      nn = current - 1;
      if ( DEBUG ) cout << "current: " << current << "\tneighbor: " << nn << endl;
      if ( lattice[ nn ] == oldspin ) {
	if ( gsl_rng_uniform_pos( gBaseRand ) < padd ) {
	  stack[sp++] = nn;
	  lattice[ nn ] = newspin;
	  *mag_ptr += (static_cast<double>(newspin*2));
	}
      }
    }

    //  counts the number of hierarchical levels, assigns to hier variable
    temp = current;
    if( current == 0 ) 
      temp = SIZE-1;
    while( temp % 2 == 0 ) { temp /= 2; hier++; }
    
    temp = current;
    if( current == 0 ) 
      temp = SIZE-1;
    while( temp % 2 == 0 ) 
      temp /= 2;
    temp++;
    while( temp % 2 == 0 ) { 
      temp /= 2; 
      count++; 
    }

    // debug specific spin sites and their neighbors
    /* 
       if ( DEBUG )
       if ( current == 3 ) {
       cout << "******************** Here's your problem ********************\n";
       cout << "current: " << current << endl;
       cout << "hier: " << hier << endl;
       cout << "count: " << count << endl;
       cout << "******************** Here's your problem ********************\n";
       }
    */

    // highest hierarchical neighbor
    switch ( current ) { 
    case 0:
      nn = SIZE-1;
      if ( DEBUG ) cout << "case 0: " << current << "\tneighbor: " << SIZE-1 << endl;
      if ( lattice[ nn ] == oldspin ) {
	if ( gsl_rng_uniform_pos( gBaseRand ) < padd ) {
	  stack[sp++] = nn;
	  lattice[ nn ] = newspin;
	  *mag_ptr += (static_cast<double>(newspin*2));
	}	   
      }
      break;
    case SIZE-1:
      nn = 0;
      if ( DEBUG ) cout << "case SIZE-1: " << current << "\tneighbor: " << 0 << endl;
      if ( lattice[ nn ] == oldspin ) {
	if ( gsl_rng_uniform_pos( gBaseRand ) < padd ) {
	  stack[sp++] = nn;
	  lattice[ nn ] = newspin;
	  *mag_ptr += (static_cast<double>(newspin*2));
	}	   
      }
      break;
    case (SIZE-1)/2:
      // do nothing
      break;
    default:
      if ( count == 1 ) {
	nn = (current + ipow( 2, hier )) % (SIZE-1);
	if ( DEBUG ) cout << "current: " << current << "\tneighbor: " << (current + ipow( 2, hier )) % (SIZE-1) << endl;
      }
      else {
	nn = (current - ipow( 2, hier )) % (SIZE-1);
	if ( DEBUG ) cout << "current: " << current << "\tneighbor: " << (current - ipow( 2, hier )) % (SIZE-1) << endl;
      }
      if ( lattice[ nn ] == oldspin ) {
	if ( gsl_rng_uniform_pos( gBaseRand ) < padd ) {
	  stack[sp++] = nn;
	  lattice[ nn ] = newspin;
	  *mag_ptr += (static_cast<double>(newspin*2));
	}	   
      }
      break;
    }

    temp = current;
    hier--;

    switch ( current ) {
    case 0:
      while ( hier > 1 ) {
	nn = ipow( 2, hier-1 );
	if ( DEBUG ) cout << "case 0: " << current << "\tneighbor: " << ipow( 2, hier-1 ) << endl;
	if ( lattice[ nn ] == oldspin ) {
	  if ( gsl_rng_uniform_pos( gBaseRand ) < y_padd ) {
	    stack[sp++] = nn;
	    lattice[ nn ] = newspin;
	    *mag_ptr += (static_cast<double>(newspin*2));
	  }
	}
	hier--;
      }
      break;
    case (SIZE-1):
      while ( hier > 1 ) {
	nn = current - ipow( 2, hier-1 );
	if ( DEBUG ) cout << "case SIZE-1: " << current << "\tneighbor: " << current - ipow( 2, hier-1 ) << endl;
	if ( lattice[ nn ] == oldspin ) {
	  if ( gsl_rng_uniform_pos( gBaseRand ) < y_padd ) {
	    stack[sp++] = nn;
	    lattice[ nn ] = newspin;
	    *mag_ptr += (static_cast<double>(newspin*2));
	  }
	}
	hier--;
      }
      break;
    case (SIZE-1)/2:
      while ( hier > 0 ) {
	nn = current + ipow( 2, hier );
	if ( DEBUG ) cout << "case (SIZE-1)/2: " << current << "\tneighbor: " << ((SIZE-1)/2) + ipow( 2, hier ) << endl;
	if ( lattice[ nn ] == oldspin ) {
	  if ( gsl_rng_uniform_pos( gBaseRand ) < y_padd ) {
	    stack[sp++] = nn;
	    lattice[ nn ] = newspin;
	    *mag_ptr += (static_cast<double>(newspin*2));
	  }
	}
	nn = current - ipow( 2, hier );
	if ( DEBUG ) cout << "case (SIZE-1)/2: " << current << "\tneighbor: " << ((SIZE-1)/2) - ipow( 2, hier ) << endl;
	if ( lattice[ nn ] == oldspin ) {
	  if ( gsl_rng_uniform_pos( gBaseRand ) < y_padd ) {
	    stack[sp++] = nn;
	    lattice[ nn ] = newspin;
	    *mag_ptr += (static_cast<double>(newspin*2));
	  }
	}
	hier--;
      }
      break;
    default: 
      while ( hier > 0 ) { 
	nn = (temp + ipow( 2, hier )) % (SIZE-1);
	if ( DEBUG ) cout << "current: " << current << "\tneighbor: " << temp + ipow( 2, hier) << endl;
	if ( lattice[ nn ] == oldspin ) {
	  if ( gsl_rng_uniform_pos( gBaseRand ) < y_padd ) {
	    stack[sp++] = nn;
	    lattice[ nn ] = newspin;
	    *mag_ptr += (static_cast<double>(newspin*2));
	  }
	}
	nn = (temp - ipow( 2, hier )) % (SIZE-1);
	if ( DEBUG ) cout << "current: " << current << "\tneighbor: " << temp - ipow( 2, hier) << endl;
	if ( lattice[ nn ] == oldspin ) {
	  if ( gsl_rng_uniform_pos( gBaseRand ) < y_padd ) {
	    stack[sp++] = nn;
	    lattice[ nn ] = newspin;
	    *mag_ptr += (static_cast<double>(newspin*2));
	  }
	}
	hier--;
      }
      count = 0;
      hier = 1;
      break;
    }
  }
}

double HK( int lattice[ SIZE ], int largestCluster[ DATA ][ SIZE ], int clusterCount, double T )
{
  int label = 0;
  int nn, divisor, index, count, big;
  int tag[ SIZE ] = { 0 };
  int tag2[ SIZE ] = { 0 };
  int labels[ SIZE ] = { 0 };
  int *lptr;
  
  // connects the backbone
  for( int i = 0; i < SIZE-1; i++ )
    {
      if ( lattice[ i ] == lattice[ i+1 ] ) tag[ i+1 ] = tag[ i ];     // keep the same label
      else { label++; tag[ i+1 ] = label; }                            // increase label by 1 and mark the new cluster in tag[ ]
    }

  for( int i = 0; i < SIZE-1; i++ )
    { 
      divisor = LIMIT;    // used to determine connectivity of spins     
      while( divisor >= 1 ) {
	if ( i % divisor == 0 ) {
	  if ( i + divisor < SIZE ) { 
	    nn = (i + divisor); 
	    if ( lattice[ i ] == lattice[ nn ] && tag[ i ] != tag[ nn ] ) {
	      tag[ i ] = min( tag[ i ], tag[ nn ] );
	      tag[ nn ]   = min( tag[ i ], tag[ nn ] );
	    }
	  }
	  if ( i - divisor >= 0 ) { 
	    nn = (i - divisor); 
	    if ( lattice[ i ] == lattice[ nn ] && tag[ i ] != tag[ nn ] ) {
	      tag[ i ] = min( tag[ i ], tag[ nn ] );
	      tag[ nn ]   = min( tag[ i ], tag[ nn ] );
	    }
	  }
	}
	else if ( (2*i) % divisor == 0 && i % divisor != 0 && (i / divisor) % 2 == 0 ) {
	  if ( i + divisor < SIZE ) { 
	    nn = (i + divisor);  
	    if ( lattice[ i ] == lattice[ nn ] && tag[ i ] != tag[ nn ] ) {
	      tag[ i ] = min( tag[ i ], tag[ nn ] );
	      tag[ nn ]   = min( tag[ i ], tag[ nn ] );
	    }
	  }
	}
	else if ( (2*i) % divisor == 0 && i % divisor != 0 && (i / divisor) % 2 != 0 ) {
	  if ( i - divisor >= 0 ) { nn = (i - divisor);  
	    if ( lattice[ i ] == lattice[ nn ] && tag[ i ] != tag[ nn ] ) {
	      tag[ i ] = min( tag[ i ], tag[ nn ] );
	      tag[ nn ]   = min( tag[ i ], tag[ nn ] );
	    }
	  }
	}	
	divisor /= 2;
      }
      if ( i == 0 ) { 
	nn = SIZE-1; 
	if ( lattice[ i ] == lattice[ nn ] && tag[ i ] != tag[ nn ] ) {
	  tag[ i ] = min( tag[ i ], tag[ nn ] );
	  tag[ nn ]   = min( tag[ i ], tag[ nn ] );
	}	
      }
      if ( i == SIZE-1 ) { 
	nn = 0; 
	if ( lattice[ i ] == lattice[ nn ] && tag[ i ] != tag[ nn ] ) {
	  tag[ i ] = min( tag[ i ], tag[ nn ] );
	  tag[ nn ]   = min( tag[ i ], tag[ nn ] );
	  cout << tag[ i ] << '\t' << tag[ nn ] << '\t' << min( tag[ i ], tag[ nn ] ) << "*" << endl;
	}
      }
    }
  for ( int i = 0; i < SIZE; i++ ) tag[ i ] = tag[ i ] + 1;
  
  //  cout << "****************************************" << '\n';
  //  print( lattice );
  //  print( tag );

  for ( int i = 0; i < SIZE; i++ ) tag2[ i ] = tag[ i ];       // copies tag to tag2

  // finds the unique elements in the tag[ ] array and 
  // stores them in the labels[ ] array
  sort( tag, tag+SIZE );
  labels[ 0 ] = tag[ 0 ];
  lptr = &labels[ 1 ];
  for ( int i = 1; i < SIZE; i++ ) {
    if ( tag[ i ] == *(lptr-1) ) ;
    else { *lptr = tag[ i ]; lptr++; }
  }
  
  //  print( tag );
  //  print( labels );
  
  big = 0;
  lptr = &labels[ 0 ];
  while ( *lptr != 0 ) {
    count = 0;
    for ( int i = 0; i < SIZE; i++ )  
      if ( tag[ i ] == *lptr ) count++;
    if ( count > big ) { big = count; index = *lptr; }
    lptr++;
  }
  //cout  << "biggest cluster size: " << big << "\nindex:" << index << "\n";
  
  for ( int i = 0; i < SIZE; i++ )
    if ( tag2[ i ] == index ) largestCluster[ clusterCount ][ i ]++;

  return big;
}

int min( int x, int y )
{
  if ( x < y ) return x;
  else if ( x > y ) return y;
  else return x;
}

int bondCountHN5( int size ) 
{
  int k = 0;
  int tempSize = size-1;
  while( tempSize % 2 == 0 ) {
    tempSize /= 2; 
    k++;
  }
  
  int bondSum = ((( 5 * ipow( 2, k )) - 8) / 2) + 2;

  return bondSum;
}

int spinEnergy( int *lptr, int index )
{
  int eSum = 0;        // running tally of energy
  int *nnptr;          // points to nearest neighbors
  int *origin = lptr;  // stationary pointer at lattice[0]
  lptr += index;       // pointer starts at beginning of array and moves to index
  nnptr = lptr;        // pointer starts at lptr
  int temp = 0;
  int count = 0;
  int hier = 1;

  
  // check the backbone neighbor forward
  if ( index < (SIZE-1) ) {
    nnptr++;
    if ( DEBUG ) cout << "index: " << index << "\tneighbor: " << index+1 << endl;
    eSum += (*lptr) * (*nnptr);
  }

  // reset nnptr
  nnptr = lptr;

  // check the backbone neighbor backward
  if ( index != 0 ) {
    nnptr--;
    if ( DEBUG ) cout << "index: " << index << "\tneighbor: " << index-1 << endl;
    eSum += (*lptr) * (*nnptr);
  }

  // reset nnptr
  nnptr = lptr;

  //  counts the number of hierarchical levels; assigns to hier variable
  temp = index;
  if( index == 0 ) 
    temp = SIZE-1;
  while( temp % 2 == 0 ) { temp /= 2; hier++; }
  
  temp = index;
  if( index == 0 ) 
    temp = SIZE-1;
  while( temp % 2 == 0 ) 
    temp /= 2;
  temp++;
  while( temp % 2 == 0 ) { 
    temp /= 2; 
    count++; 
  }
  
  // reset nnptr to origin
  nnptr = origin;

  // highest hierarchical neighbor
  switch ( index ) {
  case 0:
    nnptr += SIZE-1;
    if ( DEBUG ) cout << "index: " << index << "\tneighbor: " << SIZE-1 << endl;
    eSum += (*lptr) * (*nnptr);
    break;
  case SIZE-1:
    eSum += (*lptr) * (*nnptr);  // do nothing with nnptr
    if ( DEBUG ) cout << "index: " << index << "\tneighbor: " << 0 << endl;
    break;
  case (SIZE-1)/2:
    // do nothing
    break;
  default:
    if (count==1) {
      nnptr += (index + ipow( 2, hier )) % (SIZE-1);
      if ( DEBUG ) cout << "index: " << index << "\tneighbor: " << (index + ipow( 2, hier )) % (SIZE-1) << endl;
    }
    else {
      nnptr += (index - ipow( 2, hier )) % (SIZE-1) ;
      if ( DEBUG ) cout << "index: " << index << "\tneighbor: " << (index - ipow( 2, hier )) % (SIZE-1) << endl;
    }
    eSum += (*lptr) * (*nnptr);
    break;
  }
    
  // reset nnptr to origin
  nnptr = origin;

  ( index == 0 ) ? temp = SIZE-1 : temp = index;
  hier--; // decrement hier

  switch ( index ) { 
  case 0:
    while ( hier > 1 ) {
      nnptr = origin; // reset nnptr to origin
      nnptr += ipow(2, hier-1);
      if ( DEBUG ) cout << "index: " << index << "\tneighbor: " << ipow( 2, hier-1 ) << endl;      
      eSum += (*lptr) * (*nnptr);
      hier--;
    }
    break;
  case (SIZE-1):
    while ( hier > 1 ) {
      nnptr = lptr; // set nnptr to end of lattice
      nnptr -= ipow( 2, hier-1);
      if ( DEBUG ) cout << "index: " << index << "\tneighbor: " << (SIZE-1) - ipow( 2, hier-1 ) << endl;      
      eSum += (*lptr) * (*nnptr);
      hier--;
    }
    break;
  case (SIZE-1)/2:
    while ( hier > 0 ) {
      nnptr = lptr;
      nnptr += ipow( 2, hier );
      if ( DEBUG ) cout << "index: " << index << "\tneighbor: " << ((SIZE-1)/2) + ipow( 2, hier ) << endl;
      eSum += (*lptr) * (*nnptr);
      nnptr = lptr;
      nnptr -= ipow( 2, hier );
      if ( DEBUG ) cout << "index: " << index << "\tneighbor: " << ((SIZE-1)/2) - ipow( 2, hier ) << endl;
      eSum += (*lptr) * (*nnptr);
      hier--;
    }
    break;
  default:
    while( hier > 0 ) {
      nnptr += (temp + ipow( 2, hier )) % (SIZE-1);
      if ( DEBUG ) cout << "index: " << index << "\tneighbor: " << temp + ipow( 2, hier) << endl;
      eSum += (*lptr) * (*nnptr);
      
      //reset nnptr to origin
      nnptr = origin;
      
      nnptr += (temp - ipow( 2, hier)) % (SIZE-1);
      if ( DEBUG ) cout << "index: " << index << "\tneighbor: " << temp - ipow( 2, hier) << endl;
      eSum += (*lptr) * (*nnptr);
      
      // reset nnptr to origin
      nnptr = origin;
      
      hier--;
    }
    count = 0;
    hier = 1;
    break;
  }
  return eSum;
}
