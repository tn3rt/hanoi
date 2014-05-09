#define SIZE 64
#define LIMIT SIZE/4
#define LATTICE_SWEEPS 10000
#define EQUILIBRIUM 10000

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
gsl_rng *gBaseRand;       /* global rand number generator */

int main() 
{
  
  //---------------FUNCTION PROTOTYPES---------------//
  void randomize( int[ SIZE ] );
  void unionize( int[ SIZE ] );
  void printLattice( int[ SIZE ] );
  void sweep( int[ SIZE ], double * );
  double calc_magnetization( int [ SIZE ] );
  double power( double, int );

  //---------------VARIABLES-------------------------//
  double magnet[ LATTICE_SWEEPS ] = { 0 };
  double *magnet_ptr;
  double *magnetization_ptr;
  double magnetization;
  int pause;

  unsigned long randSeed;  
  gBaseRand = gsl_rng_alloc(gsl_rng_taus);
  srand(time(NULL));                    /* initialization for rand() */
  randSeed = rand();                    /* returns a non-negative integer */
  gsl_rng_set (gBaseRand, randSeed);    /* seed the PRNG */



  //  init_gen_rand( time( 0 ) );  
  int lattice[ SIZE ] = { 0 };  

  // start at high temperature state
  //  randomize( lattice );
  unionize( lattice );
  magnetization = calc_magnetization( lattice );
  magnetization_ptr = &magnetization;
  magnet_ptr = &magnet[ 0 ];

  for ( int i = 0; i < EQUILIBRIUM; i++ ) {
    sweep( lattice, magnetization_ptr );
    if ( i % 1000 == 0 )
      cout << "Equilibrating: " << i << " of " << EQUILIBRIUM << endl;
  }

  for ( int i = 0; i < LATTICE_SWEEPS; i++ ) {
    for ( int j = 0; j < SIZE; j++ )
      sweep( lattice, magnetization_ptr );

    // records the 1st, 2nd, and 4th moments
    *magnet_ptr = magnetization;
    magnet_ptr++;

    if ( i % 1000 == 0 )
      cout << i << endl;
  }

  ofstream outFile( "Equilibrium-HN5-4.35.txt" );
  if (!outFile) {
    cout << "Cannot open file.\n";
    return 1;
  }
  for ( int i = 0; i < LATTICE_SWEEPS; i++ )
    outFile << i << "\t" << abs( magnet[ i ] ) << endl;

  return 0;
}


//-----------------FUNCTIONS------------------------//
void randomize( int lattice[ SIZE ] )
{
  double random_num;
  int spin;
  
  for ( int i = 0; i < SIZE; i++ ) {
    //    random_num = genrand_res53();
    if ( ( random_num * 2 - 1 ) > 0 )
      spin = 1;
    else if ( ( random_num * 2 - 1 ) < 0 )
      spin = -1;
    lattice[ i ] = spin ;
  }
}
  
void unionize( int lattice[ SIZE ] )
{
  for ( int i = 0; i < SIZE; i++ )
    lattice[ i ] = 1;
}

double power( double x, int y )    // same as pow( x, y ) but w/ integer only exponents
{                                  // problem with condor_compile (cluster software) accepting
  double c;                        // the pow( x, y ) function in <cmath> library
  c = 1; 
  
  for ( int i = 0; i < y; i++ )
    c *= x;
  
  return c;
}
 

void printLattice( int lattice[ SIZE ] )
{
  for ( int i = 0; i < SIZE; i++ ) {
    if ( lattice[ i ] > 0 )
      cout << " " << lattice[ i ] << " ";
    else 
      cout << lattice[ i ] << " ";
  }
  cout << endl;
}

double calc_magnetization( int lattice[ SIZE ] ) 
{ 
  int mag = 0;
  for ( int i = 0; i < SIZE; i++ ) {
      if ( lattice[ i ] == 1 )
	mag++;
      else if ( lattice[ i ] == -1 )
	mag--;
  }
  return static_cast< double >( mag ) / SIZE;
}

 
void sweep( int lattice[ SIZE ], double *magnetization_ptr )
{
  int *lattice_ptr, *bond_ptr;
  int random_number;
  int divisor = LIMIT;
  int count = LIMIT / 2;

  double temperature = 4.35;
  double energy = 0;
  double flipped_energy;
  double random;
  double probability;

  for ( int i = 0; i < SIZE; i++ ) {

    random_number = gsl_rng_uniform_int ( gBaseRand, SIZE-1 );
    lattice_ptr = &lattice[ random_number ];

    if ( random_number == 0 )
      random_number = SIZE;

    //    cout << "random: " << random_number << "\tvalue: " << *lattice_ptr << endl;

    while( divisor >= 1 ) {
      
      if ( random_number % divisor == 0 ) {
	bond_ptr = &lattice[ (random_number + divisor) % SIZE ];
	energy += (-1) * ( *lattice_ptr * *bond_ptr );
	//	cout << (random_number + divisor) % SIZE << "-->" << *bond_ptr << "\tenergy: " << energy << endl;
	bond_ptr = &lattice[ (random_number - divisor) % SIZE ];
	energy += (-1) * ( *lattice_ptr * *bond_ptr );
	//	cout << (random_number - divisor) % SIZE << "-->" << *bond_ptr << "\tenergy: " << energy << endl;
      }

      else if ( (2*random_number) % divisor == 0 && random_number % divisor != 0 && (random_number / divisor) % 2 == 0 ) {
	bond_ptr = &lattice[ (random_number + divisor) % SIZE ];
	energy += (-1) * ( *lattice_ptr * *bond_ptr );
	//	cout << (random_number + divisor) % SIZE << "-->" << *bond_ptr << "\tenergy: " << energy << endl;
      }

      else if ( (2*random_number) % divisor == 0 && random_number % divisor != 0 && (random_number / divisor) % 2 != 0 ) {
	bond_ptr = &lattice[ (random_number - divisor) % SIZE ];
	energy += (-1) * ( *lattice_ptr * *bond_ptr );
	//	cout << (random_number - divisor) % SIZE << "-->" << *bond_ptr << "\tenergy: " << energy << endl;
      }
      divisor /= 2;
    } 
    // end of while loop
    flipped_energy = (-1) * energy;

    if ( flipped_energy - energy < 0 ) {
      *lattice_ptr = (-1) * (*lattice_ptr);
      *magnetization_ptr += 2 * ( static_cast< double >( *lattice_ptr ) / SIZE );
    }
    else if ( flipped_energy - energy > 0 ) {
      random = gsl_rng_uniform_pos ( gBaseRand );
      probability = exp( (-1) * ( flipped_energy - energy ) / temperature );
      if ( random < probability ) {
	*lattice_ptr = (-1) * (*lattice_ptr);
	*magnetization_ptr += 2 * ( static_cast< double >( *lattice_ptr ) / SIZE );
      }
      //      cout << "random: " << random << "\tprobability: " << probability << endl;
    }

    divisor = LIMIT;
    energy = 0;
  }
  //  cout << endl;
}
