// A header file with some tools for solving Diophantine equations
// Often appear when asking about rational solutions to polynomial equations
//
// Utils for recreational computuing and math
// Created by Ryan Hamilton on Aug. 26th, 2025

#ifndef diophantine_tools_h
#define diophantine_tools_h

#include <ostream>
#include "../cosmetic/typestyle.h"

typedef std::pair<uLong_64t, uLong_64t> pell_pair_t;

//-----------------------------------------------Forward declare (all methods in file):

//---------------------Tools for Lucas sequences

// Iterate a lucas sequence mod n
uLong_64t lllucasMod(Long_64t n,
                     Long_64t lucas_P,
                     Long_64t lucas_Q,
                     bool for_prime_test = false);

//---------------------Tools for solving Pell's/Bhaskara's Equation

// Compose data to the pell_pair_t datatype
pell_pair_t* makePellPair(uLong_64t noncoef, uLong_64t coef);

// Iterate via the method proposed by Bhaskara in ancient India
pell_pair_t* bhaskaraIterate(pell_pair_t* input,
                             pell_pair_t* base_sol,
                             int pell_root_D);


//-----------------------------------------------Tools for Lucas sequences




// Compute the n'th member of the Lucas sequence
// with input coefficients {P,Q}, defined by recurrence
// x_n = P*x_{n-1} - Q*x_{n-2} (mod n)
//
// if requested (useful for primality tests) returns L_{n+1} % n
uLong_64t lllucasMod(Long_64t n,
                     Long_64t lucas_P,
                     Long_64t lucas_Q,
                     bool for_prime_test) {
  uLong_64t seq_stop = n + for_prime_test;
  int get_log2 = 0; // get index of first nonzero bit
  while (seq_stop >> get_log2 != 1) ++get_log2;
  
  uLong_64t sqrt_benchmark = n / 2;
  
  // initial state: k=1
  Long_64t lucas_D = (lucas_P*lucas_P - 4*lucas_Q); // descriminant
  Long_64t lucas_U = 1;
  Long_64t lucas_V = lucas_P;
  Long_64t cpow_Q = lucas_Q;
  while (--get_log2 >= 0) {
    uLong_64t cbits = seq_stop >> get_log2; // truncate to one additional bit
    
    // make sure U,V are as small as possible before squaring them
//    if (lucas_U > INT_MAX)
    // index doubling (k->2k), using binary strategy
    lucas_U = (lucas_U*lucas_V) % n;
    lucas_V = (lucas_V*lucas_V - 2*cpow_Q) % n;
    cpow_Q  = (cpow_Q * cpow_Q) % n;
//    std::cout << "L_{x" << 2*(cbits >> 1) << "} : U = " << lucas_U << ", V = " << lucas_V << std::endl;
    
    // iterate once if next bit is 1 rather than 0. (2k->2k+1)
    if (cbits % 2 == 1) {
      Long_64t lucas_newU = lucas_P * lucas_U + lucas_V;
      Long_64t lucas_newV = lucas_D * lucas_U + lucas_P * lucas_V;
      
      // divison by 2 mod n: if result is odd, add n
      // This ensures that the correct mod is obtained.
      if (lucas_newU % 2 == 0) lucas_U = (lucas_newU / 2) % n;
      else                     lucas_U = ((lucas_newU + n) / 2) % n;
      if (lucas_newV % 2 == 0) lucas_V = (lucas_newV / 2) % n;
      else                     lucas_V = ((lucas_newV + n) / 2) % n;
      cpow_Q = (cpow_Q * lucas_Q) % n;
//      std::cout << "L_{" << cbits << "} : U = " << lucas_U << ", V = " << lucas_V << std::endl;
    }

  }return lucas_U;
}


//-----------------------------------------------Tools for solving Pell's/Bhaskara's Equation

// Form a pell pair (x,y) which solves some Pell-type equation, x^2 - D y^2 = N
pell_pair_t* makePellPair(uLong_64t noncoef, uLong_64t coef) {
  pell_pair_t* pair = new pell_pair_t;
  pair->first = noncoef;
  pair->second = coef;
  return pair;
}


// Perform two iterations of the recursive sequence generation
// for solutions to the generalized Pell/Bhaskara equation, from a known solution.
// The method also requires as an input the coefficient D on one of the squares.
// The iteration algorithm is due to Bhaskara in ancient India.
pell_pair_t* bhaskaraIterate(pell_pair_t* input,
                             pell_pair_t* base_sol,
                             int pell_root_D) {
  pell_pair_t* sol = new pell_pair_t;
  
  // The next iteration is the coefficients
  sol->first = (input->first * base_sol->first
                + pell_root_D * input->second * base_sol->second);
  sol->second = (input->first * base_sol->second
                 + input->second * base_sol->first);
  
  return sol;
}



#endif /* diophantine_tools_h */
