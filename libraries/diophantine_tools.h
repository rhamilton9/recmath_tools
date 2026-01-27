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

// Generate the periodic continued fraction representation of sqrt{N}
// Useful for finding initial solutions to Pell's Equation
std::vector<int> getRepFracForSqrt(int N,
                                   bool print = false);

// Use the continued fraction to find the convergent index that yields
// the smallest solution to x^2 - D y^2 = 1, found from the repetition period.

// TODO implement pell_pair for MPInt type (needs multiplication).


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
}// End of diophantine_tools::lllucasMod


//-----------------------------------------------Tools for solving Pell's/Bhaskara's Equation

// Form a pell pair (x,y) which solves some Pell-type equation, x^2 - D y^2 = N
pell_pair_t* makePellPair(uLong_64t noncoef, uLong_64t coef) {
  pell_pair_t* pair = new pell_pair_t;
  pair->first = noncoef;
  pair->second = coef;
  return pair;
}// End of diophantine_tools::makePellPair


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
}// End of diophantine_tools::bhaskaraIterate


// Compute the periodic continued fraction representation for the square root of N
// This output will be a vector of the form {A_0, A_1, A_2, ...} which represents the
// periodic continued fraction shown graphically as
//                                    1
//    sqrt{N} = A_0 + --------------------------------------
//                                           1
//                        A_1 + ----------------------------
//                                                  1
//                                    A_2 + ----------------
//                                                  ...
//
//            =(DEF) [A_0, A_1, A_2, ...] (notation)
//
// The algorithm finds the representation which is periodic, meaning it can be expressed as
//
//    sqrt{N} = [A_0; (A_1, A_2, ... )],
//
// implying the part in parentheses is repeated periodically.
// This representation was proven to always exist by Lagrange.
//
// This form is especially useful for finding initial solutions to Pell's Equation.
std::vector<int> getRepFracForSqrt(int N, bool print) {
  bool debug = true;
  std::vector<int> repstring;
  Long_64t root_floor = std::floor(std::sqrt(N));
  if (root_floor * root_floor == N) { // No repetition for perfect squares--terminate
    if (print) {
      std::cout << "sqrt{" << "\033[31m" << N << "\033[39m" << "} = [";
      std::cout << "\033[31m" << root_floor << "\033[39m" << "]." << std::endl;
    }return repstring;
  }// End of perfect square check
  
  // Set up for iterative loop
  Long_64t A = std::floor(std::sqrt(N));
  Long_64t D = 1;
  if (print) {
    std::cout << "sqrt{" << "\033[31m" << N << "\033[39m" << "} = [";
    std::cout << "\033[31m" << A << "\033[39m" << ";(";
  }
  do { // while (D != 1); check last since we begin with D = 1.
    
    // Debug
    if (print && D != 1) std::cout << ',';
    if (debug && (N - A*A) % D != 0 || D == 0) {
      std::cout << "ERROR :: nonwhole divisor!" << std::endl;
      break;
    }
    
    // Perform new D computation first for efficiency
    D = (N - A * A) / D;
    
    // Take the floor and assign it as the repeated index
    repstring.push_back(A / D); // Int type does floor on its own
    if (A != D || A == 1) ++repstring.back(); // Floor -> Ceiling
    
    
    // Next A from current values
    A = repstring.back()*D - A;
    if (D <= 0) {std::cout << "Failure in D." << std::endl; break;}
    while (root_floor - A >= D) { // Assure the remaining fraction is less than 1
      A += D;
      ++repstring.back();
    }// End of fraction reduction
    
    
    // Print about current iteration
    if (print) std::cout << "\033[34m" << repstring.back() << "\033[39m";
  } while (D != 1); // Iteration denominator is 1, return to sqrt(n) and repeat again.
  
  
  if (print)std::cout << ")]" << std::endl;
  return repstring;
}// End of diophantine_tools::getRepFracForSqrt


#endif /* diophantine_tools_h */
