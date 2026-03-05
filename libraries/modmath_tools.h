// A header file with some tools for modular arithmetic.
// These are especially often useful for primality tests,
// divisibility and similar relations between numbers.
//
// Utils for recreational computuing and math
// Created by Ryan Hamilton on Jun. 6th, 2025.

#ifndef modmath_tools_h
#define modmath_tools_h

#include <ostream>
#include "../cosmetic/typestyle.h"

typedef std::pair<Int_t,    std::pair<Int_t, Int_t> >       GCD_Int_t;
typedef std::pair<Long32_t, std::pair<Long32_t, Long32_t> > GCD_Long32_t;
typedef std::pair<Long64_t, std::pair<Long64_t, Long64_t> > GCD_Long64_t;

//-----------------------------------------------Forward declare (all methods in file):

//------------------------Tools for performing modular arithmetic

// An attempt to make a safe multiplication algorithm for large modulus
// It doesn't quite work though, the karatsuba algorithm
//    https://en.wikipedia.org/wiki/Karatsuba_algorithm
// just isn't the right algorithm for it...
uLong64_t karatsubaProdMod(uLong64_t a, uLong64_t b, uLong64_t modbase);

// An attempt using a more advanced algorithm, using interesting algebra properties
//    https://en.wikipedia.org/wiki/Montgomery_modular_multiplication
// 
uLong64_t montgomeryProdMod(uLong64_t a, uLong64_t b, uLong64_t modbase);

// An efficient modular exponentiation algorithm
uLong32_t fastPowMod(uLong64_t expbase, uLong32_t power, uLong32_t modbase);

uLong64_t lfastPowMod(uLong64_t expbase, uLong64_t power, uLong64_t modbase);

//------------------------Number Theoretic Modular Tools

// Legendre symbol, assumes denominator is odd prime
Int_t legendreSymbol(Int_t top, Int_t bot_oddprime);

Int_t llegendreSymbol(uLong32_t top, uLong32_t bot_oddprime);

Int_t lllegendreSymbol(uLong64_t top, uLong64_t bot_oddprime);

// Jacobi symbol, requires denominator to be odd positive integer
Int_t jacobiSymbol(Int_t top, Int_t bot);

Int_t ljacobiSymbol(Long32_t top, Long32_t bot);

Int_t lljacobiSymbol(Long64_t top, Long64_t bot);

//------------------------GCD/LCM

// Algorithms to compute the GCD and Bezout coefficients for various types
GCD_Int_t* extendedEuclidean(Int_t a, Int_t b);

GCD_Long32_t* lextendedEuclidean(Long32_t a, Long32_t b);

GCD_Long64_t* llextendedEuclidean(Long64_t a, Long64_t b);

// Compute LCM using GCD
Int_t LCM(Int_t a, Int_t b);

uLong32_t lLCM(uLong32_t a, uLong32_t b);

uLong64_t llLCM(uLong64_t a, uLong64_t b);

// Compute LCM for a vector of integers
Int_t LCM(std::vector<Int_t>& vec);

uLong32_t LCM(std::vector<uLong32_t>& vec);

uLong64_t LCM(std::vector<uLong64_t>& vec);

//------------------------Chinese Remainder Theorem

// Algorithms to compute the minimal solution guaranteed by the chinese remainder theorem
Int_t chineseRemainder(Int_t rem1, Int_t rem2,
                         Int_t mod1, Int_t mod2);

Long32_t lchineseRemainder(Long32_t rem1, Long32_t rem2,
                           Long32_t mod1, Long32_t mod2);

Long64_t llchineseRemainder(Long64_t rem1, Long64_t rem2,
                            Long64_t mod1, Long64_t mod2);



//-----------------------------------------------Tools for performing modular arithmetic

// Implements the Karatsuba algorithm for multiplication
// This allows for a safe multiplication of numbers larger than 2^32
// which would otherwise overflow and lead to errors.
uLong64_t karatsubaProdMod(uLong64_t a, uLong64_t b, uLong64_t modbase) {
  // mod x^32
  uLong64_t a_small = a % 4294967296L;
  uLong64_t b_small = b % 4294967296L;
  // shift by 32 bits: equivalent to division by x^32
  a = a >> 32; b = b >> 32;
  
  uLong64_t prod_large = a * b;
  uLong64_t prod_small = a_small * b_small;
  uLong64_t prod_mixed = (a + a_small)*(b + b_small) - prod_large - prod_small;
  
  /* debug testing:
  std::cout << "bits of a: (" <<  a << "*2^32 + " << a_small << ")." << std::endl;
  std::cout << "bits of b: (" <<  b << "*2^32 + " << b_small << ")." << std::endl;
  
  std::cout << "prod large: " << prod_large << std::endl;
  std::cout << "prod small: " << prod_small << std::endl;
  std::cout << "prod mixed: " << prod_mixed << std::endl;
  
  std::cout << "result: " << (((prod_large * (ULONG_MAX % modbase)) % modbase + prod_large +
                               (prod_mixed * (4294967296L % modbase)) % modbase ) % modbase +
                              prod_small) % modbase << std::endl;
  */
  
  
  prod_large = prod_large % modbase;
  prod_small = prod_small % modbase;
  prod_mixed = prod_mixed % modbase;
  
  
  // same problem at a deeper level: modulus is too large.
  return (((prod_large * (ULONG_MAX % modbase)) % modbase + prod_large +
           (prod_mixed * (4294967296L % modbase)) % modbase ) % modbase +
          prod_small) % modbase;
}

// Compute modular powers by a "binary" computation strategy
// There may be a faster way in general (See Problem 122)
// but for now this is the fastest method I'm aware of.
uLong32_t fastPowMod(uLong64_t expbase, uLong32_t power, uLong32_t modbase) {
  std::vector<std::pair<int, uLong32_t> >* powmod = new std::vector<std::pair<int, uLong32_t> >();
  powmod->push_back(std::make_pair(1, (expbase % modbase)));
  
  // repeatedly square until that exceeds the power
  while (2*powmod->back().first <= power)
    powmod->push_back(std::make_pair(2*powmod->back().first,
                                     (powmod->back().second*powmod->back().second) % modbase));
  
  // Fill the remaining exponential with powers of 2 (last binary digits)
  int binary_powertracker = powmod->size() - 1;
  while (powmod->back().first < power) {// search for largest power that fits
    while (powmod->back().first + powmod->at(binary_powertracker).first > power) --binary_powertracker;

    // append to list
    powmod->push_back(std::make_pair(powmod->back().first + powmod->at(binary_powertracker).first,
                                     (powmod->back().second*powmod->at(binary_powertracker).second) % modbase));
  }return powmod->back().second;
}

// For all 64-bit integers
// Take extra care with multiplication to assure no overflow
uLong64_t lfastPowMod(uLong64_t expbase, uLong64_t power, uLong64_t modbase) {
  std::vector<std::pair<uLong64_t, uLong64_t> >* powmod = new std::vector<std::pair<uLong64_t, uLong64_t> >();
  powmod->push_back(std::make_pair(1, (expbase % modbase)));
  
  // repeatedly square until that exceeds the power
  // use safe multiplication
  while (2*powmod->back().first <= power)
    powmod->push_back(std::make_pair(2*powmod->back().first,
                                     (powmod->back().second*powmod->back().second) % modbase) );
//                                     karatsubaProdMod(powmod->back().second,
//                                                      powmod->back().second,
//                                                      modbase) ) ); // karatsuba ver.
  
  //debug
//  std::cout << expbase << "^" << powmod->back().first << " % " << modbase << " = " << powmod->back().second << std::endl;
  
  // Fill the remaining exponential with powers of 2 (last binary digits)
  int binary_powertracker = powmod->size() - 1;
  while (powmod->back().first < power) {// search for largest power that fits
    while (powmod->back().first + powmod->at(binary_powertracker).first > power) --binary_powertracker;
    
    // append to list with safe multiplication
    powmod->push_back(std::make_pair(powmod->back().first + powmod->at(binary_powertracker).first,
                                     (powmod->back().second*powmod->at(binary_powertracker).second) % modbase) );
//                                     karatsubaProdMod(powmod->back().second,
//                                                      powmod->at(binary_powertracker).second,
//                                                      modbase) ) ); // karatsuba ver.
    
    // debug
//    std::cout << expbase << "^" << powmod->back().first << " % " << modbase << " = " << powmod->back().second << std::endl;
    
  }return powmod->back().second;
}


//-----------------------------------------------Number Theoretic Modular Tools

// Find the Legendre symbol (top/bot) = {-1, 0, 1}
// Assumes that the lower number is an odd prime,
// which is required for the Legendre symbol.
Int_t legendreSymbol(Int_t top, Int_t bot_oddprime) {
  uLong32_t cyclemod = fastPowMod(top % bot_oddprime, (bot_oddprime - 1)/2, bot_oddprime);
  if (cyclemod == bot_oddprime - 1) return -1;
  return static_cast<Int_t>(cyclemod);
}

Int_t llegendreSymbol(uLong32_t top, uLong32_t bot_oddprime) {
  uLong32_t cyclemod = fastPowMod(top % bot_oddprime, (bot_oddprime - 1)/2, bot_oddprime);
  if (cyclemod == bot_oddprime - 1) return -1;
  return static_cast<Int_t>(cyclemod);
}

Int_t lllegendreSymbol(uLong64_t top, uLong64_t bot_oddprime) {
  uLong64_t cyclemod = fastPowMod(top % bot_oddprime, (bot_oddprime - 1)/2, bot_oddprime);
  if (cyclemod == bot_oddprime - 1) return -1;
  return static_cast<Int_t>(cyclemod);
}

// Compute the Jacobi symbol (top/bot) = {-1, 0, 1}
// Unlike the Legendre symbol, the bottom integer need
// not be an odd prime (but must still be an odd integer.
// The two do have some similar rules, allowing us
// to compute the Jacobi symbol using a GCD-like
// algorithm approach (i.e. like Extended-Euclidean)
//
// Many ideas (bitwise operations especially) are taken from:
// https://en.wikipedia.org/wiki/Jacobi_symbol#Implementation_in_C++
Int_t jacobiSymbol(Int_t top, Int_t bot) {
  if (bot < 0 || bot % 2 == 0) // Throw error on bad input
    std::cerr << "Error in <modmath_tools>::jacobiSymbol: Bad input." << std::endl;
  
  // Since jacobi symbol is equivalent mod the lower index,
  // we can reduce top to the smallest positive modulus
  top = (top % bot + bot) % bot;
  
  // For more efficiency, we can use bitwise operations on jacobi.
  // This requires starting with zero (as a logical bit)
  Int_t jacobi = 0;
  while (top != 0) {
    while (top % 4 == 0) top /= 4; // squares never contribute
    
    // If there is a remaining factor of 2, see if it adds a minus sign
    if (top % 2 == 0) {
      // efficiently compute mod8 condition
      // (2,n) = -1 only if n % 8 = {3,5}
      jacobi ^= bot; //bitwise XOR
      top /= 2;
    }
    
    // Flip sign if top, bot are 3 mod 4.
    // As a logic operation, we require top and bot have the same
    // final 2 bits and flips only so long as those bits are ...11 (...3)
    jacobi ^= top & bot & 2;
    
    // perform Euclidean algorithm (swap top/bot, mod and iterate)
    Int_t rem = bot % top;
    bot = top;
    top = rem;
  }
  
  // loop exit test cases
  if (bot != 1) return 0; // not coprime (GCD != 1)
  if ((jacobi ^ (jacobi >> 1)) & 2) return -1; // quadratic resiprocity
  return 1; // square or otherwise fails resiprocity
}

// Jacobi symbol for 32-bit long input
Int_t ljacobiSymbol(Long32_t top, Long32_t bot) {
  if (bot < 0 || bot % 2 == 0) // Throw error on bad input
    std::cerr << "Error in <modmath_tools>::jacobiSymbol: Bad input." << std::endl;
  
  // Since jacobi symbol is equivalent mod the lower index,
  // we can reduce top to the smallest positive modulus
  top = (top % bot + bot) % bot;
  
  // For more efficiency, we can use bitwise operations on jacobi.
  // This requires starting with zero (as a logical bit)
  Long32_t jacobi = 0;
  while (top != 0) {
    while (top % 4 == 0) top /= 4; // squares never contribute
    
    // If there is a remaining factor of 2, see if it adds a minus sign
    if (top % 2 == 0) {
      // efficiently compute mod8 condition
      // (2,n) = -1 only if n % 8 = {3,5}
      jacobi ^= bot; //bitwise XOR
      top /= 2;
    }
    
    // Flip sign if top, bot are 3 mod 4.
    // As a logic operation, we require top and bot have the same
    // final 2 bits and flips only so long as those bits are ...11 (...3)
    jacobi ^= top & bot & 2;
    
    // perform Euclidean algorithm (swap top/bot, mod and iterate)
    Long32_t rem = bot % top;
    bot = top;
    top = rem;
  }
  
  // loop exit test cases
  if (bot != 1) return 0; // not coprime (GCD != 1)
  if ((jacobi ^ (jacobi >> 1)) & 2) return -1; // quadratic resiprocity
  return 1; // square or otherwise fails resiprocity
}

// Jacobi symbol for 64-bit long input
Int_t lljacobiSymbol(Long64_t top, Long64_t bot) {
  if (bot < 0 || bot % 2 == 0) // Throw error on bad input
    std::cerr << "Error in <modmath_tools>::jacobiSymbol: Bad input denom = " << bot << std::endl;
  
  // Since jacobi symbol is equivalent mod the lower index,
  // we can reduce top to the smallest positive modulus
  top = (top % bot + bot) % bot;
  
  // For more efficiency, we can use bitwise operations on jacobi.
  // This requires starting with zero (as a logical bit)
  Long64_t jacobi = 0;
  while (top != 0) {
    while (top % 4 == 0) top /= 4; // squares never contribute
    
    // If there is a remaining factor of 2, see if it adds a minus sign
    if (top % 2 == 0) {
      // efficiently compute mod8 condition
      // (2,n) = -1 only if n % 8 = {3,5}
      jacobi ^= bot; //bitwise XOR
      top /= 2;
    }
    
    // Flip sign if top, bot are 3 mod 4.
    // As a logic operation, we require top and bot have the same
    // final 2 bits and flips only so long as those bits are ...11 (...3)
    jacobi ^= top & bot & 2;
    
    // perform Euclidean algorithm (swap top/bot, mod and iterate)
    Long64_t rem = bot % top;
    bot = top;
    top = rem;
  }
  
  // loop exit test cases
  if (bot != 1) return 0; // not coprime (GCD != 1)
  if ((jacobi ^ (jacobi >> 1)) & 2) return -1; // quadratic resiprocity
  return 1; // square or otherwise fails resiprocity
}

// Print a table of jacobi symbols
// A test method to ensure the Jacobi symbols are working,
// or to print a reference table of the values.
bool printJacobi() {
  // Some terminal print colors
  char cyan[10] = "\033[36m";
  char magn[10] = "\033[35m";
  char yell[10] = "\033[33m";
  char defl[10] = "\033[39m";
  
  int max_numerator = 30;
  int max_denominator = 59;
  
  std::cout << "   ";
  for (int a = 1; a <= max_numerator ; ++a) {
    char thisint[5];
    snprintf(thisint, 5, "%3d", a);
    std::cout << thisint;
  }std::cout << std::endl;
  
  for (int n = 1; n <= max_denominator; n += 2) {
    char thisint[5];
    snprintf(thisint, 5, "%3d", n);
    std::cout << thisint;
    
    for (int a = 1; a <= max_numerator; ++a) {
      int jacobi_test = ljacobiSymbol(a, n);
      switch (jacobi_test) {
        case -1:
          std::cout << cyan << " -1" << defl;
          break;
        case 0:
          std::cout << yell << "  0" << defl;
          break;
        case 1:
          std::cout << magn << "  1" << defl;
      }
    }std::cout << std::endl;
  }
  return 1;
}


//-----------------------------------------------GCD/LCM

/// The Extended Euclidean algorithm is an algorithm dated to
/// Euclid's Elements, modified to efficiently also include computation
/// of the Bezout coefficients (the uniquely minimal coefficients (c, d)
/// satisfying gcd(a, b) = c a + d b) efficiently. These coefficients are
/// useful to the Chinese Remainder Theorem and elsewhere for computing
/// the multiplicative modular inverse of a number in some base.
///
/// The algorithm works by successively subtracting the largest multiple of
/// one number that fits within another, which does not modify the GCD, until
/// zero is obtained, where the previous iteration then gives the GCD. By also
/// performing arithmetic on other initialized coefs in the same way, we can also
/// compute the minimal Bezout coefficients.
///
/// Written for Problem 134 on 6/28/25.
GCD_Int_t* extendedEuclidean(Int_t a, Int_t b) {
  
  // Initialize arrays
  int remdr_vec[2];
  int bez_a_vec[2];
  int bez_b_vec[2];
  if (a > b) {
    remdr_vec[0] = a; remdr_vec[1] = b;
    bez_a_vec[0] = 1; bez_a_vec[1] = 0;
    bez_b_vec[0] = 0; bez_b_vec[1] = 1;
  }
  else {
    remdr_vec[0] = b; remdr_vec[1] = a;
    bez_a_vec[0] = 0; bez_a_vec[1] = 1;
    bez_b_vec[0] = 1; bez_b_vec[1] = 0;
  }
  
  // Peform the iterative subtraction (remainder/quotient gathering)
  while (remdr_vec[1] != 0) {
    // perform the division and keep the remainder (subtraction)
    std::div_t cdiv = std::div(remdr_vec[0], remdr_vec[1]);
    
    // reassign the remainder on this iteration
    remdr_vec[0] = remdr_vec[1];
    remdr_vec[1] = cdiv.rem;
    
    // iterate s and t, which track the Bezout coefficients (m-mod inverses)
    int bez_a_temp = bez_a_vec[0] - cdiv.quot*bez_a_vec[1];
    bez_a_vec[0] = bez_a_vec[1];
    bez_a_vec[1] = bez_a_temp;
    
    int bez_b_temp = bez_b_vec[0] - cdiv.quot*bez_b_vec[1];
    bez_b_vec[0] = bez_b_vec[1];
    bez_b_vec[1] = bez_b_temp;
  }// iterate until rvec is zero, meaning GCD is done.
  
  // assign the results to the GCD.
  GCD_Int_t* sol = new GCD_Int_t;
  sol->first = remdr_vec[0];
  sol->second.first = bez_a_vec[0];
  sol->second.second = bez_b_vec[0];
  
  return sol;
}

// Version with Long32_t
GCD_Long32_t* lextendedEuclidean(Long32_t a, Long32_t b) {
  
  // Initialize arrays
  Long32_t remdr_vec[2];
  Long32_t bez_a_vec[2];
  Long32_t bez_b_vec[2];
  if (a > b) {
    remdr_vec[0] = a; remdr_vec[1] = b;
    bez_a_vec[0] = 1; bez_a_vec[1] = 0;
    bez_b_vec[0] = 0; bez_b_vec[1] = 1;
  }
  else {
    remdr_vec[0] = b; remdr_vec[1] = a;
    bez_a_vec[0] = 0; bez_a_vec[1] = 1;
    bez_b_vec[0] = 1; bez_b_vec[1] = 0;
  }
  
  // Peform the iterative subtraction (remainder/quotient gathering
  while (remdr_vec[1] != 0) {
    // perform the division and keep the remainder (subtraction)
    std::ldiv_t cdiv = std::ldiv(remdr_vec[0], remdr_vec[1]);
    
    // reassign the remainder on this iteration
    remdr_vec[0] = remdr_vec[1];
    remdr_vec[1] = cdiv.rem;
    
    // iterate s and t, which track the Bezout coefficients (m-mod inverses)
    Long32_t bez_a_temp = bez_a_vec[0] - cdiv.quot*bez_a_vec[1];
    bez_a_vec[0] = bez_a_vec[1];
    bez_a_vec[1] = bez_a_temp;
    
    Long32_t bez_b_temp = bez_b_vec[0] - cdiv.quot*bez_b_vec[1];
    bez_b_vec[0] = bez_b_vec[1];
    bez_b_vec[1] = bez_b_temp;
  }// iterate until rvec is zero, meaning GCD is done.
  
  // assign the results to the GCD.
  GCD_Long32_t* sol = new GCD_Long32_t;
  sol->first = remdr_vec[0];
  sol->second.first = bez_a_vec[0];
  sol->second.second = bez_b_vec[0];
  
  return sol;
}

// Version with Long64_t
GCD_Long64_t* llextendedEuclidean(Long64_t a, Long64_t b) {
  
  // Initialize arrays
  Long64_t remdr_vec[2];
  Long64_t bez_a_vec[2];
  Long64_t bez_b_vec[2];
  if (a > b) {
    remdr_vec[0] = a; remdr_vec[1] = b;
    bez_a_vec[0] = 1; bez_a_vec[1] = 0;
    bez_b_vec[0] = 0; bez_b_vec[1] = 1;
  }
  else {
    remdr_vec[0] = b; remdr_vec[1] = a;
    bez_a_vec[0] = 0; bez_a_vec[1] = 1;
    bez_b_vec[0] = 1; bez_b_vec[1] = 0;
  }
  
  // Peform the iterative subtraction (remainder/quotient gathering
  while (remdr_vec[1] != 0) {
    // perform the division and keep the remainder (subtraction)
    std::lldiv_t cdiv = std::lldiv(remdr_vec[0], remdr_vec[1]);
    
    // reassign the remainder on this iteration
    remdr_vec[0] = remdr_vec[1];
    remdr_vec[1] = cdiv.rem;
    
    // iterate s and t, which track the Bezout coefficients (m-mod inverses)
    Long64_t bez_a_temp = bez_a_vec[0] - cdiv.quot*bez_a_vec[1];
    bez_a_vec[0] = bez_a_vec[1];
    bez_a_vec[1] = bez_a_temp;
    
    Long64_t bez_b_temp = bez_b_vec[0] - cdiv.quot*bez_b_vec[1];
    bez_b_vec[0] = bez_b_vec[1];
    bez_b_vec[1] = bez_b_temp;
  }// iterate until rvec is zero, meaning GCD is done.
  
  // assign the results to the GCD.
  GCD_Long64_t* sol = new GCD_Long64_t;
  sol->first = remdr_vec[0];
  sol->second.first = bez_a_vec[0];
  sol->second.second = bez_b_vec[0];
  
  return sol;
}

// Compute the Least Common Multiple (LCM) of two integers.
// It's fastest to use the extended Euclidean algorithm to first
// find the GCD, then apply the identity LCM(a,b) = |a*b|/GCD(a,b).
//
// Written for Problem 129 on 7/27/2025
Int_t LCM(Int_t a, Int_t b) {
  GCD_Int_t* gcd = extendedEuclidean(a,b);
  return std::abs(a/gcd->first)*b;
}

// LCM for uLong32_t
uLong32_t lLCM(uLong32_t a, uLong32_t b) {
  GCD_Long32_t* gcd = lextendedEuclidean(a,b);
  return (a/gcd->first)*b;
}

// LCM for uLong64_t
uLong64_t llLCM(uLong64_t a, uLong64_t b) {
  GCD_Long64_t* gcd = llextendedEuclidean(a,b);
  return (a/gcd->first)*b;
}

// LCM for many integers
// There should be faster ways to do this, but for now
// it's straightforward to iterate the LCM
Int_t LCM(std::vector<Int_t>& vec) {
  if (vec.size() == 0) return 0;
  if (vec.size() == 1) return vec.back();
  
  Int_t a = vec[0];
  Int_t b = vec[1];
  for (int i = 1; i < vec.size(); ++i) {
    GCD_Int_t* gcd = extendedEuclidean(a,b);
    a = std::abs(a/gcd->first)*b; // LCM
    if (i != vec.size() - 1) b = vec[i+1];
  }return a;
}

uLong32_t LCM(std::vector<uLong32_t>& vec) {
  if (vec.size() == 0) return 0;
  if (vec.size() == 1) return vec.back();
  
  uLong32_t a = vec[0];
  uLong32_t b = vec[1];
  for (int i = 1; i < vec.size(); ++i) {
    GCD_Long32_t* gcd = lextendedEuclidean(a,b);
    a = (a/gcd->first)*b; // LCM
    if (i != vec.size() - 1) b = vec[i+1];
  }return a;
}

uLong64_t LCM(std::vector<uLong64_t>& vec) {
  if (vec.size() == 0) return 0;
  if (vec.size() == 1) return vec.back();
  
  uLong64_t a = vec[0];
  uLong64_t b = vec[1];
  for (int i = 1; i < vec.size(); ++i) {
    GCD_Long64_t* gcd = llextendedEuclidean(a,b);
    a = (a/gcd->first)*b; // LCM
    if (i != vec.size() - 1) b = vec[i+1];
  }return a;
}

//-----------------------------------------------Chinese Remainder Theorem

/// An algorithm to efficiently compute the Chinese Remainder.
/// A simple description of the theorem:
/// Given an integer x has remainder r1 mod M1 and remainder r2 mod M2,
/// the integer x is uniquely specified, mod the product M1 M2.
/// The minimal value (within 0, M1 M2 - 1) is given by the minimal Bezout coefficients.
Int_t chineseRemainder(Int_t rem1, Int_t rem2,
                         Int_t mod1, Int_t mod2) {
  // Efficiently obtain GCD and Bezout coefficients
  GCD_Int_t* bezout = extendedEuclidean(mod1, mod2);
  
  // Check that mods are coprime
  if (bezout->first != 1)
    std::cout << "Warning in Chinese Remainder: Input mods not coprime." << std::endl;
  
  // The bezout coefficients provide the solution
  Int_t checkmod = (rem1 * bezout->second.second * mod2
                      + rem2 * bezout->second.first * mod1);
  
  // Check this is the minimal solution
  if (checkmod > mod1*mod2) {
    std::cout << "Warning in Chinese Remainder: result not minimal." << std::endl;
    checkmod = checkmod % (mod1*mod2);
  } else if (checkmod < 0) {
    checkmod += (1 - checkmod/(mod1*mod2))*mod1*mod2;
  } return checkmod;
}

// Version with Long32_t
Long32_t lchineseRemainder(Long32_t rem1, Long32_t rem2,
                           Long32_t mod1, Long32_t mod2) {
  // Efficiently obtain GCD and Bezout coefficients
  GCD_Long32_t* bezout = lextendedEuclidean(mod1, mod2);
  
  // Check that mods are coprime
  if (bezout->first != 1)
    std::cout << "Warning in Chinese Remainder: Input mods not coprime." << std::endl;
  
  // The bezout coefficients provide the solution
  Long32_t checkmod = (rem1 * bezout->second.second * mod2
                       + rem2 * bezout->second.first * mod1);
  
  // Check this is the minimal solution
  if (checkmod > mod1*mod2) {
    std::cout << "Warning in Chinese Remainder: result not minimal." << std::endl;
    checkmod = checkmod % (mod1*mod2);
  } else if (checkmod < 0) {
    checkmod += (1 - checkmod/(mod1*mod2))*mod1*mod2;
  } return checkmod;
}

// Version with Long64_t
Long64_t llchineseRemainder(Long64_t rem1, Long64_t rem2,
                            Long64_t mod1, Long64_t mod2) {
  // Efficiently obtain GCD and Bezout coefficients
  GCD_Long64_t* bezout = llextendedEuclidean(mod1, mod2);
  
  // Check that mods are coprime
  if (bezout->first != 1)
    std::cout << "Warning in Chinese Remainder: Input mods not coprime." << std::endl;
  
  // The bezout coefficients provide the solution
  Long64_t checkmod = (rem1 * bezout->second.second * mod2
                       + rem2 * bezout->second.first * mod1);
  
  // Check this is the minimal solution
  if (checkmod > mod1*mod2) {
//    std::cout << "Warning in Chinese Remainder: result not minimal." << std::endl;
    checkmod = checkmod % (mod1*mod2);
  } else if (checkmod < 0) {
    checkmod += (1 - checkmod/(mod1*mod2))*mod1*mod2;
  } return checkmod;
}



#endif /* modmath_tools_h */
