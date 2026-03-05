// A header file with some tools related to finding,
// enumerating, sieving, and otherwise working with prime numbers.
//
// Utils for recreational computuing and math
// Created by Ryan Hamilton on Dec. 9th, 2024.

#include <ostream>
#include "primesieve.h"

#ifndef prime_tools_h
#define prime_tools_h

//-----------------------------------------------Structs and global variables

bool debug_divisors = false;

struct divisorNode {
  int divisor;
  divisorNode* less;
  divisorNode* more;
};


divisorNode* tree_base;


//-----------------------------------------------Forward declare (all methods in file):

//------------------------Prime counters

// Sieve-based Prime Counter
uLong64_t primePi(uLong64_t x, uLong64_t y = 0);

//------------------------Number Theory Functions and Related

unsigned long divisorSum(int to_divide,
                         bool is_proper = false);

std::vector<int>* gatherDivisors(int to_divide,
                                 bool is_proper = false);

//------------------------Primality tests

bool pseudoprimeTestMR(uLong64_t to_test,
                       uLong64_t base);

bool primeTestMR(uLong64_t to_test,
                 bool check_small_primes = true);

bool primeTestBPSW(uLong64_t to_test,
                   bool check_small_primes = true);


//------------------------Test methods

bool checkPrimeTest();

//-----------------------------------------------Sieve-based Prime Counters


// Compute the number of primes less than a given number
// Without explicitly listing them. Implements a method
// by authors MEISSEL, LEHMER, LAGARIAS, MILLER, ODLYZKO
// Following optimizations suggested by M. DELEGLISE AND J. RIVAT
// See Ref.
//    MATHEMATICS OF COMPUTATION
//    Volume 65, Number 213
//    January 1996, Pages 235–245
// Written on 7/11/2025
uLong64_t primePi(uLong64_t x, uLong64_t y) {
  // default: 0 -> y = floor(x^1/3)
  if (y == 0) y = static_cast<uLong64_t>(std::floor(std::pow(x, 1./3)));
  
  // May be worthwhile to experiment with what value of y is really optimal,
  // and whether we can find an expression for such optimal values...
  
  // instantiate sieve array
  std::vector<int> primes;
  int countPrimes = 0;
  bool* isPrime = static_cast<bool*>(std::malloc(sizeof(int)*(y+1)));
  isPrime[0] = false; isPrime[1] = false;
  for (int inst = 2; inst <= y; ++inst) isPrime[inst] = true;
  
  // Sieve pi(z) up to z = y explicitly and store for later
  int sieve_pi[y+1];
  sieve_pi[0] = 0; sieve_pi[1] = 0;
  for (int i = 2; i <= y; ++i) {
    if (isPrime[i]) {
      ++countPrimes;
      primes.push_back(i);
      for (int j = 2*i; j <= y; j+=i)
        isPrime[j] = false;
    }
  }// not yet complete.
  
  return 0;
}// End of prime_tools::primePi

//-----------------------------------------------Number Theory Functions and Related

// Compute the Euler's Totient phi(n) from the number's prime factors
//     phi(n) = n * prod_{prime p | n} (1 - 1/p)
// This product formula is due to Euler.
uLong32_t eulerTotient(std::vector<int> prime_powers) {
  extern std::vector<int> extern_primes;
  if (extern_primes.size() < prime_powers.size()) return 0;
  
  uLong32_t totient = 1;
  for (int p = 0; p < prime_powers.size(); ++p) if (prime_powers[p] > 0) {
    totient *= extern_primes[p] - 1;
    for (int ppow = 2; ppow <= prime_powers[p]; ++ppow) totient *= extern_primes[p];
  }return totient;
}// End of prime_tools::eulerTotient



// Returns the sum of the divisors of a number.
// Requesting proper divisors includes all divisors except the number itself.
// To run quickly, the method expects a grid of sieved prime factors via
// primesieve::eratosthenes_withFactors();
// The method utilizes the product distributivity of the sum-of-divisors function.
// Should function effectively to about O(10^7) beyond which memory may become an issue.
//
// Note that this method was originally developed for Problem 21.
unsigned long divisorSum(int to_divide, bool is_proper) {
  extern factormap extern_factorlist;
  
  // Make sure there the factormap includes the desired number
  if (extern_factorlist.size() < to_divide) eratosthenes_withFactors(to_divide);
  std::vector<std::pair<int, short> > relevant_factors = extern_factorlist.at(to_divide);
  
  // Utilize the multiplicative distributivity of the divisorSum function
  // And a partial geometric series in the primes.
  unsigned long sum = 1;
  for (int p = 0; p < relevant_factors.size(); ++p) {
    // computing the power on the top of the pgs
    // Would be good to implement the efficient exponentiation of PE-122 once I do it...
    unsigned long pprod = 1;
    for (int w = 0; w <= relevant_factors[p].second; ++w) pprod *= relevant_factors[p].first;
    sum *= (pprod - 1);
  } for (int p = 0; p < relevant_factors.size(); ++p) {
    sum /= (relevant_factors[p].first - 1);
  }
  
  // Remove the number itself for proper divisors.
  if (is_proper) sum -= to_divide;
  return sum;
}// End of prime_tools::divisorSum


// Add a new element to the recursive divisor tree
//
// Helper method for prime_tools::gatherDivisors
void appendToDivisorTree(int new_divisor,
                         divisorNode* node = tree_base) {
  // Initializer/base case
  if (node == NULL) {
    if (debug_divisors) std::cout << "attempting to set the node!" << std::endl;
    if (node == tree_base) {// doesn't auto handle tree_base
      tree_base = new divisorNode;
      tree_base->divisor = new_divisor;
      return;
    }
    node = new divisorNode;
    node->divisor = new_divisor;
    return;
  }if (new_divisor == node->divisor) {
    if (debug_divisors) std::cout << "duplicate node " << new_divisor << std::endl;
    return; // no duplicates
  }
  
  // Case: smaller, go to less
  if (new_divisor < node->divisor) {
    if (debug_divisors) std::cout << "less than " << node->divisor << std::endl;
    if (node->less == NULL) { // append
      if(debug_divisors) std::cout << "appended new node " << new_divisor << std::endl;
      node->less = new divisorNode;
      node->less->divisor = new_divisor;
      return;
    }
    
    // continue searching
    appendToDivisorTree(new_divisor, node->less);
    return;
  }
  
  if (debug_divisors) std::cout << "more than " << node->divisor << std::endl;
  // Case larger, go to more
  if (node->more == NULL) { // append
    if (debug_divisors) std::cout << "appended new node " << new_divisor << std::endl;
    node->more = new divisorNode;
    node->more->divisor = new_divisor;
    return;
  }
  
  // continue searching
  appendToDivisorTree(new_divisor, node->more);
  return;
}// End of prime_tools::appendToDivisorTree

// Recursively add divisors to the tree.
// Necessary due to the large combinatoric nature of permutations
//
// Helper method for prime_tools::gatherDivisors
void makeAppend(std::vector<int> list, int mult = 1) {
  if (list.size() == 0) return;
  for (int i = 0; i < list.size(); ++i) {
    appendToDivisorTree(mult*list[i]);
    if (list.size() == 1) continue;
    makeAppend(*cutBelow(list, i), mult*list[i]);
  }return;
}// End of prime_tools::makeAppend

// Collapse the tree into a vector, deleting nodes as we go
//
// Helper method for prime_tools::gatherDivisors
void collapseDivisorTree(std::vector<int>* current_vector,
                         divisorNode* node = tree_base) {
  if (node == NULL) return; // base case
  
  // Continue to less first, recursively
  if (node->less != NULL) collapseDivisorTree(current_vector, node->less);
  
  // Append if none less
  current_vector->push_back(node->divisor);
  
  // continue to more first, after all less gathered
  if (node->more != NULL) collapseDivisorTree(current_vector, node->more);
  
  // all less and more added, delete this node.
  delete node;
  return;
}//End of prime_tools::collapseDivisorTree


// Return an ordered list of the divisors of an integer.
std::vector<int>* gatherDivisors(int to_divide, bool is_proper) {
  if (tree_base != NULL) tree_base = NULL; // Should add a safe deletion...
                                           // Requires the number to be factored
  extern factormap extern_factorlist;
  //  extern std::vector<int> extern_primes;
  if (extern_factorlist.size() < to_divide) return new std::vector<int>;
  
  // Form a binary search tree to assure there are no repetitions
  appendToDivisorTree(1);
  
  // Gather prime factors including multiplicities
  std::vector<int> prime_factors_multiplicity;
  for (int p = 0; p < extern_factorlist[to_divide].size(); ++p) {
    // Gather multiplicities for the main loop
    for (int m = 0; m < extern_factorlist[to_divide].at(p).second; ++m) {
      if (debug_divisors) std::cout << extern_factorlist[to_divide].at(p).first << ' ' << std::endl;
      prime_factors_multiplicity.push_back(extern_factorlist[to_divide].at(p).first);
    }
  }
  
  // Append permutations to the tree
  makeAppend(prime_factors_multiplicity);
  
  // collapse the tree
  std::vector<int>* divisors = new std::vector<int>();
  collapseDivisorTree(divisors);
  
  if (is_proper) divisors->pop_back();
  
  // debug: check divisors are in order
  if (debug_divisors) {
    std::cout << "Divisors of " << to_divide << ": { ";
    for (int i = 0; i < divisors->size(); ++i) {
      std::cout << divisors->at(i);
      if (i != divisors->size() - 1) std::cout << ',';
    }std::cout << " }" << std::endl;
  }
  
  // return
  return divisors;
}//End of prime_tools::gatherDivisors



//-----------------------------------------------Primality tests

// Miller-Rabin pseudoprimality/probable prime test
//
// This is a test based on Fermat's little theorem,
// but can sometimes produce composites which pass the test.
// The specific values of these pseudoprimes depend on the base,
// so these composites are called "pseudoprimes to base b"
//
// Assumes input is an odd integer
// It is possible to perform a rigorous primality test as a coincidence
// on the 7 bases {2, 325, 9375, 28178, 450775, 9780504, 1795265022},
// which works at least to 2^64
//
// TODO Fix for inputs larger than 2^32. Current issue: modulus is too large and yields
// overflow in multiplications which can then be larger 2^64. Solution: MPInt.
bool pseudoprimeTestMR(uLong64_t to_test, uLong64_t base) {
  if (to_test % 2 == 0) return false; // input was even
                                      // Find values satisfying n = d*2^s + 1 (n assumed odd)
  uLong64_t comp_modminus = to_test - 1;
  uLong64_t div_iterated = 2*comp_modminus;
  while (div_iterated % 2 == 0) {
    div_iterated /= 2;
    uLong64_t cmod = lfastPowMod(base, div_iterated, to_test);
    //debugging
    //std::cout << base << "^" << div_iterated << " % " << to_test << " = " << cmod << std::endl;
    if (cmod == comp_modminus) return true; // only possible for a pseudoprime
    if (cmod != 1) return false; // failed the test
  }return true; // passed the test
}//End of prime_tools::pseudoprimeTestMR

// An implementation of a coincidence primality test for the MR test
// using the coincidence of bases as defined above. It's somewhat slower
// than BPSW but helps me work around a current challenging bug in that
// method (due to squaring numbers exceeding 2^64 limit for 64-bit int)
bool primeTestMR(uLong64_t to_test, bool check_small_primes) {
  // check divisibility for a few small primes before launching into the test
  if (check_small_primes) {
    int small_primes[12] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int p = 0; p < sizeof(small_primes)/sizeof(int); ++p)
      if (to_test % small_primes[p] == 0) {
        if (to_test == small_primes[p]) return true;
        return false;
      }
  }// End of small prime divisibility test
  
  // A coincidence on the following primes is sufficient for 64-bit integers
  if (!pseudoprimeTestMR(to_test, 2)) return false;
  if (!pseudoprimeTestMR(to_test, 325)) return false;
  if (!pseudoprimeTestMR(to_test, 9375)) return false;
  if (!pseudoprimeTestMR(to_test, 28178)) return false;
  if (!pseudoprimeTestMR(to_test, 450775)) return false;
  if (!pseudoprimeTestMR(to_test, 9780504)) return false;
  if (!pseudoprimeTestMR(to_test, 1795265022)) return false;
  // currently broken, due to overflow in 64 bit long
  //  if (!pseudoprimeTestMR(to_test, 19)) return false;
  //  if (!pseudoprimeTestMR(to_test, 23)) return false;
  //  if (!pseudoprimeTestMR(to_test, 29)) return false;
  //  if (!pseudoprimeTestMR(to_test, 31)) return false;
  //  if (!pseudoprimeTestMR(to_test, 37)) return false;
  
  // passed coincidence
  return true;
}//End of prime_tools::primeTestMR


// Baillie-Pomerance-Selfridge-Wagstaff primality test
//
// A primality test based on a coincidence of two
// pseudoprime tests: Strong MR base 2 and Lucas Pseudoprime tests.
// Neither of these methods can guarantee primality but can guarantee
// compositeness. No counterexample of the combined test has been found,
// since the composite pseudoprimes that pass each test are quite different.
//
// Still, it is hypothesized that there are counterexamples for very large input,
// larger than 2^64 (the size limit for 64bit integer) so it is safe for our purposes.
bool primeTestBPSW(uLong64_t to_test, bool check_small_primes) {
  // check divisibility for a few small primes before launching into the test
  if (check_small_primes) {
    int small_primes[7] = {2, 3, 5, 7, 11, 13, 17};
    for (int p = 0; p < sizeof(small_primes)/sizeof(int); ++p)
      if (to_test % small_primes[p] == 0) {
        if (to_test == small_primes[p]) return true;
        return false;
      }
  }// End of small prime divisibility test
  
  // Perform base 2 probable prime test
  bool is_MR2_pseudoprime = pseudoprimeTestMR(to_test, 2);
  if (!is_MR2_pseudoprime) {
    //std::cout << "Failed the MR2 test" << std::endl;
    return false; // confirmed composite
  }
  
  
  // Prepare for Lucas prime test by gathering Jacobi symbols
  // We want the first D for which (D,n) = -1, where D is the
  // smallest magnitude constant within the oscillating sequence
  // D = {5, -7, 9, -11, 13, -15, ...}
  int jacobi_D = 5;
  int n_iterations = 0;
  while (true) {
    // safety condition in case of perfect square
    if (n_iterations > 7) {
      // Confirm the input is not a perfect square. If it is,
      // there will never be a value of D from Jacobi symbols
      // that satisfies the condition (D,a^2) = -1.
      uLong64_t checksq = static_cast<uLong64_t>(std::round(std::sqrt(to_test)));
      if (checksq*checksq == to_test) return false;// perfect square!
    }
    
    if (++n_iterations > 50)
      std::cerr << "Error in <primesieve>::primeTestBPSW: Long loop in Jacobi_D. Check for strange behavoir on input :: " << to_test << std::endl;
    
    // Check the Jacobi symbols via algo in modmath_tools
    if (lljacobiSymbol(jacobi_D, to_test) == -1) break;
    jacobi_D = -(jacobi_D + 2);
    if (lljacobiSymbol(jacobi_D, to_test) == -1) break;
    jacobi_D = - (jacobi_D - 2);
  }
  // Perform Lucas prime test
  uLong64_t lucas_prime_result = lllucasMod(to_test, 1, (1-jacobi_D)/4, true);
  if (lucas_prime_result != 0) {
    //std::cout << "Failed the Lucas test. D = " << jacobi_D << " was used." << std::endl;
    return false; // Failed the test
    
  }
  
  // Passed all tests
  return true;
}//End of prime_tools::primeTestBPSW



//-----------------------------------------------Test Methods

// A method to check that the primality test matches the results found by sieveing
bool checkPrimeTest() {
  int limit_tocheck = 10000000;
  std::vector<int>* primes = eratosthenes_toInt(limit_tocheck);
  
  int c_prime_index = 0;
  for (int i = 2; i < primes->back(); ++i) {
    if (i % 1000000 == 0) std::cout << "Correctly catagorized up to " << i << std::endl;
    bool isprime_fromtest = primeTestBPSW(i, true); // Replace with test of choice
    if (primes->at(c_prime_index) == i) {
      if (!isprime_fromtest) {
        std::cout << "Failed the test! Did not recognize ";
        std::cout << primes->at(c_prime_index) << " as prime!" << std::endl;
        return false;
      }
      ++c_prime_index;
      continue;
    } if (isprime_fromtest) {
      std::cout << "Failed the test! Incorrectly recognized ";
      std::cout << i << " as prime!" << std::endl;
      return false;
    }
  }
  std::cout << "Test successful! All numbers below " << limit_tocheck << " successfully categorized." << std::endl;
  return true;
}//End of prime_tools::checkPrimeTest



#endif /* prime_tools_h */
