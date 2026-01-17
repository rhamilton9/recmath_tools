// A header file with some tools related to finding,
// enumerating, sieving, and otherwise working with prime numbers.
//
// Utils for recreational computuing and math
// Created by Ryan Hamilton on Dec. 9th, 2024.

#include <ostream>
#include "../cosmetic/typestyle.h"
#include "vector_tools.h"
#include "modmath_tools.h"
#include "diophantine_tools.h"

#ifndef primesieve_h
#define primesieve_h

typedef std::vector<std::pair<int, short> > factorlist;
typedef std::vector<std::vector<std::pair<int, short> > > factormap;

bool debug_divisors = false;

struct divisorNode {
  int divisor;
  divisorNode* less;
  divisorNode* more;
};

// Global Variables
const int sievebatch = 100000001;

uLong_32t c_factored_to = 0;
std::vector<int> extern_primes;
factormap extern_factorlist;

divisorNode* tree_base;

//-----------------------------------------------Forward declare (all methods in file):

//------------------------Prime sieves

// Eratosthenes-style sieves
std::vector<int>* eratosthenes_toInt(uLong_32t n,
                                     bool is_quiet = true);

int eratosthenes_toPrime(const int n,
                         const int nprime_array = 1000,
                         const bool is_quiet = true);

// Factorizing seives
void eratosthenes_withFactors(const int n,
                              bool is_quiet = true);

std::vector<uLong_32t>* eratosthenes_factorialFactors(const int n,
                                                      bool is_quiet = true);

//------------------------Prime counters

// Sieve-based Prime Counter
uLong_64t primePi(uLong_64t x, uLong_64t y = 0);

//------------------------Derivative utility

factorlist factorInteger(uLong_64t to_factor,
                        bool is_quiet = true);

unsigned long divisorSum(int to_divide,
                         bool is_proper = false);


std::vector<int>* gatherDivisors(int to_divide,
                                 bool is_proper = false);

//------------------------Primality tests

bool pseudoprimeTestMR(uLong_64t to_test, 
                       uLong_64t base);

bool primeTestMR(uLong_64t to_test,
                 bool check_small_primes = true);

bool primeTestBPSW(uLong_64t to_test,
                   bool check_small_primes = true);


//------------------------Test methods

// Test methods
void printFactorList();

void printFactors(int i);

bool checkPrimeTest();

//-----------------------------------------------Eratosthenes-style sieves

/// The classic Sieve of Eratosthenes,
/// https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
///
/// This utility was first developed for Problems 3 and 10.
/// Modified to improve the sieve range (due to memory limitations) for Problem 800.
/// Modified to modestly improve efficiency when working on Problem 193.
std::vector<int>* eratosthenes_toInt(uLong_32t n, bool is_quiet) {
  // No need to mark multiples of primes larger than the sqrt of a number
  uLong_32t prop_lim = static_cast<uLong_32t>(std::ceil(std::sqrt(n)));
  
  // Instantiate arrays
  std::vector<int>* primes = new std::vector<int>();
  bool* isPrime = static_cast<bool*>(std::malloc(sizeof(bool)*(n + 1)));
  isPrime[0] = false; isPrime[1] = false;
  for (int inst = 2; inst < n+1; ++inst) isPrime[inst] = true;
  
  // Sieve loop: Mark multiples of each new prime
  for (int i = 2; i <= n; ++i) {
    // If a number has not already been visited (multiple of previous prime)
    // then it is a new prime. Mark all larger multiples of this prime.
    // (This is the key idea of the Sieve of Eratosthenes)
    if (isPrime[i]) {
      primes->push_back(i);
      if (!is_quiet) std::cout << i << "\tis prime #" << primes->size() << std::endl;
      
      // No need to propagate to grid, prime is too large
      if (i > prop_lim) continue;
      
      for (int j = i*i; j <= n; j+=i)
        isPrime[j] = false;
    }
  }// End of sieve loop
  
  // Mark recorded primes in case they're useful elsewhere.
  // Will be much more natural once I convert this to a class.
  extern std::vector<int> extern_primes;
  extern_primes = *primes;
  extern uLong_32t c_factored_to;
  c_factored_to = n;
  
  return primes;
}

/// A dynamically allocated Sieve of Eratosthenes
/// which allows sieveng to a desired n'th prime.
///
/// This utility was first developed for Problem 7.
int eratosthenes_toPrime(const int n,
                         const int nprime_array,
                         const bool is_quiet) {
  int primes[n];
  int prime_lastmult[n];
  int countPrimes = -1;
  int n_loop = 0;
  
  
  while (countPrimes+1 < n) {
    // re-instantiate array for this iteration
    
    bool isPrime[nprime_array];
    for (int i = 0; i < nprime_array; ++i) isPrime[i] = true;
    
    if (n_loop == 0) {
      // base prime is 2
      
      isPrime[0] = false;
      isPrime[1] = false;
    } else {
      
      // Cover primes from previous iterations
      for (int p = 0; p <= countPrimes; ++p) {
        if (!is_quiet && n < 1000) std::cout << primes[p] << "\tptime #" << p+1 << "\tlast multiple (under " << n_loop*nprime_array << ") = " << prime_lastmult[p] << std::endl;
        //if (!is_quiet && n < 1000) std::cout << Form("%i\tprime #%i\tlast multiple (under %i) = %i", primes[p], p+1, n_loop*nprime_array, prime_lastmult[p]) << std::endl;
        for (int j = prime_lastmult[p]+primes[p]; j < nprime_array*(n_loop+1); j+=primes[p]) {
          prime_lastmult[p] = j;
          isPrime[j-nprime_array*n_loop] = false;
        }
      }
      if (!is_quiet && n < 1000) std::cout << "--------------------------------------------------------------------End of loop" << n_loop << std::endl;
    }
    
    
    for (int i = nprime_array*n_loop; i < nprime_array*(n_loop+1); ++i) { // label numbers in this array for which this prime is a divisor
      if (isPrime[i-nprime_array*n_loop]) {
        primes[++countPrimes] = i;
        for (int j = i; j < nprime_array*(n_loop+1); j+=i) {
          prime_lastmult[countPrimes] = j;
          isPrime[j-nprime_array*n_loop] = false;
        }
        if (countPrimes == n-1) break;
      }
    }
    ++n_loop;
  }
  
  if (!is_quiet) for (int p = 0; p < n; ++p) {
    //std::cout << Form("%i\tprime #%i\tlast multiple (under %i) = %i", primes[i], i+1, n_loop*nprime_array, prime_lastmult[i]) << std::endl;
    std::cout << primes[p] << "\tptime #" << p+1 << "\tlast multiple (under " << n_loop*nprime_array << ") = " << prime_lastmult[p] << std::endl;
  }
  
  if (n%100 == 11 || n%100 == 12 || n%100 == 13)
    std::cout << "The " << n << "th prime is " << primes[n-1] << std::endl;
  else switch (n%10) {
    case 1:
      std::cout << "The " << n << "st prime is %i" << primes[n-1] << std::endl;
      break;
    case 2:
      std::cout << "The " << n << "nd prime is %i" << primes[n-1] << std::endl;
      break;
    case 3:
      std::cout << "The " << n << "rd prime is %i" << primes[n-1] << std::endl;
      break;
    default:
      std::cout << "The " << n << "th prime is %i" << primes[n-1] << std::endl;
      break;
  }
  std::cout << "Requested prime found after " << n_loop << " loops with base " << nprime_array << "." << std::endl;
  return primes[n-1];
}

//-----------------------------------------------Factorizing Sieves

/// Sieve to an integer while additionally recording the
/// prime factors of each integer below that sieve limit.
/// 
/// P12: First developed.
/// P21: Increased data storage efficiency using factormaps.
/// P926: (5/23/25) Improved the sieve range with static arrays. Did not end up using in final 926 ver.
void eratosthenes_withFactors(const int n, bool is_quiet) {
  extern std::vector<int> extern_primes;
  extern factormap extern_factorlist;
  
  // no need to run if we already have the needed primes
  if (extern_factorlist.size() >= n) return;
  extern_factorlist.clear();
  
  int countPrimes = 0;
  static bool isPrime[sievebatch];
  isPrime[0] = false;
  isPrime[1] = false;
  for (int inst = 0; inst < n+1; ++inst) {
    isPrime[inst] = true;
    extern_factorlist.push_back(std::vector<std::pair<int, short> >());
  }
  
  for (int i = 2; i <= n; ++i) {
    if (isPrime[i]) {
      extern_primes.push_back(i);
      
      // Label composite integers in the array and record the degree of each factor
      
      for (int c = i; c <= n; c+=i) {
        isPrime[c] = false;
        int c_mult = c;
        int count_factor_multiplicity = 0;
        while (c_mult % i == 0) {
          ++count_factor_multiplicity;
          c_mult /= i;
        }// Recursively record the powers of the current prime
        extern_factorlist[c].push_back(std::make_pair(i, count_factor_multiplicity));
      }
      
      ++countPrimes;
      if (!is_quiet) std::cout << i << "\tis prime #" << countPrimes << std::endl;
    }
  }return; // Assignments are extern
}

/// Sieve to an integer while additionally recording ALL
/// prime factors of all integer below that sieve limit.
///
/// Useful primarily for efficiently prime factorizing large factorials.
/// Information about the factors of individual integers is not recorded.
///
/// P926: First developed.
/// P231: Utility added to prevent re-sieving if used multiple times
///  also now returns a vector pointer instead of extern assignment.
std::vector<uLong_32t>* eratosthenes_factorialFactors(const int n, bool is_quiet) {
  if (n > sievebatch) { // edge case against static memory assignment
    std::cerr << "Error in <eratosthenes_toInt>: requested sieve bound is larger than local memory limit! Implement dynamic arrays, you coward!!" << std::endl;
    return new std::vector<uLong_32t>();
  }
  
  std::vector<uLong_32t> *factorial_factors = new std::vector<uLong_32t>();
  
  extern std::vector<int> extern_primes;
  
  if (n <= c_factored_to) { // already have the primes.
    
    // record total multiplicity of each prime in the factorial
    for (int p = 0; extern_primes[p] <= n && p < extern_primes.size(); ++p) {
      
      uLong_32t total_factors = 0;
      for (int pow = 1; pow <= std::round(std::log(n)/std::log(extern_primes[p])); ++pow) {
        total_factors += n/static_cast<int>(std::pow(extern_primes[p], pow));
        // might be nice to replace this with the "efficient exponentiation" of PE-122 once I do it
      }
      
      factorial_factors->push_back(total_factors);
      
      // Print some information about the total factors
      if (!is_quiet) std::cout << n << "! has " << factorial_factors->back() << " factors of " << extern_primes[p] << std::endl;
    }//end prime loop
    
  } else { // need to append as we find the primes.
    std::cout << "Sieving primes in <eratosthenes_factorialFactors>!" << std::endl;
    std::cout << "c_factored_to = " << c_factored_to << ", n = " << n << std::endl;
    
    // Instantiate the eratosthenes array
    static bool isPrime[sievebatch];
    isPrime[0] = false;
    isPrime[1] = false;
    for (int inst = 0; inst < n+1; ++inst) isPrime[inst] = true;
    
    extern_primes.clear();
    
    // perform the Sieve of Eratosthenes and mark total factors as we go
    for (int i = 2; i <= n; ++i) if (isPrime[i]) {
      extern_primes.push_back(i);
      
      // Label composite integers in the array and record the degree of each factor
      for (uLong_32t c = i; c <= n; c+=i) {
        isPrime[c] = false;
      }
      
      // record total multiplicity of this prime in the factorial
      uLong_32t total_factors = 0;
      for (int p = 1; p <= std::round(std::log(n)/std::log(i)); ++p) {
        total_factors += n/static_cast<int>(std::pow(i, p));
        // might be nice to replace this with the "efficient exponentiation" of PE-122 once I do it
      }
      
      factorial_factors->push_back(total_factors);
      
      // Print some information about the total factors
      if (!is_quiet) std::cout << n << "! has " << factorial_factors->back() << " factors of " << i << std::endl;
    }// end sieve loop
    c_factored_to = static_cast<uLong_32t>(n);
  }return factorial_factors;
}


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
uLong_64t primePi(uLong_64t x, uLong_64t y) {
  // default: 0 -> y = floor(x^1/3)
  if (y == 0) y = static_cast<uLong_64t>(std::floor(std::pow(x, 1./3)));
  
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
}

//-----------------------------------------------Derivative Utility


// returns a vector of prime factors of the input integer
factorlist factorInteger(uLong_64t to_factor, bool is_quiet) {
  extern std::vector<int> extern_primes;
  // Find primes to sqrt(to_factor). By the fundamnetal theorem of arithmetic,
  // a number can only have one prime factor larger than the square root of the number.
  int upper_bound = static_cast<int>(std::round(std::sqrt(to_factor)));
  if (extern_primes.size() == 0 || extern_primes.back() < upper_bound) eratosthenes_toInt(upper_bound, true);
  
  // Trial division by primes up to sqrt(to_factor)
  uLong_64t to_divide = to_factor;
  factorlist factors;
  for (int p = 0; p < extern_primes.size(); ++p) {
    short count_factors = 0;
    while (to_divide % extern_primes.at(p) == 0) {
      if (!is_quiet)
        std::cout << to_divide << " is divisible by " << extern_primes.at(p) << ", giving ";
      
      to_divide /= extern_primes.at(p);
      if (!is_quiet) std::cout << to_divide << std::endl;
      ++count_factors;
    }
    
    if (count_factors > 0) factors.push_back(std::make_pair(extern_primes[p], count_factors));
    if (extern_primes[p] > to_divide || extern_primes[p] > upper_bound) break;
  }
  
  // remaining factor must be larger than sqrt{n} and is therefore prime
  if (to_divide > 1) {
    factors.push_back(std::make_pair(to_divide, 1));
    if (!is_quiet) std::cout << "remaining factor " << to_divide << " is the largest prime factor." << std::endl;
  }
  
  return factors;
}

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
}


// Add a new element to the tree
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
}

// Recursively add divisors to the tree.
// Necessary due to the large combinatoric nature of permutations
void makeAppend(std::vector<int> list, int mult = 1) {
  if (list.size() == 0) return;
  for (int i = 0; i < list.size(); ++i) {
    appendToDivisorTree(mult*list[i]);
    if (list.size() == 1) continue;
    makeAppend(*cutBelow(list, i), mult*list[i]);
  }return;
}

// Collapse the tree into a vector, deleting nodes as we go
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
}


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
}

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
bool pseudoprimeTestMR(uLong_64t to_test, uLong_64t base) {
  if (to_test % 2 == 0) return false; // input was even
  // Find values satisfying n = d*2^s + 1 (n assumed odd)
  uLong_64t comp_modminus = to_test - 1;
  uLong_64t div_iterated = 2*comp_modminus;
  while (div_iterated % 2 == 0) {
    div_iterated /= 2;
    uLong_64t cmod = lfastPowMod(base, div_iterated, to_test);
    //debugging
    //std::cout << base << "^" << div_iterated << " % " << to_test << " = " << cmod << std::endl;
    if (cmod == comp_modminus) return true; // only possible for a pseudoprime
    if (cmod != 1) return false; // failed the test
  }return true; // passed the test
}

// An implementation of a coincidence primality test for the MR test
// using the coincidence of bases as defined above. It's somewhat slower
// than BPSW but helps me work around a current challenging bug in that
// method (due to squaring numbers exceeding 2^64 limit for 64-bit int)
bool primeTestMR(uLong_64t to_test, bool check_small_primes) {
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
  // currently broken, likely due to overflow in integers
//  if (!pseudoprimeTestMR(to_test, 19)) return false;
//  if (!pseudoprimeTestMR(to_test, 23)) return false;
//  if (!pseudoprimeTestMR(to_test, 29)) return false;
//  if (!pseudoprimeTestMR(to_test, 31)) return false;
//  if (!pseudoprimeTestMR(to_test, 37)) return false;
  
  // passed coincidence
  return true;
}


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
bool primeTestBPSW(uLong_64t to_test, bool check_small_primes) {
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
      uLong_64t checksq = static_cast<uLong_64t>(std::round(std::sqrt(to_test)));
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
  uLong_64t lucas_prime_result = lllucasMod(to_test, 1, (1-jacobi_D)/4, true);
  if (lucas_prime_result != 0) {
    //std::cout << "Failed the Lucas test. D = " << jacobi_D << " was used." << std::endl;
    return false; // Failed the test
    
  }
  
  // Passed all tests
  return true;
}


//-----------------------------------------------Print/Test Utility

// A simple debug print of the factormap.
void printFactorList() {
  extern factormap extern_factorlist;
  for (int i = 0; i < extern_factorlist.size(); ++i) {
    std::cout << i << ": ";
    for (int f = 0; f < extern_factorlist[i].size(); ++f) {
      std::pair<int, short> c_pair = extern_factorlist[i][f];
      std::cout << c_pair.first << "^" << c_pair.second;
      if (f != extern_factorlist[i].size()-1) std::cout << ", ";
    }std::cout << std::endl;
  }return;
}

void printFactors(int i) {
  extern factormap extern_factorlist;
  if (i > extern_factorlist.size()-1) return;
  
  std::cout << i << ": ";
  for (int f = 0; f < extern_factorlist[i].size(); ++f) {
    std::pair<int, short> c_pair = extern_factorlist[i][f];
    std::cout << c_pair.first << "^" << c_pair.second;
    if (f != extern_factorlist[i].size()-1) std::cout << ", ";
  }std::cout << std::endl;
}

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
}



/* // To eventually have a class
class PrimeSieve {
  
  ColorVal color_list;
  
public:
  
  TextModifier(ColorVal color_input) : color_list(color_input) {}
  
  friend std::ostream&
  operator<<(std::ostream& os, const TextModifier& mod) {
    return os << "\033[" << mod.color_list << "m";
  }
};

}
*/

#endif /* primesieve_h */
