// Functions for sieving primes
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



// Global Variables
const int sievebatch = 100000001;
uLong32_t c_factored_to = 0;
std::vector<int> extern_primes;
factormap extern_factorlist;


//-----------------------------------------------Forward declare (all methods in file):

//------------------------Prime sieves

// Eratosthenes-style sieves
std::vector<int>* eratosthenes_toInt(uLong32_t n,
                                     bool is_quiet = true);

int eratosthenes_toPrime(const int n,
                         const int nprime_array = 1000,
                         const bool is_quiet = true);

// Factorizing seives
void eratosthenes_withFactors(const int n,
                              bool is_quiet = true);

std::vector<uLong32_t>* eratosthenes_factorialFactors(const int n,
                                                      bool is_quiet = true);

factorlist factorInteger(uLong64_t to_factor,
                         bool is_quiet = true);


//------------------------Test methods

// Test methods
void printFactorList();

void printFactors(int i);


//-----------------------------------------------Eratosthenes-style sieves

/// The classic Sieve of Eratosthenes,
/// https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
///
/// This utility was first developed for Problems 3 and 10.
/// Modified to improve the sieve range (due to memory limitations) for Problem 800.
/// Modified to modestly improve efficiency when working on Problem 193.
std::vector<int>* eratosthenes_toInt(uLong32_t n, bool is_quiet) {
  // No need to mark multiples of primes larger than the sqrt of a number
  uLong32_t prop_lim = static_cast<uLong32_t>(std::ceil(std::sqrt(n)));
  
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
  extern uLong32_t c_factored_to;
  c_factored_to = n;
  
  return primes;
}//End of primesieve::eratosthenes_toInt

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
}//End of primesieve::eratosthenes_toPrime

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
}// End of primesieve::eratosthenes_withFactors

/// Sieve to an integer while additionally recording ALL
/// prime factors of all integer below that sieve limit.
///
/// Useful primarily for efficiently prime factorizing large factorials.
/// Information about the factors of individual integers is not recorded.
///
/// P926: First developed.
/// P231: Utility added to prevent re-sieving if used multiple times
///  also now returns a vector pointer instead of extern assignment.
std::vector<uLong32_t>* eratosthenes_factorialFactors(const int n, bool is_quiet) {
  if (n > sievebatch) { // edge case against static memory assignment
    std::cerr << "Error in <eratosthenes_toInt>: requested sieve bound is larger than local memory limit! Implement dynamic arrays, you coward!!" << std::endl;
    return new std::vector<uLong32_t>();
  }
  
  std::vector<uLong32_t> *factorial_factors = new std::vector<uLong32_t>();
  
  extern std::vector<int> extern_primes;
  
  if (n <= c_factored_to) { // already have the primes.
    
    // record total multiplicity of each prime in the factorial
    for (int p = 0; extern_primes[p] <= n && p < extern_primes.size(); ++p) {
      
      uLong32_t total_factors = 0;
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
      for (uLong32_t c = i; c <= n; c+=i) {
        isPrime[c] = false;
      }
      
      // record total multiplicity of this prime in the factorial
      uLong32_t total_factors = 0;
      for (int p = 1; p <= std::round(std::log(n)/std::log(i)); ++p) {
        total_factors += n/static_cast<int>(std::pow(i, p));
        // might be nice to replace this with the "efficient exponentiation" of PE-122 once I do it
      }
      
      factorial_factors->push_back(total_factors);
      
      // Print some information about the total factors
      if (!is_quiet) std::cout << n << "! has " << factorial_factors->back() << " factors of " << i << std::endl;
    }// end sieve loop
    c_factored_to = static_cast<uLong32_t>(n);
  }return factorial_factors;
}// End of primesieve::eratosthenes_factorialFactors


// returns a vector of prime factors of the input integer
factorlist factorInteger(uLong64_t to_factor, bool is_quiet) {
  extern std::vector<int> extern_primes;
  // Find primes to sqrt(to_factor). By the fundamnetal theorem of arithmetic,
  // a number can only have one prime factor larger than the square root of the number.
  int upper_bound = static_cast<int>(std::round(std::sqrt(to_factor)));
  if (extern_primes.size() == 0 || extern_primes.back() < upper_bound) eratosthenes_toInt(upper_bound, true);
  
  // Trial division by primes up to sqrt(to_factor)
  uLong64_t to_divide = to_factor;
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
}// End of primesieve::factorInteger


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
}// End of primesieve::printFactorList

void printFactors(int i) {
  extern factormap extern_factorlist;
  if (i > extern_factorlist.size()-1) return;
  
  std::cout << i << ": ";
  for (int f = 0; f < extern_factorlist[i].size(); ++f) {
    std::pair<int, short> c_pair = extern_factorlist[i][f];
    std::cout << c_pair.first << "^" << c_pair.second;
    if (f != extern_factorlist[i].size()-1) std::cout << ", ";
  }std::cout << std::endl;
}// End of primesieve::printFactors





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
