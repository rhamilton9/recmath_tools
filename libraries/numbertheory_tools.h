// A library of number theoretic functions useful for computing
// properties of the integers and their related numbers.
//
// Utils for recreational computuing and math
// Created by Ryan Hamilton on Jan. 14, 2026

#ifndef numbertheory_tools_h
#define numbertheory_tools_h

// Some necessary dependencies for some functions
#include "../cosmetic/typestyle.h"           // Some typedef for common types
#include "primesieve.h"                      // Tools for prime number arithmetic

// Global variables helpful for efficient computations
std::vector<int> pentagonal_numbers; // Memo-ized generalized pentagonal numbers for Euler's Algorithm


//-----------------------------------------------Forward declare (all methods in file):





// Integer partitions
std::vector<int>* partitionEuler_initialize();
std::vector<int>* partitionEuler_iterate(std::vector<int>* partitions);



//----------------------------------------------- Divisor Sum






//----------------------------------------------- Integer Partitions

// Initialize the partition vector needed to compute integer
// partitions iteratively using Euler's Algorithm.
std::vector<int>* partitionEuler_initialize() {
  extern std::vector<int> pentagonal_numbers;
  if (pentagonal_numbers.size() == 0) {pentagonal_numbers.push_back(1); pentagonal_numbers.push_back(2);}
  
  std::vector<int>* part_vec = new std::vector<int>();
  part_vec->push_back(1);
  return part_vec;
}// End of numbertheory_tools::partitionEuler_initialize()



// Iterate using the Euler recurrence to find the next Partition
std::vector<int>* partitionEuler_iterate(std::vector<int>* partitions) {
  int n = partitions->size();
  extern std::vector<int> pentagonal_numbers;
  
  // Add all partition contributions relevant to the next partition
  int next_partition = partitions->back(); // p(n) = p(n-1) + ...
  int pentag_index = 0;
  while (n >= pentagonal_numbers[++pentag_index]) {
    
    // Make sure we have the relevant pentagonal numbers
    if (n <= pentagonal_numbers.back()) {
      int next_pentag_index = pentagonal_numbers.size()/2 + 1;
      if ((pentagonal_numbers.size() & 1) == 0) //pentagonal
        pentagonal_numbers.push_back(next_pentag_index*(3*next_pentag_index - 1)/2);
      else                                    // Generalized pentagonal
        pentagonal_numbers.push_back(next_pentag_index*(3*next_pentag_index + 1)/2);
      
    }// End of pentagonal append
    
    // Append the partition contribution
    if ((pentag_index & 3) < 2)
      next_partition += partitions->at(n - pentagonal_numbers[pentag_index]);
    else
      next_partition -= partitions->at(n - pentagonal_numbers[pentag_index]);
    
  }// End of recurrence loop
  
  if (next_partition < 0) {
    std::cout << "Error in <numbertheory_tools::partitionEuler_iterate>: Overflow! Larger type must be used." << std::endl;
    return partitions;
  }
  
  partitions->push_back(next_partition);
  return partitions;
}// End of numbertheory_tools::partitionEuler_iterate()







#endif /* numbertheory_tools_h */
