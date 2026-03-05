// A header file with some tools for handling
// and performing operations on std::vector type
//
// Utils for recreational computuing and math
// Created by Ryan Hamilton on Apr. 26th, 2025

#include <stdio.h>
#include "../cosmetic/typestyle.h"

#ifndef vector_tools_h
#define vector_tools_h

//------------------------------------------------ Element handling

// Checks if a vector contains an element {
bool vectorContains(std::vector<int> &vec, int search) {
  for (std::vector<int>::iterator it = vec.begin(); it != vec.end(); ++it) {
    if (*it == search) return true;
  }return false;
}

// Add an element to an ordered vector only if it is not already present
void pushUnique(std::vector<int> &vec1, int to_push) {
  if (vec1.size() == 0 || to_push > vec1.back()) {vec1.push_back(to_push); return;}
  for (std::vector<int>::iterator it1 = vec1.begin(); it1 != vec1.end(); ++it1) {
    if (*it1 == to_push) return;
    if (*it1 > to_push) {vec1.emplace(it1, to_push); return;}
  }return;
}

// Add an element to an ordered vector only if it is not already present
void pushUnique(std::vector<uLong32_t> &vec1, uLong32_t to_push) {
  if (vec1.size() == 0 || to_push > vec1.back()) {vec1.push_back(to_push); return;}
  for (std::vector<uLong32_t>::iterator it1 = vec1.begin(); it1 != vec1.end(); ++it1) {
    if (*it1 == to_push) return;
    if (*it1 > to_push) {vec1.emplace(it1, to_push); return;}
  }return;
}

// Add an element to an ordered vector only if it is not already present
void pushUnique(std::vector<uLong64_t> &vec1, uLong64_t to_push) {
  if (vec1.size() == 0 || to_push > vec1.back()) {vec1.push_back(to_push); return;}
  for (std::vector<uLong64_t>::iterator it1 = vec1.begin(); it1 != vec1.end(); ++it1) {
    if (*it1 == to_push) return;
    if (*it1 > to_push) {vec1.emplace(it1, to_push); return;}
  }return;
}

// removes an item in the vector
int pop(std::vector<int> &vec, int position) {
  if (position > vec.size() ||
      position < 0) return 0;
  int to_return = vec[position];
  for (int i = position; i < vec.size(); ++i) vec[i] = vec[i+1];
  vec.pop_back();
  return to_return;
}

// removes an item and returns a copy of the vector
std::vector<int> popCopy(std::vector<int> vec, int position) {
  if (position > vec.size() ||
      position < 0) return *(new std::vector<int>());
  std::vector<int> copy;
  for (int i = 0; i < vec.size(); ++i) {
    if (i == position) continue;
    copy.push_back(vec[i]);
  }
  return copy;
}

std::vector<int> pushCopy(std::vector<int>& vec, int append) {
  std::vector<int> copy;
  for (std::vector<int>::iterator it = vec.begin(); it != vec.end(); ++it)
    copy.push_back(*it);
  copy.push_back(append);
  return copy;
}

// return duplicate vector with all entries below a point removed
std::vector<int>* cutBelow(std::vector<int>& vec, int position) {
  std::vector<int>* copy = new std::vector<int>();
  if (position < 0) return copy;
  for (int i = position+1; i < vec.size(); ++i) {
    if (i == position) continue;
    copy->push_back(vec[i]);
  }return copy;
}



//------------------------------------------------ Vector Operations

// Merge to vectors end-to-end
void appendVectors(std::vector<int> &vec1, std::vector<int> vec2) {
  vec1.insert(vec1.end(), vec2.begin(), vec2.end());
  return;
}

// Append only non redundant elements of a vector vec2 to vec1.
// Note that vec1 is assumed non-redundand and sorted.
void appendVectorsUnique(std::vector<int> &vec1, std::vector<int> vec2) {
  for (int i = 0; i < vec2.size(); ++i) pushUnique(vec1, vec2[i]);
  return;
}


//------------------------------------------------ std::string Operations

// Convert a string to lower case roman letters
inline void str_tolower(std::string& input_string) {
  for (std::string::iterator it = input_string.begin(); 
       it != input_string.end(); ++it) *it = std::tolower(*it);
  return;
}

// Convert a string to upper case roman letters
inline void str_toupper(std::string& input_string) {
  for (std::string::iterator it = input_string.begin();
       it != input_string.end(); ++it) *it = std::toupper(*it);
  return;
}


//------------------------------------------------ Printing tools

void printVector(std::vector<int>& vec) {
  if (vec.size() == 0) {std::cout << "{ }" << std::endl; return;}
  
  std::cout << '{';
  for (std::vector<int>::iterator it = vec.begin(); it != vec.end() - 1; ++it) std::cout << *it << ", ";
  std::cout << vec.back() << '}';
  return;
}

#endif /* vector_tools_h */
