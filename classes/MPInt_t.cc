// A class defined multiprecision integer type
// Not as fast as GMP but written for my own learning
//
// Utils for recreational computuing and math
// ------------ Changelog ------------
// Created v1.0 by Ryan Hamilton on July 17th, 2025


#include "../cosmetic/typestyle.h"


#ifndef MPInt_t
#define MPInt_t

const uLong64_t B32 = 1L << 32; // base 2^32

// Class Definition for MPInt
//
// A multiprecision integer type which stores arbitrarily
// large numbers as an array of 32-bit integers. This array
// is effectively a Base-2^32 representation of the large number,
// with each array element representing a "digit" of the number.
//
// Arithmetic oprations are implemented as methods of the class.
// See below for a brief description of each method and its use.
class MPInt {
private:
  
  // *-- Class variables
  
  bool sgn;                         // true if negative, false (default) if positive
  std::vector<uLong32_t>* number;   // Vector of machine-words for base-2^32 representation
  
  // *-- Debug flags
  
  bool fDebug_string_constructor;   // Debug constructor for input decimal string
  bool fDebug_decimal_converter;    // Debug converter to decimal string
  
public:
  
  // *-- Constructors
  
  // Default: empty vector
  MPInt () {
    SetFDebugAll(false);
    this->sgn = false;
    this->number = new std::vector<uLong32_t>();
  }
  
  // Copy Constructor
  MPInt (const MPInt& other) {
    SetFDebugAll(false);
    sgn = other.GetSign();
    number = new std::vector<uLong32_t>();
    std::vector<uLong32_t>* other_num = other.GetArray();
    for (std::vector<uLong32_t>::iterator it = other_num->begin(); it != other_num->end(); ++it)
      number->push_back(*it);
  }
  
  // Input: Unsigned 64-bit integer
  MPInt (uLong64_t input_long) {
    SetFDebugAll(false);
    this->sgn = false;
    this->number = new std::vector<uLong32_t>();
    this->number->push_back(input_long & UINT_MAX); // fast bitwise modulus
    if (input_long >> 32 != 0)
      this->number->push_back(input_long >> 32);    // fast bitwise division
  }
  
  // Input: Decimal string
  MPInt (const char* input_cstring) {
    SetFDebugAll(false);
    ConstructFromDecimalCString(input_cstring);
  }
  
  // Input: Decimal std::string
  MPInt (std::string input_string) {
    SetFDebugAll(false);
    ConstructFromDecimalCString(input_string.c_str());
  }
  
  // *-- Simple Getters
  
  size_t                    size()            const {return this->number->size();};
  size_t                    GetSize()         const {return this->number->size();};
  bool                      GetSign()         const {return this->sgn;}
  std::vector<uLong32_t>*   GetArray()        const {return this->number;}
  
  
  // *-- Simple Setters
  
  void SetSign(bool new_sgn)                                  {sgn = new_sgn; return;}
  void SetSgn(bool new_sgn)                                   {sgn = new_sgn; return;}
  void SetNumberFormatted(std::vector<uLong32_t>* vec)        {this->number = vec; return;}
  void SetFDebugStringConstructor(bool fnew)                  {fDebug_string_constructor = fnew; return;}
  void SetFDebugStringConverter(bool fnew)                    {fDebug_decimal_converter = fnew; return;}
  // add other debug flags as they're written
  
  void SetFDebugAll(bool fnew)  {
    fDebug_string_constructor = fnew;
    fDebug_decimal_converter = fnew;
    return;
  }
  
  // *-- Overwrite the number, with some special inputs
  
  void SetDecimalString(const char* input_cstring) {
    ConstructFromDecimalCString(input_cstring);
    return;
  }
  
  void SetHexString(const char* input_cstring) {
    ConstructFromHexCString(input_cstring);
    return;
  }
  
  // For Problem 16
  // Construct as a power of 2
  void SetPowerOf2(uLong64_t power2) {
    this->number = new std::vector<uLong32_t>();
    while (power2 > 32) {
      this->number->push_back(0);
      power2 -= 32;
    }
    this->number->push_back(1 << power2);
    return;
  }
  
  // For Problem 20
  // Construct as a factorial
  void SetFactorial(int factorial) {
    // would be neat to compare simple multiplication and exponentiation-combination and see which is faster/when to turn over...
    return;
  }
  
  // Construct as a Mersenne number (2^n - 1)
  void SetMersenne(int mersenne_index) {
    // push back UINT_MAX a bunch of times...
    return;
  }
  
  // For problem 56
  void SetPower(int power, int base) {
    
    return;
  }
  
  void SetRepunit(int repdigits, int base = 10) {
    if (base == 2) SetMersenne(repdigits);
    // (b^n - 1) / (b - 1)
    return;
  }
  
  // *-- C++ object operators
  
  // Reassignment
  MPInt& operator=(const MPInt& other) {
    std::vector<uLong32_t>* other_num = other.GetArray();
    if (other_num == this->number) { // check for same pointer address
      std::cout << "Caught self assignment!" << std::endl;
      return *this;
    }// Otherwise safe to copy.
    
    sgn = other.GetSign();
    number->clear();
    for (std::vector<uLong32_t>::iterator it = other_num->begin(); it != other_num->end(); ++it)
      number->push_back(*it);
    return *this;
  }
  
  // Addition of two MPInt types
  // TODO currently the operation is unsigned!
  MPInt& operator+=(const MPInt& other) {
    const std::vector<uLong32_t>* other_MP = other.GetArray();
    
    // base case--both MPInts have digits
    uLong32_t carry = 0;
    int i = 0;
    while (i < this->number->size() && 
           i < other_MP->size()) {
      uLong64_t temp = static_cast<uLong64_t>(this->number->at(i)) + other_MP->at(i) + carry;
      this->number->at(i) = temp & UINT_MAX;
      carry = temp >> 32;
      ++i;
    }// End of both digits
    
    // Case: This MPInt is longer than the other
    while (i < this->number->size()) {
      uLong64_t temp = static_cast<uLong64_t>(this->number->at(i)) + carry;
      this->number->at(i) = temp & UINT_MAX;
      carry = temp >> 32;
      ++i;
    }// End this digits
    
    // Case: The other MPInt is longer than this
    while (i < other_MP->size()) {
      uLong64_t temp = static_cast<uLong64_t>(other_MP->at(i)) + carry;
      this->number->push_back(temp & UINT_MAX);
      carry = temp >> 32;
      ++i;
    }// End other digits
    
    if (carry != 0) this->number->push_back(carry);
    
    return *this;
  }// End of MPInt::operator+=
  
  // Addition to produce a new MPInt type
  // TODO make faster than clone and add (write new vector as you add them)
  MPInt operator+(const MPInt& mp_add) {
    MPInt mp_new = *(this);
    mp_new += mp_add;
    return mp_new;
  }// End of MPInt::operator+
  
  
  // *-- Custom implementations of certain arithmetic operations
  
  // A simple instance of "single-digit" long multiplication
  void SimpleMultiply(uLong32_t to_multiply) {
    uLong32_t carry = 0;
    for (std::vector<uLong32_t>::iterator number_iterator = this->number->begin();
         number_iterator != this->number->end(); ++number_iterator) {
      uLong64_t temp = static_cast<uLong64_t>(to_multiply)*(*number_iterator) + carry;
      *number_iterator = temp & UINT_MAX;
      carry = temp >> 32;
    }
    if (carry != 0) this->number->push_back(carry);
    return;
  }//End of MPInt::SimpleMultiply
  
  // A simple implementation of the Karatsuba algorithm to square this number
  // The number previously stored is lost in the process
  void Square() {
    std::vector<uLong32_t>* number_initial = this->number;
    this->number = new std::vector<uLong32_t>();
    
    
  }
  
  
  // *-- Functional Utility
  
  // Create a copy of this MPInt object and return a pointer
  MPInt* Clone() const {
    MPInt* to_return = new MPInt();
    to_return->SetSign(this->sgn);
    std::vector<uLong32_t>* copy_vec = to_return->GetArray();
    for (std::vector<uLong32_t>::iterator number_iterator = this->number->begin();
         number_iterator != this->number->end(); ++number_iterator)
      copy_vec->push_back(*number_iterator);
    return to_return;
  }
  
  // Add the decimal digits of the number together and return the sum
  uLong64_t DecimalDigitSum() {
    std::string decimal_digits = GetDigits("decimal");
    uLong64_t digit_sum = 0;
    for (std::string::iterator it = decimal_digits.begin();
         it != decimal_digits.end(); ++it) {
      digit_sum += *it - '0';
    }return digit_sum;
  }
  
  
  // *-- Display Utility
  
  // Construct the digits of the number. Available options:
  //    - decimal/dec/10        Decimal (Base 10, default)
  //    - hexadecimal/hex/16    Hexadecimal (Base 16)
  //    - internal/B32/2^32     As stored (Base 2^32)
  std::string GetDigits(std::string option = "decimal") {
    str_tolower(option);
    if (option.compare("decimal") == 0 ||
        option.compare("dec") == 0 ||
        option.compare("10") == 0) { // Decimal digits: (Number form)
      return ConvertToDecimalString();
    } 
    if (option.compare("hexadecimal") == 0 ||
          option.compare("hex") == 0 ||
          option.compare("16") == 0) { // Hexadecimal digits: (Number form)
      // TODO convert to hex string
      
      return std::string();
    } if (option.compare("internal") == 0 ||
          option.compare("b32") == 0 ||
          option.compare("2^32") == 0) {// Internal B32 digits: (Array form)
      
      std::string outstring;
      if (this->sgn) outstring.append("(-) ");
      outstring.push_back('{');
      for (std::vector<uLong32_t>::iterator it = number->begin(); it != number->end(); ++it) {
        outstring.append(std::to_string(*it));
        if (it == number->end() - 1) break;
        outstring.append(", ");
      }
      outstring.push_back('}');
      return outstring;
    }
    
    // option not recognized
    std::cerr << "Error in <MPInt::GetDigits>: option flag not recognized. Check inputs." << std::endl;
    return std::string();
  }// End of MPInt::GetDigits
  
  // Print the output of GetDigits
  void PrintDigits(std::string option = "decimal", bool endline = true) {
    std::cout << GetDigits(option);
    if (endline) std::cout << std::endl;
    return;
  }
  
  // *-- Destructor
  
  // Use "delete" for dynamically allocated (i.e. with new keyword)
  ~MPInt() {
    delete this->number;
  }
  
  // *-- Internal Constructors
private:
  
  // Convert from input decimal string to base-2^32 vector MPInt format
  void ConstructFromDecimalCString(const char* input_cstring) {
    this->number = new std::vector<uLong32_t>();
    size_t string_size = std::strlen(input_cstring);
    
    // check that the string has content
    if (string_size == 0) {
      this->sgn = false;
      this->number->push_back(0);
      return;
    }
    
    // handle signed: first character
    this->sgn = input_cstring[0] == '-';
    
    // start with first digit
    this->number->push_back(input_cstring[string_size - 1] - '0');
    std::vector<uLong64_t> cpow_10;
    cpow_10.push_back(1);
    for (int p = string_size - 2; p >= 0 + this->sgn; --p) {
      
      // continue to the next decimal digit: multiply by 10
      uLong64_t carry = 0;
      for (int t = 0; t < cpow_10.size(); ++t) {
        uLong64_t placeholder = cpow_10.at(t) * 10 + carry;
        cpow_10.at(t) = placeholder & UINT_MAX;
        carry = placeholder >> 32;
        if (carry != 0 && t == cpow_10.size() - 1) {cpow_10.push_back(carry); break;}
      }// End of construction for next 10 power
      
      // append the current decimal digit
      for (int t = 0; t < cpow_10.size(); ++t) {
        // in the same digit
        if (t > this->number->size() - 1) this->number->push_back(0); // check if this b32 digit exisits
        
        // Append the current part of the current base10 digit
        uLong64_t c_digitpart = (static_cast<uLong64_t>(input_cstring[p] - '0')
                                 * cpow_10.at(t) ) + this->number->at(t);
        this->number->at(t) = c_digitpart & UINT_MAX;
        
        // overflow, must be here to avoid literal int overflow
        uLong32_t overflow = c_digitpart >> 32;
        int index_overflow = t + 1;
        while (overflow > 0) {// Perform the gradeschool addition algorithm
          if (index_overflow > this->number->size() - 1) {// new digit: no risk of overflow
            this->number->push_back(overflow);
            break;
          }
          
          // already content here, check for another round of overflow
          uLong64_t new_digitpart = overflow + this->number->at(index_overflow);
          this->number->at(index_overflow) = new_digitpart & UINT_MAX;
          overflow = new_digitpart >> 32;
          ++index_overflow;
        }// End of overflow
      }// End of current digit
      
      
      // debug printing
      if (fDebug_string_constructor) {
        std::cout << "current digit : " << input_cstring[p] << ", Construction at this point : ";
        PrintDigits("internal");
        
        // Print about the current power of 10 in the B32 representation
        std::cout << "\n Power 10^" << string_size - p << ": {";
        for (std::vector<uLong64_t>::iterator it = cpow_10.begin(); it != cpow_10.end() - 1; ++it) {
          std::cout << *it << ", ";
        }std::cout << cpow_10.back() << '}' << std::endl;
      }// End of debug print
    }// End of deximal digit construction
    return;
  }// End MPInt::ConstructFromDecimalCString
  
  
  
  // Convert from input hexadecimal string to base-2^32 vector MPInt format
  void ConstructFromHexCString(const char* input_cstring) {
    // TODO
    return;
  }// End MPInt::ConstructFromHexCString
  
  
  
  // convert to decimal (base 10) string
  std::string ConvertToDecimalString() {
    if (this->number->size() == 1) return std::to_string(this->number->at(0));
    std::vector<uLong64_t> decimal_digits;
    std::vector<uLong64_t> cpow_B32;
    decimal_digits.push_back(this->number->at(0));
    cpow_B32.push_back(1);
    
    // Loop over the internal B32 digits, convert to new base
    for (int i32 = 1; i32 < this->number->size(); ++i32) {
      
      // propagate forward power first, since 10 is a small base
      for (std::vector<uLong64_t>::iterator it = cpow_B32.begin(); it != cpow_B32.end(); ++it) *it *= B32;
      
      
      // carry forward for proper digits
      uLong64_t carry = 0;
      for (int t = 0; t < cpow_B32.size(); ++t) {
        std::lldiv_t divbase = lldiv(cpow_B32.at(t) + carry, 10);
        cpow_B32.at(t) = divbase.rem;
        carry = divbase.quot;
        if (carry != 0 && t == cpow_B32.size() - 1) cpow_B32.push_back(0);
      }// Completed construction of next 2^32 power in relevant base
      
      
      // propagate new digit
      for (int t = 0; t < cpow_B32.size(); ++t) {
        if (decimal_digits.size() <= t) decimal_digits.push_back(0);
        std::lldiv_t divbase = lldiv(cpow_B32.at(t) * this->number->at(i32)
                                     + decimal_digits.at(t) + carry, 10);
        decimal_digits.at(t) = divbase.rem;
        carry = divbase.quot;
      }
      
      // overflow, ok to keep out of loop since 10 is small
      int overflow_index = cpow_B32.size();
      while (carry != 0) {
        if (decimal_digits.size() <= overflow_index) decimal_digits.push_back(0);
        std::lldiv_t divbase = lldiv(carry + decimal_digits.at(overflow_index), 10);
        decimal_digits.at(overflow_index) = divbase.rem;
        carry = divbase.quot;
        ++overflow_index;
      }// Completed conversion of current B32 digit
      
      
      // debug printing
      if (fDebug_decimal_converter) {
        std::cout << "Current digit : " << this->number->at(i32) << std::endl;
        std::cout << "Decimal construction with this digit : ";
        for (std::vector<uLong64_t>::reverse_iterator it = decimal_digits.rbegin();
             it != decimal_digits.rend(); ++it) {
          std::cout << *it;
        }std::cout << std::endl;
        
        // Print about the current power of B32
        std::cout << "Power B32^" << i32 + 1 << ": ";
        for (std::vector<uLong64_t>::reverse_iterator it = cpow_B32.rbegin();
             it != cpow_B32.rend(); ++it) {
          std::cout << *it;
        }std::cout << std::endl << std::endl;
      }// End of debug print
    }// End new base construction
    
    
    // convert to ASCII characters
    std::string decimal;
    if (this->sgn) decimal.push_back('-');
    for (std::vector<uLong64_t>::reverse_iterator it = decimal_digits.rbegin();
         it != decimal_digits.rend(); ++it) {
      decimal.push_back(*it + '0');
    }return decimal;
  }// End MPInt::ConvertToDecimalString
  
  
  
  // convert to hexadecimal (base 16) string
  std::string ConvertToHexString() {
    std::string hexadecimal;
    
    
    
    return hexadecimal;
  }// End MPInt::ConvertToHexString
  
  
  
};

#endif /* MPInt_t */

