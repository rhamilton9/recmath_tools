// An import for the miscallaneous tools included
// in this directory and other standard libraries.
//
// Utils for recreational computuing and math
//
// Created by Ryan Hamilton on Aug. 26th, 2025
// ------------ Changelog ------------
// Jan. 17th, 2026     Updated diectory structure

#ifndef lib_rec_math_h
#define lib_rec_math_h

// Standard libraries
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <chrono>

// Library of my own functions made throughout this journey
#include "cosmetic/color.h"                // Terminal colored output
#include "cosmetic/typestyle.h"            // Some typedef for common types

#include "libraries/primesieve.h"          // Tools for prime number arithmetic
#include "libraries/vector_tools.h"        // Tools for std::vector manipulation
#include "libraries/modmath_tools.h"       // Tools for modular arithmetic
#include "libraries/diophantine_tools.h"   // Tools for solving diophantine equations
#include "libraries/numbertheory_tools.h"  // Implementations of some number theoretic functions

#include "classes/MPInt_t.cc"              // Multi-precision integer type

// Commonly used color modifiers, from color.h
Color::TextModifier t_red(Color::TEXT_RED);
Color::TextModifier t_grn(Color::TEXT_GREEN);
Color::TextModifier t_blu(Color::TEXT_BLUE);
Color::TextModifier t_def(Color::TEXT_DEFAULT);


#endif /* lib_rec_math_h */
