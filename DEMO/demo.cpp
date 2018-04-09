//
//  Deformabel Simplicial Complex (DSC) method
//  Copyright (C) 2013  Technical University of Denmark
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  See licence.txt for a copy of the GNU General Public License.

#include "user_interface.h"
#include <iostream>
#include "InputParser.h"

using namespace std;
double noise_level;

//// Flag: Use original DSC in build table function
//// We have improved the build_table function. It has better performance
////  but sometimes run incorrectly. At this moment, just use the old function.
bool arg_b_build_table_origin = true;


int main(int argc, char** argv)
{
    InputParser p;
    
    if(argc != 2)
        p = InputParser(argc, argv);
    else
        p = InputParser(argv[1]);
    
    noise_level = atof( p.getCmdOption("-noise", "0").c_str());

    UI ui(p);
    
    return 0;
}
