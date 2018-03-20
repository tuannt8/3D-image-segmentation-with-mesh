//
//  tet_dis_coord.hpp
//  DSC_segment
//
//  Created by Tuan Nguyen Trung on 2/9/16.
//  Copyright © 2016 Asger Nyman Christiansen. All rights reserved.
//

#ifndef tet_dis_coord_hpp
#define tet_dis_coord_hpp

#include <stdio.h>
#include <vector>

// tetrahedron sample points; upto
extern std::vector<std::vector<std::vector<double>>> tet_dis_coord;
extern std::vector<size_t> dis_coord_size;

// face sample points; up to 
extern std::vector<std::vector<std::vector<double>>> tri_dis_coord;
extern std::vector<size_t> tri_coord_size;

void set_size();

// Get position from barry centric coordinate
//  a is the barry centric coordinate, b is the triangle or tetrahedron vertices
#define get_coord(a,b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3])
#define get_coord_tri(a,b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

#endif /* tet_dis_coord_hpp */
