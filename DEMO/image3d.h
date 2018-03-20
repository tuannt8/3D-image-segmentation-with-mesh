//
//  image3d.hpp
//  DSC
//
//  Created by Tuan Nguyen Trung on 12/7/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#ifndef image3d_hpp
#define image3d_hpp

#include <stdio.h>
#include <string>
#define cimg_use_png
#define cimg_display 1
#include <CImg.h>
#include <vector>
#include <cstdint>
#include "util.h"

#define BOUND_LABEL 999
#define BOUND_INTENSITY 0

typedef cimg_library::CImg<double> cimg_byte;


class image3d
{
    enum{X,Y,Z};
    
public:
    /** Constructor */
    image3d();
    ~image3d();
    
    /**
     Load images from directory.
     */
    void load(std::string path);
    
    // return average intensity
    // ouput intensity sum and volume optional
    double get_tetra_intensity(std::vector<vec3> tet_points, double * total_inten = nullptr, double * volume = nullptr);
    
    // Return variation to intensity c
    double get_variation(std::vector<vec3> tet_points, double c);
    
    double sum_area(int x, int y, int z);
    double sum_line_z(int x, int y, int z1, int z2);

    // Interpolation
    double get_value_f(vec3 pt) const;
    double get_value_f(double x, double y, double z) const{return get_value_f(vec3(x,y,z));};
    double get_value_f(int x, int y, int z) const{return get_value_f(vec3(x,y,z));};
    
    /**
     Get - set voxel; direct
     */
    double * get_layer(const int idx);
    double get_value (const int x, const int y, const int z) const;
    inline int index(int x, int y, int z) const{
        return z*_dim[X]*_dim[Y] + y*_dim[X] + x;
    }
    
    int* dimension(){return _dim;}
    const vec3 dimension_v() const{return vec3(_dim[X], _dim[Y], _dim[Z]);}
private:
    // Currently hold all images.
    int _dim[3]; // x - y - z
    int _layer_size; // Size of image in 1 layer
    std::vector<double> _voxels;
    std::vector<double> _sum_table; // This sumtable is used for line integration, not 3D rectangle integration.
private:
    void get_integral_recur(std::vector<vec3> const & tet_points, int loops, double * total, int deep);
    std::vector<std::vector<vec3>> subdivide_tet(std::vector<vec3> const & tet_points);
    
    void build_sum_table();
public:
    void generate_sample_point(int n);
    void generate_sample_point_tri(int n);
};

#endif /* image3d_hpp */

