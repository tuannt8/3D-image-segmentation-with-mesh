//
//  image3d.cpp
//  DSC
//
//  Created by Tuan Nguyen Trung on 12/7/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#include "image3d.h"
#include <stdio.h>
#include <iostream>
// #include <boost/algorithm/string/case_conv.hpp>
// #include <boost/filesystem.hpp>
#include "CImg.h"
#include "tet_dis_coord.hpp"
#include <istream>
#include <streambuf>

#include <dirent.h>

using namespace std;
extern double noise_level;


image3d::image3d()
{

}

image3d::~image3d()
{

}
void image3d::init_raw(int dimension[3])
{
    _dim[X] = dimension[X];
    _dim[Y] = dimension[Y];
    _dim[Z] = dimension[Z];
    
    _voxels.resize(_dim[X]*_dim[Y]*_dim[Z]);
}
void image3d::load_raw(int dimension[3], std::ifstream& f)
{
    _dim[X] = dimension[X];
    _dim[Y] = dimension[Y];
    _dim[Z] = dimension[Z];
    
    _voxels.resize(_dim[X]*_dim[Y]*_dim[Z]);
    
    for (int i = 0; i < _dim[X]; i++)
    {
        for (int j = 0; j < _dim[Y]; j++)
        {
            for (int k = 0; k < _dim[Z]; k++)
            {
                f >> _voxels[index(i,j,k)];
            }
        }
    }
}

extern int num_images;
void image3d::load(std::string path)
{
    int count = num_images;
    if(num_images == 0)
    {
    DIR *dir;
    struct dirent *ent;

    if ((dir = opendir (path.c_str())) != NULL) {
        /* print all the files and directories within directory */
        while ((ent = readdir (dir)) != NULL) {
            std::string filename = ent->d_name;
            if(filename.size() < 4
               || strcasecmp(filename.substr(0,1).c_str(), ".") == 0)
                continue;
            
            auto ext = filename.substr(filename.size()-4, filename.size()-1);
            if(strcasecmp(ext.c_str(), ".PNG") == 0){
                count++;
            }
        }
        closedir (dir);
    } else {
        /* could not open directory */
        throw "Image directory not found";
    }
    
    std::cout << "Found " << count << " images " << std::endl;
    }
    
    std::vector<string> files;
    for (int i = 0; i< count; i++)
    {
        ostringstream is;
        is << path << "/im_" << i << ".png";
        files.push_back(is.str());
    }
    
    
     try{
         cimg_byte im;
         im.load(files[0].c_str());
    
         _dim[Z] = (int)files.size();
         _dim[X] = im.width();
         _dim[Y] = im.height();
    
         _voxels.resize(_dim[X]*_dim[Y]*_dim[Z]);
    
         unsigned int idx = 0;
         for (int i = 0; i < files.size(); i++)
         {
             cimg_byte im;
             im.load(files[i].c_str());

             if (noise_level > 0)
             {
                 im.noise(noise_level, 2);
             }
             
    
             for (int j = 0; j < im.height(); j++)
                 for(int i = 0; i < im.width(); i++)
                 {
                     _voxels[idx++] = (double)im(i,j)/255.0;
//                     assert(_voxels[idx-1] < 1.01 and _voxels[idx-1] >= 0);
                 }

         }
    
         set_size();
     }
     catch (exception e)
     {
         std::cout << "CImg: " << e.what() << endl;
         exit(9);
     }
     build_sum_table();
}

void image3d::build_sum_table()
{
    _sum_table.resize(_voxels.size());
// This code is for rectagle integration
//    for (int z = 0; z < _dim[2]; z++)
//    {
//        for (int y = 0; y < _dim[1]; y++)
//        {
//            for (int x = 0; x < _dim[0]; x++)
//            {
//                auto vv = _voxels[index(x, y, z)];
//                double v = _voxels[index(x, y, z)]
//                    + sum_area(x, y, z-1)
//                    + sum_area(x-1, y, z) - sum_area(x-1, y, z-1)
//                    + sum_area(x, y-1, z) - sum_area(x, y-1, z - 1)
//                    - (sum_area(x-1, y-1, z) - sum_area(x-1, y-1, z-1));
//                _sum_table[index(x, y, z)] = v;
//                assert(v >= 0);
//            }
//        }
//    }
    
    // We just need line integration
    for (int x = 0; x < _dim[0]; x++)
    {
        for(int y = 0; y < _dim[1]; y++)
        {
            _sum_table[index(x, y, 0)] = _voxels[index(x, y, 0)];
            for (int z = 1; z < _dim[2]; z++)
            {
                _sum_table[index(x, y, z)] = _sum_table[index(x, y, z-1)] + _voxels[index(x, y, z)];
            }
        }
    }
}

double image3d::sum_line_z(int x, int y, int z1, int z2)
{
//    double l2 = sum_area(x, y, z2) - sum_area(x-1, y, z2) - sum_area(x, y - 1, z2) + sum_area(x-1, y - 1, z2);
//    double l1 = sum_area(x, y, z1) - sum_area(x-1, y, z1) - sum_area(x, y - 1, z1) + sum_area(x-1, y - 1, z1);
//
//    return l2 - l1;
    
    assert(x >= 0);
    assert(x < _dim[0]);
    assert( y >= 0);
    assert( y < _dim[1]);
    assert(z1 >= 0);
    assert(z2 >= z1);
    
    if (z2 >= _dim[2])
    {
        z2 = _dim[2] -1;
        
        if (z1 >= _dim[2])
        {
            z1 = _dim[2] - 1;
        }
    }

    
    return _sum_table[index(x, y, z2)] - _sum_table[index(x, y, z1)];
}

double image3d::sum_area(int x, int y, int z)
{
    if(x < 0 or y < 0 or z < 0)
    {
        return 0;
    }

    if (!(x < _dim[0] and y < _dim[1] and z < _dim[2]))
    {
        return _sum_table.back();
    }

    // Make sure the function does not access out of range
#ifdef DEBUG
    assert(x < _dim[0] and y < _dim[1] and z < _dim[2]);
#endif
    return _sum_table[index(x, y, z)];
}

double image3d::get_value_f(vec3 pt) const
{

    CGLA::Vec3i pti(floor(pt[0]), floor(pt[1]), floor(pt[2]));

    vec3 ptif(pti[0], pti[1], pti[2]);
    vec3 relative_coord = pt - ptif;

    // TODO: write interpolate function
    double c00 = get_value(pti[0], pti[1], pti[2]) * (1 - relative_coord[0])
                + get_value(pti[0] + 1, pti[1], pti[2]) * relative_coord[0];
    double c01 = get_value(pti[0], pti[1], pti[2]+1) * (1 - relative_coord[0])
                + get_value(pti[0] + 1, pti[1], pti[2] + 1) * relative_coord[0];
    double c10 = get_value(pti[0], pti[1]+1, pti[2]) * (1 - relative_coord[0])
                + get_value(pti[0] + 1, pti[1] + 1, pti[2]) * relative_coord[0];
    double c11 = get_value(pti[0], pti[1] + 1, pti[2]  +1) * (1 - relative_coord[0])
                + get_value(pti[0] + 1, pti[1]+1, pti[2]+1) * relative_coord[0];


    double c0 = c00*(1-relative_coord[1]) + c10*relative_coord[1];
    double c1 = c01*(1-relative_coord[1]) + c11*relative_coord[1];

    double f =  c0*(1 - relative_coord[2]) + c1*relative_coord[2];
    assert(f < 1.01);
    return f;
}




void image3d::generate_sample_point(int n)
{
    double d = 1.0 / (double)n;
    vector<vec3> pts;
    pts.push_back(vec3(0,0,0));
    pts.push_back(vec3(0,1,0));
    pts.push_back(vec3(1,0,0));
    pts.push_back(vec3(0,0,1));

    vector<vector<vec3>> tets;

    for(int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < n; k++)
            {
                vec3 A1 = vec3(i*d, j*d, k*d);
                auto A2  = A1 + vec3(d, 0, 0);
                auto A3  = A1 + vec3(d, d, 0);
                auto A4  = A1 + vec3(0, d, 0);

                auto A5  = A1 + vec3(0, 0, d);
                auto A6  = A1 + vec3(d, 0, d);
                auto A7  = A1 + vec3(d, d, d);
                auto A8  = A1 + vec3(0, d, d);

                tets.push_back({A1,A4,A2,A5});
                tets.push_back({A5,A4,A2,A6});
                tets.push_back({A5,A6,A4,A8});
                tets.push_back({A2,A3,A4,A6});
                tets.push_back({A3,A4,A6,A8});
                tets.push_back({A3,A6,A7,A8});

            }
        }
    }
    std::cout<<"  //n = " << n << endl;
    cout << "  {" << endl;
    for(auto t : tets)
    {
        auto center = (t[0] + t[1] + t[2] + t[3]) / 4;
        auto f = center[0] + center[1] + center[2];
        if (f < 1)
        {
            auto coord = Util::barycentric_coords<double>(center, pts[0], pts[1], pts[2], pts[3]);
            printf("    {%f, %f, %f, %f},\n", coord[0], coord[1], coord[2], coord[3]);
        }
    }
    cout << "  }," << endl;
}

void image3d::generate_sample_point_tri(int n) {
    vector<vec3> tris;
    tris.push_back(vec3(1, 0, 0));
    tris.push_back(vec3(0, 1, 0));
    tris.push_back(vec3(0, 0, 1));

    int res = n;
    double step = 1.0/(double)res;

    printf(" { \n");
    printf(" // n = %d \n", n);

    for (int i = 0; i < res; i++) {
        double s1 = i*step;
        for (int j = res - i; j > 0; j--) {
            double s2 = j*step;

            // Triangle 1
            {
                double ep1 = s1 + step/3.0;
                double ep2 = s2 - step*2.0/3.;
                double ep3 = 1 - ep1 - ep2;
                auto pt = tris[0]*ep1 + tris[1]*ep2 + tris[2]*ep3;
                printf("  {%f, %f, %f},\n", pt[0], pt[1], pt[2]);
            }
            // Triangle 2
            {
                double ep1 = s1 + step*2/3;
                double ep2 = s2 - step*4/3;
                if (ep2 > 0) {
                    double ep3 = 1 - ep1 - ep2;
                    auto pt = tris[0]*ep1 + tris[1]*ep2 + tris[2]*ep3;
                    printf("  {%f, %f, %f},\n", pt[0], pt[1], pt[2]);
                }
            }

        }
    }

    printf("},\n");
}


double image3d::get_tetra_intensity(std::vector<vec3> tet_points, double * total_inten, double * volume)
{
    double v = Util::volume<double>(tet_points[0], tet_points[1], tet_points[2], tet_points[3]);

    long dis = std::upper_bound(dis_coord_size.begin(), dis_coord_size.end(), v) - dis_coord_size.begin() - 1;

    if(dis < 0)dis = 0;
    double total = 0;
    auto const & a = tet_dis_coord[dis];
    
    

    for (auto tb : a)
    {
        auto pt = get_coord(tet_points, tb);
        total += get_value(pt[0], pt[1], pt[2]);
    }
    

    total = total * v / a.size();


    if(total_inten)
        *total_inten = total;

    if (volume)
    {
        *volume = v;
    }

    return total / v;
}

double image3d::get_energy(std::vector<vec3> tet_points, double c)
{
    double v = Util::volume<double>(tet_points[0], tet_points[1], tet_points[2], tet_points[3]);
    
    
    long dis = std::upper_bound(dis_coord_size.begin(), dis_coord_size.end(), v) - dis_coord_size.begin() - 1;
    if(dis < 0) dis = 0;
    
    double total = 0;
    auto const a = tet_dis_coord[dis];
    
    for (auto tb : a)
    {
        auto pt = get_coord(tet_points, tb);
        total += (get_value(pt[0], pt[1], pt[2]) - c)*(get_value(pt[0], pt[1], pt[2]) - c);
    }
    
    total = total * v / a.size();
    
    return total;
}

double image3d::get_variation(std::vector<vec3> tet_points, double c)
{
    double v = Util::volume<double>(tet_points[0], tet_points[1], tet_points[2], tet_points[3]);


    long dis = std::upper_bound(dis_coord_size.begin(), dis_coord_size.end(), v) - dis_coord_size.begin() - 1;
    if(dis < 0) dis = 0;

    double total = 0;
    auto const a = tet_dis_coord[dis];
    
    for (auto tb : a)
    {
        auto pt = get_coord(tet_points, tb);
        total += std::abs(get_value(pt[0], pt[1], pt[2]) - c);
    }

    return total;
//    total = total * v / a.size();
//    return total / v;
}

void image3d::get_integral_recur(std::vector<vec3> const & tet_points, int loops, double * total, int deep)
{
    if (deep >= loops)
    {
        auto midp = (tet_points[0] + tet_points[1] + tet_points[2] + tet_points[3]) / 4.0;
        *total += get_value_f(midp);
        assert(*total < 1000);
        return;
    }
    else{
        auto subdivisions = subdivide_tet(tet_points);
        for (auto e : subdivisions)
        {
            get_integral_recur(e, loops, total, deep + 1);
        }
    }

}

std::vector<std::vector<vec3>> image3d::subdivide_tet(std::vector<vec3> const & tet_points)
{
    std::vector<std::vector<vec3>> list;

    auto A0 = tet_points[0];
    auto A1 = tet_points[1];
    auto A2 = tet_points[2];
    auto A3 = tet_points[3];

    auto B0 = (A0 + A1)/2;
    auto B1 = (A0 + A3)/2;
    auto B2 = (A2 + A3)/2;
    auto B3 = (A1 + A3)/2;
    auto B4 = (A0 + A2)/2;
    auto B5 = (A1 + A2)/2;

    list.push_back({A0, B0, B4, B1});
    list.push_back({B1, B3, B2, A3});
    list.push_back({B4, B5, A2, B2});
    list.push_back({B0, A1, B5, B3});
    list.push_back({B4, B0, B2, B1});
    list.push_back({B0, B3, B2, B1});
    list.push_back({B0, B2, B5, B3});
    list.push_back({B0, B2, B4, B5});

    return list;
}

double * image3d::get_layer(const int idx)
{
    return &_voxels[idx*_dim[X]*_dim[Y]];
}

double image3d::get_value(const int x, const int y, const int z) const
{
    if(x >= 0 && y >= 0 && z >= 0
           && x < _dim[0] && y < _dim[1] && z < _dim[2])
    {
        double f = _voxels[index(x,y,z)];
//        assert(f < 1.01 and f >= 0);
        return f;
    }

    return BOUND_INTENSITY;
}

