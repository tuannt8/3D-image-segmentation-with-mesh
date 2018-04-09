//
//  probability_image.hpp
//  DEMO
//
//  Created by Tuan Nguyen Trung on 13/12/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#ifndef probability_image_hpp
#define probability_image_hpp

#include <stdio.h>
#include "image3d.h"
#include <memory>

class probability_image
{
public:
    probability_image();
    ~probability_image();
    
    
    void load(std::string path);
    
    void draw();
    
    std::vector<double> get_avg_prob(std::vector<vec3> tet_points);
public:
    std::vector<std::shared_ptr<image3d>> m_prob_map;
    vec3i m_dimension;
    int m_num_phases;
    
public:
    void cache_binary(std::string path);
};

#endif /* probability_image_hpp */
