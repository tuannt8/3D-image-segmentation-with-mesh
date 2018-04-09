//
//  probability_image.cpp
//  DEMO
//
//  Created by Tuan Nguyen Trung on 13/12/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#include "probability_image.hpp"
#include <iostream>
#include <fstream>

#include "tet_dis_coord.hpp"

using namespace std;

probability_image::probability_image()
{
    
}

probability_image::~probability_image()
{
    
}

void probability_image::draw()
{
    
}

void probability_image::load(std::string path)
{
    cache_binary(path);

}

void probability_image::cache_binary(std::string path)
{
    ostringstream cache_path;
    cache_path << path << "_cache";

    if (access(cache_path.str().c_str(), F_OK) != -1)
    {
        // Read general info
        try
        {
            std::ifstream f(path);
            
            if (f.is_open())
            {
                f >> m_dimension[0] >> m_dimension[1] >> m_dimension[2];
                f >> m_num_phases;
                
                f.close();
            }else{
                throw "No input file";
            }
        }
        catch (std::exception e)
        {
            std::cout << "Fail to read input file " << path << " with error: " << e.what();
            exit(1);
        }
        
        std::ifstream f_cache(cache_path.str(), std::ios::binary);
        long size = m_dimension[0] * m_dimension[1] * m_dimension[2] * sizeof(double);
        for (int i = 0; i < m_num_phases; i++)
        {
            std::shared_ptr<image3d> im3d = std::shared_ptr<image3d>(new image3d);
            im3d->init_raw(m_dimension.get());
            f_cache.read((char*)im3d->_voxels.data(), size);
            
            m_prob_map.push_back(im3d);
        }
    }
    else
    {
        // Load normally and write
        // 1. Load
        try
        {
            std::ifstream f(path);
            
            if (f.is_open())
            {
                f >> m_dimension[0] >> m_dimension[1] >> m_dimension[2];
                f >> m_num_phases;
                
                for (int i = 0; i < m_num_phases; i++)
                {
                    std::shared_ptr<image3d> im3d = std::shared_ptr<image3d>(new image3d);
                    im3d->load_raw(m_dimension.get(), f);
                    m_prob_map.push_back(im3d);
                }
                

            }
        }
        catch (std::exception e)
        {
            std::cout << "Fail to read input file " << path << " with error: " << e.what();
            exit(1);
        }
        
        // 2. Write
        long size = m_dimension[0] * m_dimension[1] * m_dimension[2] * sizeof(double);
        std::ofstream f_cache_write(cache_path.str(), ios::out  | std::ios::binary);
        for(int i = 0; i < m_num_phases; i++)
        {
            f_cache_write.write((char*)m_prob_map[i]->_voxels.data(), size);
        }
    }
    
    
    for (auto im : m_prob_map)
    {
        im->build_sum_table();
    }
    set_size();
}

vector<double> probability_image::get_avg_prob(std::vector<vec3> tet_points)
{
    vector<double> avg_iten;
    avg_iten.resize(m_num_phases);
    
    for (int i = 0; i < m_num_phases; i++)
    {
        avg_iten[i] = m_prob_map[i]->get_tetra_intensity(tet_points);
    }
    
    return avg_iten;
}
