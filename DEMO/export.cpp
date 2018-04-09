//
//  export.cpp
//  DEMO
//
//  Created by Tuan Nguyen Trung on 09/03/2018.
//  Copyright © 2018 Asger Nyman Christiansen. All rights reserved.
//

#include "export.hpp"
#include <fstream>

namespace dsc_export
{
    using namespace std;
    using namespace DSC;
    
    dsc_ptr load_dsc(std::string file_name)
    {
        std::vector<vec3> points;
        std::vector<int>  tets;
        std::vector<int>  tet_labels;
        is_mesh::import_tet_mesh(file_name, points, tets, tet_labels);
        
        auto l_dsc = dsc_ptr(new DeformableSimplicialComplex<>(points, tets, tet_labels));
        
        l_dsc->set_avg_edge_length();

        return l_dsc;
    }
    
    void export_two_phase_fluid(std::string path)
    {
        auto dsc = load_dsc(path);
        
        string directory = path.substr(0, path.find_last_of("\\/"));
        string name = path.substr(path.find_last_of("\\/") + 1);
        string path_0 = directory + "/mesh0/" + name + ".obj";
        export_shared_bound(dsc, {vec3i(1,0,0), vec3i(1,2,0)}, path_0);
        string path_1 = directory + "/mesh1/" + name + ".obj";
//        export_shared_bound(dsc, {vec3i(2,0,0)}, path_1);
        export_shared_bound_no_bound(dsc, {vec3i(2,0,0)}, path_1, vec3(0.0), vec3(0.15, 0.15, 0.15));
    }
    void export_surface(dsc_ptr dsc, std::string path)
    {
        vector<bool> phase(11,0);
        for(auto tit = dsc->tetrahedra_begin(); tit != dsc->tetrahedra_end(); tit++)
        {
            phase[min(10, dsc->get_label(tit.key()))] = true;
        }
        
        for (int i = 0; i < 10; i++)
        {
            if(phase[i])
            {
                string path_i = path +  "/phase_" + std::to_string(i) + ".obj";
                export_surface(dsc, i, path_i);
            }
        }
    }
    void export_surface(std::string path)
    {
        auto dsc = load_dsc(path);

        vector<bool> phase(11,0);
        for(auto tit = dsc->tetrahedra_begin(); tit != dsc->tetrahedra_end(); tit++)
        {
            phase[min(10, dsc->get_label(tit.key()))] = true;
        }
        
        for (int i = 0; i < 10; i++)
        {
            if(phase[i])
            {
                string path_i = path.substr(0, path.length() - 4) + "_" + std::to_string(i) + ".obj";
                export_surface(dsc, i, path_i);
            }
        }
        
    }
    
    inline bool is_bound_v(vec3 p, vec3 ld, vec3 ru, double gap)
    {
//        return p[2] < gap;
        
        bool is_bound = false;
        for (int i = 0; i < 3; i++)
        {
            is_bound = is_bound || (p[i] < ld[i] + gap || p[i] > ru[i] - gap);
        }
        return is_bound;
    }
    
    void export_shared_bound_no_bound(dsc_ptr dsc, std::vector<vec3i> phases, std::string path, vec3 ld, vec3 ru)
    {
        cout << path << endl;
        
        std::vector<int> points_maps(dsc->get_no_nodes(), -1); // points_maps[dsc_idx] = new_idx
        std::vector<vec3> points;
        std::vector<vec3i> faces;
        
        int cur_point_idx = 0;
        for (auto fit = dsc->faces_begin(); fit != dsc->faces_end(); fit++)
        {
            if (!fit->is_interface())
                continue;
            
            // Ignore bound
            bool is_bound = true;
            double epsilon = dsc->get_avg_edge_length()*0.7;
            for(auto p : dsc->get_pos(dsc->get_nodes(fit.key())))
            {
                is_bound = is_bound && is_bound_v(p, ld, ru, epsilon);
            }
            if (is_bound)
            {
                continue;
            }
            
            auto tets = dsc->get_tets(fit.key());
            
            int labels[2] = {dsc->get_label(tets[0]), dsc->get_label(tets[1])};
            
            int idx_other;
            bool found = false;
            for (auto pp : phases)
            {
                if ( (pp[0] == labels[0] && pp[1] == labels[1])
                    ||  (pp[1] == labels[0] && pp[0] == labels[1])
                    )
                {
                    idx_other = pp[0] == labels[0]? 0 : 1;
                    found = true;
                    break;
                }
            }
            
            if (!found)
                continue;
            
            
            auto nodes = dsc->get_sorted_nodes(fit.key(), tets[idx_other]);
            vec3i face_map;
            for(int i = 0; i < 3; i++)
            {
                auto n = nodes[i];
                if(points_maps[n] == -1)
                {
                    points.push_back(dsc->get_pos(n));
                    points_maps[n] = cur_point_idx++;
                }
                face_map[i] = points_maps[n] + 1;
            }
            faces.push_back(face_map);
        }
        
        ofstream file(path);
        
        for(auto p : points)
            file << "v " << p[0] << " " << p[1] << " " << p[2] << endl;
        for (auto f : faces)
        {
            file << "f " << f[0] << " " << f[1] << " " << f[2] << endl;
        }
    }
    
    void export_shared_bound(dsc_ptr dsc, std::vector<vec3i> phases, std::string path)
    {
        cout << path << endl;
        
        std::vector<int> points_maps(dsc->get_no_nodes(), -1); // points_maps[dsc_idx] = new_idx
        std::vector<vec3> points;
        std::vector<vec3i> faces;
        
        int cur_point_idx = 0;
        for (auto fit = dsc->faces_begin(); fit != dsc->faces_end(); fit++)
        {
            if (!fit->is_interface())
                continue;
            
            auto tets = dsc->get_tets(fit.key());
            
            int labels[2] = {dsc->get_label(tets[0]), dsc->get_label(tets[1])};
            
            int idx_other;
            bool found = false;
            for (auto pp : phases)
            {
                if ( (pp[0] == labels[0] && pp[1] == labels[1])
                    ||  (pp[1] == labels[0] && pp[0] == labels[1])
                )
                {
                    idx_other = pp[0] == labels[0]? 0 : 1;
                    found = true;
                    break;
                }
            }
            
            if (!found)
                continue;
            
            
            auto nodes = dsc->get_sorted_nodes(fit.key(), tets[idx_other]);
            vec3i face_map;
            for(int i = 0; i < 3; i++)
            {
                auto n = nodes[i];
                if(points_maps[n] == -1)
                {
                    points.push_back(dsc->get_pos(n));
                    points_maps[n] = cur_point_idx++;
                }
                face_map[i] = points_maps[n] + 1;
            }
            faces.push_back(face_map);
        }
        
        ofstream file(path);
        
        for(auto p : points)
            file << "v " << p[0] << " " << p[1] << " " << p[2] << endl;
        for (auto f : faces)
        {
            file << "f " << f[0] << " " << f[1] << " " << f[2] << endl;
        }
    }
    
    void export_surface(dsc_ptr dsc, int phase, std::string path)
    {
        std::vector<int> points_maps(dsc->get_no_nodes(), -1); // points_maps[dsc_idx] = new_idx
        std::vector<vec3> points;
        std::vector<vec3i> faces;
        
        int cur_point_idx = 0;
        for (auto fit = dsc->faces_begin(); fit != dsc->faces_end(); fit++)
        {
            if (!fit->is_interface())
                continue;
            
            auto tets = dsc->get_tets(fit.key());
            
            int labels[2] = {dsc->get_label(tets[0]), dsc->get_label(tets[1])};
            
            if (labels[0] != phase && labels[1] != phase)
                continue;
            
            int idx_other = labels[0] == phase ? 0 : 1;
            
            auto nodes = dsc->get_sorted_nodes(fit.key(), tets[idx_other]);
            vec3i face_map;
            for(int i = 0; i < 3; i++)
            {
                auto n = nodes[i];
                if(points_maps[n] == -1)
                {
                    points.push_back(dsc->get_pos(n));
                    points_maps[n] = cur_point_idx++;
                }
                face_map[i] = points_maps[n] + 1;
            }
            faces.push_back(face_map);
        }
        
        ofstream file(path);
        
        for(auto p : points)
            file << "v " << p[0] << " " << p[1] << " " << p[2] << endl;
        for (auto f : faces)
        {
            file << "f " << f[0] << " " << f[1] << " " << f[2] << endl;
        }
    }
}
