//
//  segment_function.cpp
//  DSC_segment
//
//  Created by Tuan Nguyen Trung on 12/9/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#include "segment_function.h"
#include "tet_dis_coord.hpp"
#include <queue>

#include "profile.h"

#ifdef _WIN32 // WINDOWS
#include <GL/glut.h>
#include <GL/glew.h>
#elif defined(__APPLE__) // IOS
#include <OpenGL/gl3.h>
#include <GLUT/glut.h>
#else // LINUX
#include <GL/glew.h>
#include <GL/glut.h>
#endif

#include "otsu_multi.h"
#include "gaussian_points.h"
#include "export.hpp"

std::vector<bool> tet_touched;


using namespace std;


std::bitset<4> X_direction("0001");
std::bitset<4> Y_direction("0010");
std::bitset<4> Z_direction("0100");
std::vector<std::bitset<4>> direction_st = {std::bitset<4>("0001"), std::bitset<4>("0010"), std::bitset<4>("0100")};

inline std::bitset<4> get_direction(vec3 a)
{
    static double norm_length = 0.8;
    
    int count = 0;
    std::bitset<4> d("0000");
    for (int i = 0; i < 3; i++)
    {
        if (std::abs(a[i]) > norm_length)
        {
            count ++;
            d = d | direction_st[i];
        }
    }
    
    return d;
}

#define algin_pos(idx) \
if(pos[idx] < threshold) \
destination[idx] = 0; \
if(pos[idx] > domain_dim[idx]-1 - threshold) \
destination[idx] = domain_dim[idx] - 1; \
if(destination[idx] < threshold) destination[idx] = 0;\
if(destination[idx] > domain_dim[idx]-1 - threshold) destination[idx] = domain_dim[idx] - 1;

void segment_function::init()
{
#ifdef INTENSITY_IMAGE
    _img.load(_directory_path);
#else
    #if defined(__APPLE__) || defined(_WIN32)
        m_prob_img.load("../Large_data/Camilla/P_map.txt");
    #else
        m_prob_img.load("../../Large_data/Camilla/P_map.txt");
    #endif
    cout << "Done loading " << _directory_path << endl;
#endif
}

void segment_function::random_initialization()
{
    for (auto tit = _dsc->tetrahedra_begin(); tit != _dsc->tetrahedra_end(); tit++)
    {
        if (_dsc->get_label(tit.key()) != BOUND_LABEL)
        {
            int new_label = (rand() % NB_PHASE);
            _dsc->set_label(tit.key(), new_label);
        }
    }
}

// Threshold the histogram in a range to 3 regions
// Multilevel Thresholding for Image Segmentation through a Fast Statistical Recursive Algorithm
// Arora
void Obstu_threshold(std::vector<int> & thres_T1, std::vector<int> & thres_T2, std::vector<int> & input, int iter)
{
    if (iter == 0)
    {
        return;
    }
    
    int range_min = thres_T1.back();
    int range_max = thres_T2.front();
    
    double thres_cofficient = 1;
    // Compute mean
    int mean = 0, count=0;
    for (int i = range_min; i < range_max; i++)
    {
        mean += input[i]*i;
        count += input[i];
    }
    mean = mean/count;
    
    // Compute deviation
    double deviation = 0;
    for (int i = range_min; i < range_max; i++)
    {
        deviation += input[i]*std::abs(i - mean);
    }
    deviation = (deviation) / count;
    
    int K1 = mean - thres_cofficient*deviation;
    int K2 = mean + thres_cofficient*deviation;
    
    assert(K1 > 0 && K2 < range_max);
    
    thres_T1.push_back(K1);
    thres_T2.insert(thres_T2.begin(), K2);
    
    Obstu_threshold(thres_T1, thres_T2, input, iter-1);
}

std::vector<int> obstu_recursive(std::vector<int> input, int nb_phase)
{
    std::vector<int> thres_T1; thres_T1.push_back(0);
    std::vector<int> thres_T2; thres_T2.push_back(255);
    int nb_iter = int((nb_phase-1)/2);
    
    Obstu_threshold(thres_T1, thres_T2, input, nb_iter);
    
    thres_T1.insert(thres_T1.end(), thres_T2.begin(), thres_T2.end());
    return thres_T1;
}

void segment_function::initialization_discrete_opt()
{
    double portion_keep_for_opt = 0.7; // Depend how sparse the segmenting
    
    // Optimize the labels of the tetrahedral
#ifdef _DSC_ORIGIN_
    int no_tets = MAX_NUM_ELEMENT_MESH;
#else
    int no_tets = _dsc->get_no_tets_buffer();
#endif
    std::vector<double> total_intensity_per_tet(no_tets, -1.0);
    std::vector<double> volume_per_tet(no_tets, -1.0);
    std::vector<double> mean_inten_per_tet(no_tets, -1.0);
//    std::vector<double> variation_inten_per_tet(no_tets, -1.0);
    std::vector<int> labels(no_tets, -1);

    double max_variation = -INFINITY, min_variation = INFINITY;
    // 1. Compute mean intensity and variation of each tetrahedron
    for (auto tit = _dsc->tetrahedra_begin(); tit != _dsc->tetrahedra_end(); tit++)
    {
        if(_dsc->get_label(tit.key()) == BOUND_LABEL)
            continue;
#ifdef DSC_CACHE
        auto tet_nodes_pos = _dsc->get_pos(*_dsc->get_nodes_cache(tit.key()));
#else
        auto tet_nodes_pos = _dsc->get_pos(_dsc->get_nodes(tit.key()));
#endif
        double total_inten, volume;
        auto mean_inten = _img.get_tetra_intensity(tet_nodes_pos, &total_inten, &volume);

//        auto variation = _img.get_variation(tet_nodes_pos, mean_inten);

        mean_inten_per_tet[tit.key()] = mean_inten;
        total_intensity_per_tet[tit.key()] = total_inten;
        volume_per_tet[tit.key()] = volume;
//        variation_inten_per_tet[tit.key()] = variation;
//
//        max_variation = max(max_variation, variation);
//        min_variation = min(min_variation, variation);
    }

//    // 2. Remove tetrahedra that have high variations of intensity. More sparse, more tetrahedra to be optimized
//    // Use bin size 100 for histogram count
//    int bin_size = 200;
//    max_variation *= 1.01;
//    min_variation *= 0.99;
//    double his_step = (max_variation - min_variation)/(double)bin_size;
//    std::vector<int> his_count(bin_size, 0);
//    long total_count = 0;
//    for (auto & v : variation_inten_per_tet)
//    {
//        if(v > 0)
//        {
//            int step = (int)((v - min_variation)/his_step);
//            his_count[step] ++;
//            total_count++;
//        }
//    }
//
//    // Find the threshold
//    int index_for_thres = 0;
//    long count_cur = 0;
//    for (; index_for_thres < his_count.size(); index_for_thres++)
//    {
//        count_cur += his_count[index_for_thres];
//        if (count_cur > total_count*portion_keep_for_opt)
//        {
//            break;
//        }
//    }
//    double thres_hold = (index_for_thres+1)*his_step; // thres hold to remove high variation tets
//
    // make new list, tets with low variation
    std::vector<int> histogram_for_thresholding(256,0);
    for (int i = 0; i < mean_inten_per_tet.size(); i++)
    {
        // This tetrahedron is considered for relabeling
        int idx = (int)(mean_inten_per_tet[i]*256); // convert from [0, 1] to [0, 256] image

        if(idx > 255) idx = 255;
        histogram_for_thresholding[idx] ++;
    }

    // 3. Optimize the label
    // 3.1. Random initialize
    int nb_phases = NB_PHASE;
    vector<int> thres_hold_array = otsu_muti(histogram_for_thresholding, nb_phases);

    // Logging
    cout << "Thresholding with: ";
    for (auto tt : thres_hold_array)
    {
        cout << tt << "; ";
    }cout << endl;

    // Initialize the label
    for (auto tit = _dsc->tetrahedra_begin(); tit != _dsc->tetrahedra_end(); tit++)
    {
        if (_dsc->get_label(tit.key()) == BOUND_LABEL)
        {
            continue;
        }

        int mean_inten_tet = (int)(mean_inten_per_tet[tit.key()]*255);
        if(mean_inten_tet >= 255) mean_inten_tet = 255;



        auto v_pos = std::lower_bound(thres_hold_array.begin(), thres_hold_array.end(), mean_inten_tet);

        int label = (int)(v_pos == thres_hold_array.end())? thres_hold_array.size() : int(v_pos - thres_hold_array.begin());

        assert(label < NB_PHASE);
        _dsc->set_label(tit.key(), label+1);
    }
}

#ifndef INTENSITY_IMAGE
void segment_function::threshold_init_probability()
{
    for (auto tit = _dsc->tetrahedra_begin(); tit != _dsc->tetrahedra_end(); tit++)
    {
        if (_dsc->get_label(tit.key()) == BOUND_LABEL)
        {
            continue;
        }
        
        auto avg_prob = m_prob_img.get_avg_prob(_dsc->get_pos(_dsc->get_nodes(tit.key())));
        auto max_element = std::max_element(avg_prob.begin(), avg_prob.end());
        
        int label = (int)(max_element - avg_prob.begin());

        _dsc->set_label(tit.key(), label);
    }
}
#endif

void segment_function::initialze_segmentation()
{
//    /**
//     Hamster sample
//     */
//    // Initilization by thresholding
//    double thres = 0.6;
//    // initialize by thresholding
//    for (auto tit = _dsc->tetrahedra_begin(); tit != _dsc->tetrahedra_end(); tit++)
//    {
//        // test
//        double total, volume;
//        auto pts = _dsc->get_pos(_dsc->get_nodes(tit.key()));
//        double avgI = _img.get_tetra_intensity(pts, &total, &volume);
//        
//        assert(avgI < 1.01);
//        //assert(total < 1000);
//        
//        if (avgI > thres)
//        {
//            _dsc->set_label(tit.key(), 1);
//        }
//    }
    
    
//    // Analyzing
//    for (auto tit = _dsc->tetrahedra_begin(); tit != _dsc->tetrahedra_end(); tit++)
//    {
//        double total, volume;
//        auto pts = _dsc->get_pos(_dsc->get_nodes(tit.key()));
//        double avgI = _img.get_tetra_intensity(pts, &total, &volume);
//    }
    
    /**
     Fuel cells
     */
//    // Initialization by thresholding
//    double thres[] = {0.31, 0.57, 0.7};
//    for (auto tit = _dsc->tetrahedra_begin(); tit != _dsc->tetrahedra_end(); tit++)
//    {
//        if(_dsc->get_label(tit.key()) == BOUND_LABEL)
//            continue;
//
//        // test
//        double total, volume;
//        auto pts = _dsc->get_pos(_dsc->get_nodes(tit.key()));
//        double avgI = _img.get_tetra_intensity(pts, &total, &volume);
//
//        double var = _img.get_variation(pts, avgI);
//        if (var > 0.1)
//        {
//            continue;
//        }
//
//        assert(avgI < 1.01);
//
//        // find closest
//        int idx = 0;
//        if (std::abs(avgI - thres[1]) < 0.1 )
//        {
//            idx = 1;
//        }else if(std::abs(avgI - thres[2]) < 0.1 )
//        {
//            idx = 2;
//        }
//
//        if (idx != 0)
//        {
//            _dsc->set_label(tit.key(), idx);
//        }
//    }
}

void segment_function::update_vertex_stability()
{
    // By default, it is stable
    _vertex_stability_map = std::vector<int>(_dsc->get_no_nodes_buffer(), 1);
    
    for (auto vid = _dsc->nodes_begin(); vid != _dsc->nodes_end(); vid++)
    {
        auto dis = get_node_displacement(vid.key());
        if (dis.length() > 0.1) // not stable
        {
            _vertex_stability_map[vid.key()] = 0;
        }
    }
}

void segment_function::force_snapp()
{
    // work around boundary
    auto node_mem_size = _dsc->get_no_nodes_buffer();
    std::vector<unsigned int> is_bound_vertex(node_mem_size,0);
    std::vector<std::bitset<4>> direction_state(node_mem_size,std::bitset<4>("0000"));
    for(auto fid = _dsc->faces_begin(); fid != _dsc->faces_end(); fid++)
    {
        if (fid->is_interface())
        {
            auto tets = _dsc->get_tets(fid.key());
            if(_dsc->get_label(tets[0]) == BOUND_LABEL
               || _dsc->get_label(tets[1]) == BOUND_LABEL)
            {
                auto nodes_on_face = _dsc->get_nodes(fid.key());
                auto norm = _dsc->get_normal(fid.key());
                
                // SHOULD ALSO BASE ON ITS POSITION
                //  in some irregular cases
                std::bitset<4> direction = get_direction(norm);
                
                for (auto n : nodes_on_face)
                {
                    is_bound_vertex[(unsigned int)n] = 1;
                    direction_state[(unsigned int)n] = direction_state[(unsigned int)n] | direction;
                }
            }
        }
    }
    
    static vec3 domain_dim = _img.dimension_v();
    static double threshold = _dsc->AVG_LENGTH*0.5*0.6;
    double max_dis = 0;
    double max_boundary_dis = 0;
    
    std::vector<vec3> node_displacement(_dsc->get_no_nodes_buffer(), vec3(0));
    for (auto nid = _dsc->nodes_begin(); nid != _dsc->nodes_end(); nid++)
    {
        if ( _dsc->exists(nid.key())
            && (nid->is_interface() or nid->is_crossing()))
        {
            vec3 destination = nid->get_pos();
            
            // Align boundary
            auto pos = nid->get_pos();
            if (is_bound_vertex[nid.key()])
            {
                auto direction = direction_state[nid.key()];
                for (int i = 0; i < 3; i++)
                {
                    //                    if ((direction & direction_st[i]).to_ulong() != 0)
                    {
                        algin_pos(i);
                    }
                }
                
                _dsc->set_pos(nid.key(), destination);
            }
        }
    }
}

void segment_function::snapp_boundary(){
    // work around boundary
    auto node_mem_size = _dsc->get_no_nodes_buffer();
    std::vector<unsigned int> is_bound_vertex(node_mem_size,0);
//    std::vector<std::bitset<4>> direction_state(node_mem_size,std::bitset<4>("0000"));
//    for(auto fid = _dsc->faces_begin(); fid != _dsc->faces_end(); fid++)
//    {
//        if (fid->is_interface())
//        {
//            auto tets = _dsc->get_tets(fid.key());
//            if(_dsc->get_label(tets[0]) == BOUND_LABEL
//               || _dsc->get_label(tets[1]) == BOUND_LABEL)
//            {
//                auto nodes_on_face = _dsc->get_nodes(fid.key());
//                auto norm = _dsc->get_normal(fid.key());
//
//                // SHOULD ALSO BASE ON ITS POSITION
//                //  in some irregular cases
//                std::bitset<4> direction = get_direction(norm);
//
//                for (auto n : nodes_on_face)
//                {
//                    is_bound_vertex[(unsigned int)n] = 1;
//                    direction_state[(unsigned int)n] = direction_state[(unsigned int)n] | direction;
//                }
//            }
//        }
//    }
    
    boundary_vertices_displacements = vector<vec3>(_dsc->get_no_nodes_buffer(), vec3(0.0));
    
    // 2. Align the boundary vertices to the boundary
    static vec3 domain_dim = _img.dimension_v();
    static double threshold = _dsc->AVG_LENGTH*0.5*0.6;
    double max_dis = 0;
    double max_boundary_dis = 0;
    
    std::vector<vec3> node_displacement(_dsc->get_no_nodes_buffer(), vec3(0));
    for (auto nid = _dsc->nodes_begin(); nid != _dsc->nodes_end(); nid++)
    {
        if ( _dsc->exists(nid.key())
            && (nid->is_interface() or nid->is_crossing()))
        {
            auto dis = get_node_displacement(nid.key());
            assert(!isnan(dis.length()));
            vec3 destination = nid->get_pos() + dis;
            
            // Align boundary
            auto pos = nid->get_pos();
            if (m_vertex_bound[nid.key()])
            {
//                auto direction = direction_state[nid.key()];
                for (int i = 0; i < 3; i++)
                {
//                    if ((direction & direction_st[i]).to_ulong() != 0)
                    {
                        algin_pos(i);
                    }
                }
                
                max_boundary_dis = max(max_boundary_dis, (destination - pos).length());
                boundary_vertices_displacements[nid.key()] = destination - pos;
            }else
                max_dis = max(max_dis, (destination - pos).length());
            
            node_displacement[nid.key()] = destination - pos;
        }
    }
    
    static double stop_thres = 0.005*_dsc->get_avg_edge_length() * 0.5;
    if (max_dis < stop_thres)
    {
        cout << "================ DONE ==================" <<endl;
    }
    
    // Adapt time step
    cout << "Max displacement; bound: " << max_dis << "; " << max_boundary_dis << endl;
    double scale = 1.0;
    static double max_displacement = _dsc->get_avg_edge_length() * 0.5 * 0.4;
    if (max_dis > max_displacement)
    {
        scale = (double)max_displacement / max_dis;
        _dt *= scale;
        max_dis = max_displacement;
        cout << "Scale time step adapt to " << _dt << "; scale = " << scale << "; max displace " << max_displacement << endl;
    }

    _previous_dis.resize(_dsc->get_no_nodes_buffer(), vec3(0));
    _cur_dis.resize(_dsc->get_no_nodes_buffer(), vec3(0));
    _dt_adapt.resize(_dsc->get_no_nodes_buffer(), 1.0);
    double real_max = 0;
    for (auto nit = _dsc->nodes_begin(); nit != _dsc->nodes_end(); nit++)
    {
        if (nit->is_interface())
        {
            auto dis = node_displacement[nit.key()]*scale;
            if (m_vertex_bound[nit.key()] && dis.length() > max_dis)
            {
                dis = dis * (max_dis / dis.length());
            }
            
            _dsc->set_destination(nit.key(), nit->get_pos() + dis*_dt_adapt[nit.key()]);
            
            _cur_dis[nit.key()] = dis;
            
            real_max = max(real_max, dis.length()*_dt_adapt[nit.key()]);
        }

    }
    
    cout << "Real max displacement " << real_max << endl;
    
    for (int i = 0; i < _dt_adapt.size(); i++)
    {
        is_mesh::NodeKey nk(i);
        if (_dsc->exists(nk) && _dsc->get(nk).is_interface()
            && (!m_vertex_bound[i] ||
            (m_vertex_bound[i] && _dsc->get(nk).is_crossing()) )
            )
        {
            static double thres_minus = cos(M_PI/180.0 * 150);
            static double thres_positive = cos(M_PI/180.0 * 30);
            
            if (_cur_dis[i].length() > 0.0001)
            {
                if(!_dsc->cache.is_clean[i])
                {
                    _dsc->cache.is_clean[i] = new bool(true);
                    _dt_adapt[i] = 1.0;
                }
                
                _cur_dis[i].normalize();
                
                auto cos_sign = Util::dot(_cur_dis[i], _previous_dis[i]);
                
                if(cos_sign < thres_minus)
                    _dt_adapt[i] = max(_dt_adapt[i]*0.9, 0.1);
                if(cos_sign > thres_positive)
                    _dt_adapt[i] = min(_dt_adapt[i]*1.1, 3.0);
                
                
                _previous_dis[i] = _cur_dis[i];
            }
        }
    }
}

vec3 segment_function::get_node_displacement(is_mesh::NodeKey nkey)
{
    return (_forces[(long)nkey]  + m_alpha*_internal_forces[(long)nkey])*_dt;
}



double segment_function::get_energy_tet_assume_label(is_mesh::TetrahedronKey tkey, int assumed_label)
{
#ifdef DSC_CACHE
    auto nodes_pos = _dsc->get_pos(*_dsc->get_nodes_cache(tkey));
#else
    auto nodes_pos = _dsc->get_pos(_dsc->get_nodes(tkey));
#endif
 
#ifdef INTENSITY_IMAGE
    auto energy = 0;
    assert(0);
#else
    auto energy = m_prob_img.m_prob_map[assumed_label]->get_energy(nodes_pos, 1);
#endif
    for (auto fid : _dsc->get_faces(tkey))
    {
        auto cobound_tets = _dsc->get_tets(fid); // We assume that this face is not DSC boundary, as there is a gap between DSC boundary and the image domain
        cobound_tets -= tkey;
        assert(cobound_tets.size() == 1);
        auto other_label = _dsc->get_label(cobound_tets[0]);
        
        if (other_label != BOUND_LABEL && other_label != assumed_label)
        {
            energy += _dsc->area(fid)*m_alpha;
        }
    }
    
    return energy;
}

double segment_function::get_energy_tetrahedron(is_mesh::TetrahedronKey tkey, int assumed_label)
{
#ifdef DSC_CACHE
    auto nodes_pos = _dsc->get_pos(*_dsc->get_nodes_cache(tkey));
#else
    auto nodes_pos = _dsc->get_pos(_dsc->get_nodes(tkey));
#endif
    auto old_label = _dsc->get_label(tkey);
    assert(old_label!= BOUND_LABEL);
    auto energy = _img.get_variation(nodes_pos, _mean_intensities[assumed_label-1]);
    for (auto fid : _dsc->get_faces(tkey))
    {
        auto cobound_tets = _dsc->get_tets(fid); // We assume that this face is not DSC boundary, as there is a gap between DSC boundary and the image domain
        auto label0 = _dsc->get_label(cobound_tets[0]);
        auto label1 = _dsc->get_label(cobound_tets[1]);

        label0 = label0==old_label? label1 : label0;

        if (label0 != assumed_label)
        {
            energy += _dsc->area(fid)*m_alpha;
        }
    }

    return energy;
}

#include "tet_dis_coord.hpp"

bool compare_tet(point_to_capture &a, point_to_capture & b)
{
    return a.tet_key < b.tet_key;
}

void segment_function::remove_stable_proximity(std::vector<std::vector<double>> & barry_coord, const is_mesh::SimplexSet<is_mesh::NodeKey> & nodes)
{
    bool is_stable[4];
    for (int i = 0; i < 4; i++)
    {
//        is_stable[i] = _vertex_stability_map[nodes[i]];
        is_stable[i] = !(_dsc->get(nodes[i]).is_interface()
                         || _dsc->get(nodes[i]).is_boundary()
                         || _dsc->get(nodes[i]).is_crossing());
    }
    
    for (auto bit = barry_coord.begin(); bit != barry_coord.end();)
    {
        bool bViolate = false;
        for (int i = 0; i < 4; i++)
        {
            if (!is_stable[i] && (*bit)[i] > 0.3)
            {
                bViolate = true;
                break;
            }
        }
        
        if (bViolate)
        {
            bit = barry_coord.erase(bit);
        }
        else
            bit++;
    }
}

void segment_function::adapt_tetrahedra_1()
{

    std::vector<point_to_capture> * subdivide_tets = new std::vector<point_to_capture>();
    
    // 1. Find potential points for subdivision
    update_average_intensity();
    
    int num_relabel = 0;

    for (auto tit = _dsc->tetrahedra_begin(); tit != _dsc->tetrahedra_end(); tit++)
    {
        if( _dsc->get_label(tit.key()) == BOUND_LABEL
           || _dsc->volume(tit.key()) < min_V)
            continue;
        
        // And only stable vertices
        
        auto tet_nodes = _dsc->get_nodes(tit.key());
        auto tet_nodes_pos = _dsc->get_pos(tet_nodes);
        int label = _dsc->get_label(tit.key());
        
        int num_dis = (int)(_dsc->volume(tit.key()) / min_V);
        
        long dis = std::upper_bound(dis_coord_size.begin(), dis_coord_size.end(), num_dis) - dis_coord_size.begin() - 1;
        
        if(dis < 0)dis = 0;
        assert(dis < num_dis);
        auto a = tet_dis_coord[dis];
        
        remove_stable_proximity(a, tet_nodes);
        
        for (auto tb : a)
        {
            auto pt = get_coord(tb, tet_nodes_pos);
            
            int opt_phase = arg_min_phase_point(pt, min_edge, label);
            if(opt_phase != label)
            {
                num_relabel ++;
                
                subdivide_tets->push_back(point_to_capture(tit.key(), pt, opt_phase));
                break;
            }
        }
    }
    
    cout << "Number tet to subdivide: " << num_relabel << endl;
    
    // 2. Local subdivision
//    std::sort(subdivide_tets.begin(), subdivide_tets.end(), compare_tet);
    while (subdivide_tets->size() > 0)
    {
        std::queue<is_mesh::TetrahedronKey> queue_tets;
        cout << "recursive " << subdivide_tets->at(0).tet_key << "-------" << subdivide_tets->size() << " remain\n";
        
        recursive_divide(subdivide_tets, subdivide_tets->at(0).tet_key, 0, queue_tets);
        
//        static int c = 0;
//        if (c++ > 2)
//        {
//            break;
//        }
    }
    
    delete subdivide_tets;
}
void segment_function::adapt_tetrahedra()
{
    tet_touched = vector<bool>(_dsc->get_no_tets_buffer()*1.1, false);
    cout << "Update intensity\n";
    update_average_intensity();
    
    min_edge = 8;
    min_V = pow(min_edge, 3)/6;
    int num_relabel = 0;
    for (auto tit = _dsc->tetrahedra_begin(); tit != _dsc->tetrahedra_end(); tit++)
    {
        if(tet_touched[tit.key()] || _dsc->get_label(tit.key()) == BOUND_LABEL)
            continue;
        
        tet_touched[tit.key()] = true;
        
        cout << "Tet: " << tit.key() << "; ";
        auto tet_nodes = _dsc->get_nodes(tit.key());
        auto tet_nodes_pos = _dsc->get_pos(tet_nodes);
        int label = _dsc->get_label(tit.key());
        
        int num_dis = (int)(_dsc->volume(tit.key()) / min_V);
        
        long dis = std::upper_bound(dis_coord_size.begin(), dis_coord_size.end(), num_dis) - dis_coord_size.begin() - 1;

        if(dis < 0)dis = 0;
        auto const a = tet_dis_coord[dis];
        
        cout << a.size() << " discretized nodes\n";
        
        bool relabel = false;
        for (auto tb : a)
        {
            auto pt = get_coord(tb, tet_nodes_pos);
            
            int opt_phase = arg_min_phase_point(pt, min_edge, label);
            if(opt_phase != label)
            {
                relabel = true;
                cout << "Start recursive subdivision\n";
                // Perform subdivide the tetrahedron and
                recursive_subdivide(tit.key(), pt, opt_phase, min_V);
                break;
            }
        }
        
        if (relabel)
        {
            num_relabel ++;
            continue;
        }
    }
}

#ifndef INTENSITY_IMAGE
int segment_function::relabel_probability()
{
    int nb_relablel = 0;
    for (auto tid = _dsc->tetrahedra_begin(); tid != _dsc->tetrahedra_end(); tid++)
    {
        if (_dsc->get_label(tid.key()) == BOUND_LABEL)
        {
            continue;
        }
        
        int old_idx = _dsc->get_label(tid.key());
        
#ifdef DSC_CACHE
        auto nodes_pos = _dsc->get_pos(*_dsc->get_nodes_cache(tid.key()));
#else
        auto nodes_pos = _dsc->get_pos(_dsc->get_nodes(tid.key()));
#endif
        auto avg_prob = m_prob_img.get_avg_prob(nodes_pos);
        auto max_element = std::max_element(avg_prob.begin(), avg_prob.end());
        int new_idx = (int)(max_element - avg_prob.begin());
        
        // Check if it is worth relabeling
        if (new_idx != old_idx)
        {
            auto old_energy = get_energy_tet_assume_label(tid.key(), old_idx);
            auto new_energy = get_energy_tet_assume_label(tid.key(), new_idx);
            
            if (new_energy < old_energy)
            {
                _dsc->set_label(tid.key(), new_idx);
                nb_relablel++;
            }
        }
    }
    cout << nb_relablel << " relabeled -------" << endl;
    return nb_relablel;
}
#endif

void segment_function::relabel_tetrahedra()
{
    cout << "Relabeling --" << endl;


    for (auto tid = _dsc->tetrahedra_begin(); tid != _dsc->tetrahedra_end(); tid++)
    {
        if (_dsc->get_label(tid.key()) == BOUND_LABEL)
        {
            continue;
        }
#ifdef DSC_CACHE
        auto nodes_pos = _dsc->get_pos(*_dsc->get_nodes_cache(tid.key()));
#else
        auto nodes_pos = _dsc->get_pos(_dsc->get_nodes(tid.key()));
#endif

        
        auto mean_inten_tetra = _img.get_tetra_intensity(nodes_pos);


//        // We check if it is worth changing the label to the phase with closest mean intensity
//        double smallest_gap = INFINITY;
//        int label_of_closest_phase = -1;
//        for (int i = 0; i < _mean_intensities.size(); i++)
//        {
//            if (smallest_gap > std::abs(mean_inten_tetra - _mean_intensities[i]))
//            {
//                smallest_gap =  std::abs(mean_inten_tetra - _mean_intensities[i]);
//                label_of_closest_phase = i+1;
//            }
//        }
//
//        if (label_of_closest_phase != _dsc->get_label(tid.key()))
        for (int label_of_closest_phase = 1; label_of_closest_phase<=NB_PHASE; label_of_closest_phase++)
        {
            // Check if we reduce the energy
            auto old_energy = get_energy_tetrahedron(tid.key(), _dsc->get_label(tid.key()));
            auto new_energy = get_energy_tetrahedron(tid.key(), label_of_closest_phase);

            if (new_energy < old_energy)                {
                _dsc->set_label(tid.key(), label_of_closest_phase);
            }
        }
    }
}


// -1 to avoid singularity around the boudary

void segment_function::work_around_on_boundary_vertices()
{
//    _previous_dis.resize(_dsc->get_no_nodes_buffer(), vec3(0));
//    _cur_dis.resize(_dsc->get_no_nodes_buffer(), vec3(0));
//    _dt_adapt.resize(_dsc->get_no_nodes_buffer(), _dt);
//
//    // 1. Find boundary vertices
//
//#ifdef _DSC_ORIGIN_
//    int node_mem_size = MAX_NUM_ELEMENT_MESH;
//#else
//    auto node_mem_size = _dsc->get_no_nodes_buffer();
//#endif
//    std::vector<unsigned int> is_bound_vertex(node_mem_size,0);
//    std::vector<std::bitset<4>> direction_state(node_mem_size,std::bitset<4>("0000"));
//    for(auto fid = _dsc->faces_begin(); fid != _dsc->faces_end(); fid++)
//    {
//        if (fid->is_interface() && !fid->is_boundary())
//        {
//            auto tets = _dsc->get_tets(fid.key());
//            if(_dsc->get_label(tets[0]) == BOUND_LABEL
//               || _dsc->get_label(tets[1]) == BOUND_LABEL)
//            {
//                auto nodes_on_face = _dsc->get_nodes(fid.key());
//                auto norm = _dsc->get_normal(fid.key());
//
//                // SHOULD ALSO BASED ON ITS POSITION
//                //  in some irregular cases
//                std::bitset<4> direction = get_direction(norm);
//
//                for (auto n : nodes_on_face)
//                {
//                    is_bound_vertex[(unsigned int)n] = 1;
//                    direction_state[(unsigned int)n] = direction_state[(unsigned int)n] | direction;
//                }
//            }
//        }
//    }
//
//    // For debuging
//    d_direction_state = direction_state;
//    d_is_image_boundary = is_bound_vertex;
//
//    // 2. Align the boundary vertices to the boundary
//#ifdef INTENSITY_IMAGE
//    auto domain_dim = _img.dimension();
//#else
//    auto domain_dim = m_prob_img.m_dimension;
//#endif
//
//    double max_displacement_real = -INFINITY;
//
//    // For debuging
//    boundary_vertices_displacements = std::vector<vec3>(node_mem_size, vec3(0));
//
//    for (auto nid = _dsc->nodes_begin(); nid != _dsc->nodes_end(); nid++)
//    {
//
//        if ( (nid->is_interface() or nid->is_crossing())
//            && _dsc->exists(nid.key())
//            and !nid->is_boundary())
//        {
//            auto dis = get_node_displacement(nid.key());
//
//            assert(!isnan(dis.length()));
//
//            vec3 destination = nid->get_pos() + dis;
//            auto threshold = 0.5*_dsc->AVG_LENGTH;
//
//            // Align boundary
//            if (is_bound_vertex[nid.key()])
//            {
//                auto direct = direction_state[nid.key()];
//                auto pos = nid->get_pos();
//                if ((direct & X_direction).to_ulong() != 0) // constraint on x
//                {
//                    algin_pos(0);
//                }
//                if ((direct & Y_direction).to_ulong() != 0){ // constraint on y
//                    algin_pos(1);
//                }
//                if ((direct & Z_direction).to_ulong() != 0){ // constraint on z
//                    algin_pos(2);
//                }
//
//                boundary_vertices_displacements[nid.key()] = destination - nid->get_pos();
//
//                dis = destination - nid->get_pos();
//            }
//
//            _cur_dis[nid.key()] = dis;
//
//            _dsc->set_destination(nid.key(), destination);
//
//            if (max_displacement_real < dis.length())
//            {
//                max_displacement_real = dis.length();
//            }
//        }
//    }
//
//    // Update time step
//    // Adaptive time step
//    for (int i = 0; i < _previous_dis.size(); i++)
//    {
//        if(_dsc->cache.is_clean[i])
//        {
//            auto vp = _previous_dis[i];
//            auto vc = _cur_dis[i];
//            auto vpl = vp.length();
//            auto vcl = vc.length();
//
//            if (vpl > EPSILON
//                && vcl > EPSILON)
//            {
//                auto cosA = Util::dot(vp, vc)/vpl/vcl;
//                if ( cosA < -0.2) // 110o
//                {
//                    _dt_adapt[i] /= 2.0;
//                }
//                if (cosA > 0.2)
//                {
//                    _dt_adapt[i] *= 1.05;
//                }
//            }
//
//            if(vcl < EPSILON)
//                _dt_adapt[i] = _dt;
//        }
//        else{
//            _dsc->cache.is_clean[i] = new bool(true);
//            _dt_adapt[i] = _dt;
//        }
//
//    }
//    _previous_dis = _cur_dis;
//
//    double min_dt = INFINITY;
//    double max_dt = -INFINITY;
//    for (auto dtt : _dt_adapt)
//    {
//        min_dt = std::min(min_dt, dtt);
//        max_dt = std::max(max_dt, dtt);
//    }
//    cout<<"Min dt: " << min_dt << "; max: " << max_dt << endl;
//
//    // Log debug
//    cout << "Max displacement: " << max_displacement_real << endl;
//    double average_internal_force = 0, average_external_force = 0;
//    int count = 0;
//    for (auto nid = _dsc->nodes_begin(); nid != _dsc->nodes_end(); nid++)
//    {
//
//        if ( (nid->is_interface() or nid->is_crossing())
//            && _dsc->exists(nid.key())
//            and !nid->is_boundary()
//            && !is_bound_vertex[nid.key()])
//        {
//            if(_forces[nid.key()].length() > 0.001)
//            {
//                average_external_force += _forces[nid.key()].length();
//                average_internal_force += _internal_forces[nid.key()].length();
//
//                count ++;
//            }
//        }
//    }
//
//    cout << "average external force: " << average_external_force/count << "; average internal force: " << average_internal_force/count << endl;
}

void segment_function::compute_mesh_quality_control_force()
{
    // Adaptive force based on curvature
    std::vector<vec3> adaptive_force(_dsc->get_no_nodes_buffer(), vec3(0));
    
    std::vector<std::vector<vec3>> curvature_force(_dsc->get_no_nodes_buffer(), std::vector<vec3>());
    std::vector<std::vector<vec3>> area_force(_dsc->get_no_nodes_buffer(), std::vector<vec3>());
    
    for(auto nit = _dsc->nodes_begin(); nit != _dsc->nodes_end(); nit++)
    {
        if (nit->is_interface() && !nit->is_boundary())
        {
            auto const & nodes_in_hats = _node_in_hat[nit.key()];
            auto const & curvature_in_hats = _mean_curvature_of_each_hat[nit.key()];
            auto const & tets_in_hats = _tets_in_hat[nit.key()];
            
            assert(nodes_in_hats.size() == curvature_in_hats.size());
            assert(nodes_in_hats.size() == tets_in_hats.size());
            
            for(int i = 0; i != nodes_in_hats.size(); i++)
            {
                    vec3 node_curvature = curvature_in_hats[i];
                    auto nodes_in_current_hat = nodes_in_hats[i];
                    auto tets_in_hat = tets_in_hats[i];
                

                
                    for (auto n : nodes_in_current_hat)
                    {
                        auto const & curvatire_in_hat_of_n = _mean_curvature_of_each_hat[n];

                        vec3 curvature_n;
                        if (curvatire_in_hat_of_n.size() == 1 // bondary vertex
                            || curvatire_in_hat_of_n.size() == 2)// This vertex separates two phases. Normal case.
                        {
                            curvature_n = curvatire_in_hat_of_n[0];
                        }
                        else{
                            // Find the coressponding hat
                            auto const & tets_in_current_hats = _tets_in_hat[n];
                            auto mean_curvature_on_hats = _mean_curvature_of_each_hat[n];
                            assert(tets_in_current_hats.size() == mean_curvature_on_hats.size());
                            
                            int idx = -1;
                            for (int itet = 0; itet < tets_in_current_hats.size(); itet++)
                            {
                                auto tets_in_current_hat = tets_in_current_hats[itet];
                    
                                if((tets_in_hat & tets_in_current_hat).size() > 0) // OPTIMIZE it later
                                {
                                    idx = itet;
                                    break;
                                }
                            }
                            
                            if(idx != -1)
                            {
                                curvature_n = _mean_curvature_of_each_hat[n][idx];
                            }
                            else{
                                break;
                            }
                        }
                        
                    
                        if (node_curvature.length() > curvature_n.length())
                        {
                            // Pull this vertex to the current vertex
                            vec3 nit_pos = nit->get_pos();
                            vec3 n_pos = _dsc->get_pos(n);
                            vec3 other_node_0_pos = _dsc->get_pos(nodes_in_current_hat[(i+1)%nodes_in_current_hat.size()]);
                            vec3 other_node_1_pos = _dsc->get_pos(nodes_in_current_hat[(i-1 + nodes_in_current_hat.size())%nodes_in_current_hat.size()]);
                            
                            vec3 n_nit = nit_pos - n_pos;
                            real h0 = Util::cross(n_nit, other_node_0_pos - n_pos).length();
                            real h1 = Util::cross(n_nit, other_node_1_pos - n_pos).length();
                            
                            vec3 curvature_f = n_nit*(node_curvature.length() - curvature_n.length());
                            vec3 area_f = -n_nit*(h0 + h1);
                            
                            
                            curvature_force[n].push_back(curvature_f);
                            area_force[n].push_back(area_f);
                        }
                    }
                }
        }
    }
    
    _curvature_force = curvature_force;
    _area_force = area_force;
}

void segment_function::compute_surface_curvature()
{
    // All the 'hat' that top at the vertices
    std::vector<std::vector<is_mesh::SimplexSet<is_mesh::FaceKey>>> interface_faces_around_node(_dsc->get_no_nodes_buffer(), std::vector<is_mesh::SimplexSet<is_mesh::FaceKey>>());
    std::vector<std::vector<is_mesh::SimplexSet<is_mesh::TetrahedronKey>>> tet_in_hat(_dsc->get_no_nodes_buffer(), std::vector<is_mesh::SimplexSet<is_mesh::TetrahedronKey>>());
    
    _mean_curvature_label.resize(_dsc->get_no_nodes_buffer());
    
    // Find all connected region around nodes: Find all the hats
    for(auto nit = _dsc->nodes_begin(); nit != _dsc->nodes_end(); nit++)
    {
        if(!(nit->is_interface() && !nit->is_boundary()))
            continue;
        
        
        auto tets = _dsc->get_tets(nit.key());
        
        while (tets.size() != 0)
        {
            std::queue<is_mesh::FaceKey> growing_queue;
            is_mesh::SimplexSet<is_mesh::TetrahedronKey> region;
            region.push_back(tets.back());
            tets -= tets.back();
            
            int label = _dsc->get_label(region.back());
            
            {
                auto fs = _dsc->get_faces(region.back());
                for(auto f : fs)
                    growing_queue.push(f);
            }
            
            while (!growing_queue.empty())
            {
                auto f = growing_queue.front();
                growing_queue.pop();
                
                auto tets_f = _dsc->get_tets(f);
                //find
                for (auto t : tets_f)
                {
                    if (_dsc->get_label(t) == label)
                    {
                        int idx = tets.index(t);
                        if (idx != -1)
                        {
                            region.push_back(t);
                            tets -= t;
                            
                            // add to region
                            auto fs = _dsc->get_faces(t);
                            for(auto f : fs)
                                growing_queue.push(f);
                            
                        }
                    }
                }
            }
            
            if (label != BOUND_LABEL)
            {
                tet_in_hat[nit.key()].push_back(region);
                _mean_curvature_label[nit.key()].push_back(label);
                
                is_mesh::SimplexSet<is_mesh::FaceKey> faceset = _dsc->get_faces(region) & _dsc->get_faces(nit.key()); // optimize later
                is_mesh::SimplexSet<is_mesh::FaceKey> interface_faceset;
                for (auto f : faceset)
                {
                    if (_dsc->get(f).is_interface())
                    {
                        interface_faceset.push_back(f);
                    }
                }
                
                interface_faces_around_node[nit.key()].push_back(interface_faceset);
            }
        }
    }
    _tets_in_hat = tet_in_hat;
    
    //-----------------
    // Now compute the mean curvature
    //    std::vector<double> node_cur(_dsc->get_no_nodes_buffer(), 0);
    std::vector<std::vector<vec3>> mean_curvature_of_each_hat(_dsc->get_no_nodes_buffer(), std::vector<vec3>());
    
    auto dsc = _dsc;
    std::vector<std::vector<is_mesh::SimplexSet<is_mesh::NodeKey>>> node_in_hat(_dsc->get_no_nodes_buffer(), std::vector<is_mesh::SimplexSet<is_mesh::NodeKey>>());
    for(auto n = dsc->nodes_begin(); n!= dsc->nodes_end(); n++)
    {
        if(n->is_interface() && !n->is_boundary())
        {
            // 1. Build nodes around
            
            auto face_around = interface_faces_around_node[n.key()];
            int removed = 0;
            
            for (int i = 0; i < face_around.size(); i++)
            {
                auto fs = face_around[i];
                
                auto node_around = dsc->extract_node_around(n.key(), fs);
                
                if (!node_around) // There are many case it can be null.
                {
                    // should also update the list of relations around?
                    interface_faces_around_node[n.key()].erase(interface_faces_around_node[n.key()].begin() + i - (removed));
                    _tets_in_hat[n.key()].erase(_tets_in_hat[n.key()].begin() + i - removed);
                    removed++;
                    continue;
                }
                
                node_in_hat[n.key()].push_back(*node_around);
                
                // 2. Compute the curvature
                auto p = dsc->get_pos(n.key());
                // 2a. Mixed area
                //            vec3 normal(0.0);
                double area_mixed = 0;
                for (int i = 0; i < node_around->size(); i++)
                {
                    vec3 v0 = p;
                    vec3 v1 = dsc->get_pos((*node_around)[i]);
                    vec3 v2 = dsc->get_pos( (*node_around)[(i+1)% (node_around->size()) ] );
                    
                    double f_area = Util::area<real>(v0, v1, v2);
                    
                    double a0 = acos(dot(v1-v0, v2-v0)/(length(v1-v0)*length(v2-v0)));
                    double a1 = acos(dot(v2-v1, v0-v1)/(length(v2-v1)*length(v0-v1)));
                    double a2 = acos(dot(v0-v2, v1-v2)/(length(v0-v2)*length(v1-v2)));
                    
                    if(a0>(M_PI/2.0) && a1>(M_PI/2.0) && a2>(M_PI/2.0)) // f is non-obtuse
                    {
                        // Add Voronoi formula (see Section 3.3)
                        area_mixed += (1.0/8) *
                        ((1.0/tan(a1)) * sqr_length(v2-v0) +
                         (1.0/tan(a2)) * sqr_length(v1-v0));
                    }
                    else // Voronoi inappropriate
                    {
                        // Add either area(f)/4 or area(f)/2
                        area_mixed += f_area/3;
                    }
                }
                // 2b. unnormalize curvature normal
                auto curv_normal = vec3(0);
                for (int i = 0; i < node_around->size(); i++)
                {
                    auto nbr = dsc->get_pos((*node_around)[i]);
                    auto right = dsc->get_pos((*node_around)[(i+1)%node_around->size()]);
                    auto left = dsc->get_pos((*node_around)[( i-1 + node_around->size() ) %node_around->size()]);
                    
                    double d_left = Util::dot(Util::normalize(nbr-left), Util::normalize(p-left));
                    double d_right = Util::dot(Util::normalize(nbr-right), Util::normalize(p-right));
                    double cos_left = std::min(1.0, std::max(-1.0, d_left));
                    double cos_right = std::min(1.0, std::max(-1.0, d_right));
                    double sin_left = sqrt(1 - cos_left*cos_left);
                    double sin_right = sqrt(1 - cos_right*cos_right);
                    
                    double w = (sin_left*cos_right + sin_right*cos_left)/(1e-300 + sin_left*sin_right);
                    
                    curv_normal += w * (nbr-p);
                }
                
                // 2c. The curvature
                auto mean_curvature_norm = curv_normal / (4*area_mixed);
                mean_curvature_of_each_hat[n.key()].push_back(mean_curvature_norm);
                
                
                delete node_around;
                
            }
            
        }
    }
    _node_in_hat = node_in_hat;
    _mean_curvature_of_each_hat = mean_curvature_of_each_hat;

}

#ifndef INTENSITY_IMAGE
void segment_function::compute_external_prob_force()
{
    // Buffer to keep the forces
    std::vector<vec3> forces = std::vector<vec3>(_dsc->get_no_nodes_buffer(), vec3(0.0));
    
    // Loop on interface faces
    for(auto fid = _dsc->faces_begin(); fid != _dsc->faces_end(); fid++)
    {
        if (fid->is_interface() && !is_boundary(fid.key()))
        {
            auto tets = _dsc->get_tets(fid.key());
            
            auto verts = _dsc->get_nodes(fid.key());
            auto pts = _dsc->get_pos(verts);
            
            auto l0 = _dsc->get_label(tets[0]);
            auto l1 = _dsc->get_label(tets[1]);
            
            // get normal, from label l0 to label l1
            vec3 Norm = _dsc->get_normal(fid.key());
            auto l01 = _dsc->barycenter(tets[1]) - _dsc->barycenter(tets[0]);
            Norm = Norm*dot(Norm, l01);// modify normal direction
            Norm.normalize();
            
            // Discretize the face
            double area = Util::area<double>(pts[0], pts[1], pts[2]);
            
            for (auto const & c : gaussian_points::get_point(area))
            {
                auto p = get_coord_tri(pts, c.coord);
                auto prob0 = m_prob_img.m_prob_map[l0]->get_value_f(p);
                auto prob1 = m_prob_img.m_prob_map[l1]->get_value_f(p);
                
                auto f = Norm* ((2 - prob0 - prob1)*(prob0-prob1) * c.weight * area); // Normalized already
                
                // distribute
                for(int ii =0; ii <3; ii++)
                    forces[verts[ii]] += f*c.coord[ii];
            }
        }
    }
    
    _forces = forces;
}
#endif

// Discrete Differential-Geometry Operators for Triangulated 2-Manifolds
void segment_function::compute_internal_force()
{
    std::vector<vec3> mean_curvature_norm_a(_dsc->get_no_nodes_buffer(), vec3(0));
    
    for(auto n = _dsc->nodes_begin(); n!= _dsc->nodes_end(); n++)
    {
        if(n->is_interface() && !n->is_boundary())
        {
            auto all_cur = _mean_curvature_of_each_hat[n.key()];
            if(all_cur.size() > 0)
            {
                for(auto cc : all_cur)
                    mean_curvature_norm_a[n.key()] += cc;
                
                mean_curvature_norm_a[n.key()] /= all_cur.size();
            }
        }
    }
    
    _internal_forces = mean_curvature_norm_a;
    
    // should we normalize it?
    double max_f = -INFINITY;
    for (auto f : _internal_forces)
    {
        if (max_f < f.length())
        {
            max_f = f.length();
        }
    }
    
    double scale = 5./max_f;
    for (auto & f : _internal_forces)
    {
        f *= scale;
    }
    cout << "Normalize curvature; scale = " << scale << endl;
}

void segment_function::compute_external_force()
{
    auto c = _mean_intensities;

    // Array to store temporary forces
    // Use fixxed array for better performance. Suppose that we have less than 10000 vertices
    std::vector<vec3> forces = std::vector<vec3>(_dsc->get_no_nodes_buffer(), vec3(0.0));

    // Loop on interface faces
    for(auto fid = _dsc->faces_begin(); fid != _dsc->faces_end(); fid++)
    {
        if (fid->is_interface()
//            && !is_boundary(fid.key())
            )
        {
            auto tets = _dsc->get_tets(fid.key());
            
            auto verts = _dsc->get_nodes(fid.key());
            auto pts = _dsc->get_pos(verts);
            
            auto phase0 = _dsc->get_label(tets[0])-1;
            auto phase1 = _dsc->get_label(tets[1])-1;
            
            double c0 = phase_intensity(phase0), c1 = phase_intensity(phase1);
            
            // get normal, from label l0 to label l1
            vec3 Norm = _dsc->get_normal(fid.key());
            auto l01 = _dsc->barycenter(tets[1]) - _dsc->barycenter(tets[0]);
            Norm = Norm*dot(Norm, l01);// modify normal direction
            Norm.normalize();
            
            // Discretize the face
            double area = Util::area<double>(pts[0], pts[1], pts[2]);
            
            for (auto const & c : gaussian_points::get_point(area))
            {
                auto p = get_coord_tri(pts, c.coord);
                auto I = _img.get_value_f(p);
                
                auto f = Norm*( (2*I - c0 - c1)*(c0 - c1) * c.weight * area);
                
                for(int ii =0; ii <3; ii++)
                    forces[verts[ii]] += f*c.coord[ii];
            }
        }
    }

    _forces = forces;
}

void segment_function::face_split()
{
    
    
//    // By analyzing the external force distribution
//    double boundary_intensity = -1;
//    auto c = _mean_intensities;
//
//    std::vector<double> avg_f = std::vector<double>(_dsc->get_no_faces_buffer(), 0.0);
//    std::vector<double> avg_f_abs = std::vector<double>(_dsc->get_no_faces_buffer(), 0.0);
//
//    is_mesh::SimplexSet<is_mesh::EdgeKey> to_split_edges;
//    // Loop on interface faces
//    for(auto fid = _dsc->faces_begin(); fid != _dsc->faces_end(); fid++)
//    {
//        if (fid->is_interface() && !fid->is_boundary())
//        {
//            // Get general information
//            auto tets = _dsc->get_tets(fid.key());
//            auto verts = _dsc->get_nodes(fid.key());
//            auto pts = _dsc->get_pos(verts);
//
//            double c0 = _dsc->get_label(tets[0]) == BOUND_LABEL? boundary_intensity : c[_dsc->get_label(tets[0])];
//            double c1 = _dsc->get_label(tets[1]) == BOUND_LABEL? boundary_intensity : c[_dsc->get_label(tets[1])];
//
////            // get normal
////            vec3 Norm = _dsc->get_normal(fid.key());
////            auto l01 = _dsc->barycenter(tets[1]) - _dsc->barycenter(tets[0]);
////            Norm = Norm*dot(Norm, l01);// modify normal direction
////            Norm.normalize();
//
//            // Discretize the face
//            double area = Util::area<double>(pts[0], pts[1], pts[2]);
//
//            size_t tri_sample_index = std::ceil( sqrt(area) ) - 1;
//            if (tri_sample_index >= tri_coord_size.size())
//            {
//                tri_sample_index = tri_coord_size.size() - 1;
//            }
//            if(tri_sample_index < 1)tri_sample_index = 1;
//
//            auto a = tri_dis_coord[tri_sample_index - 1];
//
//            // Analyze each sampling point
//            for (auto coord : a)
//            {
//                auto p = get_coord_tri(pts, coord);
//                auto g = _img.get_value_f(p);
//
//                auto f = - ((2*g - c0 - c1) * (c1-c0) ); // Normalized already
//
//                avg_f[fid.key()] += f;
//                avg_f_abs[fid.key()] += std::abs(f);
//            }
//
//            avg_f[fid.key()] /= a.size();
//            avg_f_abs[fid.key()] /= a.size();
//
//            // Now analyze the faces
//            // Somehow should be related to 2\pi\alpha
//            // IMPORTANT: Parameter
//            double aa = (c1 - c0)*(c1-c0);
//            if (std::abs(avg_f[fid.key()]) < 0.02*aa
//                && avg_f_abs[fid.key()] > 0.5*aa)
//            {
//                // shall be split
//                auto e = _dsc->longest_edge(_dsc->get_edges(fid.key()));
//                to_split_edges += e;
//            }
//        }
//    }
//
//    // Now recursively subdivide these edge
//    while(to_split_edges.size() > 0)
//        recursive_divide_edges(to_split_edges[0], to_split_edges);
}
void segment_function::recursive_divide_edges(is_mesh::EdgeKey cur_edge, is_mesh::SimplexSet<is_mesh::EdgeKey> & edges)
{
    auto neighbor_tets = _dsc->get_tets(cur_edge);
    
    // Make the edge safe to be split
    while (neighbor_tets.size() > 0)
    {
        auto t0 = neighbor_tets[0];
        auto longest_e_0 = _dsc->longest_edge(_dsc->get_edges(t0));
        if (cur_edge == longest_e_0)
        {
            neighbor_tets -= t0;
        }
        else
        {
            recursive_divide_edges(longest_e_0, edges);
            neighbor_tets =  _dsc->get_tets(cur_edge);
        }
    }
    
    // Split the edge
    auto nkey = _dsc->split(cur_edge);
    edges -= cur_edge;
}

/**
 * Bounding box of point list
 */
void bounding_box(const std::vector<vec3> & pts, vec3 & ld, vec3 & ru)
{
    ld = vec3(INFINITY);
    ru = vec3(-INFINITY);
    for(auto & v : pts)
    {
        ld[0] = std::min(v[0], ld[0]);
        ld[1] = std::min(v[1], ld[1]);
        ld[2] = std::min(v[2], ld[2]);
        
        ru[0] = std::max(v[0], ru[0]);
        ru[1] = std::max(v[1], ru[1]);
        ru[2] = std::max(v[2], ru[2]);
    }
}

bool sort_intersect(intersect_pt p1, intersect_pt p2)
{
    return p1.z < p2.z;
}

void segment_function::update_average_intensity()
{
    // Using sum table. Much faster than normal loop
#ifdef LOG_DEBUG
    cout << "Computing average intensity with" << NB_PHASE << " phases " << endl;
#endif

    int nb_phase = NB_PHASE;

    // 1. Init the buffer for intersection
    auto dim = _img.dimension();

    std::vector<ray_z> init_rayz(dim[0] * dim[1]);
    for (int y = 0; y < dim[1]; y++)
    {
        for (int x = 0; x < dim[0]; x ++)
        {
            int idx = y*dim[0] + x;
            init_rayz[idx].x = x;
            init_rayz[idx].y = y;
        }
    }

    vector<std::vector<ray_z>> ray_intersect(nb_phase, init_rayz);

    // 2. Find intersection with interface
    for(auto fid = _dsc->faces_begin(); fid != _dsc->faces_end(); fid++)
    {
        if (fid->is_interface() and !fid->is_boundary())
        {
            auto tet = _dsc->get_tets(fid.key());
            auto phase0 = _dsc->get_label(tet[0])-1;
            auto phase1 = _dsc->get_label(tet[1])-1;

            // check all z-ray that intersect this triangle
            auto pts3 = _dsc->get_pos(_dsc->get_nodes(fid.key()));
            auto pts = pts3;
            pts[0][2] = 0; pts[1][2] = 0; pts[2][2] = 0;

            auto n = _dsc->get_normal(fid.key(), tet[0]);
            bool in_1 = Util::dot(n, vec3(0,0,1)) > 0;

            vec3 ld, ru;
            bounding_box(pts, ld, ru);
            for (int x = std::floor(ld[0]); x < std::round(ru[0]); x++)
            {
                for (int y = std::floor(ld[1]); y < std::round(ru[1]); y++)
                {
                    if(y < 0 || x < 0
                       || y >= dim[1] || x >= dim[0])
                        continue;

                    try
                    {
                        bool bError;
                        auto bc = Util::barycentric_coords<double>(vec3(x+0.5, y+0.5, 0), pts[0], pts[1], pts[2], &bError);

                        if (bError)
                        {
                            continue;
                        }

                        if (bc[0] > -EPSILON && bc[1] > -EPSILON && bc[2] > -EPSILON)
                        { // inside
                            auto p = pts3[0]*bc[0] + pts3[1]*bc[1] + pts3[2]*bc[2];

                            auto zz = std::nearbyint(p[2]);

                            if(phase0 != BOUND_LABEL-1)
                            {
                                ray_intersect[phase0][y*dim[0] + x].intersects.push_back(intersect_pt(zz, !in_1));

                            }
                            if(phase1 != BOUND_LABEL-1)
                            {
                                ray_intersect[phase1][y*dim[0] + x].intersects.push_back(intersect_pt(zz, in_1));

                            }
                        }
                    }
                    catch (std::exception e)
                    {

                    }
                }
            }
        }
    }

    // 3. Compute integral
    _d_rayz.clear();

    _mean_intensities.resize(nb_phase);
    std::fill(_mean_intensities.begin(), _mean_intensities.end(), 0.0);
    vector<double> area(nb_phase, 0.0);
    for (int i = 0; i < nb_phase; i++)
    {
        int count = 0;
        for (auto r : ray_intersect[i])
        {
            std::vector<intersect_pt> intersect_ps = r.intersects;
            if (intersect_ps.size() > 1)
            {
                std::sort(intersect_ps.begin(), intersect_ps.end(), sort_intersect);

                // remove identical intersections
                for (auto p = intersect_ps.begin()+1; p != intersect_ps.end(); p++)
                {
                    auto pre = p -1;
                    if (p->z == pre->z and
                        p->b_in == pre->b_in)
                    {
                        p = intersect_ps.erase(p);

                        if(p == intersect_ps.end())
                        {
                            break;
                        }
                    }
                }
                // Now count
                if (intersect_ps.size() % 2 == 0) // Should be even
                {
                    count++;
                    auto newR = r;
                    newR.intersects = intersect_ps;
                    _d_rayz.push_back(newR);

                    for (int j = 0; j < intersect_ps.size()/2; j++)
                    {
                        int z1 = intersect_ps[2*j].z;
                        int z2 = intersect_ps[2*j + 1].z;

                        if(z1 < 0) z1 = 0;
                        if(z2 < z1)
                            z2 = z1;

                        area[i] += z2 - z1;
                        _mean_intensities[i] += _img.sum_line_z(r.x, r.y, z1, z2);
                    }
                }
                else{}//???
            }
        }
    }

    _total_intensities = _mean_intensities;
    _phase_volume = area;

    for (int i  = 0; i < nb_phase; i++)
    {
        if(area[i]==0)
        {
            std::cout << "Zero; " << _mean_intensities[i] << std::endl;
        }
        _mean_intensities[i] /= area[i];
    }

    cout << "Mean intensity: ";
    for (int i = 0; i < nb_phase; i++)
    {
        cout << _mean_intensities[i] << " ; ";
    }
    cout << endl;
}

#pragma mark MAIN FUNCTION

void segment_function::estimate_time_step()
{
    update_vertex_boundary();
    
    update_average_intensity();
    
    compute_internal_force_2(); //
    compute_external_force();
    
    double max_dis = 0;
    for (auto nit = _dsc->nodes_begin(); nit != _dsc->nodes_end(); nit++)
    {
        if ( _dsc->exists(nit.key())
            && (nit->is_interface() or nit->is_crossing()))
        {
            max_dis = max(max_dis, get_node_displacement(nit.key()).length());
        }
    }
    
    // At first the force is small because of the dihedral angle
    //  when the mesh converges, sum of external forces become larger
    //  hence we estimate the time step with smaller max-displacement
    double max_displace = _dsc->get_avg_edge_length()*0.05;
    _dt = max_displace / max_dis;
    
    cout << "Estimated time step: " << _dt << endl;
}

void segment_function::adapt_surface()
{
    
    
#ifdef DSC_CACHE
    long nb_collasped = 0;
    long nb_collasped_force = 0;
    for (auto nit = _dsc->nodes_begin(); nit != _dsc->nodes_end(); nit++)
    {
        if (!(_dsc->exists(nit.key()) && !nit->is_crossing()
              && (nit->is_interface() || nit->is_boundary())))
        {
            continue;
        }
        
        // Precompute face normal
        is_mesh::SimplexSet<is_mesh::FaceKey> neighbor_faces;
        std::vector<vec3> norm_faces;
        map<is_mesh::EdgeKey,vector<vec3>> edge_cobound_norm;
        for(auto f : *_dsc->get_faces_cache(nit.key()))
        {
            if (_dsc->get(f).is_interface() || _dsc->get(f).is_boundary())
            {
                neighbor_faces += f;
                norm_faces.push_back(_dsc->get_normal(f));
                
                for (auto e : _dsc->get_edges(f))
                {
                    edge_cobound_norm[e].push_back(norm_faces.front());
                }
            }
        }
        
        // Surface flateness
        static double cos_flat = cos(88.*M_PI/180.); // Threshold of flat surface
        double cos_max_angle = 1;
        for (int i = 0; i < norm_faces.size(); i++)
        {
            for (int j = 0; j < norm_faces.size(); j++)
            {
                double cc = Util::dot(norm_faces[i], norm_faces[j]);
                cos_max_angle = min(cos_max_angle, cc);
            }
        }
        bool is_flat = cos_max_angle < cos_flat;
        
        is_mesh::SimplexSet<is_mesh::EdgeKey> ring_edge, cobound_edge;
        for(auto e : _dsc->get_edges(neighbor_faces))
        {
            auto nodes = _dsc->get_nodes(e);
            if (nodes[0] == nit.key() || nodes[1] == nit.key())
            {
                cobound_edge += e;
            }else{
                ring_edge += e;
            }
        }
        
        is_mesh::SimplexSet<is_mesh::EdgeKey> colapse_candidate_edge;
        if (is_flat)
        {
            colapse_candidate_edge = cobound_edge;
        }
        else if(nit->is_boundary() || m_vertex_bound[nit.key()])
        {
            for (auto e : cobound_edge)
            {
                auto e_norm = edge_cobound_norm[e];
                assert(e_norm.size() == 2);
                if (Util::dot(e_norm[0], e_norm[1]) < 0.1)
                {
                    // edge
                    colapse_candidate_edge.push_back(e);
                }
            }

            if (colapse_candidate_edge.size() != 2)
            {
                continue;
            }
        }
        else
        {
            continue;
        }

        assert(colapse_candidate_edge.size() > 1);
        
        // find shortest edges
        auto shortest_edge =  _dsc->shortest_edge(colapse_candidate_edge);
        auto colapse_node =_dsc->get_nodes(shortest_edge) - nit.key();
        assert(colapse_node.size()==1); // Hold only the opposite node of nit
        auto new_pos = _dsc->get_pos(colapse_node[0]);
        auto nid0 = colapse_node[0];

        // Check quality
        double min_tet_quality = INFINITY;
        auto face_link = _dsc->get_link(nit.key());
        for(auto fid : *face_link)
        {
            auto f_pts = *_dsc->get_nodes_cache(fid);

            // The nid0 belong to the triangle
            if(nid0 == f_pts[0] || nid0 == f_pts[1]  || nid0 == f_pts[2])
                continue;
            
            auto f_pos = _dsc->get_pos(f_pts);
            min_tet_quality = std::min(min_tet_quality,
                                       std::abs(Util::quality<real>(f_pos[0], f_pos[1], f_pos[2], new_pos)));
        }
        
        
        if (min_tet_quality < _dsc->pars.MIN_TET_QUALITY)
        {
            _dsc->is_edge_adapted(shortest_edge, true);
        }
        else
        {
            // 3. If improve mesh quality, check mumford-shah energy
            //   Should also consider the MS energy?
            if(!(m_vertex_bound[nit.key()] || nit->is_boundary()))
            {
                if(_dsc->collapse(shortest_edge, true))
                    nb_collasped++;
                else{
                    _dsc->collapse(shortest_edge, false);
                    nb_collasped_force ++;
                }
            }
            else{
//            _dsc->collapse_cache(shortest_edge, nid0, 0.);


            if(_dsc->collapse(shortest_edge, true))
                nb_collasped++;
            else{
                if(_dsc->collapse(shortest_edge, false))
                    nb_collasped_force ++;
            }
            }
        }
    }
    
    
    
    
    cout << nb_collasped << " collapsed; " << nb_collasped_force << " forced" << endl;
    _dsc->garbage_collect();
#else
    assert(0); // Implement non cache code
#endif
}

struct face_analyze{
    
};

void segment_function::thickenning_surface_edge_length()
{
    static double thres_long_edge = _dsc->get_avg_edge_length();
    std::vector<is_mesh::EdgeKey> edges;
    for (auto eit = _dsc->edges_begin(); eit != _dsc->edges_end(); eit++)
    {
        if (eit->is_interface() && !is_boundary(eit.key()) && _dsc->length(eit.key()) > thres_long_edge)
        {
            edges.push_back(eit.key());
        }
    }
    int i = 0;
    for(auto &e : edges)
    {
        if (_dsc->exists(e) && _dsc->length(e) > thres_long_edge)
        {
            _dsc->split(e);
            i++;
        }
    }
    _dsc->garbage_collect();
}
void segment_function::thickenning_surface()
{
    auto c = _mean_intensities;
    update_vertex_boundary();
    
    // Array to store temporary forces
    // Use fixxed array for better performance. Suppose that we have less than 10000 vertices
    std::vector<vec3> forces = std::vector<vec3>(_dsc->get_no_nodes_buffer(), vec3(0.0));
    vector<bool> edge_to_split(_dsc->get_no_edges_buffer(), false);
    
    static double thres_force_split = _dsc->get_avg_edge_length()*0.5*0.01;
    static double stable = _dsc->get_avg_edge_length()*0.5*0.1;

    
    // Loop on interface faces
    for(auto fid = _dsc->faces_begin(); fid != _dsc->faces_end(); fid++)
    {
        if (fid->is_interface()
            && !is_boundary(fid.key())
            )
        {
            auto tets = _dsc->get_tets(fid.key());
            
            auto verts = _dsc->get_nodes(fid.key());
            auto pts = _dsc->get_pos(verts);
            
            auto phase0 = _dsc->get_label(tets[0])-1;
            auto phase1 = _dsc->get_label(tets[1])-1;
            
            double c0 = phase_intensity(phase0), c1 = phase_intensity(phase1);
            
            // get normal, from label l0 to label l1
            vec3 Norm = _dsc->get_normal(fid.key());
            auto l01 = _dsc->barycenter(tets[1]) - _dsc->barycenter(tets[0]);
            Norm = Norm*dot(Norm, l01);// modify normal direction
            Norm.normalize();
            
            // Discretize the face
            double area = Util::area<double>(pts[0], pts[1], pts[2]);
            
            vector<double> face_forces;
            for (auto const & c : gaussian_points::get_point(area))
            {
                auto p = get_coord_tri(pts, c.coord);
                auto I = _img.get_value_f(p);
                
                auto magnitude =  (2*I - c0 - c1)*(c0 - c1) * c.weight * area;
                auto f = Norm*magnitude;
                
                for(int ii =0; ii <3; ii++)
                    forces[verts[ii]] += f*c.coord[ii];
                
                face_forces.push_back(magnitude);
            }
            
            // Analyze the surface
            double mean = 0;
            for (auto ff : face_forces )
            {
                mean += ff;
            }
            mean /= face_forces.size();
            double deviation = 0; // average deviation
            for (auto ff : face_forces)
            {
                deviation += std::abs(ff-mean);
            }
            deviation /= face_forces.size();
            
            // compare it with force?
            if (mean > thres_force_split || deviation > thres_force_split)
            {
                auto longest_e = _dsc->longest_edge(_dsc->get_edges(fid.key()));
                edge_to_split[longest_e] = true;
            }
        }
    }
    
    _forces = forces;
    compute_internal_force_2();
    
    vector<bool> mark_to_split(_dsc->get_no_edges_buffer(), false);
    for(auto eit = _dsc->edges_begin(); eit != _dsc->edges_end(); eit++)
    {
        if (edge_to_split[eit.key()])
        {
            auto nodes = _dsc->get_nodes(eit.key());
            if (get_node_displacement(nodes[0]).length() < stable
                && get_node_displacement(nodes[1]).length() < stable)
            {
                mark_to_split[eit.key()] = true;
            }
        }
    }
    
    int count=0;
    for(int i = 0; i < mark_to_split.size(); i++)
    {
        if (mark_to_split[i]
            && _dsc->length(is_mesh::EdgeKey(i)) > _dsc->get_avg_edge_length())
        {
            is_mesh::EdgeKey ek(i);
            if (_dsc->exists(ek))
            {
                _dsc->split(ek);
                count++;
            }
        }
    }
    
    _dsc->garbage_collect();
    _dsc->deform();
}

#define align_bound_gap(idx) \
if(pos[idx] < 0) pos[idx] = origin[idx]; \
if(pos[idx] > im_bound[idx]) pos[idx] = bound[idx];
void segment_function::pad_boundary(double scale)
{
    double gap = _dsc->get_avg_edge_length()*scale;
    vec3 origin = vec3(0.) - vec3(gap);
#ifdef INTENSITY_IMAGE
    vec3 im_bound = (_img.dimension_v());
#else
    vec3 im_bound = vec3(m_prob_img.m_dimension);
#endif
    vec3 bound = im_bound + vec3(gap);
    
    std::vector<std::bitset<4>> direction_state(_dsc->get_no_nodes_buffer(),std::bitset<4>("0000"));
    for (auto fit = _dsc->faces_begin(); fit != _dsc->faces_end(); fit++)
    {
        if (fit->is_boundary())
        {
            auto norm = _dsc->get_normal(fit.key());
            
            auto direct = get_direction(norm);
            for(auto n : _dsc->get_nodes(fit.key()))
            {
                direction_state[n] = direction_state[n] | direct;
            }
        }
    }
 
    for (auto nit = _dsc->nodes_begin(); nit != _dsc->nodes_end(); nit++)
    {
        if (nit->is_boundary())
        {
            auto direct = direction_state[nit.key()];
            auto pos = nit->get_pos();
            if ((direct & X_direction).to_ulong() != 0) // constraint on x
            {
                align_bound_gap(0);
            }
            if ((direct & Y_direction).to_ulong() != 0){ // constraint on y
                align_bound_gap(1);
            }
            if ((direct & Z_direction).to_ulong() != 0){ // constraint on z
                align_bound_gap(2);
            }
            
            _dsc->set_pos(nit.key(), pos);
        }
    }
}

void segment_function::segment_probability()
{

}

void segment_function::update_vertex_boundary()
{
    // Make sure there are gap
    for (auto nit = _dsc->nodes_begin(); nit != _dsc->nodes_end(); nit++)
    {
        if(nit->is_boundary())
        {
            for(auto t : _dsc->get_tets(nit.key()))
                _dsc->set_label(t, BOUND_LABEL);
        }
    }
    
    // Mark boundary vertex
    m_vertex_bound = vector<bool>(_dsc->get_no_nodes_buffer(), false);
    for (auto fit = _dsc->faces_begin(); fit != _dsc->faces_end(); fit++)
    {
        if(fit->is_interface())
        {
            auto tets = _dsc->get_tets(fit.key());
            if (_dsc->get_label(tets[0]) == BOUND_LABEL
                || _dsc->get_label(tets[1]) == BOUND_LABEL)
            {
                for(auto n : _dsc->get_nodes(fit.key()))
                    m_vertex_bound[n] = true;
            }
        }
    }
}

bool segment_function::is_boundary(is_mesh::EdgeKey ek)
{
    for(auto & n : _dsc->get_nodes(ek))
    {
        if (!m_vertex_bound[n])
        {
            return false;
        }
    }
    return true;
}

bool segment_function::is_boundary(is_mesh::FaceKey fk)
{
    for(auto & n : _dsc->get_nodes(fk))
    {
        if (!m_vertex_bound[n])
        {
            return false;
        }
    }
    return true;
}
void segment_function::compute_internal_force_2()
{
    std::vector<vec3> intern_f_a(_dsc->get_no_nodes_buffer(), vec3(0));
    for (auto fit = _dsc->faces_begin(); fit != _dsc->faces_end(); fit++)
    {
        if (fit->is_interface()
            && !is_boundary(fit.key())
            )
        {
            auto nodes = _dsc->get_nodes(fit.key());
            auto node_pos = _dsc->get_pos(nodes);
            for (int i = 0; i < 3; i++)
            {
                auto l0 = node_pos[i];
                auto l1 = node_pos[(i+1)%3];
                auto l2 = node_pos[(i+2)%3];
                
                auto p = Util::project_point_line(l0, l1, l2-l1);
                auto h = Util::normalize(l0-p);
                
                intern_f_a[nodes[i]] += -h*(l2-l1).length();
            }
        }
    }
    
    for (auto nit = _dsc->nodes_begin(); nit != _dsc->nodes_end(); nit++)
    {
        if (nit->is_interface() && !nit->is_crossing() && intern_f_a[nit.key()].length() > EPSILON)
        {
            auto norm = _dsc->get_normal(nit.key());
            intern_f_a[nit.key()] = norm * Util::dot(norm, intern_f_a[nit.key()]);
        }
    }
    
    _internal_forces = intern_f_a;
}


void segment_function::segment()
{
    static int iter = 0;
    cout << "Iter " << iter << " ============================" << endl;
    
    profile t("Compute force");
    
    update_vertex_boundary();
    
    update_average_intensity();
    
    compute_internal_force_2(); //
    compute_external_force();
    
    t.change("Work around boundary");
    snapp_boundary();
    
    t.change("DSC deform");
    _dsc->deform();
    
    if ( (iter % 20) == 0)
    {
        t.change("Relabel");
        update_vertex_boundary();
        update_average_intensity();
        relabel_tetrahedra();
        
        t.change("Adapt");
        update_vertex_boundary();
        _dsc->adapt();
        
//        update_vertex_boundary();
//        thickenning_surface_edge_length();
        
        adapt_surface();
        
        // fix the boundary
//        fix_snapping_boundary();
    }
    
    t.change("Measure energy");
    update_average_intensity();
    compute_energy();
    
    t.done();
    
    cout << "------------------------------" <<endl;
    profile::close();
    cout << "------------------------------" <<endl;
    iter ++;
}

void segment_function::fix_snapping_boundary()
{
    update_vertex_boundary();
    auto dim = _img.dimension_v();
    double thres_out = _dsc->get_avg_edge_length();
    double thres_in = _dsc->get_avg_edge_length()*0.5*0.6;
    bool deform = false;
    for (auto nit = _dsc->nodes_begin(); nit != _dsc->nodes_end(); nit++)
    {
        if (m_vertex_bound[nit.key()])
        {
            auto pos = nit->get_pos();
            bool snapped = false;
            for (int i = 0; i < 3; i++)
            {
                if (abs(pos[i]) < thres_in
                    || abs(pos[i] - (dim[i]-1)) < thres_in)
                {
                    snapped = true;
                }
            }
            if (!snapped)
            {
                for (int i = 0; i < 3; i++)
                {
                    if (pos[i] < thres_out)
                        pos[i] = 0;
                    
                    if(pos[i] > dim[i]-1-thres_out)
                        pos[i] = dim[i] - 1;
                }
                
                _dsc->set_destination(nit.key(), pos);
                deform = true;
            }
        }
    }
    
    if (deform)
    {
        _dsc->deform();
    }
}

void segment_function::compute_energy()
{
    // Using sum table. Much faster than normal loop
#ifdef LOG_DEBUG
    cout << "Computing average intensity with" << NB_PHASE << " phases " << endl;
#endif
    
    int nb_phase = NB_PHASE;
    
    // 1. Init the buffer for intersection
    auto dim = _img.dimension();
    
    std::vector<ray_z> init_rayz(dim[0] * dim[1]);
    for (int y = 0; y < dim[1]; y++)
    {
        for (int x = 0; x < dim[0]; x ++)
        {
            int idx = y*dim[0] + x;
            init_rayz[idx].x = x;
            init_rayz[idx].y = y;
        }
    }
    
    vector<std::vector<ray_z>> ray_intersect(nb_phase, init_rayz);
    
    // 2. Find intersection with interface
    for(auto fid = _dsc->faces_begin(); fid != _dsc->faces_end(); fid++)
    {
        if (fid->is_interface() and !fid->is_boundary())
        {
            auto tet = _dsc->get_tets(fid.key());
            auto phase0 = _dsc->get_label(tet[0])-1;
            auto phase1 = _dsc->get_label(tet[1])-1;
            
            // check all z-ray that intersect this triangle
            auto pts3 = _dsc->get_pos(_dsc->get_nodes(fid.key()));
            auto pts = pts3;
            pts[0][2] = 0; pts[1][2] = 0; pts[2][2] = 0;
            
            auto n = _dsc->get_normal(fid.key(), tet[0]);
            bool in_1 = Util::dot(n, vec3(0,0,1)) > 0;
            
            vec3 ld, ru;
            bounding_box(pts, ld, ru);
            for (int x = std::floor(ld[0]); x < std::round(ru[0]); x++)
            {
                for (int y = std::floor(ld[1]); y < std::round(ru[1]); y++)
                {
                    if(y < 0 || x < 0
                       || y >= dim[1] || x >= dim[0])
                        continue;
                    
                    try
                    {
                        bool bError;
                        auto bc = Util::barycentric_coords<double>(vec3(x+0.5, y+0.5, 0), pts[0], pts[1], pts[2], &bError);
                        
                        if (bError)
                        {
                            continue;
                        }
                        
                        if (bc[0] > -EPSILON && bc[1] > -EPSILON && bc[2] > -EPSILON)
                        { // inside
                            auto p = pts3[0]*bc[0] + pts3[1]*bc[1] + pts3[2]*bc[2];
                            
                            auto zz = std::nearbyint(p[2]);
                            
                            if(phase0 != BOUND_LABEL-1)
                            {
                                ray_intersect[phase0][y*dim[0] + x].intersects.push_back(intersect_pt(zz, !in_1));
                                
                            }
                            if(phase1 != BOUND_LABEL-1)
                            {
                                ray_intersect[phase1][y*dim[0] + x].intersects.push_back(intersect_pt(zz, in_1));
                                
                            }
                        }
                    }
                    catch (std::exception e)
                    {
                        
                    }
                }
            }
        }
    }
    
    // 3. Compute integral
    _d_rayz.clear();
    
    double difference = 0;
    vector<double> area(nb_phase, 0.0);
    for (int i = 0; i < nb_phase; i++)
    {
        int count = 0;
        double cc = _mean_intensities[i];
        
        for (auto r : ray_intersect[i])
        {
            std::vector<intersect_pt> intersect_ps = r.intersects;
            if (intersect_ps.size() > 1)
            {
                std::sort(intersect_ps.begin(), intersect_ps.end(), sort_intersect);
                
                // remove identical intersections
                for (auto p = intersect_ps.begin()+1; p != intersect_ps.end(); p++)
                {
                    auto pre = p -1;
                    if (p->z == pre->z and
                        p->b_in == pre->b_in)
                    {
                        p = intersect_ps.erase(p);
                        
                        if(p == intersect_ps.end())
                        {
                            break;
                        }
                    }
                }
                // Now count
                if (intersect_ps.size() % 2 == 0) // Should be even
                {
                    count++;
                    auto newR = r;
                    newR.intersects = intersect_ps;
                    _d_rayz.push_back(newR);
                    
                    for (int j = 0; j < intersect_ps.size()/2; j++)
                    {
                        int z1 = intersect_ps[2*j].z;
                        int z2 = intersect_ps[2*j + 1].z;
                        
                        if(z1 < 0) z1 = 0;
                        if(z2 < z1)
                            z2 = z1;
                        
                        area[i] += z2 - z1;
                        for (int k = z1+1; k<z2; k++)
                        {
                            difference += std::pow(_img.get_value(r.x, r.y, k) - cc, 2);
                        }
                    }
                }
                else{}//???
            }
        }
    }
    
    update_vertex_boundary();
    double areaa = 0;
    for (auto fit = _dsc->faces_begin(); fit != _dsc->faces_end(); fit++)
    {
        if (fit->is_interface() && !is_boundary(fit.key()))
        {
            areaa += m_alpha * _dsc->area(fit.key());
        }
    }

    cout << "=== Energy; area = " << difference << " - " << areaa << endl;
}

void segment_function::export_surface_mesh()
{
//    dsc_export::export_surface(_dsc, "./LOG");
}

int segment_function::arg_min_phase_point(vec3 pt, double radius, int current_label) {

}

void segment_function::recursive_subdivide(is_mesh::TetrahedronKey tkey, vec3 pt, int new_label, double min_volume)
{
    
    // Find the longest edge
    auto eids = _dsc->get_edges(tkey);
    auto eid = _dsc->longest_edge(eids);
    
    // Split the longest edge
    auto vid = _dsc->split(eid);
    
    // The the new tetrahedron that contain pt
    //  Optimize later !!!
    auto tet_neighbor = _dsc->get_tets(vid);
    for (auto tkey_new : tet_neighbor)
    {
        tet_touched[tkey_new] = true;
    }
    
    for (auto tkey_new : tet_neighbor)
    {
        auto tet_pos = _dsc->get_pos(_dsc->get_nodes(tkey_new));
        auto bary_coord = Util::barycentric_coords<real>(pt, tet_pos[0], tet_pos[1], tet_pos[2], tet_pos[3]);
        
        if(bary_coord[0] >= 0 && bary_coord[1] >= 0 && bary_coord[2] >= 0 && bary_coord[3] >= 0)
        {
            // Inside
            if (_dsc->volume(tkey_new) > 1.5*min_volume)
            {
                cout << "Successed. relabel to " << new_label << " \n";
                recursive_subdivide(tkey_new, pt, new_label, min_volume);
            }
            else{
                cout << "Go smaller, volume = " << _dsc->volume(tkey_new) << endl;
                _dsc->set_label(tkey_new, new_label);
            }
            break;
        }
    }
}

void segment_function::recursive_divide(std::vector<point_to_capture>* subdivide_tets, is_mesh::TetrahedronKey tkey,  int depth, std::queue<is_mesh::TetrahedronKey> & debug_tet_queue)
{

    
    debug_tet_queue.push(tkey);
    
    auto longest_e = _dsc->longest_edge(_dsc->get_edges(tkey));
//    cout << " Recursive depth " << depth << " tet " << tkey << "; to split " << longest_e << endl;
    
//    if(depth < 4)
    {
        auto neighbor_tets = _dsc->get_tets(longest_e);
        
        while (neighbor_tets.size() > 0)
        {
            auto t0 = neighbor_tets[0];
            auto longest_e_0 = _dsc->longest_edge(_dsc->get_edges(t0));
            if (longest_e == longest_e_0)
            {
                neighbor_tets -= t0;
            }
            else
            {
                recursive_divide(subdivide_tets, t0, depth+1, debug_tet_queue);
                neighbor_tets =  _dsc->get_tets(longest_e);
            }
        }
    }
    
    // Now longest_e is compatible division
    devide_element(subdivide_tets, longest_e);
}

bool is_point_inside(const vec3 & p, const vec3 & a, const vec3 & b, const vec3 & c, const vec3 & d)
{
    auto bary_coord = Util::barycentric_coords<real>(p, a, b, c, d);
    return (bary_coord[0] >= 0 && bary_coord[1] >= 0 && bary_coord[2] >= 0 && bary_coord[3] >= 0);
}

void segment_function::devide_element(std::vector<point_to_capture>* subdivide_tets, is_mesh::EdgeKey ekey)
{
//    std::cout << "Adapt edge " << ekey << endl;
    // Find capturing point affected
    auto old_tets = _dsc->get_tets(ekey);
    std::vector<point_to_capture> current_subdivide_tets;
    
//    cout << "subdivide_tets size (before): " << subdivide_tets->size() << endl;
    for (auto t : old_tets)
    {
        point_to_capture temp;
        temp.tet_key = t;
        // Optimize later
        auto tp = std::find(subdivide_tets->begin(), subdivide_tets->end(), temp);
        if (tp != subdivide_tets->end())
        {
            temp = *tp;
            current_subdivide_tets.push_back(temp);
            
            subdivide_tets->erase(tp);
//            cout << "Remove " << tp->tet_key << endl;
//            subdivide_tets.erase(std::remove(subdivide_tets.begin(), subdivide_tets.end(), temp), subdivide_tets.end());
        }
    }
//    cout << "subdivide_tets size (after): " << subdivide_tets->size() << endl;
//    cout << "related tets #: " << current_subdivide_tets.size() << endl;
    // Split the edge
    auto nkey = _dsc->split(ekey);
    
    // Update the capturing point
    auto new_tets = _dsc->get_tets(nkey);
    std::vector<point_to_capture> remain_subdivide_tets;
    for (auto t : new_tets)
    {
        auto pts = _dsc->get_pos(_dsc->get_nodes(t));
        for(int i = 0; i < current_subdivide_tets.size(); i++)
        {
            auto old_p = current_subdivide_tets[i];
            
            if (is_point_inside(old_p.pt, pts[0], pts[1], pts[2], pts[3]))
            {
//                cout << "Volume " << t << ": " << _dsc->volume(t) << "min V " << min_V << endl;
                
                // check if we could relabel
                // TUAN: Consider if needed?
                auto old_energy = get_energy_tetrahedron(t, _dsc->get_label(t));
                auto new_energy = get_energy_tetrahedron(t, old_p.new_label);
                
                if(new_energy < old_energy || _dsc->volume(t) < min_V)
                {
                    _dsc->set_label(t, old_p.new_label);
                }
                else{
                    auto new_p = old_p;
                    new_p.tet_key = t;
                    remain_subdivide_tets.push_back(new_p);
                }
                
                current_subdivide_tets.erase(current_subdivide_tets.begin() + i);
                continue;
            }
        }
    }
    
    subdivide_tets->insert(subdivide_tets->end(), remain_subdivide_tets.begin(), remain_subdivide_tets.end());
    
//    cout << "subdivide_tets size (final): " << subdivide_tets->size() << endl;
}


