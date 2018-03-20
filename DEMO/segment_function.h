//
//  segment_function.hpp
//  DSC_segment
//
//  Created by Tuan Nguyen Trung on 12/9/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#ifndef segment_function_hpp
#define segment_function_hpp

#include <stdio.h>
#include "image3d.h"
#include "DSC.h"

#include "define.h"
#include <queue>

#define STABLE_DISPLACEMENT 0.1

extern std::bitset<4> X_direction, Y_direction, Z_direction;


struct intersect_pt
{
    intersect_pt(int z_, bool b_in_z_){
        z = z_;
        b_in = b_in_z_;
    }
    int z;
    bool b_in;
};

struct ray_z
{
    int x, y;
    std::vector<intersect_pt> intersects;
};

//template <typename T>
//class vector_attribute: public std::vector<T>
//{
//public:
//    using std::vector<T>::size;
//    using std::vector<T>::resize;
//    using std::vector<T>::at;
//
//    T & operator[] (int idx){
//        if(idx >= size())
//        {
//            resize(idx+1, _default);
//        }
//        return at(idx);
//    };
//
//    vector_attribute(T defalt){
//        _default = defalt;
//    }
//
//private:
//    T _default; // default value of element when the array shrink
//};

class point_to_capture{
public:
    point_to_capture(){}
    point_to_capture(is_mesh::TetrahedronKey t, vec3 p, int l){
        tet_key = t;
        pt = p;
        new_label = l;
    };
    ~point_to_capture(){}
    
    is_mesh::TetrahedronKey tet_key;
    vec3 pt;
    int new_label;
    
    bool operator== (const point_to_capture & b )
    {
        return tet_key == b.tet_key;
    }
};

class segment_function
{
    typedef DSC::DeformableSimplicialComplex<> dsc_class;
    
public:
    segment_function(){};
    ~segment_function(){};
    
    void init(); // Load image
    void initialze_segmentation(); // Label initialization
    void random_initialization(); // Label initialization
    void initialization_discrete_opt(); // Optimize the initialization automatically
    
    void segment();
    
public:// Configuration parametters
    int NB_PHASE;
    double VARIATION_THRES_FOR_RELABELING;
    double ALPHA = 0.1;
    double QALPHA = 0.2;
    double _dt = 1;
    std::string _directory_path;
    
    // Face plit
    double ratio_signed_and_mag_mean = 0.2;
    
public:// Variables
    image3d _img; // Store crossection -> voxel
    dsc_class *_dsc; // Shared dsc
    
public:
    std::vector<double> _mean_intensities;
    std::vector<double> _total_intensities; // To update mean intensity during relabeling
    std::vector<double> _phase_volume;
    void update_average_intensity();
    void update_average_intensity1();
    
    std::vector<int> _vertex_stability_map;
    std::vector<vec3> _forces;
    std::vector<vec3> _internal_forces;
    std::vector<vec3> _quality_control_forces; // curvature based force
    std::vector<vec3> _quality_angle_forces; // angle based
    
    // Adaptive time step
    std::vector<vec3> _previous_dis;
    std::vector<vec3> _cur_dis;
    std::vector<double> _dt_adapt;
    
    void remove_stable_proximity(std::vector<std::vector<double>> & barry_coord, const is_mesh::SimplexSet<is_mesh::NodeKey> & nodes);
    vec3 get_node_displacement(is_mesh::NodeKey nkey);
    
    void compute_surface_curvature();
    void compute_external_force();
    void compute_internal_force();
    void compute_internal_force_simple();
    void compute_mesh_quality_control_force();
    
    void work_around_on_boundary_vertices();
    void update_vertex_stability();
    
    
    double get_energy_tetrahedron(is_mesh::TetrahedronKey tkey, int );
    void relabel_tetrahedra();
    
    
    double min_edge, min_V;
    void set_min_edge_length(double l){min_edge = l; min_V = pow(min_edge, 3)/6;};
    void face_split();
    
    void adapt_interface();
    
    void adapt_tetrahedra();
    void adapt_tetrahedra_1();
    void recursive_divide(std::vector<point_to_capture>* subdivide_tets, is_mesh::TetrahedronKey tkey, int depth, std::queue<is_mesh::TetrahedronKey> & debug_tet_queue);
    void recursive_divide_edges(is_mesh::EdgeKey cur_edge, is_mesh::SimplexSet<is_mesh::EdgeKey> & edges);
    void devide_element(std::vector<point_to_capture>* subdivide_tets, is_mesh::EdgeKey ekey);
    int arg_min_phase_point(vec3 pt, double radius, int current_label);
    void recursive_subdivide(is_mesh::TetrahedronKey tkey, vec3 pt, int new_label, double min_volume);
public:
    std::vector<ray_z> _d_rayz;
    
    std::vector<std::vector<vec3>> _mean_curvature_of_each_hat;
    std::vector<std::vector<is_mesh::SimplexSet<is_mesh::TetrahedronKey>>> _tets_in_hat;
    std::vector<std::vector<is_mesh::SimplexSet<is_mesh::NodeKey>>> _node_in_hat;
    
    // For debuging
    std::vector<vec3> boundary_vertices_displacements;
    
    std::vector<unsigned int> d_is_image_boundary;
    std::vector<std::bitset<4>> d_direction_state;
    
    std::vector<std::vector<vec3>> _curvature_force;
    std::vector<std::vector<vec3>> _area_force;
};

#endif /* segment_function_hpp */
