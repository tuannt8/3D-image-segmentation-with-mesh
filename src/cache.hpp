//
//  cache.hpp
//  DSC_parallel
//
//  Created by Tuan Nguyen Trung on 8/1/16.
//  Copyright © 2016 Asger Nyman Christiansen. All rights reserved.
//

#ifndef cache_hpp
#define cache_hpp

#include <stdio.h>
#include <vector>
#include <bitset>
#include "util.h"
#include "is_mesh.h"

#define DSC_CACHE

//#define MAX_ELEMENTS 1000000

#define CLEAN_GARBAGE(a, b) if(a[b]){delete a[b]; a[b] = nullptr;}

// Should consider thread safe

class dsc_cache
{
public:
    // Link of node
    std::vector<is_mesh::SimplexSet<is_mesh::FaceKey>*> link_of_node;
    // Neighbor tets of node
    std::vector<is_mesh::SimplexSet<is_mesh::TetrahedronKey>*> tets_neighbor_node;
    std::vector<is_mesh::SimplexSet<is_mesh::FaceKey>*> faces_neighbor_node;
    std::vector<is_mesh::SimplexSet<is_mesh::NodeKey>*> nodes_neighbor_node;
    
    // Neighbor of edge
    std::vector<is_mesh::SimplexSet<is_mesh::TetrahedronKey>*> tets_share_edge;
    std::vector<bool*> edge_adapted;
    
    // face
    std::vector<is_mesh::SimplexSet<is_mesh::NodeKey>*> node_on_face;
    
    // tets
    std::vector<is_mesh::SimplexSet<is_mesh::NodeKey>*> nodes_on_tet;
    std::vector<real*> quality_tet;
    
    // Colors for parallel
    std::vector<int *> node_color;
    std::vector<int *> edge_color;
    
    // for adaptive time step
    std::vector<bool *> is_clean;
public:
    dsc_cache(){};
    
    void init(int MAX_ELEMENTS)
    {
        std::cout << "Init cache with max " << MAX_ELEMENTS << " elements\n";
        // Node
        link_of_node = std::vector<is_mesh::SimplexSet<is_mesh::FaceKey>*>(MAX_ELEMENTS, nullptr);
        tets_neighbor_node = std::vector<is_mesh::SimplexSet<is_mesh::TetrahedronKey>*>(MAX_ELEMENTS, nullptr);
        faces_neighbor_node = std::vector<is_mesh::SimplexSet<is_mesh::FaceKey>*>(MAX_ELEMENTS, nullptr);
        nodes_neighbor_node = std::vector<is_mesh::SimplexSet<is_mesh::NodeKey>*>(MAX_ELEMENTS, nullptr);
        node_color = std::vector<int *>(MAX_ELEMENTS, nullptr);
        is_clean = std::vector<bool *>(MAX_ELEMENTS, nullptr);
        
        // Edge
        tets_share_edge = std::vector<is_mesh::SimplexSet<is_mesh::TetrahedronKey>*>(MAX_ELEMENTS, nullptr);
        edge_adapted = std::vector<bool *>(MAX_ELEMENTS, nullptr);
        edge_color = std::vector<int *>(MAX_ELEMENTS, nullptr);
        
        // face
        node_on_face = std::vector<is_mesh::SimplexSet<is_mesh::NodeKey>*>(MAX_ELEMENTS, nullptr);
        
        // tet
        nodes_on_tet = std::vector<is_mesh::SimplexSet<is_mesh::NodeKey>*>(MAX_ELEMENTS, nullptr);
        quality_tet = std::vector<real*>(MAX_ELEMENTS, nullptr);
    }
    
    ~dsc_cache(){};
    
    void mark_dirty(is_mesh::NodeKey nk, bool dirty)
    {
        CLEAN_GARBAGE(tets_neighbor_node, nk);
        CLEAN_GARBAGE(faces_neighbor_node, nk);
        CLEAN_GARBAGE(nodes_neighbor_node, nk);
        CLEAN_GARBAGE(link_of_node, nk);
        CLEAN_GARBAGE(node_color, nk);
        CLEAN_GARBAGE(is_clean, nk);
    }
    
    void mark_dirty(is_mesh::EdgeKey ek, bool dirty)
    {
        CLEAN_GARBAGE(tets_share_edge, ek);
        CLEAN_GARBAGE(edge_adapted, ek);
        CLEAN_GARBAGE(edge_color, ek);
    }
    
    void mark_dirty(is_mesh::FaceKey fk, bool dirty)
    {
        CLEAN_GARBAGE(node_on_face, fk);
    }
    
    void mark_dirty(is_mesh::TetrahedronKey tk, bool dirty)
    {
        CLEAN_GARBAGE(nodes_on_tet, tk);
        CLEAN_GARBAGE(quality_tet, tk);
    }
    
    void mark_dirty_tet(is_mesh::TetrahedronKey tk)
    {
        CLEAN_GARBAGE(quality_tet, tk);
    }
};

#endif /* cache_hpp */
