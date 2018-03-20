//
//  test_cache.h
//  DSC_segment
//
//  Created by Tuan Nguyen Trung on 5/25/17.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#ifndef test_cache_h
#define test_cache_h

#include "DSC.h"
#include "util.h"

namespace test_dsc
{
    typedef DSC::DeformableSimplicialComplex<> dsc_class;
    
    void reduce_mesh_quality_by_moving_vertices(dsc_class & dsc){
        double move_portion = 10; // %
        long i = 0, j = 0;
        for (auto nit = dsc.nodes_begin(); nit != dsc.nodes_end(); nit++)
        {
            if ((rand()%100) < move_portion)
            {
                if(nit->is_boundary())
                    continue;
                
                vec3 pos = dsc.get_pos(nit.key());
                vec3 direct(rand()%100,rand()%100,rand()%100);
                vec3 destination = pos + direct;
                
                real max_l = dsc.intersection_with_link(nit.key(), destination);
                
                dsc.set_pos(nit.key(), pos + direct*max_l*0.7);
                i++;
            }
            j++;
        }
        
        std::cout << "Move " << i << "/" << j << " vertices" << std::endl;
    }
    
    void topological_edge_removal(dsc_class & dsc)
    {
        dsc.topological_edge_removal();
    }
    
    void topological_face_removal(dsc_class & dsc){
        dsc.topological_face_removal();
    }
    
    void split_edges(dsc_class &dsc){
        
    }
}



#endif /* test_cache_h */
