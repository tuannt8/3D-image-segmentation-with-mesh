//
//  export.hpp
//  DEMO
//
//  Created by Tuan Nguyen Trung on 09/03/2018.
//  Copyright © 2018 Asger Nyman Christiansen. All rights reserved.
//

#ifndef export_hpp
#define export_hpp

#include <stdio.h>
#include "DSC.h"
#include <memory.h>
#include <string>
#include "define.h"

namespace dsc_export{
    typedef std::shared_ptr<DSC::DeformableSimplicialComplex<>> dsc_ptr;
    
    // Export upto 10 surface, individually
    void export_surface(dsc_ptr dsc, std::string path);
    void export_surface(std::string path);
    
    // Two phase fluid
    void export_two_phase_fluid(std::string path);
    
    // Load DSC
    dsc_ptr load_dsc(std::string path);
    
    void export_surface(dsc_ptr dsc, int phase, std::string path);
    
    // main phase is index [0]
    void export_shared_bound(dsc_ptr dsc, std::vector<vec3i> phases, std::string path);
    void export_shared_bound_no_bound(dsc_ptr dsc, std::vector<vec3i> phases, std::string path, vec3 ld, vec3 ru);
};

#endif /* export_hpp */
