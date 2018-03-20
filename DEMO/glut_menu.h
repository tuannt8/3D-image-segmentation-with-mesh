//
//  glut_menu.hpp
//  DSC_segment
//
//  Created by Tuan Nguyen Trung on 12/9/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#ifndef glut_menu_hpp
#define glut_menu_hpp

#include <stdio.h>
#include <string>
#include <map>

class glut_menu
{
public:

    
    static bool get_state(std::string name, bool default_state = false);
    
    
    ~glut_menu(){};
    
public:
    std::map<std::string, bool> _state_map;
    int menuId = -1;
    
public:
    void update_menu();
    
    /*
     For singleton - start
     */
public:
    static glut_menu & get_instance()
    {
        static glut_menu instance;
        return instance;
    }

    glut_menu(glut_menu const&) = delete;
    void operator = (glut_menu const&) = delete;
private:
    glut_menu(){};
    /* For singleton - end */
};

#endif /* glut_menu_hpp */
