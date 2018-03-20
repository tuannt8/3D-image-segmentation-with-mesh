//
//  glut_menu.cpp
//  DSC_segment
//
//  Created by Tuan Nguyen Trung on 12/9/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#include "glut_menu.h"
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
// #include <GLUT/GLUT.h>
#include <iostream>
#include <sstream>

using namespace std;
bool glut_menu::get_state(std::string name, bool default_state)
{
    if (get_instance()._state_map.find(name) == get_instance()._state_map.end())
    {
        // Not exist, create
        get_instance()._state_map.insert(std::make_pair(name, default_state));
        
        get_instance().update_menu();
        return default_state;
    }
    
    return get_instance()._state_map[name];
}

void mymenu(int idx) {
    auto it = glut_menu::get_instance()._state_map.begin();
    for (int i = 0; i < idx; i++)
    {
        it++;
    }
    it->second = !it->second;
    glut_menu::get_instance().update_menu();
}

void glut_menu::update_menu()
{
    if (menuId != -1)
    {
        glutDestroyMenu(menuId);
        menuId = -1;
    }
    
    
    menuId = glutCreateMenu(mymenu); // single menu, no need for id
    glutAddMenuEntry("===Display===", -1);
    int idx = 0;
    for (auto const &it : _state_map)
    {
        ostringstream tt;
        if (it.second){
            tt << "o | ";
        }else{
            tt << "  | ";
        }
        tt << it.first.c_str();
        glutAddMenuEntry(tt.str().c_str(), idx++);
    }
    glutAttachMenu(GLUT_RIGHT_BUTTON);
    
    glutPostRedisplay();
}
