//
//  Deformabel Simplicial Complex (DSC) method
//  Copyright (C) 2013  Technical University of Denmark
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  See licence.txt for a copy of the GNU General Public License.

#pragma once

#include "DSC.h"
#include "velocity_function.h"
#include "log.h"
#include "draw.h"
#include "segment_function.h"
#ifdef __APPLE__
#include <GLUT/GLUT.h>
#endif
#include <condition_variable>

#include <thread>

#define VERSION "3.1"

#define PROBLEM_NAME "fuelcells_smaller"
/**
 A default application which utilizes OpenGL, GLEW and GLUT for visualization. Three sample velocity functions (rotation, smoothing and expansion) can be applies to a model specified by the model_file_name variable or as input variable. See https://github.com/asny/DSC/wiki/DEMO-instructions for details on how to use this DEMO application. See https://github.com/asny/DSC/wiki/Instructions for instructions on how to build your own application which uses the implementation of the DSC method.
 */
class UI
{
public:
    std::unique_ptr<DSC::VelocityFunc<>> vel_fun;
    std::unique_ptr<DSC::DeformableSimplicialComplex<>> dsc = nullptr;
    std::unique_ptr<Log> basic_log;
    std::unique_ptr<Painter> painter;
    
    std::string model_file_name = "armadillo";
    
    vec3 eye_pos = {70., 30., 70.};
    vec3 center_pos; // look at
    
    vec3 camera_pos = {30., 30., 70.};
    vec3 light_pos = {0., 0., 70.};
    
    int WIN_SIZE_X = 1280;
    int WIN_SIZE_Y = 720;
    
    bool CONTINUOUS = false;
    bool RECORD = false;
    bool QUIT_ON_COMPLETION = false;
    
    static UI* instance;
    
#ifdef _WIN32
    const std::string obj_path = "data\\";
    const std::string log_path = "LOG\\";
#else
    const std::string obj_path = "./data/";
    const std::string log_path = "./LOG/";
#endif
    
    segment_function _seg;
    
public:
    
    UI(int &argc, char** argv);
    UI(); //no interface
    void init_data();
    
    static UI* get_instance()
    {
        return instance;
    }
    
    void display();
    
    void animate();
    
    void reshape(int width, int height);
    
    void visible(int v);
    
    void mouse(int button, int state, int x, int y);
    
    void motion(int x, int y);
    
    void update_gl();
    void setup_light();
    
    /**
     The keyboard is used for all inputs. See https://github.com/asny/DSC/wiki/DEMO-instructions for up-to-date instructions on how to use the DEMO application.
     */
    void keyboard(unsigned char key, int x, int y);
    
private:
public:
    double _min_edge_length = 7;
    double m_edge_length = 7;
    float gl_dis_max;
    vec3 _obj_dim;
    vec3 _dsc_dim; // contain boundary gap
    
    GLfloat angle = -150;   /* in degrees */
    GLfloat angle2 = 30;   /* in degrees */
    int moving, startx, starty;
    int animation = 1;
    
    vec3 _mouse_pos;
    
    int phase_draw = 0;//Which phase to draw the interface
    void init_dsc();
    void set_dsc_boundary_layer();
    
    ///////////// Thread
    GLuint m_gl_sence;
    void update_draw_list();
    std::thread m_th;
    bool force_render = false;
    /**
     Loads the .dsc file specified by the model_file_name variable.
     */
    void load_model(const std::string& file_name);
    void save_model(std::string file_name = std::string());
    
    void export_segment(std::string file_name = std::string());
    /**
     Updates the window title.
     */
    void update_title();
    
    /**
     Starts the motion.
     */
    void start(const std::string& log_folder_name);
    
    /**
     Stops the motion and deletes the DSC object.
     */
    void stop();
    
    
    void load_config_file();
    
};
