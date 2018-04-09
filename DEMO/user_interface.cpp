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

#include "user_interface.h"
#include "rotate_function.h"
#include "average_function.h"
#include "normal_function.h"
#include "draw_helper.h"
#include "tetralizer.h"
#include "glut_menu.h"

#include "test_cache.h"

#include <math.h>       /* for cos(), sin(), and sqrt() */

#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>
#include <mutex>
#include <cctype>

#include "profile.h"

#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>


using namespace DSC;
using namespace std;

int m_iters = 0;

void display_(){
    UI::get_instance()->display();
}

void keyboard_(unsigned char key, int x, int y){
    UI::get_instance()->keyboard(key, x, y);
}

void keyboard_special_(int key, int x, int y){
    UI::get_instance()->keyboard((unsigned char)key, x, y);
}

void reshape_(int width, int height){
    UI::get_instance()->reshape(width, height);
}

void visible_(int v){
    UI::get_instance()->visible(v);
}

void animate_(){
    UI::get_instance()->animate();
}

void motion_(int x, int y)
{
    UI::get_instance()->motion(x, y);
}

void mouse_(int button, int state, int x, int y){
    UI::get_instance()->mouse(button, state, x, y);
}

vec3 GetOGLPos(int x, int y)
{
    GLint viewport[4];
    GLdouble modelview[16];
    GLdouble projection[16];
    GLfloat winX, winY, winZ;
    GLdouble posX, posY, posZ;
    
    glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
    glGetDoublev( GL_PROJECTION_MATRIX, projection );
    glGetIntegerv( GL_VIEWPORT, viewport );
    
    winX = (float)x;
    winY = (float)viewport[3] - (float)y;
    glReadPixels( x, int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );
    
    gluUnProject( winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ);
    
//    gluUnProject( x - viewport[0], y-viewport[1], 0, modelview, projection, viewport, &posX, &posY, &posZ);
    
    return vec3(posX, posY, posZ);
}

void UI::mouse(int button, int state, int x, int y)
{
    if (button == GLUT_LEFT_BUTTON) {
        if (state == GLUT_DOWN) {
            moving = 1;
            startx = x;
            starty = y;
        }
        if (state == GLUT_UP) {
            moving = 0;
        }
        
        if(state == GLUT_UP && glutGetModifiers() == GLUT_ACTIVE_SHIFT)
        {
            _mouse_pos = GetOGLPos(x, y);
            cout << _mouse_pos << endl;
        }
    }
}



void UI::motion(int x, int y)
{
    if (moving) {
        angle2= angle2 - (x - startx)*0.03;
        angle = angle + (y - starty)*0.03;

//        cout << angle*180/3.14 <<"/"<<angle2*180/3.14 << endl;
        
        startx = x;
        starty = y;
        glutPostRedisplay();
    }
    

}

UI* UI::instance = NULL;

void UI::setup_light()
{
    vec3 center = _obj_dim/ 2.0;
    vec3 eye = center + vec3(gl_dis_max*2.0*cos(angle)*cos(angle2),
                             gl_dis_max*2.0*cos(angle)*sin(angle2),
                             gl_dis_max*2.0*sin(angle));
    

    GLfloat light_position[] = { (GLfloat)eye[0], (GLfloat)eye[1], (GLfloat)eye[2], 0.0 };
    glShadeModel (GL_SMOOTH);
    
    GLfloat amb = 0.2, diff = 1., spec = 1.;
    
    GLfloat light_ambient[] = { amb,amb,amb, 1.0 };
    GLfloat light_diffuse[] = {diff, diff, diff, 1.0 };
    GLfloat light_specular[] = {spec, spec, spec, 1.0 };
    
//    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
//    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
//    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
    
//    GLfloat g_amb = .0;
//    GLfloat global_ambient[] = {g_amb, g_amb, g_amb, 0.1};
//    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
    
    
//    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
//    GLfloat mat_shininess[] = { 5.0 };
//    
//    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
//    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    
    glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE ) ;
    glEnable(GL_COLOR_MATERIAL);
    

}

int num_images;

//extern string config_file;

inline string get_option(std::map<std::string, std::string> & option, string key, bool bAcceptVoid = false)
{
    if (option.find(key) == option.end())
    {
        ostringstream err;
        err << "Can not read key: " << key;
        if (bAcceptVoid)
        {
            return "";
        }
        else
            throw err.str().c_str();
    }else{
        return option[key];
    }
}

#define get_opt(v,k) \
if(options.find(k) != options.end()){\
    if(typeid(v) == typeid(string)){ \
        v = options[k];\
    }\
    if(typeid(v) == typeid(int)){ \
        v = stoi(options[k]);\
    }\
    if(typeid(v) == typeid(double)){ \
        v = stof(options[k]);\
    }\
}else{\
cout<<"Could not load key: " << k << endl;\
}

//void UI::load_config_file()
//{
//}

void UI::update_draw_list()
{
}

UI::UI()
{
    init_data();
}

//    extern std::string config_file;

void UI::init_data()
{
    // Load cross sections
    _seg.init();
    
    bool test = false;
#ifdef __APPLE__
    test = true;
#endif
    
#ifdef INTENSITY_IMAGE
    _obj_dim = _seg._img.dimension_v();
#else
    _obj_dim = vec3(_seg.m_prob_img.m_dimension);
#endif
    _dsc_dim = _obj_dim + vec3(2*m_edge_length);
    
    cout << "Image dimension " << _obj_dim[0] << " " << _obj_dim[1] << " " <<  _obj_dim[2] << " " ;
    
    gl_dis_max = fmax(_obj_dim[0], fmax(_obj_dim[1], _obj_dim[2]));
    
    // Generate DSC
    if(!test)
        init_dsc();
    else
    {
//        load_model("./important_log/cement.dsc");
        load_model("./LOG/test/iter_340.dsc");
//        load_model("./LOG/bundle/bundle.dsc");

//        load_model(output_path + "/iter_100.dsc");
//        dsc->cache.init(dsc->get_no_faces_buffer()*4);
    }
    
    dsc->set_min_edge_length(m_edge_length);
    
    if(!test)
        set_dsc_boundary_layer();
    
    _seg._dsc = &*dsc;
    if(!test)
        _seg.initialization_discrete_opt();
    else
    {
//        for(auto tit = dsc->tetrahedra_begin(); tit != dsc->tetrahedra_end(); tit++)
//        {
//            auto l = dsc->get_label(tit.key());
//            if (l == 999)
//            {
//                dsc->set_label(tit.key(), 0);
//            }else{
//                dsc->set_label(tit.key(), l+1);
//            }
//        }

        _seg.update_average_intensity();
        _seg.compute_energy();
    }

    std::cout << "Mesh initialized: " << dsc->get_no_nodes() << " nodes; "
    << dsc->get_no_tets() << " tets" << endl;

}


UI::UI(InputParser p)
{
    bool bDisplay = !p.cmdOptionExists("-no-display");
    _seg.num_iter = std::stoi(p.getCmdOption("-no-iter", "500"));
    
    _seg.NB_PHASE = stoi(p.getCmdOption("-nb-phase", "5"));
    _seg.m_alpha = stof(p.getCmdOption("-alpha", "0.01"));
    _seg.m_max_dis = stof(p.getCmdOption("-max-dis", "0.3"));
    _seg._directory_path = p.getCmdOption("-data-path", "../Large_data/square_round");
#ifdef __linux__
    _seg._directory_path = std::string("../") + _seg._directory_path;
#endif
    m_edge_length = stof(p.getCmdOption("-edge-length", "20"));
    
    if (p.cmdOptionExists("-log-path"))
    {
        output_path = p.getCmdOption("-log-path", "./LOG");
    }else{
        std::ostringstream out;
#ifdef __APPLE__
        out << "./LOG/seg_" << std::setprecision(2) << m_edge_length << "_" << std::setprecision(3) << _seg.m_alpha;
#else
        out << "./seg_" << std::setprecision(2) << m_edge_length << "_"<< std::setprecision(3) << _seg.m_alpha;
#endif
        output_path = out.str();
    }
    cout << "Log to: " << output_path << endl;
    
    
    for (int i = 0; i < 10; i++)
    {
        
    }
    // create output path
    mkdir(output_path.c_str(), S_IRWXU);
    
    p.print();
    
    if(p.cmdOptionExists("-h"))
        exit(0);
    
    if (bDisplay)
    {
        instance = this;
        
        int a=0;
        char * b;
        glutInit(&a, &b);
        glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_STENCIL | GLUT_MULTISAMPLE);
        
        glutCreateWindow("3D segmentation");
        
        glutDisplayFunc(display_);
        glutKeyboardFunc(keyboard_);
        glutSpecialFunc(keyboard_special_);
        glutSetKeyRepeat(GLUT_KEY_REPEAT_ON);
        glutVisibilityFunc(visible_);
        glutReshapeFunc(reshape_);
        glutMotionFunc(motion_);
        glutMouseFunc(mouse_);
        
        
        glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
        glEnable(GL_DEPTH_TEST);
        glLineWidth(1.0);
        
        setup_light();
        
        glutReshapeWindow(WIN_SIZE_X, WIN_SIZE_Y);
        check_gl_error();
    }
    
    // Init data
    init_data();
    
    if (bDisplay)
    {
        // Update texture draw
#ifdef INTENSITY_IMAGE
        draw_helper::update_texture(_seg._img, 0,0,0);
#else
        draw_helper::update_texture(*_seg.m_prob_img.m_prob_map[0], 0,0,0);
#endif
        glutMainLoop();
    }
    else{
        for(int i = 0; i < _seg.num_iter; i++)
        {
            if(i%20==0)
                save_model(output_path + "/iter_" + std::to_string(i) + ".dsc");
            
            _seg.segment();
        }
    }
}

UI::UI(int &argc, char** argv)
{
    std::cout << "Version " << VERSION << endl;
    
    instance = this;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_STENCIL | GLUT_MULTISAMPLE);
    
    glutCreateWindow("3D segmentation");
    
    glutDisplayFunc(display_);
    glutKeyboardFunc(keyboard_);
    glutSpecialFunc(keyboard_special_);
    glutSetKeyRepeat(GLUT_KEY_REPEAT_ON);
    glutVisibilityFunc(visible_);
    glutReshapeFunc(reshape_);
    glutMotionFunc(motion_);
    glutMouseFunc(mouse_);


    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glEnable(GL_DEPTH_TEST);
    glLineWidth(1.0);
    
    setup_light();
    
	glutReshapeWindow(WIN_SIZE_X, WIN_SIZE_Y);
    check_gl_error();
    
    // Init data
    init_data();
    
    // Update texture draw
#ifdef INTENSITY_IMAGE
    draw_helper::update_texture(_seg._img, 0,0,0);
#else
    draw_helper::update_texture(*_seg.m_prob_img.m_prob_map[0], 0,0,0);
#endif
}

// Label the gap between DSC boundary and image boundary to BOUND_LABEL (999)
void UI::set_dsc_boundary_layer()
{
    std::vector<bool> is_tet_bound(dsc->get_no_tets_buffer(), false);
    for (auto nit =dsc->nodes_begin(); nit != dsc->nodes_end(); nit++)
    {
        if(nit->is_boundary())
        {
            auto tets = dsc->get_tets(nit.key());
            for (auto t : tets)
            {
                is_tet_bound[t] = true;
            }
        }
    }
    for(auto tit = dsc->tetrahedra_begin(); tit != dsc->tetrahedra_end(); tit++)
    {
        if (!is_tet_bound[tit.key()])
        {
            dsc->set_label(tit.key(), 1);
        }
    }
}

void UI::pad_boundary(double scale)
{

}

void UI::load_model(const std::string& file_name)
{
    std::cout << "\nLoading " << file_name << std::endl;
    dsc = nullptr;
    std::vector<vec3> points;
    std::vector<int>  tets;
    std::vector<int>  tet_labels;
    is_mesh::import_tet_mesh(file_name, points, tets, tet_labels);
    
    dsc = std::unique_ptr<DeformableSimplicialComplex<>>(new DeformableSimplicialComplex<>(points, tets, tet_labels));
    
    std::cout << "Loading done" << std::endl << std::endl;
}

void UI::export_segment(std::string file_name)
{
    if(file_name.empty())
    {
        file_name = std::string("segment.dsc");
    }
    
    std::ofstream file(file_name.data());
    
    std::map<int, int> vertices_index;
    int idx = 0;
    for (auto vit = dsc->nodes_begin(); vit != dsc->nodes_end(); vit++)
    {
        if (!vit->is_boundary())
        {
            vertices_index.insert(std::make_pair((int)vit.key(), idx++));
            auto pos = dsc->get_pos(vit.key());
            
            file << "v " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
        }
    }
    
    for (auto tit = dsc->tetrahedra_begin(); tit != dsc->tetrahedra_end(); tit++)
    {
        if (dsc->get_label(tit.key()) != BOUND_LABEL)
        {
            auto nodes = dsc->get_nodes(tit.key());
            file << "t " << vertices_index[(int)nodes[0]] << " "
                    << vertices_index[(int)nodes[1]] << " "
                    << vertices_index[(int)nodes[2]] << " "
                    << vertices_index[(int)nodes[3]] << " "
            << dsc->get_label(tit.key()) << std::endl;
        }
    }
}

void UI::save_model( std::string file_name){
    if(file_name.empty())
    {
        // Find a file name to avoid overwriting
        for (int i = 0; i < 100; i++)
        {
            std::string temp_path = log_path + std::string(PROBLEM_NAME) + std::to_string(i) + std::string(".dsc");
            if (!std::ifstream(temp_path))
            {
                file_name = temp_path;
                break;
            }
        }
        if (file_name.empty())
        {
            cout<<"There are over 100 files in the log folder. Clean temporary saved file to save more file" << endl;
            return;
        }
    }
    
    cout << "Save to " << file_name << endl;
    
    std::vector<vec3> points;
    std::vector<int> faces;
    std::vector<int> tets;
    std::vector<int> tet_labels;
    dsc->extract_tet_mesh(points, tets, tet_labels);
    is_mesh::export_tet_mesh(file_name, points, tets, tet_labels);
}

#define index_cube(x,y,z) ((z)*NX*NY + (y)*NX + (x))
void UI::init_dsc()
{
    std::vector<vec3> points;
    std::vector<int> tets;
    std::vector<int> tet_labels;
    

    double delta = m_edge_length;
    
    cout << "delta " << delta << endl;
    cout << _dsc_dim[0] << " " << _dsc_dim[1] << " " <<  _dsc_dim[2] << " " ;
    
    int NX = round(_dsc_dim[0] / delta) + 1; // number of vertices
    int NY = round(_dsc_dim[1] / delta) + 1;
    int NZ = round(_dsc_dim[2] / delta) + 1;
    
        cout << "Compute point" << NX << " " << NY << " " << NZ << "\n";

    double deltax = _dsc_dim[0]/(double)(NX-1);
    double deltay = _dsc_dim[1]/(double)(NY - 1);
    double deltaz = _dsc_dim[2]/(double)(NZ - 1);
    
    
    // points
    for (int iz = 0; iz < NZ; iz++)
    {
        for (int iy = 0; iy < NY; iy++)
        {
            for (int ix = 0; ix < NX; ix++)
            {
                points.push_back(vec3(ix*deltax, iy*deltay, iz*deltaz) - vec3(m_edge_length));
            }
        }
    }

    cout << "Compute tets\n";
    
    // tets
    for (int iz = 0; iz < NZ - 1; iz++)
    {
        for (int iy = 0; iy < NY - 1; iy++)
        {
            for (int ix = 0; ix < NX - 1; ix++)
            {
                // 8 vertices
                int vertices[] = {
                    index_cube(ix, iy, iz),
                    index_cube(ix+1, iy, iz),
                    index_cube(ix+1, iy+1, iz),
                    index_cube(ix, iy+1, iz),
                    index_cube(ix, iy, iz + 1),
                    index_cube(ix+1, iy, iz + 1),
                    index_cube(ix+1, iy+1, iz + 1),
                    index_cube(ix, iy+1, iz + 1)
                };
                
                int tetras[] = {
                    0, 4, 5, 7,
                    0, 7, 5, 1,
                    0, 1, 3, 7,
                    1, 5, 6, 7,
                    1, 6, 7, 3,
                    1, 2, 6, 3
                };
                
                for(int i = 0; i < 6*4; i++)
                {
                    tets.push_back(vertices[tetras[i]]);
                }
            }
        }
    }
    
    long nbTet = tets.size()/4;
    tet_labels = std::vector<int>(nbTet, 0);
    
    cout << "Init DSC from point\n";
    
    dsc = std::unique_ptr<DeformableSimplicialComplex<>>(new DeformableSimplicialComplex<>(points, tets, tet_labels));
    dsc->set_avg_edge_length(delta);
}

void UI::update_title()
{
    std::ostringstream oss;
    oss << "3D DSC\t" << vel_fun->get_name() << ", Time step " << vel_fun->get_time_step();
    oss << " (Nu = " << vel_fun->get_velocity() << ", Delta = " << dsc->get_avg_edge_length() << ", Alpha = " << vel_fun->get_accuracy() << ")";
    std::string str(oss.str());
    glutSetWindowTitle(str.c_str());
}

void UI::update_gl()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective( /* field of view in degree */ 40.0,
                   /* aspect ratio */ 1.0,
                   /* Z near */ 10.0, /* Z far */ gl_dis_max*5.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    vec3 center = _obj_dim/ 2.0;
    vec3 eye = center + vec3(gl_dis_max*3.0*cos(angle)*cos(angle2),
                             gl_dis_max*3.0*cos(angle)*sin(angle2),
                             gl_dis_max*3.0*sin(angle));
    vec3 head = vec3(-sin(angle)*cos(angle2),
                     -sin(angle)*sin(angle2),
                     cos(angle));
    gluLookAt(eye[0], eye[1], eye[2],
              center[0], center[1], center[2],
              head[0], head[1], head[2]);
    
    int size = std::min(WIN_SIZE_Y, WIN_SIZE_X);
    glViewport((WIN_SIZE_X-size)/2.0, (WIN_SIZE_Y-size)/2.0, size, size);
    
    glClearColor(1.0, 1.0, 1.0, 1.0);
    
    eye_pos = eye;
    center_pos = center;
    
    //
//    glEnable(GL_BLEND);
//    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//    
//    glEnable(GL_MULTISAMPLE_ARB);
//    glEnable(GL_MULTISAMPLE);
//    glEnable(GL_POINT_SMOOTH);
//    glEnable(GL_LINE_SMOOTH);
//    glEnable(GL_POLYGON_SMOOTH);
//
//    glHint( GL_POLYGON_SMOOTH_HINT,GL_NICEST );
//    glHint( GL_POINT_SMOOTH,GL_NICEST );
//    glHint( GL_LINE_SMOOTH,GL_NICEST );
}

void UI::display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    update_gl();
    setup_light();
    
    static double total_time = 100;
    static auto init_time = std::chrono::system_clock::now();
    std::chrono::duration<real> t = std::chrono::system_clock::now() - init_time;
    total_time += t.count();
    init_time = std::chrono::system_clock::now();

    //
    if (glut_menu::get_state("Ray line", 0))
    {
        glDisable(GL_LIGHTING);
        glPointSize(2.0);
        glBegin(GL_POINTS);
        glColor3f(1, 0, 0);
        for (auto r : _seg._d_rayz)
        {
            for (int i = 0; i < r.intersects.size()/2; i++)
            {
                for (int j = r.intersects[2*i].z; j < r.intersects[2*i + 1].z; j++)
                {
                    glVertex3f(r.x, r.y, j);
                }
            }
        }
        glEnd();
        glEnable(GL_LIGHTING);
    }
    if (glut_menu::get_state("boundary displace", 0))
    {
        auto ff = _seg.boundary_vertices_displacements;
        glBegin(GL_LINES);
        glColor3d(1, 0, 0);
        for (auto nit = dsc->nodes_begin(); nit != dsc->nodes_end(); nit++)
        {
            if (ff[nit.key()].length() > EPSILON)
            {
                auto pos = nit->get_pos();
                glVertex3dv(pos.get());
                glVertex3dv((pos+ff[nit.key()]*4).get());
            }
        }
        glEnd();
    }
    
    if (glut_menu::get_state("tet cross", 0))
    {
#ifdef INTENSITY_IMAGE
        draw_helper::draw_tet_cross(*dsc, _seg.NB_PHASE, _seg._img.dimension_v());
#else
        draw_helper::draw_tet_cross(*dsc, _seg.NB_PHASE, vec3(_seg.m_prob_img.m_dimension));
#endif
    }
    
    if (glut_menu::get_state("Cross section", 0))
    {
        draw_helper::draw_cross(*dsc, _seg.NB_PHASE,  _seg._img.dimension_v());
    }
    
    if (glut_menu::get_state("Surface curvature", 0))
    {
        draw_helper::draw_curvature(*dsc, _seg._mean_curvature_of_each_hat, phase_draw+1, _seg._mean_curvature_label);
    }

    if (glut_menu::get_state("Transparent surface", 0))
    {
        draw_helper::draw_transparent_surface(*dsc, _seg.NB_PHASE);
    }

    if (glut_menu::get_state("Triple interface", 0))
    {
        draw_helper::draw_triple_interface(*dsc);
    }

    if (glut_menu::get_state("Ray cross section", 0))
    {
        glDisable(GL_LIGHTING);
        glPointSize(2.0);
        glBegin(GL_POINTS);
        glColor3f(1, 0, 0);
        auto zz = draw_helper::get_instance()._cur_cross_poss[2];
        for (auto r : _seg._d_rayz)
        {
            for (int i = 0; i < r.intersects.size()/2; i++)
            {
                if (zz > r.intersects[2*i].z and zz < r.intersects[2*i + 1].z)
                {
                    glVertex3f(r.x, r.y, zz);
                }
            }

        }
        glEnd();
        glEnable(GL_LIGHTING);
    }
    
    if (glut_menu::get_state("Draw DEMO", 1))
    {
        int idx[]={1,2};
        
        glDisable(GL_LIGHTING);
        glColor3f(1, 0, 0);
        draw_helper::dsc_draw_triple_edge(*dsc);
        
        glColor3f(0.5, 0.3, 1.0);
        glLineWidth(2.0);
        draw_helper::dsc_draw_one_interface_edge(*dsc, idx[0]);
        draw_helper::dsc_draw_one_interface_edge(*dsc, idx[1]);
        
        //        glEnable(GL_CULL_FACE);
        glEnable(GL_LIGHTING);
        //        glEnable(GL_BLEND);
        //        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glColor4f(0.7, 0.7, 0.7, 1.0);
        draw_helper::dsc_draw_one_interface(*dsc, idx[0]);
        
        glColor4f(1, 0.7, 0.7, 1.0);
        draw_helper::dsc_draw_one_interface(*dsc, idx[1]);
        

        
        glDisable(GL_BLEND);
    }

    if (glut_menu::get_state("Draw DSC single interface edge", 0))
    {
//        glDisable(GL_CULL_FACE);
        glDisable(GL_LIGHTING);
        glColor3f(0.5, 0.3, 1.0);
        draw_helper::dsc_draw_one_interface_edge(*dsc, phase_draw+1);
    }

    if (glut_menu::get_state("Draw DSC single interface", 0))
    {
//        glEnable(GL_CULL_FACE);
        glEnable(GL_LIGHTING);
//        glEnable(GL_BLEND);
//        glBlendFunc(GL_SRC_ALPHA, GL_SRC_ALPHA);
//        glBlendEquation(GL_FUNC_ADD);
        glColor4f(0.7, 0.7, 0.7, 0.3);
        draw_helper::dsc_draw_one_interface(*dsc, phase_draw+1);
        glDisable(GL_BLEND);
    }


    if (glut_menu::get_state("Draw DSC edges", 0))
    {
        glColor3f(0, 0, 0);
        draw_helper::dsc_draw_edge(*dsc);
    }

    if (glut_menu::get_state("Draw DSC domain", 0))
    {
        glColor3f(0.3, 0.3, 0.3);
        draw_helper::dsc_draw_domain(*dsc);
    }

    if (glut_menu::get_state("Draw DSC interface edge", 0))
    {
        glDisable(GL_LIGHTING);
        glColor3f(0.0, 0.0, 1.0);
        draw_helper::dsc_draw_interface_edge(*dsc);
    }

    if (glut_menu::get_state("Draw DSC interface", 0))
    {
        glDisable(GL_CULL_FACE);
        glEnable(GL_LIGHTING);
        draw_helper::dsc_draw_interface(*dsc);
    }


    if (glut_menu::get_state("Draw DSC face normal", 0))
    {
        draw_helper::dsc_draw_face_norm(*dsc);
    }


    if (glut_menu::get_state("Draw Image slide", 1))
    {
        draw_helper::draw_image_slice( _seg._img);
    }

    if (glut_menu::get_state("Draw internal force", 0))
    {
        glColor3f(1, 0, 0);
        draw_helper::dsc_draw_node_arrow(*dsc, _seg._internal_forces);
    }

    if (glut_menu::get_state("area force", 0))
    {
        glColor3f(0, 0, 1);
        draw_helper::dsc_draw_node_multi_arrow(*dsc, _seg._area_force, 0.002);
    }

//    if (glut_menu::get_state("Mean curvature", 0))
//    {
//        glColor3f(0, 1, 1);
//        draw_helper::dsc_draw_node_multi_arrow(*dsc, _seg._curvature_force, 5);
//    }

    if (glut_menu::get_state("Draw boundary destination", 0))
    {
        draw_helper::draw_boundary_destination(_seg, &*dsc);
    }

    if (glut_menu::get_state("Interface Vertices indices", 0))
    {
        draw_helper::draw_dsc_interface_vertices_indices(*dsc, phase_draw+1);
    }

    if (glut_menu::get_state("Tets indices", 0))
    {
        draw_helper::draw_dsc_tet_indices(*dsc);
    }

    if (glut_menu::get_state("Debug", 0))
    {
        // Draw incorect interface
        static is_mesh::SimplexSet<is_mesh::FaceKey> wrong_face;
        if(wrong_face.size()== 0)
        {
            for(auto tit = dsc->tetrahedra_begin(); tit != dsc->tetrahedra_end(); tit++)
            {
                assert(dsc->exists(tit.key()));
                // Check faces:
                auto faces = dsc->get_faces(tit.key());
                assert(faces.size() == 4);
                for (auto f : faces) {
                    assert(dsc->exists(f));
                    auto cotets = dsc->get_tets(f);
                    int labels[2] = {dsc->get_label(cotets[0]), dsc->get_label(cotets[1])};
                    if((dsc->get(f).is_boundary() && cotets.size() == 2) || (!dsc->get(f).is_boundary() && cotets.size() == 1)
                       )
                    {
                        wrong_face += f;
                        goto OUT;
                    }
                }
            }
        OUT:
            cout << endl << wrong_face.size() << " invalid faces" << endl;
        }
        glColor3f(1, 0, 0);
        glBegin(GL_TRIANGLES);
        for(auto f : wrong_face)
        {
            for( auto pos : dsc->get_pos(dsc->get_nodes(f)))
                glVertex3dv(pos.get());
            
        }
        glEnd();
    }
  
    glutSwapBuffers();
    
    std::ostringstream os;
    os << m_iters;
    glutSetWindowTitle(os.str().c_str());

    check_gl_error();
    
    if(CONTINUOUS)
    {
        if (m_iters %20 ==0)
        {
            save_model(output_path + "/iter_" + std::to_string(m_iters) + ".dsc");
        }
        
        if (m_iters > _seg.num_iter)
        {
            save_model(output_path +  "/output.dsc");
            exit(0);
        }
     
        _seg.segment();
        
        m_iters++;
    }
}

void UI::reshape(int width, int height)
{
    WIN_SIZE_X = width;
    WIN_SIZE_Y = height;
    
    update_gl();
    
    glutReshapeWindow(WIN_SIZE_X, WIN_SIZE_Y);
}

void UI::animate()
{
}

void UI::keyboard(unsigned char key, int x, int y) {
    switch(key) {
        case GLUT_KEY_UP:
            draw_helper::update_texture(_seg._img, 0,0,1);
            break;
        case GLUT_KEY_DOWN:
            draw_helper::update_texture(_seg._img, 0,0,-1);
            break;
        case ' ':
            CONTINUOUS = !CONTINUOUS;
            break;
        case '\t':
            draw_helper::save_painting(WIN_SIZE_X, WIN_SIZE_Y);
            break;
        case 's':
            save_model();
            break;
        case 'i':
            dsc->validity_check();
            break;
        case 'l':
        {
        }
            break;
        case 'p':// Display time counter
            profile::close();
            break;
        case 'v':// Change surface type
            phase_draw = (phase_draw+1) % _seg.NB_PHASE;
            break;
        case 'u':
        {
            _seg.update_average_intensity();
            _seg.thickenning_surface();
        }
            break;
        case 't':
            dsc->deform();
            break;
        case 'a':
            std::cout << "Mesh initialized: " << dsc->get_no_nodes() << " nodes; "
            << dsc->get_no_tets() << " tets" << endl;
            break;
        case 'z':
            dsc->adapt();
            break;
        case 'd':
            _seg.update_vertex_boundary();
            _seg.adapt_surface();
//            _seg.force_snapp();
            break;
        default:
            break;
    }
    
    
}

void UI::take_screen_shot(string name)
{
#ifdef __APPLE__
    if (name.length() == 0)
    {
        // date and time
        auto now = std::chrono::system_clock::now();
        auto in_time_t = std::chrono::system_clock::to_time_t(now);
        std::stringstream ss;
        ss << "./LOG/"<< std::put_time(std::localtime(&in_time_t), "%Y%m%d_%X") << ".png";
    }
    int success = SOIL_save_screenshot(name.c_str(), SOIL_SAVE_TYPE_PNG, 0, 0, WIN_SIZE_X, WIN_SIZE_Y);
#endif
}


void idle(void)
{
    glutPostRedisplay();
}

void UI::visible(int v)
{
    if (v == GLUT_VISIBLE) {
        if (animation)
            glutIdleFunc(idle);
    } else {
        if (!animation)
            glutIdleFunc(NULL);
    }
}

void UI::stop()
{
    if(RECORD && basic_log)
    {
        basic_log->write_message("MOTION STOPPED");
        basic_log->write_log(*dsc);
        basic_log->write_log(*vel_fun);
        basic_log->write_timings(*vel_fun);
        
        std::vector<vec3> points;
        std::vector<int> faces;
        std::vector<int> tets;
        std::vector<int> tet_labels;
        dsc->extract_tet_mesh(points, tets, tet_labels);
        is_mesh::export_tet_mesh(basic_log->get_path() + std::string("/mesh.dsc"), points, tets, tet_labels);
        points.clear();
        dsc->extract_surface_mesh(points, faces);
        is_mesh::export_surface_mesh(basic_log->get_path() + std::string("/mesh.obj"), points, faces);
        basic_log = nullptr;
    }
    
    CONTINUOUS = false;
    update_title();
    glutPostRedisplay();
}

void UI::start(const std::string& log_folder_name)
{
    glutPostRedisplay();
}
