//
//  draw_helper.cpp
//  DSC_segment
//
//  Created by Tuan Nguyen Trung on 12/9/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#include "draw_helper.h"
#ifdef __APPLE__
#include <GLUT/GLUT.h>
#endif
#include "DSC.h"
#include <iostream>
#include <string>
#include <SOIL/SOIL.h>
#include <fstream>
#include "define.h"

void draw_helper::draw_image_slice(const image3d & im)
{
    glDisable(GL_LIGHTING);
    // Draw bounding box
    auto dim = im.dimension_v();
    
    glPushMatrix();
    
    glTranslatef(dim[0]/2.0, dim[1]/2.0, dim[2]/2.0);
    glScalef(dim[0], dim[1], dim[2]);
    
    glColor3f(1, 0, 0);
    glutWireCube(1.0);
    
    glPopMatrix();
    
    // slice
    vec3 ld(0,0, get_instance()._cur_cross_poss[2]);
    vec3 ru(dim[0], dim[1], get_instance()._cur_cross_poss[2]);
    
    glEnable(GL_TEXTURE_2D);
    glColor3f(1, 1, 1);
    glBegin(GL_QUADS);

    glTexCoord3f(0., 0., 0.);
    glVertex3f(ld[0], ld[1], ld[2]);

    glTexCoord3f(1., 0., 0.);
    glVertex3f(ru[0], ld[1], ld[2]);

    glTexCoord3f(1., 1., 0.);
    glVertex3f(ru[0], ru[1], ld[2]);

    glTexCoord3f(0., 1., 0.);
    glVertex3f(ld[0], ru[1], ld[2]);

    glEnd();
    glDisable(GL_TEXTURE_2D);
    
//    // border
//    glColor3f(1, 0, 0);
//    glBegin(GL_LINE_LOOP);
//    glVertex3f(ld[0], ld[1], ld[2]);
//    glVertex3f(ru[0], ld[1], ld[2]);
//    glVertex3f(ru[0], ru[1], ld[2]);
//    glVertex3f(ld[0], ru[1], ld[2]);
//    glEnd();
}

void RenderString(vec3 pos, const std::string &string)
{
    glColor3d(1.0, 0.0, 0.0);
    glRasterPos3d(pos[0], pos[1], pos[2]);
    for (int n=0; n<string.size(); ++n) {
        glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, string[n]);
    }
}

void draw_helper::draw_dsc_tet_indices(dsc_class& dsc)
{
    for (auto tit = dsc.tetrahedra_begin(); tit != dsc.tetrahedra_end(); tit++)
    {
        auto pts = dsc.get_pos(dsc.get_nodes(tit.key()));
        auto midPt = (pts[0] + pts[1] + pts[2] + pts[3])*0.25;
        
        RenderString(midPt, std::to_string((int)tit.key()));
    }
}


void draw_helper::draw_tet_list(dsc_class& dsc, std::vector<int> tet_list)
{
    for (auto tk : tet_list)
    {
        is_mesh::TetrahedronKey tkey(tk);
        
        auto pts = dsc.get_pos(dsc.get_nodes(tkey));
        auto midPt = (pts[0] + pts[1] + pts[2] + pts[3])*0.25;
        
        glColor3f(0, 0, 1);
        RenderString(midPt, std::to_string((int)tkey));
        
        glColor3f(1, 0, 0);
        glBegin(GL_LINES);
        for (auto e : dsc.get_edges(tkey))
        {
            auto pt_e = dsc.get_pos(dsc.get_nodes(e));
            glVertex3dv(pt_e[0].get());
            glVertex3dv(pt_e[1].get());
            
        }
        glEnd();
        
    }
    
}
void draw_helper::draw_dsc_interface_vertices_indices( dsc_class &dsc, int phase)
{
    std::vector<bool> bIs_interface_nodes(dsc.get_no_nodes_buffer(), false);
    for (auto f = dsc.faces_begin(); f != dsc.faces_end(); f++)
    {
        if (f->is_interface() && !f->is_boundary())
        {
            auto tets = dsc.get_tets(f.key());
            if (!(dsc.get_label(tets[0]) == phase
                  or dsc.get_label(tets[1]) == phase))
            {
                continue;
            }
            
#ifdef DSC_CACHE
            auto nodes=*dsc.get_nodes_cache(f.key());
#else
            auto nodes=dsc.get_nodes(f.key());
#endif
            for (auto n : nodes)
            {
                bIs_interface_nodes[n] = true;
            }
        }
    }
    
    for (auto nit = dsc.nodes_begin(); nit != dsc.nodes_end(); nit++)
    {
        if (bIs_interface_nodes[nit.key()])
        {
            auto pos = nit->get_pos();
            RenderString(pos, std::to_string((int)nit.key()));
        }
    }
}

void draw_helper::dsc_draw_node_multi_arrow(dsc_class & dsc, std::vector<std::vector<vec3>> arrows, double scale)
{
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    
    for (auto nit = dsc.nodes_begin(); nit != dsc.nodes_end(); nit++)
    {
        auto arrow = arrows[nit.key()];
        for(auto a : arrow)
        {
            glVertex3dv(nit->get_pos().get());
            glVertex3dv((nit->get_pos() + a*scale).get());
        }
    }
    glEnd();
}

void draw_helper::dsc_draw_node_arrow(dsc_class & dsc, std::vector<vec3> arrow)
{
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    
    for (auto nit = dsc.nodes_begin(); nit != dsc.nodes_end(); nit++)
    {
        if (arrow[nit.key()].length() > 0.1)
        {
            glVertex3dv(nit->get_pos().get());
            glVertex3dv((nit->get_pos() + arrow[nit.key()]).get());
        }
    }
    glEnd();
}

void draw_helper::draw_boundary_destination(segment_function &_seg, dsc_class *dsc)
{
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    glColor3f(1, 0, 0);
    for (auto nit = dsc->nodes_begin(); nit != dsc->nodes_end(); nit++)
    {
        if (_seg.d_is_image_boundary[nit.key()])
        {
            glVertex3dv(nit->get_pos().get());
            glVertex3dv((nit->get_pos() + _seg.boundary_vertices_displacements[nit.key()]).get());
        }
    }
    glEnd();
}

void draw_helper::draw_boundary_direction(segment_function &_seg, dsc_class *dsc)
{
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    for (auto nit = dsc->nodes_begin(); nit != dsc->nodes_end(); nit++)
    {
        if (_seg.d_is_image_boundary[nit.key()])
        {
            auto state = _seg.d_direction_state[nit.key()];
            
            auto a = (state & X_direction).to_ulong();
            
            if ((state & X_direction).to_ulong() != 0)
            {
                glColor3f(1, 0, 0);
                glVertex3dv(nit->get_pos().get());
                glVertex3dv((nit->get_pos() + vec3(4,0,0)).get());
            }
            if ((state & Y_direction).to_ulong() != 0)
            {
                glColor3f(0, 1, 0);
                glVertex3dv(nit->get_pos().get());
                glVertex3dv((nit->get_pos() + vec3(0,4,0)).get());
            }
            if ((state & Z_direction).to_ulong() != 0)
            {
                glColor3f(0, 0, 1);
                glVertex3dv(nit->get_pos().get());
                glVertex3dv((nit->get_pos() + vec3(0,0,4)).get());
            }
        }
    }
    glEnd();
}

void draw_helper::draw_coord(float length)
{
    glBegin(GL_LINES);
    glColor3f(1, 0, 0);
    glVertex3f(length, 0., 0.);
    glVertex3f(0., 0., 0.);
    
    glColor3f(0, 1, 0);
    glVertex3f(0., length, 0.);
    glVertex3f(0., 0., 0.);
    
    glColor3f(0, 0, 1.);
    glVertex3f(0., 0., length);
    glVertex3f(0., 0., 0.);
    glEnd();
}

void draw_helper::dsc_draw_edge(dsc_class &dsc)
{
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    for (auto eit = dsc.edges_begin(); eit != dsc.edges_end(); eit++)
    {
        auto pos = dsc.get_pos(dsc.get_nodes(eit.key()));
        glVertex3dv(pos[0].get());
        glVertex3dv(pos[1].get());
    }
    glEnd();
    glEnable(GL_LIGHTING);
}

void draw_helper::dsc_draw_face_norm(dsc_class & dsc)
{
    glColor3f(0, 0, 1);
    for (auto fid = dsc.faces_begin(); fid != dsc.faces_end(); fid++)
    {
        if (fid->is_interface())
        {
            auto pts = dsc.get_pos(dsc.get_nodes(fid.key()));
            auto center = (pts[0] + pts[1] + pts[2]) / 3.0;

            vec3 Norm = dsc.get_normal(fid.key());

            
            glBegin(GL_LINES);
            glVertex3dv(center.get());
            glVertex3dv((center + Norm*5).get());
            glEnd();
        }
    }
}

void draw_helper::dsc_draw_interface_edge(dsc_class & dsc)
{

    
    glBegin(GL_LINES);
    for (auto eit = dsc.edges_begin(); eit != dsc.edges_end(); eit++)
    {
        if (eit->is_interface())
        {
            glColor3f(0, 0, 1);
            if (eit->is_crossing())
            {
                glColor3f(1, 0, 0);
            }
            auto pos = dsc.get_pos(dsc.get_nodes(eit.key()));
            glVertex3dv(pos[0].get());
            glVertex3dv(pos[1].get());
        }
    }
    
    glEnd();
}

void draw_helper::save_painting(int WIDTH, int HEIGHT, std::string folder)
{
    std::ostringstream s;
    s << folder << "/scr";
    int i = 0;
    while (1)
    {
        std::ostringstream name;
        name << s.str() << "_" << i << ".png";
        std::ifstream file(name.str().c_str());
        if (!file)
        {
            // could not open
            s << "_" << i << ".png";
            break;
        }
        file.close();
        i++;
    }
    
    int success = SOIL_save_screenshot(s.str().c_str(), SOIL_SAVE_TYPE_PNG, 0, 0, WIDTH, HEIGHT); 
    if(!success)
    {
        std::cout << "ERROR: Failed to take screen shot: " << s.str().c_str() << std::endl;
        return;
    }
    else{
        std::cout <<"Screen shot: " << s.str().c_str() << std::endl;
    }
}

enum interface_type
{
    type_0_2 = 0,
    type_0_1,
    type_1_2
};

std::vector<vec3> draw_helper::node_normal_vector;
void draw_helper::update_normal_vector_interface(dsc_class & dsc, int phase, vec3 eye_pos)
{
#ifdef _DSC_ORIGIN_
    std::vector<int> neighbor_faces_count(MAX_NUM_ELEMENT_MESH, 0);
#else
    node_normal_vector = std::vector<vec3>(dsc.get_no_nodes_buffer(), vec3(0.0));
    std::vector<int> neighbor_faces_count(dsc.get_no_nodes_buffer(), 0);
#endif
    
    for (auto f = dsc.faces_begin(); f != dsc.faces_end(); f++)
    {
        if (f->is_interface() && !f->is_boundary())
        {
            auto tets = dsc.get_tets(f.key());
            
            if (!(dsc.get_label(tets[0]) == phase
                  or dsc.get_label(tets[1]) == phase))
            {
                continue;
            }
#ifdef DSC_CACHE
            auto nodes = *dsc.get_nodes_cache(f.key());
#else
            auto nodes = dsc.get_nodes(f.key());
#endif
            
            
            // normalize the normal to the eye
            is_mesh::TetrahedronKey other_tet = (dsc.get_label(tets[0]) == phase)? tets[0] : tets[1];
            
            auto norm = dsc.get_normal(f.key(), other_tet);
            
            
            for(auto n : nodes)
            {
                node_normal_vector[n] += norm;
                neighbor_faces_count[n]++;
            }
        }
    }
    
    for(int i = 0; i < node_normal_vector.size(); i++)
    {
        if(neighbor_faces_count[i] > 0)
            node_normal_vector[i] = Util::normalize(node_normal_vector[i]);
    }
}

void draw_helper::dsc_draw_debug_node(is_mesh::NodeKey nk)
{
    
}

void draw_helper::dsc_draw_one_interface_edge(dsc_class & dsc, int phase)
{
    for (auto f = dsc.faces_begin(); f != dsc.faces_end(); f++)
    {
        if (f->is_interface() && !f->is_boundary())
        {
            auto tets = dsc.get_tets(f.key());
            if (!(dsc.get_label(tets[0]) == phase
                  or dsc.get_label(tets[1]) == phase)
                //                || dsc.get_label(tets[0]) == BOUND_LABEL
                //                || dsc.get_label(tets[1]) == BOUND_LABEL
                )
            {
                continue;
            }
#ifdef DSC_CACHE
            auto nodes=*dsc.get_nodes_cache(f.key());
#else
            auto nodes=dsc.get_nodes(f.key());
#endif
            auto pts = dsc.get_pos(nodes);
            
            // normalize the normal to the eye
            is_mesh::TetrahedronKey other_tet = (dsc.get_label(tets[0]) == phase)? tets[0] : tets[1];
            
            
            auto norm = dsc.get_normal(f.key(), other_tet);
            
            // Draw triangle
            glBegin(GL_LINES);
            for (int i =0; i < 3; i++)
            {
                auto v0 = pts[i];
                auto v1 = pts[(i+1)%3];
                
                glVertex3dv(v0.get());
                glVertex3dv(v1.get());
            }
            glEnd();
            
            
        }
    }
}

void draw_helper:: dsc_draw_one_interface_no_share(dsc_class & dsc, int phase)
{
    for (auto f = dsc.faces_begin(); f != dsc.faces_end(); f++)
    {
        if (f->is_interface() && !f->is_boundary())
        {
            auto tets = dsc.get_tets(f.key());
            int label[2] = {dsc.get_label(tets[0]), dsc.get_label(tets[1])};
            if (!(label[0] == phase
                  or label[1] == phase)
                )
            {
                continue;
            }
            
            int other_label = label[0] == phase? label[1] : label[0];
//            if (other_label < phase)
//            {
//                continue;
//            }
            
#ifdef DSC_CACHE
            auto nodes=*dsc.get_nodes_cache(f.key());
#else
            auto nodes=dsc.get_nodes(f.key());
#endif
            auto pts = dsc.get_pos(nodes);
            
            // normalize the normal to the eye
            is_mesh::TetrahedronKey other_tet = (dsc.get_label(tets[0]) == phase)? tets[1] : tets[0];
            
            
            auto norm = -dsc.get_normal(f.key(), other_tet);
            
            
            
            // Draw triangle
            glEnable(GL_LIGHTING);
            glBegin(GL_TRIANGLES);
            for (int i =0; i < 3; i++)
            {
                auto v = pts[i];
                auto n = nodes[i];
                vec3 real_norm = norm;
                
                glNormal3dv(real_norm.get());
                glVertex3dv(v.get());
            }
            glEnd();
            
        }
    }
}

void draw_helper::dsc_draw_one_interface(dsc_class & dsc, int phase)
{
    
    glBegin(GL_TRIANGLES);
    for (auto f = dsc.faces_begin(); f != dsc.faces_end(); f++)
    {
        if (f->is_interface() && !f->is_boundary())
        {
            auto tets = dsc.get_tets(f.key());
            if (!(dsc.get_label(tets[0]) == phase
                  or dsc.get_label(tets[1]) == phase)
                )
            {
                continue;
            }
#ifdef DSC_CACHE
            auto nodes=*dsc.get_nodes_cache(f.key());
#else
            auto nodes=dsc.get_nodes(f.key());
#endif
            auto pts = dsc.get_pos(nodes);
            
            is_mesh::TetrahedronKey other_tet = (dsc.get_label(tets[0]) == phase)? tets[0] : tets[1];
            
            auto norm = dsc.get_normal(f.key(), other_tet);
            
            // Draw triangle
            glNormal3dv(norm.get());
            for (int i =0; i < 3; i++)
            {
                auto v = pts[i];
                auto n = nodes[i];
                vec3 real_norm = norm;
//                if((int)n < node_normal_vector.size() && node_normal_vector[n].length() > 0.1)
//                    real_norm = node_normal_vector[n];
                
//                glNormal3dv(real_norm.get());
                glVertex3dv(v.get());
            }
        }
    }
    
    glEnd();
}

void draw_helper::draw_cross(dsc_class & dsc, int nb_phase, vec3 domain_size)
{
    vector<vector<double>> colors = {{1,0,0}, {0,1,0}, {0,0,1}, {1,1,0}, {1,0,1}};
    // draw all phase near boundary
    double y_lim = domain_size[1] / 2.1;
    double y_lim_2 = y_lim - dsc.get_avg_edge_length()*1.5;
    double epsilon = dsc.get_avg_edge_length();
    

    for (auto fit = dsc.faces_begin(); fit != dsc.faces_end(); fit++)
    {
        auto node_pos = dsc.get_pos(dsc.get_nodes(fit.key()));
        auto mid_point = (node_pos[0] + node_pos[1] + node_pos[2])/3;
        auto tets = dsc.get_tets(fit.key());
        
        bool should_draw = false;
        bool between = false;
        if (!fit->is_boundary() &&
            mid_point[1] < y_lim)
        {
            for(int i = 0; i < 3; i++)
            {
                if (mid_point[i] < epsilon
                    || domain_size[i] - mid_point[i]  < epsilon)
                {
                    should_draw = true;
                }
            }
            if (mid_point[1] > y_lim_2)
            {
                should_draw = true;
                between = true;
            }
        }
        if (should_draw)
        {
            auto norm = dsc.get_normal(fit.key());
            if(between)
                norm = norm*(Util::dot(norm, vec3(0,1,0)) > 0? 1 : -1);
            else
                norm = -norm;
            
            int l = std::min(dsc.get_label(tets[0]), dsc.get_label(tets[1]));
            
            glBegin(GL_TRIANGLES);
            glEnable(GL_LIGHTING);
            if(l != BOUND_LABEL)
            {
                glColor3dv(&colors[l][0]);
                glNormal3dv(norm.get());
                for(auto v : node_pos)
                    glVertex3dv(v.get());
            }
            glEnd();
            
            if (between)
            {
                glBegin(GL_LINES);
                glDisable(GL_LIGHTING);
                glColor3f(1, 1, 1);
                for (int i = 0; i < 3; i++)
                {
                    glVertex3dv(node_pos[i].get());
                    glVertex3dv(node_pos[((i+1)%3)].get());
                }
                glEnd();
            }
        }
    }

}

void draw_helper::draw_tet_cross(dsc_class & dsc, int nb_phase, vec3 domain_size)
{
    vector<vector<double>> colors = { {0.5, 0.3, 1.0}, {1, 0.7, 0.7}, {0.7, 0.7, 0.7}, {0.7,0.8,1.0},  {1,1,0}, {1,0,1}};
    // draw all phase near boundary
    double y_lim = domain_size[1] / 3;
    double y_lim_2 = y_lim + dsc.get_avg_edge_length()*2;
    vec3 eye(0,0,1);
    
    for(auto tit = dsc.tetrahedra_begin(); tit != dsc.tetrahedra_end(); tit++)
    {
        if (tit->label() == 0)
        {
            continue;
        }
        
        auto node_pos = dsc.get_pos(dsc.get_nodes(tit.key()));
        bool b_draw = false;
        for(auto & p : node_pos)
            if (p[2] > y_lim && p[2] < y_lim_2)
            {
                b_draw = true;
            }
        
        if (b_draw)
        {
            glEnable(GL_LIGHTING);
            glColor3dv(colors[tit->label()].data());
            glBegin(GL_TRIANGLES);
            for(auto f : dsc.get_faces(tit.key()))
            {
                auto cotets = dsc.get_tets(f);
                if (dsc.get_label(cotets[0]) == 0 || dsc.get_label(cotets[1]) == 0)
                {
                    continue;
                }
                
                auto norm = dsc.get_normal(f);
                norm = norm*Util::dot(norm, eye);
                glNormal3dv(norm.get());
                for(auto v : dsc.get_pos(dsc.get_nodes(f)))
                    glVertex3dv(v.get());
            }
            glEnd();
            
            glDisable(GL_LIGHTING);
            glBegin(GL_LINES);
            glLineWidth(2.0);
            glColor3dv(colors[0].data());
            for (int i = 0; i < 4; i++)
            {
                glVertex3dv(node_pos[i].get());
                glVertex3dv(node_pos[(i+1)%4].get());
            }
            glEnd();
        }
            
    }
    
    for (auto fit = dsc.faces_begin(); fit != dsc.faces_end(); fit++)
    {
        if (!fit->is_interface())
        {
            continue;
        }
        
        auto cotets = dsc.get_tets(fit.key());
        int labels[2] = {dsc.get_label(cotets[0]), dsc.get_label(cotets[1])};
        if (labels[0] != 0 && labels[1] != 0)
        {
            continue;
        }

        auto pts = dsc.get_pos(dsc.get_nodes(fit.key()));
        bool b_draw = false;
        for(auto & p : pts)
            if ( p[2] < y_lim_2)
            {
                b_draw = true;
            }
        
        if (b_draw)
        {
            int l = labels[0] == 0? labels[1] : labels[0];
            glColor3dv(colors[l].data());
            
            glEnable(GL_LIGHTING);
            glBegin(GL_TRIANGLES);
            auto norm = dsc.get_normal(fit.key());
            glNormal3dv(norm.get());
            for(auto p : pts)
                glVertex3dv(p.get());
            glEnd();
            
            glDisable(GL_LIGHTING);
            glBegin(GL_LINES);
            glColor3dv(colors[0].data());
            for (int i = 0; i < 3; i++)
            {
                glVertex3dv(pts[i].get());
                glVertex3dv(pts[(i+1)%3].get());
            }
            glEnd();
        }
    }
}


void draw_helper::draw_curvature(dsc_class & dsc, std::vector<std::vector<vec3>> const & mean_curvature_hats, int phase, std::vector<std::vector<int>> const & correspond_label)
{
    if(mean_curvature_hats.size() == 0)
        return;
    
    // Find the vertices color base on curvature
    static std::vector<double> vertices_color; //(dsc.get_no_nodes_buffer(), 0.0);
    double max = 0.0;

    static bool bool_t_inited = false;
//    if (!bool_t_inited)
    {
        bool_t_inited = true;
     
        vertices_color = std::vector<double>(dsc.get_no_nodes_buffer(), 0.0);
        for (auto f = dsc.faces_begin(); f != dsc.faces_end(); f++)
        {
            if (f->is_interface() && !f->is_boundary())
            {
                auto tets = dsc.get_tets(f.key());
                if (!(dsc.get_label(tets[0]) == phase
                      or dsc.get_label(tets[1]) == phase)
                    )
                {
                    continue;
                }
            }
            
            for (auto & node : dsc.get_nodes(f.key()))
            {
                if (std::abs(vertices_color[node]) > 0)
                {
                    continue;
                }
                
                auto mean_curvature = mean_curvature_hats[node];
                auto corres_label = correspond_label[node];
                
                auto find_pos = std::find(corres_label.begin(), corres_label.end(), phase);
                
                if (find_pos == corres_label.end())
                {
                    // error
                    vertices_color[node] = -1;
                }else{
                    vertices_color[node] = mean_curvature[find_pos - corres_label.begin()].length();
                }
                
                max = std::max(max, vertices_color[node]);
            }
        }
        // Normalize the curvature???
        for (auto & cc : vertices_color)
        {
            cc /= max;
        }
    }
    

    
    
    for (auto f = dsc.faces_begin(); f != dsc.faces_end(); f++)
    {
        if (f->is_interface() && !f->is_boundary())
        {
            auto tets = dsc.get_tets(f.key());
            if (!(dsc.get_label(tets[0]) == phase
                  or dsc.get_label(tets[1]) == phase)
                //                || dsc.get_label(tets[0]) == BOUND_LABEL
                //                || dsc.get_label(tets[1]) == BOUND_LABEL
                )
            {
                continue;
            }
#ifdef DSC_CACHE
            auto nodes=*dsc.get_nodes_cache(f.key());
#else
            auto nodes=dsc.get_nodes(f.key());
#endif
            auto pts = dsc.get_pos(nodes);
            
            // normalize the normal to the eye
            is_mesh::TetrahedronKey other_tet = (dsc.get_label(tets[0]) == phase)? tets[0] : tets[1];
            
            
            auto norm = dsc.get_normal(f.key(), other_tet);
            
            
            
            
            // Draw triangle
            glBegin(GL_TRIANGLES);
            for (int i =0; i < 3; i++)
            {
                auto v = pts[i];
                auto n = nodes[i];
                vec3 real_norm = norm;
                if((int)n < node_normal_vector.size() && node_normal_vector[n].length() > 0.1)
                    real_norm = node_normal_vector[n];
                
                glColor3f(vertices_color[n], 0, 1 - vertices_color[n]);
                
                glNormal3dv(real_norm.get());
                glVertex3dv(v.get());
            }
            glEnd();
            
            
        }
    }
}

void draw_helper::draw_transparent_surface(dsc_class & dsc, int nb_phase)
{
    std::vector<vec3> color = {vec3(0,1,1), vec3(1,0,0), vec3(0,1,0), vec3(0,0,1), vec3(1,1,0), vec3(1,0,1)};
 
//    glEnable(GL_BLEND);
//    glBlendFunc (GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
//    glDisable(GL_LIGHTING);
    for(int i = 1; i<= nb_phase; i++)
    {
        auto c = color[i];
        glColor4f(c[0], c[1], c[2], 0.3);
        dsc_draw_one_interface_no_share(dsc, i);
    }
    
//    glColor3f(1.0, 0.0, 0.0);
//    dsc_draw_one_interface(dsc, 1);
//    
//    glColor4f(0.0, 0.0, 1.0, 0.4);
//    dsc_draw_one_interface(dsc, 2);
    
    glDisable(GL_BLEND);
}

void draw_helper::draw_triple_interface(dsc_class &dsc)
{
    std::vector<double> is_triple_edges(dsc.get_no_edges_buffer(), 0);
    for (auto fit = dsc.faces_begin(); fit != dsc.faces_end(); fit++)
    {
        if(fit->is_interface())
        {
            auto tets = dsc.get_tets(fit.key());
            double w = 1;
            if(dsc.get_label(tets[0]) == BOUND_LABEL || dsc.get_label(tets[1]) == BOUND_LABEL)
            {
                w = 0.4;
            }
            
            for(auto e : dsc.get_edges(fit.key()))
            {
                is_triple_edges[e] += w;
            }
        }
        
    }
    
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    glColor3f(1, 0, 0);
    for (auto eit = dsc.edges_begin(); eit != dsc.edges_end(); eit++)
    {
        if (is_triple_edges[eit.key()] > 2)
        {
            auto pts = dsc.get_pos(dsc.get_nodes(eit.key()));
            glVertex3dv(pts[0].get());
            glVertex3dv(pts[1].get());
        }
    }
    glEnd();
}

void draw_helper::dsc_draw_interface(dsc_class & dsc)
{
    for (auto f = dsc.faces_begin(); f != dsc.faces_end(); f++)
    {
        if (f->is_interface() && !f->is_boundary())
        {

            auto pts = dsc.get_pos(dsc.get_nodes(f.key()));
            //auto norm = Util::normal_direction(pts[0], pts[1], pts[2]);
            auto norm = -dsc.get_normal(f.key());
            
            glColor3f(0.7, 0.0, 0);
            glBegin(GL_TRIANGLES);
            for (auto v : pts)
            {
                glNormal3dv(norm.get());
                glVertex3dv(v.get());
            }
            glEnd();
            
//            
//            glDisable(GL_LIGHTING);
//            glColor3f(0, 0, 0);
//            glBegin(GL_LINES);
//            
//            auto edges = dsc.get_edges(f.key());
//            
//            for (int i = 0; i < 3; i++)
//            {
//                
//             //   glNormal3dv(norm.get());
//                glVertex3dv(pts[i].get());
//                
//             //   glNormal3dv(norm.get());
//                glVertex3dv(pts[(i+1)%3].get());
//            }
//            glEnd();
//            glEnable(GL_LIGHTING);
            
        }
    }
}

#define P_NONE  0x0000
#define P_ZERO  0x0001
#define P_ONE   0x0010
#define P_TWO   0x0100

#define P_ALL   0x0111

void draw_helper::dsc_draw_triple_edge(dsc_class & dsc)
{
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    for (auto eid = dsc.edges_begin(); eid != dsc.edges_end(); eid++)
    {
        if (eid->is_crossing())
        {
            auto pts = dsc.get_pos(dsc.get_nodes(eid.key()));
            glVertex3dv(pts[0].get());
            glVertex3dv(pts[1].get());
        }
    }
    glEnd();
}

void draw_helper::dsc_draw_domain(dsc_class & dsc)
{
    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT);
    glBegin(GL_TRIANGLES);
    for (auto fit = dsc.faces_begin(); fit != dsc.faces_end(); fit++)
    {
        if (fit->is_boundary())
        {
            auto verts = dsc.get_pos(dsc.get_sorted_nodes(fit.key()));
            vec3 normal = Util::normal_direction(verts[0], verts[1], verts[2]);
            for (auto v : verts)
            {
                glNormal3dv(normal.get());
                glVertex3dv(v.get());
            }
            
        }
    }
    
    
    glEnd();
    glDisable(GL_CULL_FACE);
}

void draw_helper::update_texture(const image3d & im,
                    int const &off_x, int const & off_y, int const & off_z)
{
    
    
    auto dim = im.dimension_v();
    
    get_instance()._cur_cross_poss += CGLA::Vec3i(off_x, off_y, off_z);
    for (int i = 0; i < 3; i++)
    {
        if (get_instance()._cur_cross_poss[i] < 0)
        {
            get_instance()._cur_cross_poss[i] = 0;
        }
        if (get_instance()._cur_cross_poss[i] > dim[i] - 1)
        {
            get_instance()._cur_cross_poss[i] = dim[i] - 1;
        }
    }

    int z = get_instance()._cur_cross_poss[2];
    int width = dim[0];
    int height = dim[1];
    uint8_t * data = (uint8_t *)malloc(dim[0] * dim[1] * 3 * sizeof(uint8_t));
    
    uint8_t *ptr = data;
    for (int j = 0; j < dim[1]; j++)
    {
        for (int i = 0; i < dim[0]; i++)
        {
            uint8_t v = (uint8_t)( im.get_value_f(i,j,z) * 255 );
            *(ptr++) = v;
            *(ptr++) = v;
            *(ptr++) = v;
        }
    }
    
    
    
    static GLuint tex_ID = 0;
    if (tex_ID != 0)
    {
        glDeleteTextures(1, &tex_ID);
        tex_ID = 0;
    }

    
    glGenTextures(1, &tex_ID);
    
    // "Bind" the newly created texture : all future texture functions will modify this texture
    glBindTexture(GL_TEXTURE_2D, tex_ID);
    
    glPixelStorei ( GL_UNPACK_ALIGNMENT,   1 );
    // Give the image to OpenGL
    glTexImage2D(GL_TEXTURE_2D, 0,GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
    
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    
    glBindTexture(GL_TEXTURE_2D, tex_ID);
    
    delete [] data;
}
