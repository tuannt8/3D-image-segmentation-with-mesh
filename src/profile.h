//
//  profile.hpp
//  DSC
//
//  Created by Tuan Nguyen Trung on 2/4/16.
//  Copyright © 2016 Asger Nyman Christiansen. All rights reserved.
//

#ifndef profile_hpp
#define profile_hpp

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <chrono>
#include <map>
#include <string>

#define P_TIME_NOW (std::chrono::system_clock::now())
typedef std::chrono::duration<double> p_duration_t;
typedef std::chrono::system_clock::time_point p_time_point;

struct profile_att
{
    int count = 0;
    double total_time = 0;
    p_time_point m_start;
};

class profile
{
public:
    static void init();
    static void close();
    
    profile(std::string name);
    
    void change(std::string name);
    void done();
    
    profile();
    ~profile();
private:

    static std::map<std::string, profile_att> m_objects;
    static profile_att * get_object(const std::string &  name );
private:
    
    std::string m_name;
};

#endif /* profile_hpp */
