//
//  config.h
//  DSC_segment_probability
//
//  Created by Tuan Nguyen Trung on 14/12/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#ifndef config_h
#define config_h

#include <stdio.h>
#include <string>


class config
{
private:
    std::map<std::string, std::string> m_option;
    
public:
    static config * get_instance()
    {
        static config instance;
        return & instance;
    }
    
    static void set_arg(int num_arg, char *arg[])
    {
        auto instance = get_instance();
        for (int i = 1; i < num_arg; i++)
        {
            
        }
    }
    static void set_config_file_path()
    {
        
    }
    
    static std::string get_setting(std::string key)
    {
        
    }
}

#endif /* config_h */
