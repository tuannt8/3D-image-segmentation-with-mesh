//
//  InputParser.h
//  DEMO
//
//  Created by Tuan Nguyen Trung on 01/02/2018.
//

// Parser input arguments of the main function

#ifndef InputParser_h
#define InputParser_h

#include <vector>
#include <string>
#include <fstream>

class InputParser{
public:
    InputParser(){};
    InputParser(std::string path)
    {
        std::ifstream input;
        input.open(path);
        if (input.is_open())
        {
            std::string word;
            while(input >> word)
            {
                this->tokens.push_back(word);
            }
        }
        input.close();
    }
    
    InputParser (int &argc, char **argv){
        for (int i=1; i < argc; ++i)
            this->tokens.push_back(std::string(argv[i]));
    }
    /// @author iain
    const std::string& getCmdOption(const std::string &option){
        std::vector<std::string>::const_iterator itr;
        itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
        if (itr != this->tokens.end() && ++itr != this->tokens.end()){
            argument_used.push_back(option + " [" + *itr +  "]");
            return *itr;
        }
        static const std::string empty_string("");
        return empty_string;
    }
    
    std::string getCmdOption(const std::string &option, const std::string &default_){
        std::vector<std::string>::const_iterator itr;
        itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
        if (itr != this->tokens.end() && ++itr != this->tokens.end()){
            argument_used.push_back(option + " [" + *itr +  "]");
            return *itr;
        }
        argument_used.push_back(option + " [" + default_ +  "]");
        return default_;
    }
    /// @author iain
    bool cmdOptionExists(const std::string &option){
        argument_used.push_back(option);
        
        return std::find(this->tokens.begin(), this->tokens.end(), option)
        != this->tokens.end();
    }
    
    void print(){
        std::sort(argument_used.begin(), argument_used.end());
        argument_used.erase(std::unique(argument_used.begin(), argument_used.end()), argument_used.end());
        
        std::cout << "----------------Command lines available--------------\n" ;
        for(auto v : argument_used)
        {
            std::cout << v << "  ;  ";
        }
        std::cout << "\n--------------------------------------------\n";
    }
private:
    std::vector <std::string> tokens;
    std::vector <std::string> argument_used;
};

#endif /* InputParser_h */
