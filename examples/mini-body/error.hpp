#pragma once

#include <exception>
#include <string>

class bad_arguments : public std::exception
{
    private:
        std::string message;
        
    public:
    bad_arguments(std::string mes): message(mes) {}
    
    const char* what()const noexcept {return message.c_str();}
    
};
