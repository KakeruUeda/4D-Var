#ifndef PREPROCESSDIRECT_H
#define PREPROCESSDIRECT_H

#include <iostream>
#include "Config.h"

namespace Direct
{
    class Preprocess
    {
        public:
            Preprocess(std::string inputFile, std::string appName);
            ~Preprocess(){}
    
            Config conf;
    };
}


#endif
