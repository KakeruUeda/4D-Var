#include "Config.h"

Config::Config(std::string inputFile)
{
    tryOpenConfigFile(inputFile);
    tryReadConfigFile();
}

void Config::tryOpenConfigFile(std::string inputFile)
{
    try
    {
        int error;
        if ((error = tp.read(inputFile)) != TP_NO_ERROR)
            throw std::runtime_error("Open error");
    }
    catch (const std::runtime_error& e)
    {
		std::cout << e.what() << std::endl;
		isReadingError = true;
    }
    
    return;
}

void Config::tryReadConfigFile()
{
    try
    {
        readConfigFile();
    }
    catch(const std::runtime_error& e)
    {
        std::cout << e.what() << std::endl;
        isReadingError = true;
    }

    return;
}

void Config::readConfigFile()
{
    readGridType();
    if(gridType == GridType::STRUCTURED)
    {
        readBase(); readPysicalParam();
        readGrid();
    }
    else if(gridType == GridType::UNSTRUCTURED) 
    {
    }

    return;
}
