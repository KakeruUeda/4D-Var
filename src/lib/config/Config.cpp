#include "Config.h"

Config::Config(std::string inputFile, std::string appName)
{
    setApplication(appName);
    tryOpenConfigFile(inputFile);
    tryReadConfigFile();
}

void Config::setApplication(std::string appName)
{   
    if     (appName == "STRGRID")  app = Application::STRGRID;
    else if(appName == "SNS")      app = Application::SNS;
    else if(appName == "USNS")     app = Application::USNS;
    else if(appName == "TDVAR")    app = Application::TDVAR;
    else if(appName == "FDVAR")    app = Application::FDVAR;
    else
        if(mpi.myId == 0)
            std::cout << "Unknown appName" << std::endl;

    return;
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
		if(mpi.myId == 0) 
            std::cout << e.what() << std::endl;
		if(mpi.myId == 0) 
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
        if(mpi.myId == 0)
            std::cout << e.what() << std::endl;
        if(mpi.myId == 0)
            isReadingError = true;
    }

    return;
}

void Config::readConfigFile()
{
    switch(app)
    {
        case Application::STRGRID:
            readBasicParameter(); 
            readGridParameter(); 
            readBoundaryParameter();
            break;

        case Application::USNS:
            readBasicParameter();
            break;

        default:
            if(mpi.myId == 0)
                throw std::runtime_error("Unknown Application");
            break;
    }
    return;
}
