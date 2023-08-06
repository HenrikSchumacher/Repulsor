#pragma once

#define MATHEMATICA

#include "mathlink.h"
#include <string>
#include <cstdint>
#include <ostream>
#include <sstream>

namespace mma
{
    WolframLibraryData libData;
    
    inline void print(const char *msg)
    {
        if (libData->AbortQ())
        {
            return; // trying to use the MathLink connection during an abort appears to break it
        }

        MLINK link = libData->getMathLink(libData);
        
        MLPutFunction(link, "EvaluatePacket", 1);
        
        MLPutFunction(link, "Print", 1);
        
        MLPutString(link, msg);
        
        libData->processMathLink(link);
        
        int pkt = MLNextPacket(link);
        
        if (pkt == RETURNPKT)
        {
            MLNewPacket(link);
        }
    }

    /// Call _Mathematica_'s `Print[]`, `std::string` argument version.
    inline void print(const std::string &msg)
    {
        print(msg.c_str());
    }
    
}


extern "C" DLLEXPORT mint WolframLibrary_getVersion()
{
    return WolframLibraryVersion;
}

extern "C" DLLEXPORT int WolframLibrary_initialize(WolframLibraryData libData)
{
    mma::libData = libData;
    return LIBRARY_NO_ERROR;
}

extern "C" DLLEXPORT void WolframLibrary_uninitialize(WolframLibraryData libData)
{
    return;
}
