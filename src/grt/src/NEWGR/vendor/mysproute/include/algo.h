#ifndef ALGO_H
#define ALGO_H
#include <iostream>
#include <string>

enum Algo {
    FineGrain = 0,
    NonDet,
    DetPart,
    Astar,
    RUDY,
    DetPart_Astar,
    PIN_DENSITY,
    DetPart_Astar_Data,
    DetPart_Astar_RUDY,
    DRC_MAP,
    DetPart_Astar_Local,
    DetPart_Astar_Small,
    Region,
    DetPart_Astar_Region,
    DetPart_Astar_Regular_Region,
    WRITE_GUIDE_POINT,
    OUT_OF_GUIDE,
};

Algo StrToAlgo(string algo_str) {
    if(algo_str == "FineGrain")
        return FineGrain;
    else if(algo_str == "NonDet") 
        return NonDet;
    else if(algo_str == "DetPart") 
        return DetPart;
    else if(algo_str == "Astar") 
        return Astar;
    else if(algo_str == "RUDY") 
        return RUDY;
    else if(algo_str == "DetPart_Astar") 
        return DetPart_Astar;
    else if(algo_str == "PIN_DENSITY") 
        return PIN_DENSITY;
    else if(algo_str == "DetPart_Astar_Data")
        return DetPart_Astar_Data;
    else if(algo_str == "DetPart_Astar_RUDY")
        return DetPart_Astar_RUDY;
    else if(algo_str == "DRC_MAP")
        return DRC_MAP;
    else if(algo_str == "DetPart_Astar_Local")
        return DetPart_Astar_Local;
    else if(algo_str == "DetPart_Astar_Small")
        return DetPart_Astar_Small;
    else if(algo_str == "Region")
        return Region;
    else if(algo_str == "DetPart_Astar_Region")
        return DetPart_Astar_Region;
    else if(algo_str == "DetPart_Astar_Regular_Region")
        return DetPart_Astar_Regular_Region;
    else if(algo_str == "WRITE_GUIDE_POINT")
        return WRITE_GUIDE_POINT;
    else if(algo_str == "OUT_OF_GUIDE")
        return OUT_OF_GUIDE;
    else {
        std::cout << "unknown algorithm: " << algo_str << std::endl;
        exit(1);
    }

}

#endif