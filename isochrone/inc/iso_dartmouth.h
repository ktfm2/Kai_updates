#ifndef ISODART_H
#define ISODART_H
//=============================================================================
#include "iso_base.h"
//=============================================================================
// Data for Dartmouth isochrones
namespace Dartmouth{
    const VecDoub ZList = {0.0000547,0.000172,0.000547,0.00172,0.00537,0.01885,0.02558,0.0348,0.05207};
    const VecDoub FeHList = {-2.5,-2.,-1.5,-1.,-0.5,0.07,0.21,0.36,0.56};
    const VecDoub ageList = {0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.,1.25,1.5,1.75,2.,2.25,2.5,2.75,3.,3.25,3.5,3.75,4.,4.25,4.5,4.75,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.,10.5,11.,11.5,12.,12.5,13.,13.5};
    const VecDoub ageList_1 = {0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};
    const VecDoub ageList_2 = {1.25,1.5,1.75,2.,2.25,2.5,2.75,3.,3.25,3.5,3.75,4.,4.25,4.5,4.75,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.,10.5,11.,11.5,12.,12.5,13.,13.5};
}
//=============================================================================
/**
 * @brief Class to store a single Dartmouth isochrone
 */
 class isochrone_dartmouth: public isochrone{
    public:
        isochrone_dartmouth(void);
        double get_metallicity(std::vector<std::string> input_iso, std::string dir);
        void fill(std::vector<std::string> input_iso, std::string dir, double Age=0., double thin_mag=-1.);
};

//=============================================================================
#endif
//=============================================================================
