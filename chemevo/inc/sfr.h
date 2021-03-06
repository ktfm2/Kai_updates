#ifndef SFR_H
#define SFR_H
//=============================================================================
#include <map>
#include <string>
//=============================================================================
#include "params.h"
#include "GSLInterface/GSLInterface.h"
#include "utils.h"
#include "in_out.h"
#include "grid.h" //Kai - to use for creating SFR grid
//=============================================================================
/**
 * @brief Star Formation rate base class
 */
class StarFormationRate{
private:
protected:
    double Rmin, Rmax; // Minimum and maximum radius of star formation
    double GalaxyAge;
    double presentSFR;
    double assign_sfr_value; //Kai - the value needed for the SFR
    unsigned sfr_grid_radial_point; //Kai - the radial grid point for SFR
    unsigned sfr_grid_time_point; //Kai - the time grip point for SFR
    std::unique_ptr<Grid> sfr_grid; //Kai - for SFR grid
public:
    //StarFormationRate(double Rmin, double Rmax, double Tmax)
    //    :Rmin(Rmin),Rmax(Rmax),GalaxyAge(Tmax),presentSFR(1.){}
    //StarFormationRate(ModelParameters M);
    StarFormationRate(double Rmin, double Rmax, double Tmax, double asfrv, unsigned sfrgrp, unsigned sfrgtp)
        :Rmin(Rmin),Rmax(Rmax),GalaxyAge(Tmax), assign_sfr_value(asfrv), sfr_grid_radial_point(sfrgrp), sfr_grid_time_point(sfrgtp), presentSFR(1.){}
    StarFormationRate(ModelParameters M, double asfrv, unsigned sfrgrp, unsigned sfrgtp);
//    StarFormationRate(ModelParameters M);
    double MaxAge(void){return GalaxyAge;}
    /**
     * @brief Star formation rate
     *
     * @param R radius
     * @param t time
     *
     * @return star formation rate at radius R and time t.
     */
    virtual double operator()(double R, double t)=0;
    /**
     * @brief gas consumed at radius R
     * @details integrated star formation rate at radius R
     *
     * @param R radius
     * @return gas consumed at radius R
     */
    double gas_consumed_per_unit_radius(double R);
    /**
     * @brief gas consumed over all radii and time
     * @details integrated SFR over radius and time
     * @return integrated SFR over radius and time
     */
    double gas_consumed(void);
    //Sets the SFR grid 
    virtual void set_sfr_grid(double asfrv, unsigned sfrgrp, unsigned sfrgtp){std::cout<<"Test check works. Not important what is here.";};
    //Expands the SFR grid
    virtual void sfr_add_time(unsigned nt, double t){std::cout<<"Test to check works. Not really important what is here.";};
};

//=============================================================================
/**
 * @brief Star formation rate used in Sanders & Binney (2015)
 */
class SFR_SB15:public StarFormationRate{
private:
    double t_d; // Long time-scale for decline in SFR
    double t_s; // Short time-scale for initial rise in SFR
    double Rd;  // Scale-length of population of stars formed.
    double Rb;  // Radius at which exponential profile is truncated.
    double KSA; // Kennicutt-Schmidt A
    double KSN; // Kennicutt-Schmidt N
public:
    //SFR_SB15(double t_d=8.,double t_s=0.5,double Rd=4.,
    //         double Rmin=0.3, double Rmax=20., double Tmax=12.,double Rb=1000.)
	//   :StarFormationRate(Rmin,Rmax,Tmax),t_d(t_d),t_s(t_s),Rd(Rd),Rb(Rb){}
    SFR_SB15(double t_d=8.,double t_s=0.5,double Rd=4., double asfrv=0., unsigned sfrgrp=0u, unsigned sfrgtp=0u,
             double Rmin=0.3, double Rmax=20., double Tmax=12.,double Rb=1000.)
	   :StarFormationRate(Rmin,Rmax,Tmax,asfrv,sfrgrp,sfrgtp),t_d(t_d),t_s(t_s),Rd(Rd),Rb(Rb){}
    SFR_SB15(ModelParameters M, double asfrv=0., unsigned sfrgrp=0u, unsigned sfrgtp=0u,double t_d=8.,double t_s=0.5,double Rdx=4.,double Rbx=1000.)
        :StarFormationRate(M, asfrv, sfrgrp, sfrgtp),
	 t_d(t_d),
	 t_s(t_s),
	 Rd(Rdx),
     Rb(Rbx){
        Rd=extract_param(M.parameters["fundamentals"],
                         "StarScaleLength", Rdx);
        Rb=extract_param(M.parameters["fundamentals"],
                         "TruncationRadius", Rbx);
        presentSFR=1.;
        presentSFR=(double)M.parameters["fundamentals"]["PresentSFR"]/(
                (*this)(M.parameters["fundamentals"]["SolarRadius"],
                        M.parameters["fundamentals"]["GalaxyAge"]));
        KSA=extract_param(M.parameters["fundamentals"],
                         "Kennicutt-Schmidt_A", 0.067);
        KSN=extract_param(M.parameters["fundamentals"],
                         "Kennicutt-Schmidt_Coeff", 1.4);
	sfr_grid = make_unique<Grid>(M);
        }
    void set_sfr_grid(double asfrv, unsigned sfrgrp, unsigned sfrgtp);
    void sfr_add_time(unsigned nt, double t);
    double operator()(double R, double t);
};
class SFR_ExpDecay:public StarFormationRate{
private:
    double t_d; // Long time-scale for decline in SFR
    double Rd;  // Scale-length of population of stars formed.
    double Rb;  // Radius at which exponential profile is truncated.
    double KSA; // Kennicutt-Schmidt A
    double KSN; // Kennicutt-Schmidt N
public:
    //SFR_ExpDecay(double t_d=8.,double Rd=4.,
      //       double Rmin=0.3, double Rmax=20., double Tmax=12.,double Rb=1000.)
      // :StarFormationRate(Rmin,Rmax,Tmax),t_d(t_d),Rd(Rd){}
    SFR_ExpDecay(double t_d=8.,double Rd=4., double asfrv=0., unsigned sfrgrp=0u, unsigned sfrgtp=0u,
             double Rmin=0.3, double Rmax=20., double Tmax=12.,double Rb=1000.)
       :StarFormationRate(Rmin,Rmax,Tmax,asfrv,sfrgrp,sfrgtp),t_d(t_d),Rd(Rd){}
    SFR_ExpDecay(ModelParameters M, double asfrv=0., unsigned sfrgrp=0u, unsigned sfrgtp=0u, double t_dx=8.,double Rdx=4.,double Rbx=1000.)
        :StarFormationRate(M, asfrv, sfrgrp, sfrgtp),
     t_d(t_dx),
     Rd(Rdx),
     Rb(Rbx){
        Rd=extract_param(M.parameters["fundamentals"],
                         "StarScaleLength", Rdx);
        Rb=extract_param(M.parameters["fundamentals"],
                         "TruncationRadius", Rbx);
        t_d=extract_param(M.parameters["fundamentals"],
                         "SFR_decay_scale", t_dx);
        presentSFR=1.;
        presentSFR=(double)M.parameters["fundamentals"]["PresentSFR"]/(
                (*this)(M.parameters["fundamentals"]["SolarRadius"],
                        M.parameters["fundamentals"]["GalaxyAge"]));
        KSA=extract_param(M.parameters["fundamentals"],
                         "Kennicutt-Schmidt_A", 0.067);
        KSN=extract_param(M.parameters["fundamentals"],
                         "Kennicutt-Schmidt_Coeff", 1.4);
	sfr_grid = make_unique<Grid>(M);
        }
    void set_sfr_grid(double asfrv, unsigned sfrgrp, unsigned sfrgtp);
    void sfr_add_time(unsigned nt, double t);
    double operator()(double R, double t);
};

class SFR_DoubleInfall:public StarFormationRate{
private:
    double t_d; // Long time-scale for decline in SFR
    double Rd;  // Scale-length of population of stars formed.
    double Rb;  // Radius at which exponential profile is truncated.
    double KSA; // Kennicutt-Schmidt A
    double KSN; // Kennicutt-Schmidt N
    double coeff; 
    double acoeff; //Kai 
//    double t_d_2; 
//    double a_t_d_2; //Kai 
    std::unique_ptr<GasDump> gasdumpflow;
    std::unique_ptr<AlternateGasDump> alternategasdumpflow; //Kai
    double original_SFR_scaling;
public:
    //SFR_DoubleInfall(double t_d=8.,double Rd=4.,double Rmin=0.3, double Rmax=20., 
      //       double Tmax=12., double Rb=1000.)
      // :StarFormationRate(Rmin,Rmax,Tmax),t_d(t_d),Rd(Rd){}
    SFR_DoubleInfall(double t_d=8.,double Rd=4.,double Rmin=0.3, double Rmax=20., 
             double Tmax=12., double Rb=1000, double asfrv=0., unsigned sfrgrp=0u, unsigned sfrgtp=0u)
       :StarFormationRate(Rmin,Rmax,Tmax,asfrv,sfrgrp,sfrgtp),t_d(t_d),Rd(Rd){}
    SFR_DoubleInfall(ModelParameters M, double asfrv=0., unsigned sfrgrp=0u, unsigned sfrgtp=0u, double t_dx=8.,double Rdx=4.,double Rbx=1000.)
        :StarFormationRate(M, asfrv, sfrgrp, sfrgtp), t_d(t_dx), Rd(Rdx), Rb(Rbx)
    {
        Rd=extract_param(M.parameters["fundamentals"],
                         "StarScaleLength", Rdx);
        Rb=extract_param(M.parameters["fundamentals"],
                         "TruncationRadius", Rbx);
        t_d=extract_param(M.parameters["fundamentals"],
                         "SFR_decay_scale", t_dx);
        KSA=extract_param(M.parameters["fundamentals"],
                         "Kennicutt-Schmidt_A", 0.067);
        KSN=extract_param(M.parameters["fundamentals"],
                         "Kennicutt-Schmidt_Coeff", 1.4);
//        t_d_2=extract_param(M.parameters["flows"]["gasdump"],
//                         "star_formation_boost_timescale", 1.);
//        a_t_d_2=extract_param(M.parameters["flows"]["alternategasdump"], //Kai
//                         "star_formation_boost_timescale", 1.);
    
    std::shared_ptr<SolarAbundances> S = solar_types["Asplund"](M);
    auto F = M.parameters["flows"];
    gasdumpflow = gasdump_types[M.parameters["flows"]["gasdump"]["Form"]](M,S);
    double gS = (*gasdumpflow)(M.parameters["fundamentals"]["SolarRadius"],
                               gasdumpflow->dump_time(),1.);
    alternategasdumpflow = alternategasdump_types[M.parameters["flows"]["alternategasdump"]["Form"]](M,S);
    double agS = (*alternategasdumpflow)(M.parameters["fundamentals"]["SolarRadius"],
                               alternategasdumpflow->dump_time(),1.);
    coeff = 0.;
    acoeff = 0.;
    presentSFR = 1.;
    original_SFR_scaling=1.;
    original_SFR_scaling=(double)M.parameters["fundamentals"]["PresentSFR"]/(
                (*this)(M.parameters["fundamentals"]["SolarRadius"],
                        M.parameters["fundamentals"]["GalaxyAge"]));
 
    coeff=extract_param(M.parameters["flows"]["gasdump"],
                        "star_formation_boost_coeff", 1.);
    acoeff=extract_param(M.parameters["flows"]["alternategasdump"],
                        "star_formation_boost_coeff", 1.);
    
    presentSFR=(double)M.parameters["fundamentals"]["PresentSFR"]/(
                (*this)(M.parameters["fundamentals"]["SolarRadius"],
                        M.parameters["fundamentals"]["GalaxyAge"]));
    sfr_grid = make_unique<Grid>(M);
    }
    void set_sfr_grid(double asfrv, unsigned sfrgrp, unsigned sfrgtp);
    void sfr_add_time(unsigned nt, double t);
    double operator()(double R, double t);
};
class SFR_KS:public StarFormationRate{ //Kai - SFR working straight from KS law
private:
    double KSA; // Kennicutt-Schmidt A
    double KSN; // Kennicutt-Schmidt N
public:
    //SFR_KS(double Rmin=0.3, double Rmax=20., double Tmax=12.)
//	   :StarFormationRate(Rmin,Rmax,Tmax){}
    SFR_KS(double Rmin=0.3, double Rmax=20., double Tmax=12., double asfrv=0., unsigned sfrgrp=0u, unsigned sfrgtp=0u)
	   :StarFormationRate(Rmin,Rmax,Tmax,asfrv,sfrgrp,sfrgtp){}
    SFR_KS(ModelParameters M, double asfrv=0., unsigned sfrgrp=0u, unsigned sfrgtp=0u)
        :StarFormationRate(M, asfrv, sfrgrp, sfrgtp){
        KSA=extract_param(M.parameters["fundamentals"],
                         "Kennicutt-Schmidt_A", 0.067);
        KSN=extract_param(M.parameters["fundamentals"],
                         "Kennicutt-Schmidt_Coeff", 1.4);
	sfr_grid = make_unique<Grid>(M);
    }
    void set_sfr_grid(double asfrv, unsigned sfrgrp, unsigned sfrgtp);
    void sfr_add_time(unsigned nt, double t);
    double operator()(double R, double t);
};

class SFR_Neige2020: public StarFormationRate{
private:
    double t_sfr;
    double x_io;
    double t_m;
    double t_to;
    double Rd;  // Scale-length of population of stars formed.
    double Rb;  // Radius at which exponential profile is truncated.
    double R0;  // Solar radius.
    double KSA; // Kennicutt-Schmidt A
    double KSN; // Kennicutt-Schmidt N
public:
    SFR_Neige2020(ModelParameters M, double asfrv=0., unsigned sfrgrp=0u, unsigned sfrgtp=0u): StarFormationRate(M, asfrv, sfrgrp, sfrgtp), t_m(6.){
        t_sfr = extract_param(M.parameters["fundamentals"],"SFR_decay_scale",1.02);
        x_io = extract_param(M.parameters["fundamentals"],"InsideOutScale",0.67);
        Rd = extract_param(M.parameters["fundamentals"],"StarScaleLength",2.85);
        Rb = extract_param(M.parameters["fundamentals"],"TruncationRadius",1000.);
        t_to = extract_param(M.parameters["fundamentals"],"SFR_ShortTimescale",1.);
	R0 = M.parameters["fundamentals"]["SolarRadius"];
        KSA=extract_param(M.parameters["fundamentals"],
                         "Kennicutt-Schmidt_A", 0.067);
        KSN=extract_param(M.parameters["fundamentals"],
                         "Kennicutt-Schmidt_Coeff", 1.4);
        presentSFR=1.;
        presentSFR=(double)M.parameters["fundamentals"]["PresentSFR"]/(
                (*this)(M.parameters["fundamentals"]["SolarRadius"],
                        M.parameters["fundamentals"]["GalaxyAge"]));
	sfr_grid = make_unique<Grid>(M);
    }
    void set_sfr_grid(double asfrv, unsigned sfrgrp, unsigned sfrgtp);
    void sfr_add_time(unsigned nt, double t);
    double operator()(double R, double t);
};

//=============================================================================
// Simple structure for integrating the SFR wrt time.
struct gas_con_st{
    StarFormationRate *SFR;
    double R;
};
//=============================================================================
// Map for creating new instances of SFR from a string giving the class name
// sfr_types[std::string s="<ClassName>"](ModelParameters M)
// produces a shared pointer of Class <ClassName> initialized with the model
// parameters M.
extern shared_map<StarFormationRate,ModelParameters> sfr_types;
//=============================================================================
#endif
//=============================================================================
