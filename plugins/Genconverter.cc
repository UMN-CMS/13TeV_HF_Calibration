/* 
 * File:   Genconverter.cc
 * Author: lesko
 * 
 * Created on May 14, 2014, 4:50 PM
 */

#include "Genconverter.h"
#include <math.h>
Genconverter::Genconverter() { }

Genconverter::Genconverter(const Genconverter& orig) { }
//don't think I need this

Genconverter::~Genconverter() { }

int Genconverter::Eta2IEta(double eta)
{
    int etasign = (int) (fabs(eta) / eta);
    eta = fabs(eta);
    int tower = 0;



    if(eta <= 2.964) tower = 0;
    else if(eta <= 3.139) tower = 30;
    else if(eta <= 3.314) tower = 31;
    else if(eta <= 3.489) tower = 32;
    else if(eta <= 3.664) tower = 33;
    else if(eta <= 3.839) tower = 34;
    else if(eta <= 4.013) tower = 35;
    else if(eta <= 4.191) tower = 36;
    else if(eta <= 4.363) tower = 37;
    else if(eta <= 4.538) tower = 38;
    else if(eta <= 4.716) tower = 39;
    else tower = 0;
    

    tower *= etasign;



    return tower;

}

int Genconverter::Phi2Iphi(double phi, double eta)
{
    int iphi;
    eta = fabs(eta);
    if(eta > 1.740 && eta < 4.716)
    {
        if(phi >= 0)
        {
            iphi = ((int) ((phi)*180 / M_PI) / 10)*2 + 1; //so this converts into degrees since everything is nicer there. followed by conversion it iphi
        } else
        {
            iphi = ((int) (phi*180 / M_PI + 360) / 10)*2 + 1;
        }
    }
    else if(eta > 4.716 && eta < 5.191)
    {
        if(phi >= 0)
        {
            iphi = ((int) ((phi)*180 / M_PI) / 20)*4 + 3; //so this converts into degrees since everything is nicer there. followed by conversion it iphi
        } else
        {
            iphi = ((int) (phi*180 / M_PI + 360) / 20)*4 + 3;
        }

    }
    else if(eta>0&&eta<1.653)
    {
        if(phi >= 0)
        {
            iphi = ((int) ((phi)*180 / M_PI) / 5) + 1; //so this converts into degrees since everything is nicer there. followed by conversion it iphi
        } else
        {
            iphi = ((int) (phi*180 / M_PI + 360) / 5) + 1;
        }
        
    }
    else iphi=0;
           
    return iphi;
}
