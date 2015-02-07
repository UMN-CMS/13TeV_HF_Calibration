/* 
 * File:   Genconverter.h
 * Author: lesko
 *
 * Created on May 14, 2014, 4:50 PM
 */

#ifndef GENCONVERTER_H
#define	GENCONVERTER_H

class Genconverter {
public:
    Genconverter();
    Genconverter(const Genconverter& orig);
    virtual ~Genconverter();
    int Eta2IEta(double eta);
    int Phi2Iphi(double phi, double eta);
private:

};

#endif	/* GENCONVERTER_H */

