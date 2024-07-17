#ifndef __PARAMETER__
#define __PARAMETER__

#include <cstdint>
#include "util.h"
#include "database.h"

class mlwe_parameter
{
private:
    int degree;
    int rank;
    int numInstance;

    uint64_t ptxt_modulus;
    uint64_t ctxt_modulus;
    ringMatrix crs;

public:
    mlwe_parameter();
    mlwe_parameter(const database& db, const int degree, const int rank);

    int getDegree() const {return degree;}
    int getRank() const {return rank;}
    int getNumInstance() const {return numInstance;}

    void setCRS(ringMatrix crs) {this -> crs = crs;}
    ringMatrix getCRS() const {return crs;}

    uint64_t getPtxtModulus() const {return ptxt_modulus;}
    uint64_t getCtxtModulus() const {return ctxt_modulus;}

    void print();
};


class pir_parameter
{
private:
    // uint64_t ctxt_modulus;
    matrix crs;

public:
    pir_parameter()
    {
        // ctxt_modulus = 0;
        crs.clear();
    }
    pir_parameter(database& db);  

    // int getCtxtModulus() const {return ctxt_modulus;}
    void setCRS(matrix crs) {this -> crs = crs;}
    matrix getCRS() const {return crs;}

    void print();
};


#endif