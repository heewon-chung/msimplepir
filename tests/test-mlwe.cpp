#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <random>

#include "util.h"

using namespace std;

int main()
{
    int rank = 3;
    int degree = 16; 
    int64_t ctxt_modulus = 7919;
    int64_t ptxt_modulus = 31;
    int scale = ctxt_modulus / ptxt_modulus;
    int mean = 0;
    double stddev = 6.4;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<uint64_t> msg_space(0, ptxt_modulus); // Use uint64_t for distribution range

    vector<poly> sk(rank), a(rank);

    
    poly b(degree), tmp, msg_poly;
    for(int i = 0; i < rank; i++)
    {
        randVector(sk[i], degree, ctxt_modulus);
        randVector(a[i],  degree, ctxt_modulus);

        multiply_ntt(a[i], sk[i], tmp, ctxt_modulus, 3);
        for(int j = 0; j < tmp.size(); j++)
        {
            b[j] += tmp[j];
            b[j] = (b[j] % ctxt_modulus + ctxt_modulus) % ctxt_modulus;
        }
    }
    // invert_ntt(b, b, ctxt_modulus, 3);

    for(int i = 0; i < degree; i++)
    {
        auto msg = msg_space(gen);
        msg_poly.push_back(msg);
        b[i] += msg_poly[i] * scale + generateDiscreteGaussian(mean, stddev);
        b[i] = (b[i] % ctxt_modulus + ctxt_modulus) % ctxt_modulus;
    }

    printf("\nEncrypting...\n");
    print(msg_poly);
    
    poly dec_msg = b;
    for(int i = 0; i < rank; i++)
    {
        multiply_ntt(a[i], sk[i], tmp, ctxt_modulus, 3);
        for(int j = 0; j < degree; j++)
        {
            dec_msg[j] -= tmp[j];
            dec_msg[j] = (dec_msg[j] % ctxt_modulus + ctxt_modulus) % ctxt_modulus;
        }
    }
    
    for(int i = 0; i < degree; i++)
    {
        // dec_msg[i] /= scale;
        dec_msg[i] = (2 * dec_msg[i] + scale) / (scale << 1);
    }

    printf("\nDecrypting...\n");
    print(dec_msg);

    for(int i = 0; i < degree; i++)
    {
        if(msg_poly[i] != dec_msg[i])
        {
            cerr << i << "th msg is wrong!" << endl;
        }
    }

    return 0;
}