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
    int degree = 512; 
    int scale = 10;
    int64_t ctxt_modulus = 7919;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<uint64_t> msg_space(0, 30); // Use uint64_t for distribution range
    uniform_int_distribution<uint64_t> error_distribution(0, 3); // Use uint64_t for distribution range

    vector<poly> sk(rank), a(rank);

    cout << "Encrypting...\n";
    poly b(degree), tmp;
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
        b[i] += msg * scale + error_distribution(gen);
        b[i] = (b[i] % ctxt_modulus + ctxt_modulus) % ctxt_modulus;
        cout << msg << "\t";
    }
    
    poly dec_msg = b;
    cout << "\nDecrypting...\n";
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
        dec_msg[i] /= scale;
    }

    print(dec_msg);

    return 0;
}