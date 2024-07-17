#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include "util.h"

using namespace std;

void negacyclic(int degree, int rank, int numInstance, int64_t ctxt_modulus, ringMatrix& mlwe_crs, matrix& crs)
{
    // ringMatrix mlwe_crs;
    randMatrix(mlwe_crs, numInstance, rank, degree, ctxt_modulus);

    cout << "\nInput Matrix: ";
    cout << mlwe_crs.size() << " * " << mlwe_crs[0].size() << endl;
#if __DEBUG == 1
    print(mlwe_crs);
#endif

    crs.resize(numInstance * degree);
    for(int i = 0; i < crs.size(); i++)
    {
        crs[i].resize(rank * degree);
    }

    cout << "Output Matrix: ";
    cout << crs.size() << " * " << crs[0].size() << endl;

    for(int i = 0; i < numInstance; i++)
    {
        // j: column
        int colIdx = 0;
        int rowIdx = i * degree;
        for(int j = 0; j < rank; j++)
        {
            // tmp \in R_q
            auto tmp = mlwe_crs[i][j];
            for(int k = 0; k < degree; k++)
            {
                // transpose
                for(int deg = 0; deg < degree; deg++)
                {
                    crs[rowIdx + deg][colIdx] = tmp[deg];
                }
                uint64_t lastElement = tmp.back();
                for(int deg = degree - 1; deg > 0; deg--)
                {
                    tmp[deg] = tmp[deg - 1];
                }
                tmp[0] = ctxt_modulus - lastElement;
                colIdx++;   
            }
        }
    }

#if __DEBUG == 1
    print(crs);
#endif

}


void instance_conversion(int degree, int rank, int numInstance, int64_t ctxt_modulus)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<uint64_t> msg_space(0, 10); // Use uint64_t for distribution range
    
    ringMatrix mlwe_crs;
    matrix crs;
    negacyclic(degree, rank, numInstance, ctxt_modulus, mlwe_crs, crs);

    vector<poly> sk(mlwe_crs[0].size());
    for(int i = 0; i < mlwe_crs[0].size(); i++)
    {
        randVector(sk[i], degree, ctxt_modulus);
    }

    printf("\nFor MLWE Instnace....\n");
    printf("    A: "); cout << mlwe_crs.size() << " * " << mlwe_crs[0].size() << endl;
#if __DEBUG == 1
    print(mlwe_crs);
#endif

    poly tmp, msg;
    vector<poly> ctxt(numInstance);
    for(int i = 0; i < numInstance; i++)
    {
        ctxt[i].resize(degree);
        for(int j = 0; j < rank; j++)
        {
            multiply_ntt(mlwe_crs[i][j], sk[j], tmp, ctxt_modulus, 3);
            for(int k = 0; k < degree; k++)
            {
                ctxt[i][k] += tmp[k];
            }
        }
        
        for(int j = 0; j < degree; j++)
        {
            msg.push_back(msg_space(gen));
            ctxt[i][j] += msg.back();
            ctxt[i][j] = (ctxt[i][j] % ctxt_modulus + ctxt_modulus) % ctxt_modulus;
        }
    }

    printf("    b: "); cout << ctxt.size() << " * " << ctxt[0].size() << endl;
    printf("    - Underlying message is...  "); print(msg);
#if __DEBUG == 1
    print(ctxt);
#endif

    printf("\nFor LWE Instance....\n");
    printf("    A: "); cout << crs.size() << " * " << crs[0].size() << endl;
    // print(crs);
    vector<int64_t> sk_lwe, ctxt_lwe;
    for(int i = 0; i < sk.size(); i++)
    {
        assert(sk[i].size() == degree);
        for(int j = 0; j < degree; j++)
        {
            sk_lwe.push_back(sk[i][j]);
        }
    }

    for(int i = 0; i < ctxt.size(); i++)
    {
        assert(ctxt[i].size() == degree);
        for(int j = 0; j < degree; j++)
        {
            ctxt_lwe.push_back(ctxt[i][j]);
        }
    }

    printf("    b: "); cout << ctxt_lwe.size() << endl;
#if __DEBUG == 1
    print(ctxt_lwe);
#endif

    vector<int64_t> dec_msg;
    matrixMultiply(crs, sk_lwe, ctxt_modulus, dec_msg);
    
    for(int i = 0; i < dec_msg.size(); i++)
    {
        dec_msg[i] = ctxt_lwe[i] - dec_msg[i];
        dec_msg[i] = (dec_msg[i] % ctxt_modulus + ctxt_modulus) % ctxt_modulus;
    }
    
    printf("    - Decrypting message is...  ");print(dec_msg);

    for(int i = 0; i < msg.size(); i++)
    {
        if(msg[i] != dec_msg[i])
            cerr << i << "th msg is different\n";
    }

}


int main()
{
    int degree = 512;
    int rank = 3;
    int numInstance = 5;
    int64_t ctxt_modulus = 998244353;
    ringMatrix mlwe_crs;
    matrix crs;

    // negacyclic(degree, rank, numInstance, ctxt_modulus, mlwe_crs, crs);
    instance_conversion(degree, rank, numInstance, ctxt_modulus);

    return 0;
}