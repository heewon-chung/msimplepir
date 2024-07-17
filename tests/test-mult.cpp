#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <random>

#include "database.h"
#include "util.h"

using namespace std;

int main()
{
    int numRow = 3;
    int numCol = 4;
    int numInstance = 5;

    int qryCol = 1;

    int ptxt_modulus = 101;
    int ctxt_modulus = 173;

    // DB: numRow * numCol
    database db(numRow, numCol);
    cout << "Database:\n";
    print(db.getDB());

    vector<int64_t> qry(numCol), res;
    qry[qryCol]++;

    matrixMultiply(db.getDB(), qry, ptxt_modulus, res);
    cout << "Plaintext Multiplication:\t"; print(res);


    // A: numCol * numInstance
    matrix A(numCol, vector<int64_t>(numInstance)), hint;
    randMatrix(A, numCol, numInstance, ctxt_modulus);
    // hint: numRow * numInstance     
    matrixMultiply(db.getDB(), A, ctxt_modulus, hint);

    // sk: numInstance * 1
    vector<int64_t> sk(numInstance), qryCtxt;
    randVector(sk, sk.size(), ctxt_modulus);
    // qryCtxt = A * sk
    // qryCtxt: numCol * 1
    matrixMultiply(A, sk, ctxt_modulus, qryCtxt);
    qryCtxt[qryCol]++;

    vector<int64_t> tmp, multCtxt;
    // multCtxt = DB * qryCtxt
    // multCtxt: numRow * 1
    matrixMultiply(db.getDB(), qryCtxt, ctxt_modulus, multCtxt);
    // tmp = hint * sk
    // tmp: numRow * 1
    matrixMultiply(hint, sk, ctxt_modulus, tmp);

    for(int i = 0; i < tmp.size(); i++)
    {
        multCtxt[i] -= tmp[i];
        multCtxt[i] = (multCtxt[i] % ctxt_modulus + ctxt_modulus) % ctxt_modulus;
        multCtxt[i] %= ptxt_modulus;
    }

    cout << "Ciphertext Multiplication:\t"; print(multCtxt);
}