#include "parameter.h"

mlwe_parameter::mlwe_parameter()
{
    // security = 128;
    degree = 0;
    rank = 0;
    
    ptxt_modulus = 0;
    ctxt_modulus = 0;

    ringMatrix crs;
    // randMatrix(crs, rank, degree, ctxt_modulus);
    crs.clear();

    this -> crs = crs;
}

mlwe_parameter::mlwe_parameter(const database& db)
{
    int numCol = db.getNumCol();
    int numRow = db.getNumRow();

    // p = k 2^m + 1
    // 2m-th root of unity
    // ctxt_modulus = 17;
    ctxt_modulus = 998244353;
    ptxt_modulus = 10;

    degree = 4;
    rank = 3;

    numInstance = numRow / degree;

    ringMatrix crs;
    randMatrix(crs, numInstance, rank, degree, ctxt_modulus);

    this -> crs = crs;
}


void mlwe_parameter::print()
{
    cout << "   - degree: " << degree << endl;
    cout << "   - rank: " << rank << endl;
    cout << "   - number of instances: " << numInstance << endl;
    cout << "   - MLWE CRS size: " << crs.size() << " * " << crs[0].size() << endl;
}


pir_parameter::pir_parameter(database& db)
{
    matrix ref_db = db.getDB();

    crs.resize(ref_db.size());
    for(int i = 0; i < crs.size(); i++)
    {
        crs[i].resize(ref_db[i].size());
    }
}


void pir_parameter::print()
{
    // cout << "   - ctxt modulus: " << ctxt_modulus << endl;
    cout << "   - PIR CRS size: " << crs.size() << " * " << crs[0].size() << endl;
}