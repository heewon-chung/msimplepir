#include "simplepir.h"

void MLWEtoLWE(const mlwe_parameter& mlwe_param, pir_parameter& pir_param)
{
    int rank = mlwe_param.getRank();
    int degree = mlwe_param.getDegree();
    int numInstance = mlwe_param.getNumInstance();
    uint64_t ctxt_modulus = mlwe_param.getCtxtModulus();
    ringMatrix mlwe_crs = mlwe_param.getCRS();

    assert(mlwe_crs.size() == numInstance);

    int numRow = numInstance * degree;
    int numCol = rank * degree;
    int rowIdx = 0;

    matrix crs;
    crs.resize(numRow, vector<int64_t>(numCol));

    // i: row
    for(int i = 0; i < numInstance; i++)
    {
        // j: column
        int colIdx = 0;
        int rowIdx = i * degree;
        for(int j = 0; j < rank; j++)
        {
            // tmp \in R_q
            poly tmp = mlwe_crs[i][j];
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

    pir_param.setCRS(crs);

}


void setup(const mlwe_parameter& mlwe_param, pir_parameter& pir_param, const database& db, matrix& hint_client)
{
    int db_col = db.getNumCol();
    int degree = mlwe_param.getDegree();
    uint64_t ctxtModulus = mlwe_param.getCtxtModulus();

    MLWEtoLWE(mlwe_param, pir_param);
    matrix pir_crs = pir_param.getCRS();

#if __DEBUG == 1
    cout << "=== PIR CRS ===" << endl;
    printMatrix(pir_crs);
#endif

    // hint_client = DB * CRS
    matrixMultiply(db.getDB(), pir_crs, ctxtModulus, hint_client);

#if __DEBUG == 1
    cout << "=== Hint ===" << endl;
    printMatrix(hint_client);
#endif
}

void query(const mlwe_parameter& mlwe_param, const pir_parameter& pir_param, const int col, vector<poly>& qry, vector<poly>& sk)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<uint64_t> error_distribution(0, 3); // Use uint64_t for distribution range
    
    int rank = mlwe_param.getRank();
    int degree = mlwe_param.getDegree();
    int numInstance = mlwe_param.getNumInstance();
    uint64_t ctxt_modulus = mlwe_param.getCtxtModulus();
    ringMatrix mlwe_crs = mlwe_param.getCRS();

    assert(mlwe_crs.size() == numInstance);
    assert(mlwe_crs[0].size() == rank);
    
    int blk = col / degree;
    int pos = col % degree;
    
    sk.resize(rank);
    for(int i = 0; i < sk.size(); i++)
    {
        randVector(sk[i], degree, ctxt_modulus);
    }

    // mlwe_crs * sk
    qry.resize(numInstance);
    poly tmp(degree);
    for(int i = 0; i < numInstance; i++)
    {
        // mlwe_crs \cdot sk
        qry[i].resize(degree);
        for(int j = 0; j < rank; j++)
        {
            // root 3
            multiply_ntt(mlwe_crs[i][j], sk[j], tmp, ctxt_modulus, 3);
            // add to qry
            for(int k = 0; k < tmp.size(); k++)
            {
                qry[i][k] += tmp[k];
                qry[i][k] = (qry[i][k] % ctxt_modulus + ctxt_modulus) % ctxt_modulus;
            }
        }
    }

    qry[blk][pos]++;

}

void answer(const mlwe_parameter& mlwe_param, const pir_parameter& pir_param, const database& db, const vector<poly>& qry, vector<int64_t>& ans)
{
    assert(qry.size() * qry[0].size() == db.getNumCol());
    vector<int64_t> qryLWE;
    for(int i = 0; i < qry.size(); i++)
    {
        for(int j = 0; j < qry[i].size(); j++)
        {
            qryLWE.push_back(qry[i][j]);
        }
    }

    matrixMultiply(db.getDB(), qryLWE, mlwe_param.getCtxtModulus(), ans);
}

void recover(const mlwe_parameter& mlwe_param, const vector<int64_t>& ans, const matrix& hint_client, const vector<poly>& sk, const int qryRow, int64_t& res)
{
    // vector<int64_t> sk_tmp, db_tmp;
    uint64_t ptxt_modulus = mlwe_param.getPtxtModulus();
    uint64_t ctxt_modulus = mlwe_param.getCtxtModulus();
    res = ans[qryRow];

    int idx = 0;
    int64_t tmp = 0;
    for(int i = 0; i < sk.size(); i++)
    {
        assert(sk[i].size() == mlwe_param.getDegree());
        for(int j = 0; j < sk[i].size(); j++)
        {
            tmp += hint_client[qryRow][idx] * sk[i][j];
            idx++;
        }
    }
    tmp = (tmp % ctxt_modulus + ctxt_modulus) % ctxt_modulus;
    res = res - tmp;
    res = (res % ctxt_modulus + ctxt_modulus) % ctxt_modulus;
    assert(idx == sk.size() * sk[0].size());

}