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

    matrix crs(numRow, vector<int64_t>(numCol));

    for(int i = 0; i < numInstance; ++i)
    {
        int rowIdx = i * degree;
        for(int j = 0; j < rank; ++j)
        {
            poly tmp = mlwe_crs[i][j];
            for(int k = 0; k < degree; ++k)
            {
                int colIdx = j * degree + k;
                for(int deg = 0; deg < degree; ++deg)
                {
                    crs[rowIdx + deg][colIdx] = tmp[deg];
                }
                std::rotate(tmp.rbegin(), tmp.rbegin() + 1, tmp.rend());
                tmp[0] = ctxt_modulus - tmp[0];
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
    // #pragma omp parallel for
    for(int j = 0; j < rank; j++)
    {
        randVector(sk[j], degree, ctxt_modulus);
    }

    // mlwe_crs * sk
    qry.resize(numInstance, poly(degree, 0));
    // #pragma omp parallel for
    for(int i = 0; i < numInstance; ++i)
    {
        poly tmp(degree);
        for(int j = 0; j < rank; ++j)
        {
            // Multiply mlwe_crs[i][j] with sk[j] using NTT
            multiply_ntt(mlwe_crs[i][j], sk[j], tmp, ctxt_modulus, 3);
            // Add tmp to qry[i] and ensure results are within modulus
            for(int k = 0; k < degree; ++k)
            {
                // #pragma omp atomic
                qry[i][k] = (qry[i][k] + tmp[k]) % ctxt_modulus;
            }
        }
    }

    qry[blk][pos]++;

}


void query(const pir_parameter& pir_param, const int qryCol, vector<int64_t> qry, int64_t ctxt_modulus)
{
    matrix crs = pir_param.getCRS();
    int numRow = crs.size();
    int numCol = crs[0].size();
    qry.resize(numRow);

    vector<int64_t> sk_lwe;
    randVector(sk_lwe, numCol, ctxt_modulus);

    for(int i = 0; i < numRow; i++)
    {
        for(int j = 0; j < crs[i].size(); j++)
        {
            qry[i] += crs[i][j] * sk_lwe[j];
        }
        qry[i] = (qry[i] % ctxt_modulus + ctxt_modulus) % ctxt_modulus;
    }

    qry[qryCol]++;
}


void answer(const mlwe_parameter& mlwe_param, const pir_parameter& pir_param, const database& db, const vector<poly>& qry, vector<int64_t>& ans)
{
    assert(qry.size() * qry[0].size() == db.getNumCol());
    vector<int64_t> qryLWE;
    qryLWE.reserve(qry.size() * qry[0].size());

    // Flatten the qry vector into qryLWE
    for(const auto& poly : qry)
    {
        qryLWE.insert(qryLWE.end(), poly.begin(), poly.end());
    }

    matrixMultiply(db.getDB(), qryLWE, mlwe_param.getCtxtModulus(), ans);
}

void recover(const mlwe_parameter& mlwe_param, const vector<int64_t>& ans, const matrix& hint_client, const vector<poly>& sk, const int qryRow, int64_t& res)
{
    uint64_t ctxt_modulus = mlwe_param.getCtxtModulus();
    res = ans[qryRow];

    int idx = 0;
    int64_t tmp = 0;

    for(const auto& poly : sk)
    {
        assert(poly.size() == mlwe_param.getDegree());
        for(const auto& element : poly)
        {
            tmp += hint_client[qryRow][idx] * element;
            tmp = (tmp % ctxt_modulus + ctxt_modulus) % ctxt_modulus;
            idx++;
        }
    }
    res = (res - tmp) % ctxt_modulus;
    res = (res + ctxt_modulus) % ctxt_modulus;

    assert(idx == static_cast<int>(sk.size() * sk[0].size()));
}