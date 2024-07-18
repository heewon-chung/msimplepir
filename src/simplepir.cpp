#include "simplepir.h"

void MLWEtoLWE(parameter& param)
{
    int rank = param.getRank();
    int degree = param.getDegree();
    int numInstance = param.getNumInstance();
    uint64_t ctxt_modulus = param.getCtxtModulus();
    // ringMatrix mlwe_crs = param.getCRSforMLWE();
    ringMatrix mlwe_crs;
    randMatrix(mlwe_crs, numInstance, rank, degree, ctxt_modulus); // move to setup

    assert(mlwe_crs.size() == numInstance);

    int numRow = numInstance * degree;
    int numCol = rank * degree;
    int rowIdx = 0;

    matrix lwe_crs(numRow, vector<int64_t>(numCol));

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
                    lwe_crs[rowIdx + deg][colIdx] = tmp[deg];
                }
                std::rotate(tmp.rbegin(), tmp.rbegin() + 1, tmp.rend());
                tmp[0] = ctxt_modulus - tmp[0];
            }
        }
    }

    param.setCRSforMLWE(mlwe_crs);
    param.setCRSforLWE(lwe_crs);

}


void setup(parameter& param, const database& db, matrix& hint_client)
{
    int db_col = db.getNumCol();
    int degree = param.getDegree();
    int rank = param.getRank();
    int numInstance = param.getNumInstance();
    uint64_t ctxtModulus = param.getCtxtModulus();

    MLWEtoLWE(param);
    // hint_client = DB * CRS
    matrixMultiply(db.getDB(), param.getCRSforLWE(), ctxtModulus, hint_client);

}

void query(const parameter& param, const int col, vector<poly>& qry, vector<poly>& sk)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<uint64_t> error_distribution(0, 3); // Use uint64_t for distribution range
    
    int rank = param.getRank();
    int degree = param.getDegree();
    int numInstance = param.getNumInstance();
    uint64_t ctxt_modulus = param.getCtxtModulus();
    int root = param.getRoot();
    const ringMatrix& mlwe_crs = param.getCRSforMLWE();

    assert(mlwe_crs.size() == numInstance);
    assert(mlwe_crs[0].size() == rank);
    
    int blk = col / degree;
    int pos = col % degree;
    
    sk.resize(rank);
    vector<poly> sk_ntt(rank);

    // Generate random vectors and compute their NTT in parallel
    #pragma omp parallel for
    for(int j = 0; j < rank; j++)
    {
        randVector(sk[j], degree, ctxt_modulus); // TODO sample ternary for the efficient
        sk_ntt[j] = poly(sk[j].begin(), sk[j].end());
        sk_ntt[j].resize(2 * degree);
        ntt(sk_ntt[j], ctxt_modulus, root);
    }

    // Initialize qry
    qry.resize(numInstance, poly(degree, 0));

    // Perform matrix-vector multiplication in parallel
    #pragma omp parallel for
    for(int i = 0; i < numInstance; ++i)
    {
        vector<poly> tmp(rank, poly(2 * degree));
        // Compute NTT of mlwe_crs rows
        for(int j = 0; j < rank; ++j)
        {
            tmp[j] = poly(mlwe_crs[i][j].begin(), mlwe_crs[i][j].end());
            tmp[j].resize(2 * degree);
            ntt(tmp[j], ctxt_modulus, root);
            multiply_ntt(tmp[j], sk_ntt[j], tmp[j], ctxt_modulus, root, true);
            // Add tmp to qry[i] and ensure results are within modulus
            for(int k = 0; k < degree; ++k)
            {
                #pragma omp atomic
                qry[i][k] = (qry[i][k] + tmp[j][k]) % ctxt_modulus;
            }
        }
    }

    qry[blk][pos] = (qry[blk][pos] + param.getScale()) % ctxt_modulus;
}


void query(const parameter& param, const int qryCol, vector<int64_t> qry, int64_t ctxt_modulus)
{
    matrix crs = param.getCRSforLWE();
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


void answer(const parameter& param, const database& db, const vector<poly>& qry, vector<int64_t>& ans)
{
    assert(qry.size() * qry[0].size() == db.getNumCol());
    vector<int64_t> qryLWE;
    qryLWE.reserve(qry.size() * qry[0].size());

    // Flatten the qry vector into qryLWE
    for(const auto& poly : qry)
    {
        qryLWE.insert(qryLWE.end(), poly.begin(), poly.end());
    }

    matrixMultiply(db.getDB(), qryLWE, param.getCtxtModulus(), ans);
}

void recover(const parameter& param, const vector<int64_t>& ans, const matrix& hint_client, const vector<poly>& sk, const int qryRow, int64_t& res)
{
    uint64_t ctxt_modulus = param.getCtxtModulus();
    res = ans[qryRow];

    int idx = 0;
    int64_t tmp = 0;

    for(const auto& poly : sk)
    {
        assert(poly.size() == param.getDegree());
        for(const auto& element : poly)
        {
            tmp += hint_client[qryRow][idx] * element;
            tmp = (tmp % ctxt_modulus + ctxt_modulus) % ctxt_modulus;
            idx++;
        }
    }
    res = (res - tmp) % ctxt_modulus;
    res = (res + ctxt_modulus) % ctxt_modulus;
    res /= param.getScale();

    assert(idx == static_cast<int>(sk.size() * sk[0].size()));
}