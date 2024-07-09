#include "database.h"
#include "simplepir.h"

#define __DEBUG 1

int main()
{
    int numRow = 20;
    int numCol = 20;

    int degree = 4;
    int rank = 3;

    database db(numRow, numCol);
    
#if __DEBUG == 1
    cout << "=== Database ===" << endl;
    print(db.getDB());
#endif

    int qryRow = 1;
    int qryCol = 1;

    mlwe_parameter mlwe_param(db);
    pir_parameter pir_param;

#if __DEBUG == 1
    cout << "\n Parameter Set " << endl;
    db.print();
    mlwe_param.print();
    // pir_param.print();
#endif

    cout << "\n Executing Private Information Retrieval" << endl;
    cout << "   - [Server] preprocessing......" << endl;
    matrix hint_client;
    // MLWEtoLWE(mlwe_param, pir_param);
    setup(mlwe_param, pir_param, db, hint_client);
    assert(hint_client.size() == mlwe_param.getDegree() * mlwe_param.getNumInstance());
    assert(hint_client[0].size() == mlwe_param.getDegree() * mlwe_param.getRank());

    cout << "   - [Client] generating query......" << endl;
    vector<poly> qry, sk;
    query(mlwe_param, pir_param, qryCol, qry, sk);
    assert(qry.size() == mlwe_param.getNumInstance());
    assert(qry[0].size() == mlwe_param.getDegree());

    cout << "   - [Server] answer the query......" << endl;
    vector<int64_t> ans;
    answer(mlwe_param, pir_param, db, qry, ans);
    assert(ans.size() == mlwe_param.getDegree() * mlwe_param.getNumInstance());

    cout << "   - [Client] recover......" << endl;
    int64_t res;
    recover(mlwe_param, ans, hint_client, sk, qryRow, res);

    cout << "\nAnswer:\t" << db.getDB()[qryRow][qryCol] << endl;
    cout << "Result:\t" << res << endl;

#if __DEBUG == 2
    poly tmp1, tmp2;
    int64_t ctxt_modulus = mlwe_param.getCtxtModulus();
    int numInstance = mlwe_param.getNumInstance();


    for(int i = 0; i < numInstance; i++)
    {
        tmp1 = qry[i];
        for(int j = 0; j < sk.size(); j++)
        {
            multiply_ntt(mlwe_param.getCRS()[i][j], sk[j], tmp2, ctxt_modulus, 3);
            for(int k = 0; k < degree; k++)
            {
                tmp1[k] = tmp1[k] - tmp2[k];
                tmp1[k] = (tmp1[k] % ctxt_modulus + ctxt_modulus) % ctxt_modulus;
            }
        }
        print(tmp1);
    }

    vector<int64_t> sk_lwe, tmp_ans;
    for(int i = 0; i < sk.size(); i++)
    {
        for(int j = 0; j < degree; j++)
        {
            sk_lwe.push_back(sk[i][j]);
        }
    }
    matrixMultiply(hint_client, sk_lwe, ctxt_modulus, tmp_ans);
    print(tmp_ans);
    for(int i = 0; i < ans.size(); i++)
    {
        tmp_ans[i] = ans[i] - tmp_ans[i];
        tmp_ans[i] = (tmp_ans[i] % ctxt_modulus + ctxt_modulus) % ctxt_modulus;
    }
    print(tmp_ans);
#endif

    return 0;
}
