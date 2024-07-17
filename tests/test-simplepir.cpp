#include <ctime>

#include "database.h"
#include "simplepir.h"

#define __DEBUG 1

long double get_time() 
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec * 1000000L + t.tv_nsec / 1000;  // Convert seconds to microseconds
}

int main()
{
    int numRow = 2048;
    int numCol = 2048;

    int degree = 256;
    int rank = 3;

    database db(numRow, numCol);
    
#if __DEBUG == 3
    cout << "=== Database ===" << endl;
    print(db.getDB());
#endif

    int qryRow = 1;
    int qryCol = 1;

    long double start, end, total_time, qry_time;

#if __DEBUG == 1
    cout << "\n Parameter Set " << endl;
    db.print();
    // pir_param.print();
#endif

    mlwe_parameter mlwe_param(db, degree, rank);
    pir_parameter pir_param;
    matrix hint_client;

    cout << "\n Executing Private Information Retrieval" << endl;
    cout << "   - [Server] preprocessing......\t";
    start = get_time();
    setup(mlwe_param, pir_param, db, hint_client);
    end = get_time();
    total_time = end - start;
    printf("%LF microseconds\n", total_time);
    assert(hint_client.size() == mlwe_param.getDegree() * mlwe_param.getNumInstance());
    assert(hint_client[0].size() == mlwe_param.getDegree() * mlwe_param.getRank());

    vector<poly> qry, sk;
    cout << "   - [Client] generating query......\t";
    start = get_time();
    query(mlwe_param, pir_param, qryCol, qry, sk);
    end = get_time(); 
    qry_time = end - start;
    printf("%LF microseconds\n", qry_time);
    assert(qry.size() == mlwe_param.getNumInstance());
    assert(qry[0].size() == mlwe_param.getDegree());

    cout << "   - [Server] answer the query......\t";
    vector<int64_t> ans;
    start = get_time(); 
    answer(mlwe_param, pir_param, db, qry, ans);
    end = get_time();
    total_time = end - start;
    printf("%LF microseconds\n", total_time); 
    assert(ans.size() == mlwe_param.getDegree() * mlwe_param.getNumInstance());

    cout << "   - [Client] recover the data......\t";
    int64_t res;
    start = get_time();
    recover(mlwe_param, ans, hint_client, sk, qryRow, res);
    end = get_time();
    total_time = end - start;
    printf("%LF microseconds\n", total_time);

    cout << "\nAnswer:\t" << db.getDB()[qryRow][qryCol] << endl;
    cout << "Result:\t" << res << endl;

    vector<int64_t> qry_lwe;
    start = get_time();
    query(pir_param, qryCol, qry_lwe, mlwe_param.getCtxtModulus());
    end = get_time();
    total_time = end - start;
    printf("\n\tCompared to LWE-based query, we can save %Lf query time.", qry_time / total_time * 100);

#if __DEBUG == 3
    poly tmp1, tmp2;
    int64_t ctxt_modulus = mlwe_param.getCtxtModulus();
    int numInstance = mlwe_param.getNumInstance();

    printf("\nCheck query ciphertexts\n");
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

    // printf("\n Check matrix-vector multiplication\n");
    vector<int64_t> sk_lwe, tmp_ans;
    for(int i = 0; i < sk.size(); i++)
    {
        for(int j = 0; j < degree; j++)
        {
            sk_lwe.push_back(sk[i][j]);
        }
    }
    matrixMultiply(hint_client, sk_lwe, ctxt_modulus, tmp_ans);
    // print(tmp_ans);

    printf("\n Check answer ciphertext\n");
    for(int i = 0; i < ans.size(); i++)
    {
        tmp_ans[i] = ans[i] - tmp_ans[i];
        tmp_ans[i] = (tmp_ans[i] % ctxt_modulus + ctxt_modulus) % ctxt_modulus;
    }
    print(tmp_ans);

    for(int i = 0; i < ans.size(); i++)
    {
        if(tmp_ans[i] != db.getDB()[i][qryCol]) cerr << i << "th data is wrong!" << endl;
    }
#endif

    return 0;
}
