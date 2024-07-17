#ifndef __SIMPLEPIR
#define __SIMPLEPIR
#include <cassert>

#include "parameter.h"
#include "database.h"

void MLWEtoLWE(const mlwe_parameter& mlwe_param, pir_parameter& pir_param);
void setup(const mlwe_parameter& mlwe_param, pir_parameter& pir_param, const database& db, matrix& hint_client);
void query(const mlwe_parameter& mlwe_param, const pir_parameter& pir_param, const int col, vector<poly>& qry, vector<poly>& sk);
void answer(const mlwe_parameter& mlwe_param, const pir_parameter& pir_param, const database& db, const vector<poly>& qry, vector<int64_t>& ans);
void recover(const mlwe_parameter& mlwe_param, const vector<int64_t>& ans, const matrix& hint_client, const vector<poly>& sk, const int qryRow, int64_t& res);

void query(const pir_parameter& pir_param, const int qryCol, vector<int64_t> qry, int64_t ctxt_modulus);

#endif