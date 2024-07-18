#ifndef __SIMPLEPIR
#define __SIMPLEPIR
#include <cassert>

#include "parameter.h"
#include "database.h"

void MLWEtoLWE(parameter& param);
void setup(parameter& param, const database& db, matrix& hint_client);
void query(const parameter& param, const int col, vector<poly>& qry, vector<poly>& sk);
void answer(const parameter& param, const database& db, const vector<poly>& qry, vector<int64_t>& ans);
void recover(const parameter& param, const vector<int64_t>& ans, const matrix& hint_client, const vector<poly>& sk, const int qryRow, int64_t& res);

void query(const parameter& param, const int qryCol, vector<int64_t>& qry, int64_t ctxt_modulus);

#endif