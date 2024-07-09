#ifndef __DATABASE
#define __DATABASE
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

#include "util.h"

using namespace std;

class database
{
private:
    int numRow;
    int numCol;
    // uint64_t ptxt_modulus;
    matrix db;

public:
    database();
    database(int numRow, int numCol);

    void randDB();
    matrix getDB() const {return db;}

    int getNumRow() const {return numRow;}
    int getNumCol() const {return numCol;}

    void setElem(int row, int col, int val) {db[row][col] = val;}    
    int getElem(int row, int col) {return db[row][col];}

    void print();

};
#endif