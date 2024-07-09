#include "database.h"


database::database(/* args */)
{
    numRow = 0;
    numCol = 0;
}

database::database(int numRow, int numCol)
{
    this -> numRow = numRow;
    this -> numCol = numCol;

    db.resize(numRow);
    for(int i = 0; i < numRow; i++)
    {
        db[i].resize(numCol);
    }

    randDB();
}


void database::randDB()
{
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<uint64_t> msg_space(0, 150); // Use uint64_t for distribution range

    for(int i = 0; i < numRow; i++)
    {
        for(int j = 0; j < numCol; j++)
        {
            // setElem(i, j, cnt); // TODO: change to < ptxt_modulus
            setElem(i, j, msg_space(gen));
        }
    }
}

void database::print()
{
    cout << "   - DB size: " << numRow << " * " << numCol << endl;
}