#ifndef ___TEST_
#define ___TEST_

/*------------テストを行うための関数をこのヘッダーにまとめる--------------*/

#include <math.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>

template <typename T>
void test1(std::vector<T> &data1, std::vector<T> &data2, std::string filepath, char mode)
{
    int N = data1.size();
    std::ofstream ofs(filepath);
    // fileをopenできるか確認する
    if (!ofs)
    {
        printf("@test1(data1, data2, filepath, mode) error:: \"%s\" could not open.", filepath);
    }
    else
    {
        if (mode == 'c') // csv出力
        {
            for (int i = 0; i < N; i++)
            {
                ofs << data1[i] << "," << data2[i] << std::endl;
            }
        }
        else // txt出力
        {
            for (int i = 0; i < N; i++)
            {
                ofs << data1[i] << " " << data2[i] << std::endl;
            }
        }
    }
    ofs.close();
};

#endif