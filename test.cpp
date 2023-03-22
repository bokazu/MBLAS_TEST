#include "MBLAS.h"

using namespace std;

void test_mv_copy(int test_num, int mat_dim, double *u, double **M)
{
    bool flag = true;
    for (int t = 0; t < test_num; t++)
    {
        mv_copy(mat_dim, u, M);
        for (int i = 0; i < mat_dim; i++)
        {
            for (int j = 0; j < mat_dim; j++)
            {
                if (u[i * mat_dim + j] != M[i][j])
                    flag = false;
            }
        }
    }

    if (flag)
        cout << "test_mv_copy::success!\n";
    else
        cout << "test_mv_copy::success!\n";
}

// test_num回だけ内積計算を行い、計算結果と計算時間をファイル出力する
// 内積計算が正しく行われていれば、"success"という文字列を標準出力する
void test_mm_ddot(int test_num, int mat_dim, double *u1, double *u2,
                  double **M1, double **M2, string PATH_result_of_time, string PATH_result_of_calc)
{

    bool flag = true;

    int state_num = mat_dim * mat_dim;
    double vec_ddot = 0.;
    double mat_ddot = 0.;

    // u1とM1、u2とM2が対応するようにする
    mv_copy(mat_dim, u1, M1);
    mv_copy(mat_dim, u2, M2);

    for (int i = 0; i < test_num; i++)
    {
        vec_ddot = cblas_ddot(state_num, u1, 1, u2, 1);
        mat_ddot = mm_ddot(mat_dim, M1, M2);

        if (vec_ddot != mat_ddot)
        {
            flag = false;
        }
    }
}