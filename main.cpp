#include <mkl.h>

#include <iomanip>
#include <iostream>
#include <random>
#include <vector>
#include <chrono>
#include <omp.h>

#include "MBLAS.h"
#include "test_output.hpp"

#define NUM_THREADS 20

using namespace std;

void isoB_mmprod()
{
    int mat_dim = 4;
    vector<int> row_ptr = {0, 1, 3, 5, 6};
    vector<int> col = {0, 1, 2, 1, 2, 3};
    vector<double> val = {1, 2, 3, 4, 5, 6};
    vector<vector<double>> M0(mat_dim, vector<double>(mat_dim, 1.0));
    vector<vector<double>> M1(mat_dim, vector<double>(mat_dim, 0.0));

    M0[1][1] = 2.;
    M0[1][2] = 2.;
    M0[2][1] = 3.;
    M0[2][2] = 3.;

    for (int i = 0; i < mat_dim; i++)
    {
        for (int j = 0; j < mat_dim; j++)
        {
            for (int k = row_ptr[j]; k < row_ptr[j + 1]; k++)
            {
                M1[i][j] += M0[i][col[k]] * val[k];
            }
        }
    }
    cout << "M1 = \n";
    print_mat(mat_dim, M1);
}

void mmprod()
{
    int dim_A = 4;
    int dim_B = 3;
    vector<vector<double>> A(dim_A, vector<double>(dim_A, 0.0));
    vector<vector<double>> B(dim_B, vector<double>(dim_B, 0.0));
    vector<vector<double>> psi0(dim_A, vector<double>(dim_B, 1.0));
    vector<vector<double>> psi1(dim_A, vector<double>(dim_B, 0.0));

    A[0][1] = 2.0;
    A[1][3] = 2.0;

    B[0][2] = 2.0;
    B[2][1] = 2.0;

    for (int i = 0; i < dim_A; i++)
    {
        for (int j = 0; j < dim_B; j++)
        {
            for (int k = 0; k < dim_A; k++)
            {
                for (int l = 0; l < dim_B; l++)
                {
                    psi1[i][j] += A[i][k] * psi0[k][l] * B[l][j];
                }
            }
        }
    }
    cout << "psi1 = \n";
    print_mat(dim_A, dim_B, psi1);
}
int main()
{
    omp_set_num_threads(NUM_THREADS);
    cout << "max_thread_num = " << omp_get_max_threads() << endl;
    cout << "thread_num = " << omp_get_thread_num() << endl;

    mmprod();
    // isoB_mmprod();
    int state_num = 16;
    int mat_dim = 50000;
    // /*-----------状態ベクトルのメモリ確保&random初期化-----------*/
    // double *u1 = new double[state_num];
    // double *u2 = new double[state_num];

    random_device rand;
    mt19937 mt(rand());
    uniform_real_distribution<> rand1(0, 1);
    uniform_real_distribution<> rand2(0, 1);
    // for (int i = 0; i < state_num; i++)
    //     u1[i] = rand1(mt);
    // for (int i = 0; i < state_num; i++)
    //     u2[i] = rand1(mt);
    // /*-----------------------------------------------------------*/

    /*------------------------状態行列の用意(簡単のために正方行列として用意する)---------------------*/
    double **M1 = new double *[mat_dim];
    for (int i = 0; i < mat_dim; i++)
        M1[i] = new double[mat_dim];

    double **M2 = new double *[mat_dim];
    for (int i = 0; i < mat_dim; i++)
        M2[i] = new double[mat_dim];

    for (int i = 0; i < mat_dim; i++)
    {
        for (int j = 0; j < mat_dim; j++)
        {
            M1[i][j] = rand1(mt);
            M2[i][j] = rand2(mt);
        }
    }
    /*-----------------------------------------------------------*/

    /*----------------Test1 : uの要素をMにコピーする------------*/
    // mv_copy(mat_dim, u1, M1);
    // mv_copy(mat_dim, u2, M2);

    // //結果の確認
    // cout << "\n\n Test1 : mv_copyt()\n";
    // cout << "====================================================\n";
    // cout << "u1 = ";
    // print_vec(state_num, u1);
    // cout << "M1 = \n";
    // print_mat(mat_dim, M1);

    // cout << "u2 = ";
    // print_vec(state_num, u2);
    // cout << "M2 = \n";
    // print_mat(mat_dim, M2);
    // cout << "====================================================\n";
    // /*-----------------------------------------------------------*/

    /*----------------Test2 : M1とM2の内積計算を行う------------*/
    chrono::system_clock::time_point mm_ddot_start, MP_mm_ddot_start, mm_ddot_end, MP_mm_ddot_end;
    // double vec_ddot = cblas_ddot(state_num, u1, 1, u2, 1);

    cout << "\n\n Test2 : mm_ddot()\n";
    cout << "====================================================\n";
    string outputfile1 = "./test_ddot/time_ddot.csv";
    string outputfile2 = "./test_ddot/result_ddot.csv";
    vector<double> time_mm_ddot;
    vector<double> time_MP_ddot;
    vector<double> res_mm_ddot;
    vector<double> res_MP_ddot;

    for (int i = 0; i < 100; i++)
    {
        mm_ddot_start = chrono::system_clock::now();
        double mat_ddot = mm_ddot(mat_dim, M1, M2);
        mm_ddot_end = chrono::system_clock::now();
        MP_mm_ddot_start = chrono::system_clock::now();
        double MP_mat_ddot = MP_mm_ddot(mat_dim, M1, M2);
        MP_mm_ddot_end = chrono::system_clock::now();

        auto mm_ddot_time_msec = chrono::duration_cast<chrono::milliseconds>(mm_ddot_end - mm_ddot_start).count();
        auto MP_mm_ddot_time_msec = chrono::duration_cast<chrono::milliseconds>(MP_mm_ddot_end - MP_mm_ddot_start).count();

        /*test1関数に渡す用のデータ*/
        res_mm_ddot.push_back(mat_ddot);
        res_MP_ddot.push_back(MP_mat_ddot);
        time_mm_ddot.push_back(mm_ddot_time_msec);
        time_MP_ddot.push_back(MP_mm_ddot_time_msec);
    }

    test1<double>(res_mm_ddot, res_MP_ddot, outputfile1, 'c');
    test1<double>(time_mm_ddot, time_MP_ddot, outputfile2, 'c');
    // cout << "@vec_ddot = " << scientific << setprecision(15) << vec_ddot << endl;
    // cout << "@mat_ddot = " << scientific << setprecision(15) << mat_ddot << endl;
    // cout << "@time(mm_ddot) = " << mm_ddot_time_msec << "msec" << endl;
    // cout << "@MP_mat_ddot = " << scientific << setprecision(15) << MP_mat_ddot << endl;
    // cout << "@time(MP_ddot) = " << MP_mm_ddot_time_msec << "msec" << endl;
    cout << "====================================================\n";
    /*-----------------------------------------------------------*/

    // /*----------------Test3 : M1とM2のノルム計算を行う------------*/
    // double u1_norm = cblas_dnrm2(state_num, u1, 1);
    // double M1_norm = mm_dnrm2(mat_dim, M1);

    // double u2_norm = cblas_dnrm2(state_num, u2, 1);
    // double M2_norm = mm_dnrm2(mat_dim, M2);

    // cout << "\n\n Test3 : mm_dnrm2()\n";
    // cout << "====================================================\n";
    // cout << "@u1_norm = " << scientific << setprecision(15) << u1_norm << endl;
    // cout << "@M1_norm = " << scientific << setprecision(15) << M1_norm << endl;
    // cout << "\n";
    // cout << "@u2_norm = " << scientific << setprecision(15) << u2_norm << endl;
    // cout << "@M2_norm = " << scientific << setprecision(15) << M2_norm << endl;
    // cout << "====================================================\n";
    // /*-----------------------------------------------------------*/

    // /*----------------Test4 : dscalの計算の確認------------------*/
    // double alpha = 2.0;
    // cblas_dscal(state_num, alpha, u1, 1);
    // mm_dscal(mat_dim, alpha, M1);

    // cblas_dscal(state_num, alpha, u2, 1);
    // mm_dscal(mat_dim, alpha, M2);

    // cout << "\n\n Test : mm_dscal() 4\n";
    // cout << "====================================================\n";
    // cout << "u1 = ";
    // print_vec(state_num, u1);
    // cout << "M1 = \n";
    // print_mat(mat_dim, M1);

    // cout << "u2 = ";
    // print_vec(state_num, u2);
    // cout << "M2 = \n";
    // print_mat(mat_dim, M2);
    // cout << "====================================================\n";
    // /*--------------------------------------------------------------*/

    // /*-------------------Test5 : daxpyの計算の確認-----------------*/
    // cblas_daxpy(state_num, alpha, u1, 1, u2, 1);
    // mm_daxpy(mat_dim, alpha, M1, M2);

    // cout << "\n\n Test5 : mm_daxpy()\n";
    // cout << "====================================================\n";
    // cout << "u2 = ";
    // print_vec(state_num, u2);
    // cout << "M2 = \n";
    // print_mat(mat_dim, M2);
    // cout << "====================================================\n";
    // /*---------------------------------------------------------------*/

    // /*---------------------Test5 : sdzの計算の確認-------------------*/
    // sdz(state_num, u1);
    // mm_sdz(mat_dim, M1);

    // cout << "\n\n Test5 mm_sdz()\n";
    // cout << "====================================================\n";
    // cout << "u1 = ";
    // print_vec(state_num, u1);
    // cout << "M1 = \n";
    // print_mat(mat_dim, M1);
    // cout << "====================================================\n";
    /*---------------------------------------------------------------*/

    /*----------------Test3 : M1とM2のノルム計算を行う------------*/
    // u1_norm = cblas_dnrm2(state_num, u1, 1);
    // M1_norm = mm_dnrm2(mat_dim, M1);

    // cout << "\n\n Test3\n";
    // cout << "====================================================\n";
    // cout << "@u1_norm = " << u1_norm << endl;
    // cout << "@M1_norm = " << M1_norm << endl;
    // cout << "====================================================\n";
    /*-----------------------------------------------------------*/
    /*----------------Test2 : M1の要素をM2にコピーする----------*/
    // mm_dcopy(mat_dim, M1, M2);

    // //結果の確認
    // cout << "====================================================\n";
    // cout << "M1 = \n";
    // print_mat(mat_dim, M1);
    // cout << "M2 = \n";
    // print_mat(mat_dim, M2);
    // cout << "====================================================\n";
    /*-----------------------------------------------------------*/
    /*---------------------- -メモリの開放-----------------------*/
    // delete[] u1;
    // delete[] u2;

    for (int i = 0; i < mat_dim; i++)
        delete[] M1[i];
    delete[] M1;

    for (int i = 0; i < mat_dim; i++)
        delete[] M2[i];
    delete[] M2;
}
