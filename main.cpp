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

int main()
{
    omp_set_num_threads(NUM_THREADS); //[comment : ここでsetするのではなく、ターミナルで一度設定するだけで良いらしい]
    cout << "max_thread_num = " << omp_get_max_threads() << endl;
    // cout << "thread_num = " << omp_get_thread_num() << endl;

    // mmprod();
    // isoB_mmprod();
    int state_num = 16;
    int dim_A = 50000;
    int dim_B = 50000;
    int daxpy_mat_dim = 10000; // daxpyのテストコードは他と比べて消費メモリ量が大きいので別途変数を用意した

    /*--------------------Test1--------------------*/
    // cout << "\n\n Test1 : mm_dcopy()\n";
    // cout << "====================================================\n";
    // test_mm_dcopy(10, dim_A, dim_B,"./test_dcopy/result_dcopy.csv", "./test_dcopy/time_dcopy.csv");
    // cout << "====================================================\n";
    /*---------------------------------------------*/

    /*----------------Test2 : M1とM2の内積計算を行う------------*/
    // cout << "\n\n Test2 : mm_ddot()\n";
    // cout << "====================================================\n";
    // test_mm_ddot(10, dim_A, dim_B ,"./test_ddot/result_ddot.csv", "./test_ddot/time_ddot.csv");
    // cout << "====================================================\n";
    /*-----------------------------------------------------------*/

    /*----------------Test3 : M1とM2のノルム計算を行う------------*/
    // cout << "\n\n Test3 : mm_dnrm2()\n";
    // cout << "====================================================\n";
    // test_mm_dnrm2(10, dim_A, dim_B,"./test_dnrm2/result_dnrm2.csv", "./test_dnrm2/time_dnrm2.csv");
    // cout << "====================================================\n";
    /*-----------------------------------------------------------*/

    /*----------------Test4 : dscalの計算の確認------------------*/
    // cout << "\n\n Test4 : mm_dscal()\n";
    // cout << "====================================================\n";
    // test_mm_dscal(10, dim_A, dim_B,"./test_dscal/result_dscal.csv", "./test_dscal/time_dscal.csv");
    // cout << "====================================================\n";
    /*-----------------------------------------------------------*/

    /*----------------Test5 : daxpyの計算の確認------------------*/
    // cout << "\n\n Test5 : mm_daxpy()\n";
    // cout << "====================================================\n";
    // test_mm_daxpy(10, daxpy_mat_dim, daxpy_mat_dim,"./test_daxpy/result_daxpy.csv", "./test_daxpy/time_daxpy.csv");
    // cout << "====================================================\n";
    /*-----------------------------------------------------------*/

    /*----------------Test6 : sdzの計算の確認------------------*/
    // cout << "\n\n Test6 : mm_sdz()\n";
    // cout << "====================================================\n";
    // test_mm_sdz(10, dim_A, dim_B,"./test_sdz/result_sdz.csv", "./test_sdz/time_sdz.csv");
    // cout << "====================================================\n";
    /*-----------------------------------------------------------*/

    /*---------------Test7 : 疎行列-行列積計算の確認-------------*/
    //[Memo]fileへのデータの書き込みは追記モードで行っているので、テストを行うごとにファイルを削除する必要あり
    // test_isoA_mmprod(10, 4, 4, "./test_isoA_mmprod/result_isoA_mmprod.csv", "./test_isoA_mmprod/time_isoA_mmprod.csv");
    // test_isoB_mmprod(10, 4, 4, "./test_isoB_mmprod/result_isoB_mmprod.csv", "./test_isoB_mmprod/time_isoB_mmprod.csv");

    /*---------------Test8 : S^+ V S^+行列積計算の確認-------------*/
    // test_int_rise_mmprod(10, 4, 4, "./test_rise_mmprod/result_rise_mmprod.csv", "./test_rise_mmprod/time_rise_mmprod.csv");
    // test_int_dsmn_mmprod(10, 4, 4, "./test_dsmn_mmprod/result_dsmn_mmprod.csv", "./test_dsmn_mmprod/time_dsmn_mmprod.csv");
    // test_int_zz_mmprod(10, 4, 4, "./test_zz_mmprod/result_zz_mmprod.csv", "./test_zz_mmprod/time_zz_mmprod.csv");

    /*---------------Test8 : Lanczos法における三重対角行列の主対角成分-------------*/
    // test_calc_alpha(10, dim_A, dim_B, "./test_calc_alpha/result_calc_alpha.csv", "./test_calc_alpha/time_calc_alpha.csv");

    /*---------------Test8 : Lanczos法における三重対角行列の次対角成分-------------*/
    test_calc_beta_even(10, dim_A, dim_B, "./test_calc_beta/result_calc_beta.csv", "./test_calc_beta/time_calc_beta.csv");

    // /*-----------状態ベクトルのメモリ確保&random初期化-----------*/
    // double *u1 = new double[state_num];
    // double *u2 = new double[state_num];

    // random_device rand;
    // mt19937 mt(rand());
    // uniform_real_distribution<> rand1(0, 1);
    // uniform_real_distribution<> rand2(0, 1);
    // for (int i = 0; i < state_num; i++)
    //     u1[i] = rand1(mt);
    // for (int i = 0; i < state_num; i++)
    //     u2[i] = rand1(mt);
    // /*-----------------------------------------------------------*/

    /*------------------------状態行列の用意(簡単のために正方行列として用意する)---------------------*/
    // double **M1 = new double *[mat_dim];
    // for (int i = 0; i < mat_dim; i++)
    //     M1[i] = new double[mat_dim];

    // double **M2 = new double *[mat_dim];
    // for (int i = 0; i < mat_dim; i++)
    //     M2[i] = new double[mat_dim];

    // for (int i = 0; i < mat_dim; i++)
    // {
    //     for (int j = 0; j < mat_dim; j++)
    //     {
    //         M1[i][j] = rand1(mt);
    //         M2[i][j] = rand1(mt);
    //     }
    // }
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
    // chrono::system_clock::time_point mm_dcopy_start, MP_mm_dcopy_start, mm_dcopy_end, MP_mm_dcopy_end;
    // mm_dcopy_start = chrono::system_clock::now();
    // mm_dcopy(mat_dim, M1, M2);
    // mm_dcopy_end = chrono::system_clock::now();

    // string mm_result = "success";
    // bool is_copy = true;
    // for (int i = 0; i < mat_dim; i++)
    // {
    //     if (is_copy)
    //     {
    //         for (int j = 0; j < mat_dim; j++)
    //         {
    //             if (M1[i][j] != M2[i][j])
    //             {
    //                 is_copy = false;
    //                 mm_result = "failed.";
    //                 break;
    //             }
    //         }
    //     }
    //     else
    //         break;
    // }

    // // 行列M2の初期化
    // for (int i = 0; i < mat_dim; i++)
    // {
    //     for (int j = 0; j < mat_dim; j++)
    //     {
    //         M2[i][j] = 0.;
    //     }
    // }

    // MP_mm_dcopy_start = chrono::system_clock::now();
    // MP_mm_dcopy(mat_dim, M1, M2);
    // MP_mm_dcopy_end = chrono::system_clock::now();

    // string MP_result = "success";
    // is_copy = true;
    // for (int i = 0; i < mat_dim; i++)
    // {
    //     if (is_copy)
    //     {
    //         for (int j = 0; j < mat_dim; j++)
    //         {
    //             if (M1[i][j] != M2[i][j])
    //             {
    //                 is_copy = false;
    //                 MP_result = "failed.";
    //                 break;
    //             }
    //         }
    //     }
    //     else
    //         break;
    // }
    // auto mm_dcopy_time_msec = chrono::duration_cast<chrono::milliseconds>(mm_dcopy_end - mm_dcopy_start).count();
    // auto MP_mm_dcopy_time_msec = chrono::duration_cast<chrono::milliseconds>(MP_mm_dcopy_end - MP_mm_dcopy_start).count();

    // cout << "@mat_dcopy = " << mm_result << endl;
    // cout << "@mat_MP_dcopy = " << MP_result << endl;
    // cout << "@time(mm_ddot) = " << mm_dcopy_time_msec << "msec" << endl;
    // cout << "@time(MP_ddot) = " << MP_mm_dcopy_time_msec << "msec" << endl;

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

    // for (int i = 0; i < mat_dim; i++)
    //     delete[] M1[i];
    // delete[] M1;

    // for (int i = 0; i < mat_dim; i++)
    //     delete[] M2[i];
    // delete[] M2;
}
