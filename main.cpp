#include <mkl.h>

#include <iomanip>
#include <iostream>
#include <random>

#include "MBLAS.h"

using namespace std;

int main()
{
    int state_num = 16;
    int mat_dim = 4;
    /*-----------状態ベクトルのメモリ確保&random初期化-----------*/
    double *u1 = new double[state_num];
    double *u2 = new double[state_num];

    random_device rand;
    mt19937 mt(rand());
    uniform_real_distribution<> rand1(0, 1);
    for (int i = 0; i < state_num; i++) u1[i] = rand1(mt);
    for (int i = 0; i < state_num; i++) u2[i] = rand1(mt);
    /*-----------------------------------------------------------*/

    /*------------------------状態行列の用意---------------------*/
    double **M1 = new double *[mat_dim];
    for (int i = 0; i < mat_dim; i++) M1[i] = new double[mat_dim];

    double **M2 = new double *[mat_dim];
    for (int i = 0; i < mat_dim; i++) M2[i] = new double[mat_dim];
    /*-----------------------------------------------------------*/

    /*----------------Test1 : uの要素をMにコピーする------------*/
    mv_copy(mat_dim, u1, M1);
    mv_copy(mat_dim, u2, M2);

    //結果の確認
    cout << "\n\n Test1\n";
    cout << "====================================================\n";
    cout << "u1 = ";
    print_vec(state_num, u1);
    cout << "M1 = \n";
    print_mat(mat_dim, M1);

    cout << "u2 = ";
    print_vec(state_num, u2);
    cout << "M2 = \n";
    print_mat(mat_dim, M2);
    cout << "====================================================\n";
    /*-----------------------------------------------------------*/

    /*----------------Test2 : M1とM2の内積計算を行う------------*/
    double vec_ddot = cblas_ddot(state_num, u1, 1, u2, 1);
    double mat_ddot = mm_ddot(mat_dim, M1, M2);

    cout << "\n\n Test2\n";
    cout << "====================================================\n";
    cout << "@vec_ddot = " << vec_ddot << endl;
    cout << "@mat_ddot = " << mat_ddot << endl;
    cout << "====================================================\n";
    /*-----------------------------------------------------------*/

    /*----------------Test3 : M1とM2のノルム計算を行う------------*/
    double u1_norm = cblas_dnrm2(state_num, u1, 1);
    double M1_norm = mm_dnrm2(mat_dim, M1);

    double u2_norm = cblas_dnrm2(state_num, u2, 1);
    double M2_norm = mm_dnrm2(mat_dim, M2);

    cout << "\n\n Test3\n";
    cout << "====================================================\n";
    cout << "@u1_norm = " << u1_norm << endl;
    cout << "@M1_norm = " << M1_norm << endl;
    cout << "\n";
    cout << "@u2_norm = " << u2_norm << endl;
    cout << "@M2_norm = " << M2_norm << endl;
    cout << "====================================================\n";
    /*-----------------------------------------------------------*/

    /*----------------Test4 : dscalの計算の確認------------------*/
    double alpha = 2.0;
    cblas_dscal(state_num, alpha, u1, 1);
    mm_dscal(mat_dim, alpha, M1);

    cblas_dscal(state_num, alpha, u2, 1);
    mm_dscal(mat_dim, alpha, M2);

    cout << "\n\n Test4\n";
    cout << "====================================================\n";
    cout << "u1 = ";
    print_vec(state_num, u1);
    cout << "M1 = \n";
    print_mat(mat_dim, M1);

    cout << "u2 = ";
    print_vec(state_num, u2);
    cout << "M2 = \n";
    print_mat(mat_dim, M2);
    cout << "====================================================\n";
    /*--------------------------------------------------------------*/

    /*-------------------Test5 : daxpyの計算の確認-----------------*/
    cblas_daxpy(state_num, alpha, u1, 1, u2, 1);
    mm_daxpy(mat_dim, alpha, M1, M2);

    cout << "\n\n Test5\n";
    cout << "====================================================\n";
    cout << "u2 = ";
    print_vec(state_num, u2);
    cout << "M2 = \n";
    print_mat(mat_dim, M2);
    cout << "====================================================\n";
    /*---------------------------------------------------------------*/

    /*---------------------Test5 : sdzの計算の確認-------------------*/
    sdz(state_num, u1);
    mm_sdz(mat_dim, M1);

    cout << "\n\n Test5\n";
    cout << "====================================================\n";
    cout << "u1 = ";
    print_vec(state_num, u1);
    cout << "M1 = \n";
    print_mat(mat_dim, M1);
    cout << "====================================================\n";
    /*---------------------------------------------------------------*/

    /*----------------Test3 : M1とM2のノルム計算を行う------------*/
    u1_norm = cblas_dnrm2(state_num, u1, 1);
    M1_norm = mm_dnrm2(mat_dim, M1);

    cout << "\n\n Test3\n";
    cout << "====================================================\n";
    cout << "@u1_norm = " << u1_norm << endl;
    cout << "@M1_norm = " << M1_norm << endl;
    cout << "====================================================\n";
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
    delete[] u1;

    for (int i = 0; i < mat_dim; i++) delete[] M1[i];
    delete[] M1;

    for (int i = 0; i < mat_dim; i++) delete[] M2[i];
    delete[] M2;
}