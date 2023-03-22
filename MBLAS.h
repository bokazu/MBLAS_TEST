#ifndef ___MBLAS_
#define ___MBLAS_

#include <math.h>
#include <mkl.h>

#include <iomanip>
#include <iostream>
#include <random>

/*------------------非並列処理version-----------------*/
double mm_ddot(int mat_dim, double **M1, double **M2);
void mm_dcopy(int mat_dim, double **M1, double **M2);
void mm_dscal(int mat_dim, double alpha, double **M1);
void mm_daxpy(int mat_dim, double alpha, double **M1, double **M2);
double mm_dnrm2(int mat_dim, double **M);
void mv_copy(int mat_dim, double *vec, double **M);
void print_mat(int mat_dim, double **M);
void print_vec(int mat_dim, double *vec);
void sdz(int dim, double *vec);
void mm_sdz(int mat_dim, double **M);
/*--------------------------------------------------*/

/*------------------OpenMP利用version(schedulingは非実施)---------------*/
double MP_mm_ddot(int mat_dim, double **M1, double **M2);              // tested
void MP_mm_dcopy(int mat_dim, double **M1, double **M2);               // tested
void MP_mm_dscal(int mat_dim, double alpha, double **M1);              // tested
void MP_mm_daxpy(int mat_dim, double alpha, double **M1, double **M2); // tested
double MP_mm_dnrm2(int mat_dim, double **M);                           // tested
void MP_mv_copy(int mat_dim, double *vec, double **M);
void MP_sdz(int dim, double *vec);
void MP_mm_sdz(int mat_dim, double **M);
void MP_schedule_mm_sdz(int mat_dim, double **M);
/*----------------------------------------------------------------------*/

/*---------------OpenMP利用かつschedulingを実施---------------*/
double MP_schedule_mm_ddot(int mat_dim, double **M1, double **M2);
void MP_schedule_mm_dcopy(int mat_dim, double **M1, double **M2);
void MP_schedule_mm_dscal(int mat_dim, double alpha, double **M1);
void MP_schedule_mm_daxpy(int mat_dim, double alpha, double **M1, double **M2);
double MP_schedule_mm_dnrm2(int mat_dim, double **M);
/*------------------------------------------------------------*/
/*動作確認用コード*/
void test_ddot(int mat_dim, double **M1, double **M2); // mm_ddot、MP_mm_ddot, MP_schedule_mm_ddotのテストを行う
void print_mat(int mat_dim, std::vector<std::vector<double>> M);
void print_mat(int dim_A, int dim_B, std::vector<std::vector<double>> M);
#endif