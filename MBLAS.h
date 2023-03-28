#ifndef ___MBLAS_
#define ___MBLAS_

#include <math.h>
#include <mkl.h>

#include <iomanip>
#include <iostream>
#include <random>

/*------------------非並列処理version-----------------*/
double mm_ddot(int row_dim, int col_dim, double **M1, double **M2);
void mm_dcopy(int row_dim, int col_dim, double **M1, double **M2);
void mm_dscal(int row_dim, int col_dim, double alpha, double **M1);
void mm_daxpy(int row_dim, int col_dim, double alpha, double **M1, double **M2);
double mm_dnrm2(int row_dim, int col_dim, double **M);
void mv_copy(int row_dim, int col_dim, double *vec, double **M);
void print_mat(int mat_dim, double **M);
void print_vec(int mat_dim, double *vec);
void sdz(int dim, double *vec);
void mm_sdz(int row_dim, int col_dim, double **M);
void isoA_mmprod(int dim_A, int dim_B, double **u_i, double **u_j);
void isoB_mmprod(int dim_A, int dim_B, double **u_i, double **u_j);
void int_rise_mmprod(int nnz_pA, int nnz_pB, int *prow_ind_A, int *pcol_ind_A, int *prow_ind_B, int *pcol_ind_B, double **V0, double **V1_dic1);
void int_dsmn_mmprod(int nnz_mA, int nnz_mB, int *mrow_ind_A, int *mcol_ind_A, int *mrow_ind_B, int *mcol_ind_B, double **V0, double **V1_inc1);
void int_zz_mmprod(int dim_A, int dim_B, int *sz_A, int *sz_B, double **V0, double **V1);
void calc_alpha(const int ls, const int bm_A_size, const int bm_B_size, double *alpha, double **V0, double **V1);
void calc_beta_odd_step(const int tri_mat_dim, const int ls, double *alpha, double *beta, double **V0, double **V1);
void calc_beta_even_step(const int tri_mat_dim, const int ls, double *alpha, double *beta, double **V0, double **V1);
/*--------------------------------------------------*/

/*------------------OpenMP利用version(schedulingは非実施)---------------*/
double MP_mm_ddot(int row_dim, int col_dim, double **M1, double **M2);              // tested
void MP_mm_dcopy(int row_dim, int col_dim, double **M1, double **M2);               // tested
void MP_mm_dscal(int row_dim, int col_dim, double alpha, double **M1);              // tested
void MP_mm_daxpy(int row_dim, int col_dim, double alpha, double **M1, double **M2); // tested
double MP_mm_dnrm2(int row_dim, int col_dim, double **M);                           // tested
void MP_mv_copy(int row_dim, int col_dim, double *vec, double **M);
void MP_sdz(int dim, double *vec);
void MP_mm_sdz(int row_dim, int col_dim, double **M);
void MP_schedule_mm_sdz(int row_dim, int col_dim, double **M);
void MP_isoA_mmprod(int dim_A, int dim_B, double **u_i, double **u_j);
void MP_isoB_mmprod(int dim_A, int dim_B, double **u_i, double **u_j);
void MP_int_rise_mmprod(int nnz_pA, int nnz_pB, int *prow_ind_A, int *pcol_ind_A, int *prow_ind_B, int *pcol_ind_B, double **V0, double **V1_dic1);
void MP_int_dsmn_mmprod(const int nnz_mA, const int nnz_mB, int *mrow_ind_A, int *mcol_ind_A, int *mrow_ind_B, int *mcol_ind_B, double **V0, double **V1_inc1);
void MP_int_zz_mmprod(int dim_A, int dim_B, int *sz_A, int *sz_B, double **V0, double **V1);
void MP_calc_alpha(const int ls, const int bm_A_size, const int bm_B_size, double *alpha, double **V0, double **V1);
void MP_calc_beta_odd_step(const int tri_mat_dim, const int ls, double *alpha, double *beta, double **V0, double **V1);
void MP_calc_beta_even_step(const int tri_mat_dim, const int ls, double *alpha, double *beta, double **V0, double **V1);
/*----------------------------------------------------------------------*/

/*---------------OpenMP利用かつschedulingを実施---------------*/
double MP_schedule_mm_ddot(int row_dim, int col_dim, double **M1, double **M2);
void MP_schedule_mm_dcopy(int row_dim, int col_dim, double **M1, double **M2);
void MP_schedule_mm_dscal(int row_dim, int col_dim, double alpha, double **M1);
void MP_schedule_mm_daxpy(int row_dim, int col_dim, double alpha, double **M1, double **M2);
double MP_schedule_mm_dnrm2(int row_dim, int col_dim, double **M);
/*------------------------------------------------------------*/
/*動作確認用コード*/
void test_ddot(int mat_dim, double **M1, double **M2); // mm_ddot、MP_mm_ddot, MP_schedule_mm_ddotのテストを行う
void print_mat(int mat_dim, std::vector<std::vector<double>> M);
void print_mat(int dim_A, int dim_B, std::vector<std::vector<double>> M);
#endif