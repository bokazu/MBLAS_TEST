#ifndef ___TEST_
#define ___TEST_

/*------------テストを行うための関数をこのヘッダーにまとめる--------------*/

#include <math.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include "MBLAS.h"

// mm_ddot, MP_mm_ddot, MP_schedule_mm_ddotの(i)計算結果,(ii)計算時間のテストを行う
void test_mm_ddot(int trial_num, int mat_dim, std::string PATH_result_of_calc, std::string PATH_result_of_time, char mode = 'c');
void test_mm_dcopy(int trial_num, int mat_dim, std::string PATH_result_of_calc, std::string PATH_result_of_time, char mode = 'c');
void test_mm_dnrm2(int trial_num, int mat_dim, std::string PATH_result_of_calc, std::string PATH_result_of_time, char mode = 'c');
void test_mm_dscal(int trial_num, int mat_dim, std::string PATH_result_of_calc, std::string PATH_result_of_time, char mode = 'c');
void test_mm_daxpy(int trial_num, int mat_dim, std::string PATH_result_of_calc, std::string PATH_result_of_time, char mode = 'c');
void test_mm_sdz(int trial_num, int mat_dim, std::string PATH_result_of_calc, std::string PATH_result_of_time, char mode = 'c');
void test_isoA_mmprod(int trial_num, int dim_A, int dim_B, std::string PATH_result_of_calc, std::string PATH_result_of_time, char mode = 'c');
void test_isoB_mmprod(int trial_num, int dim_A, int dim_B, std::string PATH_result_of_calc, std::string PATH_result_of_time, char mode = 'c');
void test_int_rise_mmprod(int trial_num, int dim_A, int dim_B, std::string PATH_result_of_calc, std::string PATH_result_of_time, char mode = 'c');
void test_int_dsmn_mmprod(int trial_num, int dim_A, int dim_B, std::string PATH_result_of_calc, std::string PATH_result_of_time, char mode = 'c');
void test_int_zz_mmprod(int trial_num, int dim_A, int dim_B, std::string PATH_result_of_calc, std::string PATH_result_of_time, char mode = 'c');
#endif