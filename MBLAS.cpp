#include "MBLAS.h"
#include <omp.h>

using namespace std;

// mm_ddot : 状態ベクトルの内積計算の状態行列version
double mm_ddot(int row_dim, int col_dim, double **M1, double **M2)
{
    double val = 0.;
    for (int i = 0; i < row_dim; i++)
    {
        for (int j = 0; j < col_dim; j++)
        {
            val += M1[i][j] * M2[i][j];
        }
    }
    return val;
}

// mm_dcopy : 状態行列のコピー
void mm_dcopy(int row_dim, int col_dim, double **M1, double **M2)
{
    for (int i = 0; i < row_dim; i++)
    {
        for (int j = 0; j < col_dim; j++)
        {
            M2[i][j] = M1[i][j];
        }
    }
};

// mm_dscal M1 = α*M1
void mm_dscal(int row_dim, int col_dim, double alpha, double **M1)
{
    for (int i = 0; i < row_dim; i++)
    {
        for (int j = 0; j < col_dim; j++)
        {
            M1[i][j] *= alpha;
        }
    }
}

// mm_daxpy M2 += α*M1
void mm_daxpy(int row_dim, int col_dim, double alpha, double **M1, double **M2)
{
    for (int i = 0; i < row_dim; i++)
    {
        for (int j = 0; j < col_dim; j++)
        {
            M2[i][j] += alpha * M1[i][j];
        }
    }
}

double mm_dnrm2(int row_dim, int col_dim, double **M)
{
    double val = 0.;
    for (int i = 0; i < row_dim; i++)
    {
        for (int j = 0; j < col_dim; j++)
        {
            val += M[i][j] * M[i][j];
        }
    }
    return sqrt(val);
}

// 状態行列の規格化
void mm_sdz(int row_dim, int col_dim, double **M)
{
    double a = 1. / mm_dnrm2(row_dim, col_dim, M);
    mm_dscal(row_dim, col_dim, a, M);
};

// 状態ベクトルの成分を状態行列にコピーする M = vec
void mv_copy(int mat_dim, double *vec, double **M)
{
    for (int i = 0; i < mat_dim; i++)
    {
        for (int j = 0; j < mat_dim; j++)
        {
            M[i][j] = vec[i * mat_dim + j];
        }
    }
}

void sdz(int dim, double *vec)
{
    double a = 1. / cblas_dnrm2(dim, vec, 1);
    cblas_dscal(dim, a, vec, 1);
}

// 状態行列の表示
void print_mat(int mat_dim, double **M)
{
    double mtmp;
    for (int row_num = 0; row_num < mat_dim; row_num++)
    {
        std::cout << "[";
        for (int col_num = 0; col_num < mat_dim; col_num++)
        {
            mtmp = M[row_num][col_num];
            // printf("%5.8e", mtmp);
            std::cout << std::scientific << std::setprecision(15) << mtmp;
            if (col_num < mat_dim - 1)
            {
                std::cout << "  ";
            }
        }
        if (row_num < mat_dim - 1)
        {
            std::cout << "];" << std::endl;
        }
        else
        {
            std::cout << "]";
        }
    }
    std::cout << "]" << std::endl;
}
// M x N状態行列の表示
void print_mat(int dim_A, int dim_B, std::vector<std::vector<double>> M)
{
    double mtmp;
    for (int row_num = 0; row_num < dim_A; row_num++)
    {
        std::cout << "[";
        for (int col_num = 0; col_num < dim_B; col_num++)
        {
            mtmp = M[row_num][col_num];
            // printf("%5.8e", mtmp);
            std::cout << std::scientific << std::setprecision(15) << mtmp;
            if (col_num < dim_B - 1)
            {
                std::cout << "  ";
            }
        }
        if (row_num < dim_A - 1)
        {
            std::cout << "];" << std::endl;
        }
        else
        {
            std::cout << "]";
        }
    }
    std::cout << "]" << std::endl;
}

void print_mat(int mat_dim, std::vector<std::vector<double>> M)
{
    double mtmp;
    for (int row_num = 0; row_num < mat_dim; row_num++)
    {
        std::cout << "[";
        for (int col_num = 0; col_num < mat_dim; col_num++)
        {
            mtmp = M[row_num][col_num];
            // printf("%5.8e", mtmp);
            std::cout << std::scientific << std::setprecision(15) << mtmp;
            if (col_num < mat_dim - 1)
            {
                std::cout << "  ";
            }
        }
        if (row_num < mat_dim - 1)
        {
            std::cout << "];" << std::endl;
        }
        else
        {
            std::cout << "]";
        }
    }
    std::cout << "]" << std::endl;
}

// 状態ベクトルの表示
void print_vec(int mat_dim, double *vec)
{
    double mtmp;

    std::cout << "[";
    for (int col_num = 0; col_num < mat_dim; col_num++)
    {
        mtmp = vec[col_num];
        // printf("%5.8e", mtmp);
        std::cout << std::scientific << std::setprecision(15) << mtmp;
        if (col_num < mat_dim - 1)
        {
            std::cout << "  ";
        }
    }
    std::cout << "]";

    std::cout << "]" << std::endl;
}

/*---------------OpenMPを利用した並列処理を施したversion---------------*/
// mm_ddot : 状態ベクトルの内積計算の状態行列version
double MP_mm_ddot(int row_dim, int col_dim, double **M1, double **M2)
{
    double val = 0.;
#pragma omp parallel for reduction(+ \
                                   : val)
    for (int i = 0; i < row_dim; i++)
    {
        for (int j = 0; j < col_dim; j++)
        {
            val += M1[i][j] * M2[i][j];
        }
    }
    return val;
}

double MP_schedule_mm_ddot(int row_dim, int col_dim, double **M1, double **M2)
{
    double val = 0.;
#pragma omp parallel for schedule(runtime) reduction(+ \
                                                     : val)
    for (int i = 0; i < row_dim; i++)
    {
        for (int j = 0; j < col_dim; j++)
        {
            val += M1[i][j] * M2[i][j];
        }
    }
    return val;
}

void MP_mm_dcopy(int row_dim, int col_dim, double **M1, double **M2)
{
#pragma omp parallel for
    for (int i = 0; i < row_dim; i++)
    {
        for (int j = 0; j < col_dim; j++)
        {
            M2[i][j] = M1[i][j];
        }
    }
};

void MP_schedule_mm_dcopy(int row_dim, int col_dim, double **M1, double **M2)
{
#pragma omp parallel for schedule(runtime)
    for (int i = 0; i < row_dim; i++)
    {
        for (int j = 0; j < col_dim; j++)
        {
            M2[i][j] = M1[i][j];
        }
    }
}

double MP_mm_dnrm2(int row_dim, int col_dim, double **M1)
{
    double val = 0.;
#pragma omp parallel for reduction(+ \
                                   : val)
    for (int i = 0; i < row_dim; i++)
    {
        for (int j = 0; j < col_dim; j++)
        {
            val += M1[i][j] * M1[i][j];
        }
    }
    return sqrt(val);
}

double MP_schedule_mm_dnrm2(int row_dim, int col_dim, double **M)
{
    double val = 0.;
#pragma omp parallel for schedule(runtime) reduction(+ \
                                                     : val)
    for (int i = 0; i < row_dim; i++)
    {
        for (int j = 0; j < col_dim; j++)
        {
            val += M[i][j] * M[i][j];
        }
    }
    return sqrt(val);
}

// mm_dscal M1 = α*M1
void MP_mm_dscal(int row_dim, int col_dim, double alpha, double **M1)
{
#pragma omp parallel for
    for (int i = 0; i < row_dim; i++)
    {
        for (int j = 0; j < col_dim; j++)
        {
            M1[i][j] *= alpha;
        }
    }
}

void MP_schedule_mm_dscal(int row_dim, int col_dim, double alpha, double **M1)
{
#pragma omp parallel for schedule(runtime)
    for (int i = 0; i < row_dim; i++)
    {
        for (int j = 0; j < col_dim; j++)
        {
            M1[i][j] *= alpha;
        }
    }
}

// mm_daxpy M2 += α*M1
void MP_mm_daxpy(int row_dim, int col_dim, double alpha, double **M1, double **M2)
{
    int j;
#pragma omp parallel for private(j)
    for (int i = 0; i < row_dim; i++)
    {
        for (j = 0; j < col_dim; j++)
        {
            M2[i][j] += alpha * M1[i][j];
        }
    }
}

void MP_schedule_mm_daxpy(int row_dim, int col_dim, double alpha, double **M1, double **M2)
{
#pragma omp parallel for schedule(runtime)
    for (int i = 0; i < row_dim; i++)
    {
        for (int j = 0; j < col_dim; j++)
        {
            M2[i][j] += alpha * M1[i][j];
        }
    }
}

void MP_mm_sdz(int row_dim, int col_dim, double **M)
{
    double a = 1. / MP_mm_dnrm2(row_dim, col_dim, M);
    MP_mm_dscal(row_dim, col_dim, a, M);
};

void MP_schedule_mm_sdz(int row_dim, int col_dim, double **M)
{
    double a = 1. / MP_schedule_mm_dnrm2(row_dim, col_dim, M);
    MP_schedule_mm_dscal(row_dim, col_dim, a, M);
};

/*[Memo]並列化コードのテストを行うために用意した関数。実際のHamiltonian行列対角化コードとは実装が異なることに注意*/
void isoA_mmprod(int dim_A, int dim_B, double **u_i, double **u_j)
{
    dim_A = 4;
    dim_B = 4;

    // memory確保
    int *row = new int[5];
    int *col_ptr = new int[6];
    double *val = new double[6];

    // 具体的な要素の代入
    row[0] = 0;
    row[1] = 1;
    row[2] = 3;
    row[3] = 5;
    row[4] = 6;

    col_ptr[0] = 0;
    col_ptr[1] = 1;
    col_ptr[2] = 2;
    col_ptr[3] = 1;
    col_ptr[4] = 2;
    col_ptr[5] = 3;

    val[0] = 1.;
    val[1] = 2.;
    val[2] = 3.;
    val[3] = 4.;
    val[4] = 5.;
    val[5] = 6.;

    // 行列ベクトル積の計算
    for (int row_num = 0; row_num < dim_A; row_num++)
    {
        for (int k = col_ptr[row_num]; k < col_ptr[row_num + 1]; k++)
        {
            for (int col_num = 0; col_num < dim_B; col_num++)
            {
                u_j[row_num][col_num] += val[k] * u_i[row[k]][col_num];
            }
        }
    }

    delete[] row;
    delete[] col_ptr;
    delete[] val;
}

void MP_isoA_mmprod(int dim_A, int dim_B, double **u_i, double **u_j)
{
    dim_A = 4;
    dim_B = 4;

    // memory確保
    int *row = new int[5];
    int *col_ptr = new int[6];
    double *val = new double[6];

    // 具体的な要素の代入
    row[0] = 0;
    row[1] = 1;
    row[2] = 3;
    row[3] = 5;
    row[4] = 6;

    col_ptr[0] = 0;
    col_ptr[1] = 1;
    col_ptr[2] = 2;
    col_ptr[3] = 1;
    col_ptr[4] = 2;
    col_ptr[5] = 3;

    val[0] = 1.;
    val[1] = 2.;
    val[2] = 3.;
    val[3] = 4.;
    val[4] = 5.;
    val[5] = 6.;

    // 行列ベクトル積の計算

    int k;
    int col_num;
#pragma omp parallel for private(k, col_num)
    for (int row_num = 0; row_num < dim_A; row_num++)
    {
        for (k = col_ptr[row_num]; k < col_ptr[row_num + 1]; k++)
        {
            for (col_num = 0; col_num < dim_B; col_num++)
            {
                u_j[row_num][col_num] += val[k] * u_i[row[k]][col_num];
            }
        }
    }
    delete[] row;
    delete[] col_ptr;
    delete[] val;
}

void isoB_mmprod(int dim_A, int dim_B, double **u_i, double **u_j)
{
    dim_A = 4;
    dim_B = 4;

    // memory確保
    int *row = new int[5];
    int *col_ptr = new int[6];
    double *val = new double[6];

    // 具体的な要素の代入
    row[0] = 0;
    row[1] = 1;
    row[2] = 3;
    row[3] = 5;
    row[4] = 6;

    col_ptr[0] = 0;
    col_ptr[1] = 1;
    col_ptr[2] = 2;
    col_ptr[3] = 1;
    col_ptr[4] = 2;
    col_ptr[5] = 3;

    val[0] = 1.;
    val[1] = 2.;
    val[2] = 3.;
    val[3] = 4.;
    val[4] = 5.;
    val[5] = 6.;

    for (int row_num = 0; row_num < dim_A; row_num++)
    {
        for (int col_num = 0; col_num < dim_B; col_num++)
        {
            for (int k = col_ptr[col_num]; k < col_ptr[col_num + 1]; k++)
            {
                u_j[row_num][col_num] += u_i[row_num][row[k]] * val[k];
            }
        }
    }
    /*--------------------------------------------------*/
    delete[] row;
    delete[] col_ptr;
    delete[] val;
}

void MP_isoB_mmprod(int dim_A, int dim_B, double **u_i, double **u_j)
{
    dim_A = 4;
    dim_B = 4;

    // memory確保
    int *row = new int[5];
    int *col_ptr = new int[6];
    double *val = new double[6];

    // 具体的な要素の代入
    row[0] = 0;
    row[1] = 1;
    row[2] = 3;
    row[3] = 5;
    row[4] = 6;

    col_ptr[0] = 0;
    col_ptr[1] = 1;
    col_ptr[2] = 2;
    col_ptr[3] = 1;
    col_ptr[4] = 2;
    col_ptr[5] = 3;

    val[0] = 1.;
    val[1] = 2.;
    val[2] = 3.;
    val[3] = 4.;
    val[4] = 5.;
    val[5] = 6.;

    int k, col_num;

#pragma omp parallel for private(col_num, k)
    for (int row_num = 0; row_num < dim_A; row_num++)
    {
        for (col_num = 0; col_num < dim_B; col_num++)
        {
            for (k = col_ptr[col_num]; k < col_ptr[col_num + 1]; k++)
            {
                u_j[row_num][col_num] += u_i[row_num][row[k]] * val[k];
            }
        }
    }

    delete[] row;
    delete[] col_ptr;
    delete[] val;
}

void int_rise_mmprod(int nnz_pA, int nnz_pB, int *prow_ind_A, int *pcol_ind_A, int *prow_ind_B, int *pcol_ind_B, double **V0, double **V1_dic1)
{
    for (int i = 0; i < nnz_pA; i++)
    {
        for (int j = 0; j < nnz_pB; j++)
        {
            V1_dic1[prow_ind_A[i]][pcol_ind_B[j]] += 0.5 * 1.0 * V0[pcol_ind_A[i]][prow_ind_B[j]]; // テストの簡略化のためにJ.val(0) = 1.0とした
        }
    }
};

void MP_int_rise_mmprod(int nnz_pA, int nnz_pB, int *prow_ind_A, int *pcol_ind_A, int *prow_ind_B, int *pcol_ind_B, double **V0, double **V1_dic1)
{
    int j;
#pragma omp parallel for private(j)
    for (int i = 0; i < nnz_pA; i++)
    {
        for (j = 0; j < nnz_pB; j++)
        {
            V1_dic1[prow_ind_A[i]][pcol_ind_B[j]] += 0.5 * 1.0 * V0[pcol_ind_A[i]][prow_ind_B[j]]; // テストの簡略化のためにJ.val(0) = 1.0とした
        }
    }
};

void int_dsmn_mmprod(int nnz_mA, int nnz_mB, int *mrow_ind_A, int *mcol_ind_A, int *mrow_ind_B, int *mcol_ind_B, double **V0, double **V1_inc1)
{
    for (int i = 0; i < nnz_mA; i++)
    {
        for (int j = 0; j < nnz_mB; j++)
        {
            V1_inc1[mrow_ind_A[i]][mcol_ind_B[j]] += 0.5 * 1.0 * V0[mcol_ind_A[i]][mrow_ind_B[j]]; // テストの簡略化のためにJ.val(0)=1.0とした
        }
    }
};

void MP_int_dsmn_mmprod(int nnz_mA, int nnz_mB, int *mrow_ind_A, int *mcol_ind_A, int *mrow_ind_B, int *mcol_ind_B, double **V0, double **V1_inc1)
{
    int j;
#pragma omp parallel for private(j)
    for (int i = 0; i < nnz_mA; i++)
    {
        for (j = 0; j < nnz_mB; j++)
        {
            V1_inc1[mrow_ind_A[i]][mcol_ind_B[j]] += 0.5 * 1.0 * V0[mcol_ind_A[i]][mrow_ind_B[j]]; // テストの簡略化のためにJ.val(0)=1.0とした
        }
    }
};

void int_zz_mmprod(int dim_A, int dim_B, int *sz_A, int *sz_B, double **V0, double **V1)
{
    for (int i = 0; i < dim_A; i++)
    {
        for (int j = 0; j < dim_B; j++)
        {
            V1[i][j] += 0.25 * 1.0 * sz_A[i] * sz_B[j] * V0[i][j]; // テスト簡略化のためにJ.val(0)=1.0とした
        }
    }
};

void MP_int_zz_mmprod(int dim_A, int dim_B, int *sz_A, int *sz_B, double **V0, double **V1)
{
    int j;
#pragma omp parallel for private(j)
    for (int i = 0; i < dim_A; i++)
    {
        for (j = 0; j < dim_B; j++)
        {
            V1[i][j] += 0.25 * 1.0 * sz_A[i] * sz_B[j] * V0[i][j]; // テスト簡略化のためにJ.val(0)=1.0とした
        }
    }
};

// testを行う際はls < 10とする
void calc_alpha(const int ls, const int bm_A_size, const int bm_B_size, double *alpha, double **V0, double **V1)
{
    int pair_num = 5;

    for (int No = 0; No < pair_num; No++)
        alpha[ls] += mm_ddot(bm_A_size, bm_B_size, V1, V0);
}

void MP_calc_alpha(const int ls, const int bm_A_size, const int bm_B_size, double *alpha, double **V0, double **V1)
{
    int pair_num = 5;

    double tmp = 0.;
#pragma omp parallel for reduction(+ \
                                   : tmp)
    for (int No = 0; No < pair_num; No++)
    {
        tmp += MP_mm_ddot(bm_A_size, bm_B_size, V1, V0);
    }
    alpha[ls] = tmp; //[Memo] alpha[ls] += tmpのほうがいいかな？
}

// testを行う際はls < 10, tri_mat_dim = 10とする
void calc_beta_odd_step(const int tri_mat_dim, const int ls, double *alpha, double *beta, double **V1, double **V0)
{
    int pair_num = 5;
    int bm_A_size;
    int bm_B_size;

    // 並列化を意識してloop中のif文を除外したコード
    if (ls != tri_mat_dim - 1)
    {
        for (int No = 0; No < pair_num; No++)
        {
            mm_daxpy(-alpha[ls], bm_A_size, bm_B_size, V1, V0);
            beta[ls] += mm_ddot(bm_A_size, bm_B_size, V0, V0);
        }
        beta[ls] = sqrt(beta[ls]);
    }
    else
    {
        for (int No = 0; No < pair_num; No++)
        {
            mm_daxpy(-alpha[ls], bm_A_size, bm_B_size, V1, V0);
        }
    }

    // 修正前のコード
    //  for (int No = 0; No < pair_num; No++)
    //  {
    //      mm_daxpy(-alpha[ls], bm_A_size, bm_B_size, V1, V0);
    //      if (ls != tri_mat_dim - 1)
    //          beta[ls] += mm_ddot(bm_A_size, bm_B_size, V0, V0);
    //  }
    //  if (ls != tri_mat_dim - 1)
    //  {
    //      beta[ls] = sqrt(beta[ls]);
    //  }
}

void MP_calc_beta_odd_step(const int tri_mat_dim, const int ls, double *alpha, double *beta, double **V1, double **V0)
{
    int pair_num = 5;
    int bm_A_size;
    int bm_B_size;

    double tmp = 0;
    // 並列化を意識してloop中のif文を除外したコード
    if (ls != tri_mat_dim - 1)
    {
#pragma omp parallel for reduction(+ \
                                   : tmp)
        for (int No = 0; No < pair_num; No++)
        {
            MP_mm_daxpy(-alpha[ls], bm_A_size, bm_B_size, V1, V0);
            tmp += MP_mm_ddot(bm_A_size, bm_B_size, V0, V0);
        }
        beta[ls] = sqrt(tmp);
    }
    else
    {
#pragma omp parallel for
        for (int No = 0; No < pair_num; No++)
        {
            MP_mm_daxpy(-alpha[ls], bm_A_size, bm_B_size, V1, V0);
        }
    }
}

void calc_beta_even_step(const int tri_mat_dim, const int ls, double *alpha, double *beta, double **V0, double **V1)
{
    int pair_num = 5;
    int bm_A_size = 50000;
    int bm_B_size = 50000;

    if (ls != tri_mat_dim - 1)
    {
        for (int No = 0; No < pair_num; No++)
        {
            mm_daxpy(-alpha[ls], bm_A_size, bm_B_size, V0, V1);
            beta[ls] += mm_ddot(bm_A_size, bm_B_size, V1, V1);
        }
        beta[ls] = sqrt(beta[ls]);
    }
    else
    {
        for (int No = 0; No < pair_num; No++)
        {
            mm_daxpy(-alpha[ls], bm_A_size, bm_B_size, V0, V1);
        }
    }

    // for (int No = 0; No < pair_num; No++)
    // {
    //     mm_daxpy(-alpha[ls], bm_A_size, bm_B_size, V0, V1);
    //     if (ls != tri_mat_dim - 1)
    //         beta[ls] += mm_ddot(bm_A_size, bm_B_size, V1, V1);
    // }
    // if (ls != tri_mat_dim - 1)
    // {
    //     beta[ls] = sqrt(beta[ls]);
    // }
}

void MP_calc_beta_even_step(const int tri_mat_dim, const int ls, double *alpha, double *beta, double **V0, double **V1)
{
    int pair_num = 5;
    int bm_A_size = 50000;
    int bm_B_size = 50000;
    double tmp = 0.;
    if (ls != tri_mat_dim - 1)
    {
#pragma omp parallel for reduction(+ \
                                   : tmp)
        for (int No = 0; No < pair_num; No++)
        {
            mm_daxpy(-alpha[ls], bm_A_size, bm_B_size, V0, V1);
            tmp += mm_ddot(bm_A_size, bm_B_size, V1, V1);
        }
        beta[ls] = sqrt(tmp);
    }
    else
    {
#pragma omp parallel for
        for (int No = 0; No < pair_num; No++)
        {
            mm_daxpy(-alpha[ls], bm_A_size, bm_B_size, V0, V1);
        }
    }
}