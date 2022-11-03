#include <math.h>
#include <mkl.h>

#include <iomanip>
#include <iostream>

// mm_ddot : 状態ベクトルの内積計算の状態行列version
double mm_ddot(int mat_dim, double **M1, double **M2)
{
    double val = 0.;
    for (int i = 0; i < mat_dim; i++)
    {
        for (int j = 0; j < mat_dim; j++)
        {
            val += M1[i][j] * M2[i][j];
        }
    }
    return val;
}

// mm_dcopy : 状態行列のコピー
void mm_dcopy(int mat_dim, double **M1, double **M2)
{
    for (int i = 0; i < mat_dim; i++)
    {
        for (int j = 0; j < mat_dim; j++)
        {
            M2[i][j] = M1[i][j];
        }
    }
};

// mm_dscal M1 = α*M1
void mm_dscal(int mat_dim, double alpha, double **M1)
{
    for (int i = 0; i < mat_dim; i++)
    {
        for (int j = 0; j < mat_dim; j++)
        {
            M1[i][j] *= alpha;
        }
    }
}

// mm_daxpy M2 += α*M1
void mm_daxpy(int mat_dim, double alpha, double **M1, double **M2)
{
    for (int i = 0; i < mat_dim; i++)
    {
        for (int j = 0; j < mat_dim; j++)
        {
            M2[i][j] += alpha * M1[i][j];
        }
    }
}

double mm_dnrm2(int mat_dim, double **M)
{
    double val = 0.;
    for (int i = 0; i < mat_dim; i++)
    {
        for (int j = 0; j < mat_dim; j++)
        {
            val += M[i][j] * M[i][j];
        }
    }
    return sqrt(val);
}

//状態行列の規格化
void mm_sdz(int mat_dim, double **M)
{
    double a = 1. / mm_dnrm2(mat_dim, M);
    mm_dscal(mat_dim, a, M);
};

//状態ベクトルの成分を状態行列にコピーする M = vec
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

//状態行列の表示
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
            std::cout << std::scientific << std::setprecision(4) << mtmp;
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

//状態ベクトルの表示
void print_vec(int mat_dim, double *vec)
{
    double mtmp;

    std::cout << "[";
    for (int col_num = 0; col_num < mat_dim; col_num++)
    {
        mtmp = vec[col_num];
        // printf("%5.8e", mtmp);
        std::cout << std::scientific << std::setprecision(4) << mtmp;
        if (col_num < mat_dim - 1)
        {
            std::cout << "  ";
        }
    }
    std::cout << "]";

    std::cout << "]" << std::endl;
}