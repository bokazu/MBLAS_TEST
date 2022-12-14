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