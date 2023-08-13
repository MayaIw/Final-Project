#ifndef HEADER_FILE
#define HEADER_FILE

double** sym_c(double**, int , int);
double** ddg_c(double**, int, int);
double** norm_c(double**, int , int);
double** symnmf_c(double**, double**, int, int);
double** matrix_allocation(int, int);
void free_matrix(double**);

#endif