# include <stdio.h>
# include <stdlib.h>
# include <math.h>

double sum_of_row(double **mat, int row_index, int num_of_rows){
    int j;
    double sum=0.0;
    for(j=0; j<num_of_rows; j++){
        sum += mat[row_index][j];
    }
    return sum;
}

double squared_distance(double *p, double *q, int d)
{
    int i;
    double sum_of_squares = 0;
    for(i=0; i<d; i++)
    {
        sum_of_squares += pow(p[i]-q[i],2);
    }
    return sum_of_squares;
}

double** multiply_matrices(double** mat1, int rows1, int cols1, double** mat2, int cols2) {
    double** result;
    double* result_1d;
    int i;
    int j;
    int l;
    /*memory allocation for the multiplication matrix*/
    result_1d = calloc(rows1*cols2, sizeof(double));
    if(result_1d == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    result = calloc(rows1,sizeof(double *));
    if(result == NULL){
        printf("An Error Has Occurred\n");
        free(result_1d);
        exit(1);
    }
    for(i=0; i<rows1; i++)
    {
        result[i] = result_1d+i*cols2;
    }

    /*calculation*/
    for (int i = 0; i < rows1; i++) {
        for (int j = 0; j < cols2; j++) {
            for (int l = 0; l < cols1; l++) {
                result[i][j] += mat1[i][l] * mat2[l][j];
            }
        }
    }
    return result;
}

double** mult_by_transpose(double **mat, int rows, int cols){
    double** result;
    double* result_1d;
    int i;
    int j;
    int l;
    /*memory allocation for the multiplication matrix*/
    result_1d = calloc(rows*rows, sizeof(double));
    if(result_1d == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    result = calloc(rows,sizeof(double *));
    if(result == NULL){
        printf("An Error Has Occurred\n");
        free(result_1d);
        exit(1);
    }
    for(i=0; i<rows; i++)
    {
        result[i] = result_1d+i*rows;
    }
    /*calculation*/
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < rows; j++) {
            for (int l = 0; l < cols; l++) {
                result[i][j] += mat[i][l] * mat[j][l];
            }
        }
    }
    return result;
}

double** sym(double **X, int num_of_elements, int d){
    int i=0;
    int j=0;
    double *matrix_1d;
    double **matrix;
    /*memory allocation for the similarity matrix*/
    matrix_1d = calloc(num_of_elements*num_of_elements, sizeof(double));
    if(matrix_1d == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    matrix = calloc(num_of_elements,sizeof(double *));
    if(matrix == NULL){
        printf("An Error Has Occurred\n");
        free(matrix_1d);
        exit(1);
    }
    for(i=0; i<num_of_elements; i++)
    {
        matrix[i] = matrix_1d+i*num_of_elements;
    } 

    /*computation of the similarity matrix*/
    for(i=0; i<num_of_elements; i++){
        for(j=i+1; j<num_of_elements; j++){
            matrix[i][j]=exp(-0.5*squared_distance(X[i], X[j], d));
            matrix[j][i]=matrix[i][j];
        }
    }
    return matrix;
}

double** diagonal_degree_mat(double **sym_mat, int num_of_elements){
    int i=0;
    double *matrix_1d;
    double **matrix;
    /*memory allocation for the diagonal degree matrix*/
    matrix_1d = calloc(num_of_elements*num_of_elements, sizeof(double));
    if(matrix_1d == NULL){
        printf("An Error Has Occurred\n");
        free(sym_mat);
        exit(1);
    }
    matrix = calloc(num_of_elements,sizeof(double *));
    if(matrix == NULL){
        printf("An Error Has Occurred\n");
        free(sym_mat);
        free(matrix_1d);
        exit(1);
    }
    for(i=0; i<num_of_elements; i++)
    {
        matrix[i] = matrix_1d+i*num_of_elements;
    } 

    /*computation of the diagonal degree matrix*/
    for(i=0; i<num_of_elements; i++){
        matrix[i][i]=sum_of_row(sym_mat, i, num_of_elements);
    }
    return matrix;
}

double** ddg(double **X, int num_of_elements, int d){
    double **sym_mat;
    double **dd_mat;
    sym_mat=(sym(X, num_of_elements, d));
    dd_mat=(diagonal_degree_mat(sym_mat, num_of_elements));
    return dd_mat;
}

double** normalized_mat(double **sym_mat, double **dd_mat, int num_of_elements){
    double **mult_mat; /*D^-0.5*/
    double **mult_on_rows; /*D^-0.5*A*/
    double **mult_on_cols; /*D^-0.5*A*D^-0.5*/
    int i=0;
    int j=0;
    mult_mat = dd_mat;
    mult_on_rows = sym_mat;
    for(i=0; i<num_of_elements; i++){
        mult_mat[i][i]= pow(mult_mat[i][i], -0.5);
    }

    for(i=0; i<num_of_elements; i++){
        for(j=0; j<num_of_elements; j++){
            mult_on_rows[i][j]*=mult_mat[i][i];
        }
    }

    mult_on_cols=mult_on_rows;
    for(i=0; i<num_of_elements; i++){
        for(j=0; j<num_of_elements; j++){
            mult_on_cols[j][i]*=mult_mat[i][i];
        }
    }

    return mult_on_cols;
}

double** update_H(double **H, double **W, int k, int num_of_elements){
    double **numerator_mat;
    double **denominator_mat;
    double **new_H;
    double *new_H_1d;
    int i, j;
    numerator_mat = multiply_matrices(W,num_of_elements,num_of_elements,H,k);
    denominator_mat = mult_by_transpose(H,num_of_elements,k);
    denominator_mat = multiply_matrices(denominator_mat,num_of_elements,num_of_elements,H,k);
    
    /*memory allocation for the new H*/
    new_H_1d = calloc(num_of_elements*k, sizeof(double));
    if(new_H_1d == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    new_H = calloc(num_of_elements,sizeof(double *));
    if(new_H == NULL){
        printf("An Error Has Occurred\n");
        free(new_H_1d);
        exit(1);
    }
    for(i=0; i<num_of_elements; i++)
    {
        matrix[i] = matrix_1d+i*k;
    }
    
    /*new H calculation*/ 
    for(i=0; i<num_of_elements; i++){
        for(j=0; j<k; j++){
            new_H[i][j]= H[i][j]*(0.5+0.5*numerator_mat[i][j]/denominator_mat[i][j]);
        }
    }
    free(numerator_mat);
    free(denominator_mat);
    return new_H;
} 

