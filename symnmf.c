# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include "symnmf.h"

double** matrix_allocation(int num_of_rows, int num_of_cols){
    int i=0;
    double *matrix_1d;
    double **matrix;

    matrix_1d = calloc(num_of_rows*num_of_cols, sizeof(double));
    if(matrix_1d == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    matrix = calloc(num_of_rows,sizeof(double *));
    if(matrix == NULL){
        printf("An Error Has Occurred\n");
        free(matrix_1d);
        exit(1);
    }
    for(i=0; i<num_of_rows; i++)
    {
        matrix[i] = matrix_1d+i*num_of_cols;
    } 
    return matrix;
}

void free_matrix(double **matrix){
    free(matrix[0]);
    free(matrix);
}

double sum_of_row(double **mat, int row_index, int num_of_cols){
    int j;
    double sum=0.0;
    for(j=0; j<num_of_cols; j++){
        sum += mat[row_index][j];
    }
    return sum;
}

double squared_distance(double *p, double *q, int d)
{
    int i;
    double sum_of_squares = 0.0;
    for(i=0; i<d; i++)
    {
        sum_of_squares += pow(p[i]-q[i],2);
    }
    return sum_of_squares;
}

double** multiply_matrices(double** mat1, int rows1, int cols1, double** mat2, int cols2) {
    double** result;
    int i,j,l;

    /*memory allocation for the multiplication matrix*/
    result = matrix_allocation(rows1, cols2);

    /*calculation*/
    for(i = 0; i < rows1; i++) {
        for(j = 0; j < cols2; j++) {
            for(l = 0; l < cols1; l++) {
                result[i][j] += mat1[i][l] * mat2[l][j];
            }
        }
    }
    return result;
}

double** mult_by_transpose(double **mat, int rows, int cols){
    double** result;
    int i,j,l;

    /*memory allocation for the multiplication matrix*/
    result = matrix_allocation(rows, rows);
    
    /*calculation*/
    for(i = 0; i < rows; i++) {
        for(j = 0; j < rows; j++) {
            for(l = 0; l < cols; l++) {
                result[i][j] += mat[i][l] * mat[j][l];
            }
        }
    }
    return result;
}

int check_convergence(double **new_H, double **H, int num_of_elements, int k){
    int i, j;
    double sum=0.0;
    for(i=0; i<num_of_elements; i++){
        for(j=0; j<k; j++){
            sum += pow(new_H[i][j]-H[i][j], 2);
        }
    }
    return (sum<0.0001);
}

double** sym_c(double **X, int num_of_elements, int d){
    int i=0;
    int j=0;
    double **matrix;

    /*memory allocation for the similarity matrix*/
    matrix = matrix_allocation(num_of_elements, num_of_elements);

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
    double **matrix;

    /*memory allocation for the diagonal degree matrix*/
    matrix = matrix_allocation(num_of_elements, num_of_elements);
    
    /*computation of the diagonal degree matrix*/
    for(i=0; i<num_of_elements; i++){
        matrix[i][i]=sum_of_row(sym_mat, i, num_of_elements);
    }
    return matrix;
}

double** ddg_c(double **X, int num_of_elements, int d){
    double **sym_mat;
    double **dd_mat;
    sym_mat=(sym_c(X, num_of_elements, d));
    dd_mat=(diagonal_degree_mat(sym_mat, num_of_elements));
    free_matrix(sym_mat);
    return dd_mat;
}

double** normalized_mat(double **sym_mat, double **dd_mat, int num_of_elements){
    double **mult_mat; /*will store D^-0.5*/
    double **mult_on_rows; /*will store D^-0.5*A*/
    double **mult_on_cols; /*will store D^-0.5*A*D^-0.5*/
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

double** norm_c(double **X, int num_of_elements, int d){
    double **sym_mat;
    double **dd_mat;
    double **norm_mat;
    sym_mat= sym_c(X, num_of_elements, d);
    dd_mat= diagonal_degree_mat(sym_mat, num_of_elements);
    norm_mat= normalized_mat(sym_mat, dd_mat, num_of_elements);
    free_matrix(dd_mat);
    return norm_mat;
}

double** update_H(double **H,double **H_alloc, double **W, int k, int num_of_elements){
    double **numerator_mat;
    double **first_denominator_mat;
    double **second_denominator_mat;
    int i, j;
    numerator_mat = multiply_matrices(W,num_of_elements,num_of_elements,H,k);
    first_denominator_mat = mult_by_transpose(H,num_of_elements,k);
    second_denominator_mat = multiply_matrices(first_denominator_mat,num_of_elements,num_of_elements,H,k);
    
    /*new H calculation*/ 
    for(i=0; i<num_of_elements; i++){
        for(j=0; j<k; j++){
            H_alloc[i][j]= H[i][j]*(0.5+0.5*numerator_mat[i][j]/second_denominator_mat[i][j]);
        }
    }
    free_matrix(numerator_mat);
    free_matrix(first_denominator_mat);
    free_matrix(second_denominator_mat);
    return H_alloc;
}

void printMatrix(double **matrix, int num_of_rows,  int num_of_cols){
    int i=0;
    int j=0;
    for(i=0; i<num_of_rows; i++){
        for(j=0; j<num_of_cols; j++){
            printf("%.4f", matrix[i][j]);
            if(j<num_of_cols-1){
                printf("%c", ',');
            }
        }
        printf("\n");
    }
    free_matrix(matrix);
}

double** symnmf_c(double **H, double **W, int k, int num_of_elements){
    int max_iter=300;
    int i, j, l;
    double **new_H;

    /*memory allocation for the new H*/
    new_H = matrix_allocation(num_of_elements, k);

    /*updating H until max iteration number is reached OR intil the convergence of H*/
    for(i=0; i<max_iter; i++){
        new_H = update_H(H, new_H, W, k, num_of_elements);
        if(check_convergence(new_H, H, num_of_elements, k)){
            break;
        }
        for(j=0; j<num_of_elements; j++){
            for(l=0; l<k; l++){
               H[j][l] = new_H[j][l]; 
            }
        }
    }

    return(new_H);
}

int main(int argc, char **argv){ 
    int num_of_elements=0;
    int d=1;
    double **elements;
    FILE *points;
    char *goal;
    char c, next_char;
    int num_rows, num_cols;
    char delimiter;

    if(argc!=3){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    
    points = fopen(argv[2], "r");
    if(points==NULL){
        printf("An Error Has Occurred\n");
        fclose(points);
        exit(1);
    }

    /*reading the first line, and counting the number of coordinates of a datapoint.*/
    while ((c = fgetc(points)) != EOF)
    {
        if(c==','){
            d+=1;
        }
        if(c=='\n'){
            num_of_elements += 1;
            break;
        }
    }

    num_of_elements += 1; /*counting the last line if there is no '\n' character at the end*/

    /*reading the rest pf the file, and counting the number of datapoints.*/
    while ((c = fgetc(points)) != EOF){
        if (c == '\n'){
            num_of_elements += 1;
        } 
    }
    
    /*returning to the begining of the file*/
    rewind(points);

    /*memory allocation for all the datapoints in the file*/
    elements = matrix_allocation(num_of_elements, d);

    /*put all the datapoints in one array*/
    num_rows = 0;
    while (num_rows < num_of_elements) {
        num_cols = 0;
        while (num_cols < d) {
            if (fscanf(points, "%lf", &elements[num_rows][num_cols]) == 1) {
                num_cols++;
            } else {
                break;  
            }
            delimiter = fgetc(points);
            if (delimiter == ',') {
                continue; 
            } else if (delimiter == '\n' || delimiter == EOF) {
                break;
            } else { /*Invalid delimiter*/
                printf("An Error Has Occurred\n");
                fclose(points);
                exit(1);
            }
        }
        num_rows++;
        next_char = fgetc(points);
        if (next_char == '\n' || feof(points)) {
            break;  
        }
        ungetc(next_char, points);
    }

    fclose(points);

    /*performing the given goal*/
    goal = argv[1];
    if(!strcmp(goal, "sym")){
        printMatrix(sym_c(elements, num_of_elements, d),num_of_elements, num_of_elements);
    }
    else if(!strcmp(goal, "ddg")){
        printMatrix(ddg_c(elements, num_of_elements, d),num_of_elements, num_of_elements);
    }
    else if(!strcmp(goal, "norm")){
        printMatrix(norm_c(elements, num_of_elements, d),num_of_elements, num_of_elements);
    }
    else{
        printf("An Error Has Occurred\n");
        free_matrix(elements);
        exit(1);
    }

    free_matrix(elements);

    return 0;
}