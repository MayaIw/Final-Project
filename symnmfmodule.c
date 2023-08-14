# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# define PY_SSIZE_T_CLEAN
# include <Python.h>
# include "symnmf.h"

static PyObject* sym(PyObject *self, PyObject *args){
    PyObject *elements_lst;
    PyObject *element;
    double element_coord;
    double **elements;
    double **returned_mat;
    int i, j;
    int num_of_elements;
    int d;

    if(!PyArg_ParseTuple(args, "Oii", &elements_lst,  &num_of_elements, &d)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                    PyObject* so it is used to signal that an error has occurred. */
    }

    /*memory allocation for all the points in the file*/
    elements = matrix_allocation(num_of_elements, d);

    /*reading the data points from python and passing into c matrix*/
    for(i=0; i<num_of_elements; i++){
        element = PyList_GetItem(elements_lst, i);
        for(j=0; j<d; j++){
            element_coord = PyFloat_AsDouble(PyList_GetItem(element, j)); 
            elements[i][j] = element_coord;
        }
    }

    returned_mat = sym_c(elements, num_of_elements, d);

    free_matrix(elements);

    PyObject* matrix;
    PyObject* vector;
    PyObject* python_float;
    matrix = PyList_New(num_of_elements);
    for (int i = 0; i < num_of_elements; i++)
    {   
        vector = PyList_New(num_of_elements);
        for(j=0; j<num_of_elements; j++){
            python_float = PyFloat_FromDouble(returned_mat[i][j]);
            PyList_SetItem(vector, j, python_float); 
        }
        PyList_SetItem(matrix, i, vector);
    }

    free_matrix(returned_mat);

    return Py_BuildValue("O", matrix); 
}

static PyObject* ddg(PyObject *self, PyObject *args){
    PyObject *elements_lst;
    PyObject *element;
    double element_coord;
    double **elements;
    double **returned_mat;
    int i, j;
    int num_of_elements;
    int d;

    if(!PyArg_ParseTuple(args, "Oii", &elements_lst,  &num_of_elements, &d)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                    PyObject* so it is used to signal that an error has occurred. */
    }

    /*memory allocation for all the points in the file*/
    elements = matrix_allocation(num_of_elements, d);

    /*reading the data points from python and passing into c matrix*/
    for(i=0; i<num_of_elements; i++){
        element = PyList_GetItem(elements_lst, i);
        for(j=0; j<d; j++){
            element_coord = PyFloat_AsDouble(PyList_GetItem(element, j)); 
            elements[i][j] = element_coord;
        }
    }

    returned_mat = ddg_c(elements, num_of_elements, d);

    free_matrix(elements);

    PyObject* matrix;
    PyObject* vector;
    PyObject* python_float;
    matrix = PyList_New(num_of_elements);
    for (int i = 0; i < num_of_elements; i++)
    {   
        vector = PyList_New(num_of_elements);
        for(j=0; j<num_of_elements; j++){
            python_float = PyFloat_FromDouble(returned_mat[i][j]);
            PyList_SetItem(vector, j, python_float); 
        }
        PyList_SetItem(matrix, i, vector);
    }

    free_matrix(returned_mat);

    return Py_BuildValue("O", matrix); 
}

static PyObject* norm(PyObject *self, PyObject *args){
    PyObject *elements_lst;
    PyObject *element;
    double element_coord;
    double **elements;
    double **returned_mat;
    int i, j;
    int num_of_elements;
    int d;

    if(!PyArg_ParseTuple(args, "Oii", &elements_lst,  &num_of_elements, &d)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                    PyObject* so it is used to signal that an error has occurred. */
    }

    /*memory allocation for all the points in the file*/
    elements = matrix_allocation(num_of_elements, d);

    /*reading the data points from python and passing into c matrix*/
    for(i=0; i<num_of_elements; i++){
        element = PyList_GetItem(elements_lst, i);
        for(j=0; j<d; j++){
            element_coord = PyFloat_AsDouble(PyList_GetItem(element, j)); 
            elements[i][j] = element_coord;
        }
    }

    returned_mat = norm_c(elements, num_of_elements, d);

    free_matrix(elements);

    PyObject* matrix;
    PyObject* vector;
    PyObject* python_float;
    matrix = PyList_New(num_of_elements);
    for (int i = 0; i < num_of_elements; i++)
    {   
        vector = PyList_New(num_of_elements);
        for(j=0; j<num_of_elements; j++){
            python_float = PyFloat_FromDouble(returned_mat[i][j]);
            PyList_SetItem(vector, j, python_float); 
        }
        PyList_SetItem(matrix, i, vector);
    }

    free_matrix(returned_mat);

    return Py_BuildValue("O", matrix); 
}

static PyObject* symnmf(PyObject *self, PyObject *args){
    PyObject *H_lst;
    PyObject *H;
    PyObject *W_lst;
    PyObject *W;
    PyObject *H_float;
    double H_entry;
    double **H_c;
    double W_entry;
    double **W_c;
    double **returned_mat;
    int i, j;
    int num_of_elements;
    int k;

    if(!PyArg_ParseTuple(args, "OOii", &H_lst, &W_lst, &k, &num_of_elements)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                    PyObject* so it is used to signal that an error has occurred. */
    }

    /*memory allocation for H*/
    H_c = matrix_allocation(num_of_elements, k);

    /*memory allocation for W*/
    W_c = matrix_allocation(num_of_elements, num_of_elements);

    /*reading H from python and passing into c matrix*/
    for(i=0; i<num_of_elements; i++){
        H = PyList_GetItem(H_lst, i);
        for(j=0; j<k; j++){
            H_float = PyList_GetItem(H, j);
            H_entry = PyFloat_AsDouble(H_float); 
            H_c[i][j] = H_entry;
        }
    }

    /*reading W from python and passing into c matrix*/
    for(i=0; i<num_of_elements; i++){
        W = PyList_GetItem(W_lst, i);
        for(j=0; j<k; j++){
            W_entry = PyFloat_AsDouble(PyList_GetItem(W, j)); 
            W_c[i][j] = W_entry;
        }
    }

    returned_mat = symnmf_c(H_c, W_c, k, num_of_elements);

    free_matrix(W_c);

    PyObject* matrix;
    PyObject* vector;
    PyObject* python_float;
    matrix = PyList_New(num_of_elements);
    for (int i = 0; i < num_of_elements; i++)
    {   
        vector = PyList_New(k);
        for(j=0; j<k; j++){
            python_float = PyFloat_FromDouble(returned_mat[i][j]);
            PyList_SetItem(vector, j, python_float); 
        }
        PyList_SetItem(matrix, i, vector);
    }

    free_matrix(returned_mat);

    return Py_BuildValue("O", matrix); 
}

static PyMethodDef symnmfMethods[] = {
    {"sym",                   /* the Python method name that will be used */
      (PyCFunction) sym, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parameters accepted for this function */
      PyDoc_STR("A geometric series up to n. sum_up_to_n(z^n)")}, /*  The docstring for the function */
      {"ddg",                   /* the Python method name that will be used */
      (PyCFunction) ddg, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parameters accepted for this function */
      PyDoc_STR("A geometric series up to n. sum_up_to_n(z^n)")}, /*  The docstring for the function */
      {"norm",                   /* the Python method name that will be used */
      (PyCFunction) norm, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parameters accepted for this function */
      PyDoc_STR("A geometric series up to n. sum_up_to_n(z^n)")}, /*  The docstring for the function */
      {"symnmf",                   /* the Python method name that will be used */
      (PyCFunction) symnmf, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parameters accepted for this function */
      PyDoc_STR("A geometric series up to n. sum_up_to_n(z^n)")}, /*  The docstring for the function */
    {NULL, NULL, 0, NULL}     /* The last entry must be all NULL */
};

static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "mysymnmf", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    symnmfMethods /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC PyInit_mysymnmf(void){
    PyObject *m;
    m = PyModule_Create(&symnmfmodule);
    if (!m) {
        return NULL;
    }
    return m;
}