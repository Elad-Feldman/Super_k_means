//
// Created by TomBarzilay on 22/07/2021.
//
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "kmeans.c"



/* ************  API FUNCTION CONVERT PYTHON OBJECT INTO C ARRAY AND BACK ************ */
static void check_for_py_list(PyObject * _list)
{
    assert (!PyList_Check(_list) && "this is not a list ");
}

static void free_dots_list(double ** list, Py_ssize_t n, Py_ssize_t d )
{
    Py_ssize_t i;
    for (i = 0; i < n; i++)
        free(list[i]);

    free(list);
}

static  PyObject* get_py_lst_from_c_array(double* c_dot, Py_ssize_t  d)
{
    Py_ssize_t i;
    PyObject * dot ;
    PyObject* item;
    dot = PyList_New(d); /* https://docs.python.org/3/c-api/list.html */
    for (i = 0; i < d; i++) {
        item = PyFloat_FromDouble(c_dot[i]);
        PyList_SetItem(dot,i,item);
        // TO DO FREE MEMORY
    }
    return dot;
}

static PyObject* get_py_lst_from_c_matrix(double ** c_dot_list ,  Py_ssize_t n, Py_ssize_t d)
{
    Py_ssize_t i;
    PyObject * dot_list ;
    PyObject* dot;

    dot_list = PyList_New(n);
    for (i = 0; i < n; i++)
    {
        dot = get_py_lst_from_c_array(c_dot_list[i],d);
        PyList_SetItem(dot_list,i,dot);
    }
    // TO DO FREE MEMORY
    return dot_list;

}


static  double* get_c_array_from_py_lst(PyObject * list, Py_ssize_t  d)
{
    Py_ssize_t i;
    PyObject *item;
    /* NEVER EVER USE malloc/calloc/realloc or free on PyObject */
    double *c_dot = malloc(sizeof(double) * d);
    assert(c_dot != NULL && "Problem in  load_python_list_to_c_array");
    for (i = 0; i < d; i++) {
        item = PyList_GetItem(list, i); /* DON'T FREE - cause problems */
        c_dot[i] =  PyFloat_AsDouble((item));

        if (c_dot[i]  == -1 && PyErr_Occurred()){
            puts("Something bad ...");
            free(c_dot);
            return c_dot;
        }

    }

    return c_dot;
}


static double** get_c_matrix_from_py_lst(PyObject * _list,  Py_ssize_t n, Py_ssize_t d)
{
    check_for_py_list(_list);

    Py_ssize_t i;
    PyObject  *item;
    double ** c_dot_list;
    c_dot_list = (double **) malloc( n* sizeof(double *));
    assert(c_dot_list != NULL && "Problem in  load_python_list_to_c_array");

    for (i = 0; i < n; i++)
    {
        item = PyList_GetItem(_list, i);
        check_for_py_list(item);
        c_dot_list[i] = get_c_array_from_py_lst(item,d);
    }

    return c_dot_list;

}

/* ************  API FUNCTION CONVERT PYTHON OBJECT INTO C ARRAY AND BACK ************ */


/*
 * API functions
 */


/*
 * Print list of lists of ints without changing it
 */
static PyObject* get_flag(PyObject *self, PyObject *args){

}
static PyObject* fit(PyObject *self, PyObject *args)
{
    PyObject *_dots_list,*_cluster_list,*_cluster_index_list;
    Py_ssize_t n,k, d;
    int n_c, k_c, d_c;
    double ** c_dot_list;
    double ** c_cluster_list;
    double * c_cluster_index_list;

    if(!PyArg_ParseTuple(args, "OOO", &_dots_list,&_cluster_list,&_cluster_index_list)) { return NULL;}

    /* Check that we got lists */
    check_for_py_list(_dots_list);
    check_for_py_list(_cluster_list);
    check_for_py_list(_cluster_index_list);

    n = PyList_Size(_dots_list);
    k = PyList_Size(_cluster_list);
    d = PyList_Size(PyList_GetItem(_dots_list, 0));

    /*  create c arrays from python lists */
    c_dot_list = get_c_matrix_from_py_lst(_dots_list,n,d);
    c_cluster_list = get_c_matrix_from_py_lst(_cluster_list,k,d);
    c_cluster_index_list = get_c_array_from_py_lst(_cluster_index_list,k);

    n_c =(int) n;
    k_c =(int) k;
    d_c =(int) d;

    simple_kmean(c_dot_list, c_cluster_list ,c_cluster_index_list, n_c, k_c, d_c);

    _cluster_list =  get_py_lst_from_c_matrix(c_cluster_list ,k,d);

    /* free memory */
    free_dots_list(c_dot_list,n,d);
    free_dots_list(c_cluster_list,k,d);
    free(c_cluster_index_list);

    return _cluster_list;

}

/*
 * A macro to help us with defining the methods
 * Compare with: {"f1", (PyCFunction)f1, METH_NOARGS, PyDoc_STR("No input parameters")}
*/
#define FUNC(_flag, _name, _docstring) { #_name, (PyCFunction)_name, _flag, PyDoc_STR(_docstring) }

static PyMethodDef _methods[] = {
        FUNC(METH_VARARGS, fit, "Print list of lists of ints without changing it"),
        {NULL, NULL, 0, NULL}   /* sentinel */
};

static struct PyModuleDef _moduledef = {
        PyModuleDef_HEAD_INIT,
        "spkmeans",
        NULL,
        -1,
        _methods
};

PyMODINIT_FUNC
PyInit_spkmeans(void)
{
    return PyModule_Create(&_moduledef);
}



