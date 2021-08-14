//
// Created by TomBarzilay on 22/07/2021.
//
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "spkmeans.h"



/* ************  API FUNCTION CONVERT PYTHON OBJECT INTO C ARRAY AND BACK ************ */
static void check_for_py_list(PyObject * _list)
{
    assert (!PyList_Check(_list) && "this is not a list ");
}

static void free_observations(double ** list, Py_ssize_t n, Py_ssize_t d )
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

static PyObject* get_py_lst_from_c_matrix(double ** c_observations ,  Py_ssize_t n, Py_ssize_t d)
{

    Py_ssize_t i;
    PyObject * observations;
    PyObject* dot;

    observations = PyList_New(n);
    for (i = 0; i < n; i++)
    {
        dot = get_py_lst_from_c_array(c_observations[i],d);
        PyList_SetItem(observations,i,dot);
    }
    // TO DO FREE MEMORY
    return observations;

}


static  double* get_c_array_from_py_lst(PyObject * list, Py_ssize_t  d)
{
    Py_ssize_t i;
    PyObject *item;
    /* NEVER EVER USE malloc/calloc/realloc or free on PyObject */
    double *c_dot =(double *) calloc(d, sizeof(double));
    assert(c_dot != NULL && "Problem in  load_python_list_to_c_array");
    for (i = 0; i < d; i++) {
        item = PyList_GetItem(list, i); /* DON'T FREE - cause problems */
        c_dot[i] =  PyFloat_AsDouble((item));

        if (c_dot[i]  == -1 && PyErr_Occurred()){
            puts("get_c_array_from_py_lst - ERROR");
            free(c_dot);
            return c_dot;
        }

    }

    return c_dot;
}
static  int* get_int_c_array_from_py_lst(PyObject * list, Py_ssize_t  d)
{
    Py_ssize_t i;
    PyObject *item;
    /* NEVER EVER USE malloc/calloc/realloc or free on PyObject */
    int*c_dot = malloc(sizeof(double) * d);
    assert(c_dot != NULL && "Problem in  load_python_list_to_c_array");
    for (i = 0; i < d; i++) {

        item = PyList_GetItem(list, i); /* DON'T FREE - cause problems */
        c_dot[i] = (int) PyFloat_AsDouble((item));

        if (c_dot[i]  == -1 && PyErr_Occurred()){
            puts("Something bad ... 2");
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
    double ** c_observations;
    c_observations = (double **) calloc( n, sizeof(double *));
    assert(c_observations != NULL && "Problem in  load_python_list_to_c_array");

    for (i = 0; i < n; i++)
    {
        item = PyList_GetItem(_list, i);
        check_for_py_list(item);
        c_observations[i] = get_c_array_from_py_lst(item,d);

    }
    return c_observations;

}

/* ************  API FUNCTION CONVERT PYTHON OBJECT INTO C ARRAY AND BACK ************ */


/*
 * API functions
 */


/*
 * Print list of lists of ints without changing it
 */
static PyObject* get_flag(PyObject *self, PyObject *args){
   char* goal;
    int k;
    PyObject* _observations, * _T;
    PyObject* T_K;
    double** observations;
    int n;
    int d;
    spk_results res;
    T_K = PyList_New( (Py_ssize_t) 2 );
    _T =  Py_None;
    if(!PyArg_ParseTuple(args, "siO",&goal,&k ,&_observations)) {//getting data from python
        printf("An Error Has Occured");
        return NULL;
    }
    check_for_py_list(_observations);//checking to see if the python object is actually a list
    n = (int) PyList_Size(_observations);
    d = (int)PyList_Size(PyList_GetItem(_observations, 0));

    observations = get_c_matrix_from_py_lst(_observations, n, d);
    assert_goal(goal);


     res  = activate_flag(goal,observations,k, n, d);

    // printf("C fit: n=%d,  k=%d\n",res.T_size, res.k); /* TODO DELETE */
     _T = get_py_lst_from_c_matrix(res.T ,res.T_size, res.k);
     PyList_SetItem(T_K,0, _T );
     PyList_SetItem(T_K,1,Py_BuildValue("i",res.k));

     free_matrix(observations, n);
     return T_K;
    }



static PyObject* fit(PyObject *self, PyObject *args)
{
    PyObject *_T, *_cluster_list, *_cluster_index_list;
    Py_ssize_t n, k;

    int n_c, k_c;
    double** T_c;
    double ** cluster_list_c;
    int* cluster_index_list_c ;


    if(!PyArg_ParseTuple(args, "OOO",&_T, &_cluster_list, &_cluster_index_list)) { return NULL;}

    /* Check that we got lists */
    check_for_py_list(_T);
    check_for_py_list(_cluster_list);
    check_for_py_list(_cluster_index_list);

    n = PyList_Size(_T);
    k = PyList_Size(PyList_GetItem(_T, 0));
     // printf("C fit: n=%d,  k=%d\n",n, k); /* TODO DELETE */

    /*  create c arrays from python lists */
    T_c = get_c_matrix_from_py_lst(_T,n,k);
    cluster_list_c = get_c_matrix_from_py_lst(_cluster_list,k,k);
    cluster_index_list_c = get_int_c_array_from_py_lst(_cluster_index_list,k);
    n_c = (int) n;
    k_c = (int) k;

    simple_kmean(T_c, cluster_list_c , cluster_index_list_c, n_c, k_c, k_c);
    _cluster_list =  get_py_lst_from_c_matrix(cluster_index_list_c ,k,k);

    /* free memory */
    free_matrix(cluster_list_c,k_c);
    free_matrix(T_c,n_c);
    free(cluster_index_list_c );
    return _cluster_list;

}

/*
 * A macro to help us with defining the methods
 * Compare with: {"f1", (PyCFunction)f1, METH_NOARGS, PyDoc_STR("No input parameters")}
*/
#define FUNC(_flag, _name, _docstring) { #_name, (PyCFunction)_name, _flag, PyDoc_STR(_docstring) }

static PyMethodDef _methods[] = {
        FUNC(METH_VARARGS, fit, "Print list of lists of ints without changing it"),
        {"get_flag",(PyCFunction)get_flag,METH_VARARGS,PyDoc_STR("recive the goal as an argument and do the requested operation")},
        {NULL, NULL, 0, NULL}
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



