//
// Created by TomBarzilay on 22/07/2021.
//
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "spkmeans.c"



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
    double *c_dot = malloc(sizeof(double) * d);
    assert(c_dot != NULL && "Problem in  load_python_list_to_c_array");
        printf("\n");
    for (i = 0; i < d; i++) {
        printf("here \n");
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
static  int* get_int_c_array_from_py_lst(PyObject * list, Py_ssize_t  d)
{
    Py_ssize_t i;
    PyObject *item;
    /* NEVER EVER USE malloc/calloc/realloc or free on PyObject */
    int *c_dot = malloc(sizeof(double) * d);
    assert(c_dot != NULL && "Problem in  load_python_list_to_c_array");
    for (i = 0; i < d; i++) {

        item = PyList_GetItem(list, i); /* DON'T FREE - cause problems */
        c_dot[i] =  PyLong_AsLong((item));

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
    double ** c_observations;
    c_observations = (double **) malloc( n* sizeof(double *));
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
    PyObject* _observations;
    double** observations;
    int n;
    int d;
    spk_results res;
    if(!PyArg_ParseTuple(args, "siO",&goal,&k ,&_observations)) {//getting data from python
        printf("An Error Has Occured");
        return NULL;
    }

    check_for_py_list(_observations);//checking to see if the python object is actually a list
    n = (int) PyList_Size(_observations);
    d = (int)PyList_Size(PyList_GetItem(_observations, 0));
    observations = get_c_matrix_from_py_lst(_observations,n,d);
    assert_goal(goal);
     PyObject* T_K = PyList_New((Py_ssize_t)2);
     res  = activate_flag(goal,observations,k,n,d);
     if(is_goal("spk")){
     PyList_SetItem(T_K,0, get_py_lst_from_c_matrix(res.mat,n,d));
     PyList_SetItem(T_K,1,Py_BuildValue("i",k));
     }



     // TODO free over
       free_matrix(observations,n);
       return T_K;
    }



static PyObject* fit(PyObject *self, PyObject *args)
{
    PyObject *T,*_cluster_list,*_cluster_index_list,*_observations;
    Py_ssize_t n,k, d;
    int n_c, k_c, d_c;
    double** T_c;
    double ** c_observations;
    double ** c_cluster_list;
    int * c_cluster_index_list;

    if(!PyArg_ParseTuple(args, "OOOO",&T_c ,&_cluster_list,&_cluster_index_list, &_observations)) { return NULL;}

    /* Check that we got lists */
    check_for_py_list(T);
    check_for_py_list(_observations);
    check_for_py_list(_cluster_list);
    check_for_py_list(_cluster_index_list);

    n = PyList_Size(_observations);
    k = PyList_Size(_cluster_list);
    d = PyList_Size(PyList_GetItem(_observations, 0));
    /*  create c arrays from python lists */
    T_c = get_c_matrix_from_py_lst(T,n,k);
    c_cluster_list = get_c_matrix_from_py_lst(_cluster_list,k,d);
    c_cluster_index_list = get_int_c_array_from_py_lst(_cluster_index_list,k);
    c_observations = get_c_matrix_from_py_lst(_observations,n,d);
    n_c=(int) n;
    k_c=(int) k;
    d_c=(int) d;
    simple_kmean(T_c, c_cluster_list , c_cluster_index_list, c_observations, n_c, k_c, d_c);//double ** T_mat, double ** T_cluster_list, int * cluster_index_list,double ** observations, int n, int k, int d

    _cluster_list =  get_py_lst_from_c_matrix(c_cluster_list ,k,d);

    /* free memory */
    free_observations(c_observations,n_c,d_c);
    free_observations(c_cluster_list,k_c,d_c);
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



