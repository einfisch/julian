#include <Python.h>
#include <iostream>
int
call_python(double x, double y, double z, double epsilon)
{
    PyObject *pName, *pModule, *pFunc;
    PyObject *pArgs, *pValue;
    int i;

/*
    if (argc < 3) {
        fprintf(stderr,"Usage: call pythonfile funcname [args]\n");
        return 1;
    }
*/
    const char python_filename[] = "fiber_vector";
    const char python_funcname[] = "get_fiber_values";

    Py_Initialize();
    /* python script to be loaded is in the same folder */
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append(\".\")");
    pName = PyString_FromString(python_filename);
    /* Error checking of pName left out */

    pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule != NULL) {
        pFunc = PyObject_GetAttrString(pModule, python_funcname);
        /* pFunc is a new reference */

        double argv[] = {0, 0,0, x, y, z, epsilon};
        if (pFunc && PyCallable_Check(pFunc)) {
            pArgs = PyTuple_New(4);
            for (i = 0; i < 4; ++i) {
                pValue = PyFloat_FromDouble(argv[i + 3]);
                if (!pValue) {
                    Py_DECREF(pArgs);
                    Py_DECREF(pModule);
                    fprintf(stderr, "Cannot convert argument\n");
                    return 1;
                }
                /* pValue reference stolen here: */
                PyTuple_SetItem(pArgs, i, pValue);
            }
            pValue = PyObject_CallObject(pFunc, pArgs);
            Py_DECREF(pArgs);
            if (pValue != NULL) {
                //printf("Result of call: %ld\n", PyList_AsTuple(pValue));//PyInt_AsLong
                std::cout << PyList_AsTuple(pValue);
                Py_DECREF(pValue);
            }
            else {
                Py_DECREF(pFunc);
                Py_DECREF(pModule);
                PyErr_Print();
                fprintf(stderr,"Call failed\n");
                return 1;
            }
        }
        else {
            if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "Cannot find function \"%s\"\n", python_funcname);
        }
        Py_XDECREF(pFunc);
        Py_DECREF(pModule);
    }
    else {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"%s\"\n", python_filename);
        return 1;
    }
    Py_Finalize();
    return 0;
}

int main(int argc, char* argv[]) {
    call_python(atof(argv[3]), atof(argv[4]), atof(argv[5]), atof(argv[6]));
    return 0;
}