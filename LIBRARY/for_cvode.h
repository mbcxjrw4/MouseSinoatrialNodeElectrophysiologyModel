#ifndef FOR_CVODE_H
#define FOR_CVODE_H

#include "single_cell.h"
#include "initial_value.h"
#include <malloc.h>
#include <cvode/cvode_dense.h>            // explicitly include this header file
#include <cvode/cvode.h>                  // prototypes for CVODE fcts., consts.
#include <nvector/nvector_serial.h>       // serial N_Vector types, fcts., macros
#define Ith(v,i)    NV_Ith_S(v,i)         // Ith numbers components 0..NEQ-1
#define IJth(A,i,j) DENSE_ELEM(A,i,j)     // IJth numbers rows,cols 1..NEQ-1

void initial_for_cvode(N_Vector y, MY_CELL *data);
void prepare_for_cvode(double ddt, N_Vector y, MY_CELL *data);
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

void initial_for_cvode(N_Vector y, MY_CELL *data)
{
    int NEQ;
    if (data->cell_type >= 3)
        NEQ = 41;
    else if (data->cell_type == 2 || data->cell_type == 1)
        NEQ = 33;

    double *para;
    para = (double *)malloc((NEQ) * sizeof(double));

    initial(para, data);

    for (int i = 0; i < NEQ; i++)
        Ith(y, i) = para[i];

    free(para);
}


void prepare_for_cvode(double ddt, N_Vector y, MY_CELL *data)
{
    int NEQ;
    if (data->cell_type >= 3)
        NEQ = 41;
    else if (data->cell_type == 2 || data->cell_type == 1)
        NEQ = 33;

    data->cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);         // creat CVODE object
    CVodeInit(data->cvode_mem, f, 0.0, y);                  // initialize CVODE solver
    CVodeSStolerances(data->cvode_mem, 10e-6, 1e-6);          // specify intergration tolerances
    CVDense(data->cvode_mem, NEQ);                            // attach liner solver module
    CVodeSetMaxStep(data->cvode_mem, ddt);
    CVodeSetUserData(data->cvode_mem, data);
}

/* Functions Called by the Solver ******************************************************/

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    MY_CELL *data = (MY_CELL *)user_data;

    int NEQ;
    if (data->cell_type >= 3)
        NEQ = 41;
    else if(data->cell_type == 2 || data->cell_type ==1) 
        NEQ = 33;

    double *para;
    double *rate;
    para = (double *)malloc((NEQ) * sizeof(double));
    rate = (double *)malloc((NEQ) * sizeof(double));

    for (int i = 0; i < NEQ; i++)
    {
        para[i] = Ith(y, i);
        rate[i] = 0.0;
    }

    if (data->cell_type >= 3)
        SAN(para, rate, data);
    else if(data->cell_type == 2 || data->cell_type ==1) 
        ATRIUM(para, rate, data);

    for (int i = 0; i < NEQ; i++)
        Ith(ydot, i) = rate[i];

    free(para);
    free(rate);

    return 0;
}

#endif
