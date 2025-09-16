#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "LIBRARY/for_cvode.h"
#define PARALLEL
#define CPU_NUM 10
#define NUM_CELLS 192

#ifdef PARALLEL
#include <omp.h>
#endif

#define D1 0.0001  //0.0002*6*0.001
#define D3 0.04    //diffusion for atrium
#define dx 0.015   // mm  0.075  0.04  // mm

/**** Gradient **************************************************************/

double Cm_gradient(int n)
{
    double Cm_SAN, c_cent, c_peri;

    c_cent = 0.025;
    c_peri = 0.05;

    Cm_SAN = c_cent + (c_peri - c_cent) * (1.0 / (1.0 + exp(-0.3 * ((double)(n - 117))))
                                           + 1.0 / (1.0 + exp(0.3 * ((double)(n - 59)))));
    return Cm_SAN;
}

double DIFFUSION(int n, MY_CELL **data)
{
    double dDiffusiondt;
    double diffusion, diffusion_former, diffusion_latter;
    double ddiffusion; // first derivative of diffusion.
    // set up the diffusion. Put in your diffusion gradient into diffusion[].

    diffusion = D1 + (D3 - D1) * (1.0 / (1.0 + exp(-0.25 * ((double)(n - 119))))
                                  + 1.0 / (1.0 + exp(0.25 * ((double)(n - 57)))));

    diffusion_former = D1 + (D3 - D1) * (1.0 / (1.0 + exp(-0.25 * ((double)(n - 1 - 119))))
                                         + 1.0 / (1.0 + exp(0.25 * ((double)(n - 1 - 57)))));

    diffusion_latter = D1 + (D3 - D1) * (1.0 / (1.0 + exp(-0.25 * ((double)(n + 1 - 119))))
                                         + 1.0 / (1.0 + exp(0.25 * ((double)(n + 1 - 57)))));

    /*  if(n>=(NUM_CELLS/2-10) && n<=(NUM_CELLS/2+10))
    {
        diffusion = D1 + (D2 - D1)*(1.0/(1.0+exp(-0.5*((double)(n-NUM_CELLS/2-10))))
                                   + 1.0/(1.0+exp(0.5*((double)(n-NUM_CELLS/2+10)))));

        diffusion_former = D1 + (D2 - D1)*(1.0/(1.0+exp(-0.5*((double)(n-1-NUM_CELLS/2-10))))
                                          + 1.0/(1.0+exp(0.5*((double)(n-1-NUM_CELLS/2+10)))));

        diffusion_latter = D1 + (D2 - D1)*(1.0/(1.0+exp(-0.5*((double)(n+1-NUM_CELLS/2-10))))
                                          + 1.0/(1.0+exp(0.5*((double)(n+1-NUM_CELLS/2+10)))));
    }
    else
    {
        diffusion = D1 + (D3 - D1)*(1.0/(1.0+exp(-0.6*((double)(n-NUM_CELLS/2-14))))
                                   + 1.0/(1.0+exp(0.6*((double)(n-NUM_CELLS/2+14)))));

        diffusion_former = D1 + (D3 - D1)*(1.0/(1.0+exp(-0.6*((double)(n-1-NUM_CELLS/2-14))))
                                          + 1.0/(1.0+exp(0.6*((double)(n-1-NUM_CELLS/2+14)))));

        diffusion_latter = D1 + (D3 - D1)*(1.0/(1.0+exp(-0.6*((double)(n+1-NUM_CELLS/2-14))))
                                          + 1.0/(1.0+exp(0.6*((double)(n+1-NUM_CELLS/2+14)))));
    }*/

    if (n == 0 || n == (NUM_CELLS - 1))
    {
        ddiffusion = 0.0;    // at the boundary, it does not matter what it is.
    }
    else
    {
        ddiffusion = (diffusion_latter - diffusion_former) / (2.0 * dx);
    }

    // set up the diffusion

    if (n > 0 && n < NUM_CELLS)
    {
        dDiffusiondt = diffusion * (data[n + 1]->V + data[n - 1]->V - 2.0 * data[n]->V) / (dx * dx)
                       + ddiffusion * (data[n + 1]->V - data[n - 1]->V) / (2.0 * dx);
    }
    if (n == 0)             dDiffusiondt = 2.0 * D3 / (dx * dx) * (data[n + 1]->V - data[n]->V);
    if (n == (NUM_CELLS - 1)) dDiffusiondt = 2.0 * D3 / (dx * dx) * (data[n - 1]->V - data[n]->V);

    return dDiffusiondt;
}

double DIFFUSION_OP(int n, double *array)
{
    double dDiffusiondt;
    double diffusion, diffusion_former, diffusion_latter;
    double ddiffusion; // first derivative of diffusion.
    // set up the diffusion. Put in your diffusion gradient into diffusion[].

    diffusion = D1 + (D3 - D1) * (1.0 / (1.0 + exp(-0.25 * ((double)(n - 119))))
                                  + 1.0 / (1.0 + exp(0.25 * ((double)(n - 57)))));

    diffusion_former = D1 + (D3 - D1) * (1.0 / (1.0 + exp(-0.25 * ((double)(n - 1 - 119))))
                                         + 1.0 / (1.0 + exp(0.25 * ((double)(n - 1 - 57)))));

    diffusion_latter = D1 + (D3 - D1) * (1.0 / (1.0 + exp(-0.25 * ((double)(n + 1 - 119))))
                                         + 1.0 / (1.0 + exp(0.25 * ((double)(n + 1 - 57)))));

    if (n == 0 || n == (NUM_CELLS - 1))
    {
        ddiffusion = 0.0;    // at the boundary, it does not matter what it is.
    }
    else
    {
        ddiffusion = (diffusion_latter - diffusion_former) / (2.0 * dx);
    }

    // set up the diffusion

    if (n > 0 && n < NUM_CELLS)
    {
        dDiffusiondt = diffusion * (array[n + 1] + array[n - 1] - 2.0 * array[n]) / (dx * dx)
                       + ddiffusion * (array[n + 1] - array[n - 1]) / (2.0 * dx);
    }
    if (n == 0)             dDiffusiondt = 2.0 * D3 / (dx * dx) * (array[n + 1] - array[n]);
    if (n == (NUM_CELLS - 1)) dDiffusiondt = 2.0 * D3 / (dx * dx) * (array[n - 1] - array[n]);

    return dDiffusiondt;
}

int main()
{
    omp_set_num_threads(CPU_NUM);
    FILE *output;
    FILE *cvv;
    int i, outputcounter = 0;
    double total_time, ddt, du, * new_para, * old_para, * t;
    int flag3, flag2, flag1;
    double t1, t2, t3, t4, cvsan, cvatrium;
    total_time = 1000.0;
    ddt = 0.005;
    clock_t time_begin, time_end;

    time_begin = clock();
    new_para = (double *)malloc((NUM_CELLS) * sizeof(double));
    old_para = (double *)malloc((NUM_CELLS) * sizeof(double));
    t = (double *)malloc((NUM_CELLS) * sizeof(double));

    flag3 = -1;
    flag2 = -1;
    flag1 = -1;
    t1 = -1000.0;

    realtype tout;
    N_Vector y[NUM_CELLS];

    for (int i = 0; i < NUM_CELLS; i++)
        y[i] = NULL;

    MY_CELL **cell;
    cell = (MY_CELL **)malloc((NUM_CELLS) * sizeof(MY_CELL *));
    for (i = 0; i < NUM_CELLS; i++)
    {
        cell[i] = (MY_CELL *)malloc(sizeof(MY_CELL));

        cell[i]->HCN = 0;

        if (i >= 52 && i <= 124)
            cell[i]->cell_type = 3;
        else
            cell[i]->cell_type = 2;

        if (cell[i]->cell_type >= 3)
        {
            cell[i]->capacitance = Cm_gradient(i);
            y[i] = N_VNew_Serial(39);
            initial_for_cvode(y[i], cell[i]);
            prepare_for_cvode(ddt, y[i], cell[i]);
        }
        else// if(cell[i]->cell_type == 2)
        {
            y[i] = N_VNew_Serial(33);
            initial_for_cvode(y[i], cell[i]);
            prepare_for_cvode(ddt, y[i], cell[i]);
        }

        new_para[i] = cell[i]->V;
        old_para[i] = cell[i]->V;
        t[i] = 0.0;
    }

#ifdef PARALLEL
    #pragma omp parallel \
    default(none) \
    shared(output, outputcounter, total_time, new_para, old_para, cell, t, y, ddt, \
           flag3, flag2, flag1, t1, t2, t3, t4, cvsan, cvv, cvatrium) \
    private(i, du, tout)
#endif
    {
        for (tout = ddt; tout < total_time; tout += ddt)
        {
#ifdef PARALLEL
            #pragma omp for schedule(static)
#endif
            for (i = 1; i < NUM_CELLS - 1; i++)
            {
                du = DIFFUSION_OP(i, old_para);
                new_para[i] = old_para[i] + (ddt / 2 ) * du;
                cell[i]->V = new_para[i];
            } // end of for loop

#ifdef PARALLEL
            #pragma omp for schedule(static)
#endif
            for (i = 1; i < NUM_CELLS - 1; i++)
            {
                CVode(cell[i]->cvode_mem, tout, y[i], &t[i], CV_NORMAL); // 1 time step solution.
                new_para[i] = new_para[i] + ddt * cell[i]->dvdt;
            } // end of for loop

            // measurements here: because you have old_para and new_para.
#ifdef PARALLEL
            #pragma omp single
#endif
            {
                if (tout > 600.0 && flag3 == -1)
                {
                    if (cell[88]->V > 0.0 && t1 == -1000.0)
                    {
                        t1 = tout;
                    }

                    if (cell[124]->V > 0.0 && flag1 == -1)
                    {
                        t2 = tout;
                        cvsan = 5.0 * dx / (t2 - t1);
                        flag1 = 0;
                    }

                    if (cell[135]->V > 0.0 && flag2 == -1)
                    {
                        t3 = tout;
                        flag2 = 0;
                    }
                    if (cell[175]->V > 0.0 && flag3 == -1)
                    {
                        t4 = tout;
                        cvatrium = 40.0 * dx / (t4 - t3);
                        flag3 = 0;
                    }
                    if (flag3 == 0)
                    {
                        cvv = fopen("cv.dat", "a+");
                        fprintf(cvv, "Propagation velocity in SAN is %f m/s \n", cvsan);
                        fprintf(cvv, "Propagation velocity in Atrium is %f m/s \n", cvatrium);
                        fprintf(cvv, "t1 = %f\n", t1);
                        fprintf(cvv, "t2 = %f\n", t2);
                        fprintf(cvv, "t3 = %f\n", t3);
                        fprintf(cvv, "t4 = %f\n", t4);
                        fclose(cvv);
                    }
                }
            }
            // renew the voltage
#ifdef PARALLEL
            #pragma omp for schedule(static)
#endif
            for (i = 1; i < NUM_CELLS - 1; i++)
            {
                du = DIFFUSION_OP(i, new_para);
                old_para[i] = new_para[i] + (ddt / 2 ) * du;
            } // end of for loop

#ifdef PARALLEL
            #pragma omp single
#endif
            {
                // at boundary
                cell[0]->V = cell[1]->V;
                old_para[0] = old_para[1];
                new_para[0] = new_para[1];
                cell[NUM_CELLS - 1]->V = cell[NUM_CELLS - 2]->V;
                old_para[NUM_CELLS - 1] = old_para[NUM_CELLS - 2];
                new_para[NUM_CELLS - 1] = new_para[NUM_CELLS - 2];

                outputcounter++;
                if (outputcounter % 200 == 0)
                {
                    // output the voltage.
                    output = fopen("1dresult.dat", "a+");
                    for ( i = 0; i < NUM_CELLS; i++)
                    {
                        fprintf(output, "%5.5f\t ", cell[i]->V);
                    }
                    fprintf(output, "\n");

                    fclose(output);

                    for ( i = 0; i < NUM_CELLS; i++) if (cell[i]->V != cell[i]->V)
                        {
                            printf("possible nans\n");
                            exit(1);
                        };
                } // end of output
            }
            #pragma omp barrier
        } // end of time loop
    }  //end of parallel

    for (i = 0; i < NUM_CELLS; i++)
    {
        N_VDestroy_Serial(y[i]);
        CVodeFree(&cell[i]->cvode_mem);
        free(cell[i]);
    }
    free(cell);
    free(old_para);
    free(new_para);

    time_end = clock();
    printf("CPU t1ime=%f seconds\n", ((double)(time_end - time_begin)) / (CLOCKS_PER_SEC * CPU_NUM));

    return 0;
}
