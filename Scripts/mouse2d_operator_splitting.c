#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "LIBRARY/for_cvode.h"
#define PARALLEL
#define CPU_NUM 48

#define X 192
#define Z 125
#define dx 0.015 //0.04

/*double DIFFUSION(int x, int z, double **D, MY_CELL *** data)
{
    int cell_type = data[x][z]->cell_type;
    int x_latter = data[x + 1][z]->cell_type;
    int x_former = data[x - 1][z]->cell_type;
    int z_latter = data[x][z + 1]->cell_type;
    int z_former = data[x][z - 1]->cell_type;

    double dudx2, dudz2, dudx, dudz;
    double dDiffusiondt;
    double diffusion;
    double ddiffusion_x, ddiffusion_z; // first derivative of diffusion.

    // set up the diffusion_coefficient. Put in your diffusion gradient into diffusion.

    diffusion = D[x][z];

    ddiffusion_x = (D[x + 1][z] - D[x - 1][z]) / (2.0 * dx);
    ddiffusion_z = (D[x][z + 1] - D[x][z - 1]) / (2.0 * dx);

    //set up the dudx2, dudz2, dudx

    dudx2 = (data[x - 1][z]->V + data[x + 1][z]->V - 2.0 * data[x][z]->V) / (dx * dx);
    dudz2 = (data[x][z - 1]->V + data[x][z + 1]->V - 2.0 * data[x][z]->V) / (dx * dx);
    dudx  = (data[x + 1][z]->V - data[x - 1][z]->V) / (2.0 * dx);
    dudz  = (data[x][z + 1]->V - data[x][z - 1]->V) / (2.0 * dx);

    if (x_latter == 0)
    {
        dudx2 = 2.0 * (data[x - 1][z]->V - data[x][z]->V) / (dx * dx);
        dudx = 0.0;
    }
    if (x_former == 0)
    {
        dudx2 = 2.0 * (data[x + 1][z]->V - data[x][z]->V) / (dx * dx);
        dudx = 0.0;
    }
    if (x_latter == 0 && x_former == 0)
    {
        dudx2 = 0.0;
    }

    if (z_latter == 0)
    {
        dudz2 = 2.0 * (data[x][z - 1]->V - data[x][z]->V) / (dx * dx);
        dudz = 0.0;
    }
    if (z_former == 0)
    {
        dudz2 = 2.0 * (data[x][z + 1]->V - data[x][z]->V) / (dx * dx);
        dudz = 0.0;
    }
    if (z_latter == 0 && z_former == 0)
    {
        dudz2 = 0.0;
    }

    if (cell_type >= 3)
    {
        if (x_latter == 2)
        {
            dudx2 = 2.0 * (data[x - 1][z]->V - data[x][z]->V) / (dx * dx);
            dudx = 0.0;
        }
        if (x_former == 2)
        {
            dudx2 = 2.0 * (data[x + 1][z]->V - data[x][z]->V) / (dx * dx);
            dudx = 0.0;
        }
        if (x_latter == 2 && x_former == 2)
        {
            dudx2 = 0.0;
        }

        if (z_latter == 2)
        {
            dudz2 = 2.0 * (data[x][z - 1]->V - data[x][z]->V) / (dx * dx);
            dudz = 0.0;
        }
        if (z_former == 2)
        {
            dudz2 = 2.0 * (data[x][z + 1]->V - data[x][z]->V) / (dx * dx);
            dudz = 0.0;
        }
        if (z_latter == 2 && z_former == 2)
        {
            dudz2 = 0.0;
        }
    }

    if (cell_type == 2)
    {
        if (x_latter >= 3)
        {
            dudx2 = 2.0 * (data[x - 1][z]->V - data[x][z]->V) / (dx * dx);
            dudx = 0.0;
        }
        if (x_former >= 3)
        {
            dudx2 = 2.0 * (data[x + 1][z]->V - data[x][z]->V) / (dx * dx);
            dudx = 0.0;
        }
        if (x_latter >= 3 && x_former >= 3)
        {
            dudx2 = 0.0;
        }

        if (z_latter >= 3)
        {
            dudz2 = 2.0 * (data[x][z - 1]->V - data[x][z]->V) / (dx * dx);
            dudz = 0.0;
        }
        if (z_former >= 3)
        {
            dudz2 = 2.0 * (data[x][z + 1]->V - data[x][z]->V) / (dx * dx);
            dudz = 0.0;
        }
        if (z_latter >= 3 && z_former >= 3)
        {
            dudz2 = 0.0;
        }
    }

    // set up the diffusion

    dDiffusiondt = diffusion * (dudx2 + dudz2) + ddiffusion_x * dudx + ddiffusion_z * dudz;

    return dDiffusiondt;
}*/

double DIFFUSION_OP(int x, int z, double **D, MY_CELL *** data, double **array)
{
    int cell_type = data[x][z]->cell_type;
    int x_latter = data[x + 1][z]->cell_type;
    int x_former = data[x - 1][z]->cell_type;
    int z_latter = data[x][z + 1]->cell_type;
    int z_former = data[x][z - 1]->cell_type;

    double dudx2, dudz2, dudx, dudz;
    double dDiffusiondt;
    double diffusion;
    double ddiffusion_x, ddiffusion_z; // first derivative of diffusion.

    // set up the diffusion_coefficient. Put in your diffusion gradient into diffusion.

    diffusion = D[x][z];

    ddiffusion_x = (D[x + 1][z] - D[x - 1][z]) / (2.0 * dx);
    ddiffusion_z = (D[x][z + 1] - D[x][z - 1]) / (2.0 * dx);

    //set up the dudx2, dudz2, dudx

    dudx2 = (array[x - 1][z] + array[x + 1][z] - 2.0 * array[x][z]) / (dx * dx);
    dudz2 = (array[x][z - 1] + array[x][z + 1] - 2.0 * array[x][z]) / (dx * dx);
    dudx  = (array[x + 1][z] - array[x - 1][z]) / (2.0 * dx);
    dudz  = (array[x][z + 1] - array[x][z - 1]) / (2.0 * dx);

    if (x_latter == 0)
    {
        dudx2 = 2.0 * (array[x - 1][z] - array[x][z]) / (dx * dx);
        dudx = 0.0;
    }
    if (x_former == 0)
    {
        dudx2 = 2.0 * (array[x + 1][z] - array[x][z]) / (dx * dx);
        dudx = 0.0;
    }
    if (x_latter == 0 && x_former == 0)
    {
        dudx2 = 0.0;
    }

    if (z_latter == 0)
    {
        dudz2 = 2.0 * (array[x][z - 1] - array[x][z]) / (dx * dx);
        dudz = 0.0;
    }
    if (z_former == 0)
    {
        dudz2 = 2.0 * (array[x][z + 1] - array[x][z]) / (dx * dx);
        dudz = 0.0;
    }
    if (z_latter == 0 && z_former == 0)
    {
        dudz2 = 0.0;
    }

    if (cell_type >= 3)
    {
        if (x_latter == 2)
        {
            dudx2 = 2.0 * (data[x - 1][z]->V - data[x][z]->V) / (dx * dx);
            dudx = 0.0;
        }
        if (x_former == 2)
        {
            dudx2 = 2.0 * (data[x + 1][z]->V - data[x][z]->V) / (dx * dx);
            dudx = 0.0;
        }
        if (x_latter == 2 && x_former == 2)
        {
            dudx2 = 0.0;
        }

        if (z_latter == 2)
        {
            dudz2 = 2.0 * (data[x][z - 1]->V - data[x][z]->V) / (dx * dx);
            dudz = 0.0;
        }
        if (z_former == 2)
        {
            dudz2 = 2.0 * (data[x][z + 1]->V - data[x][z]->V) / (dx * dx);
            dudz = 0.0;
        }
        if (z_latter == 2 && z_former == 2)
        {
            dudz2 = 0.0;
        }
    }

    if (cell_type == 2)
    {
        if (x_latter >= 3)
        {
            dudx2 = 2.0 * (data[x - 1][z]->V - data[x][z]->V) / (dx * dx);
            dudx = 0.0;
        }
        if (x_former >= 3)
        {
            dudx2 = 2.0 * (data[x + 1][z]->V - data[x][z]->V) / (dx * dx);
            dudx = 0.0;
        }
        if (x_latter >= 3 && x_former >= 3)
        {
            dudx2 = 0.0;
        }

        if (z_latter >= 3)
        {
            dudz2 = 2.0 * (data[x][z - 1]->V - data[x][z]->V) / (dx * dx);
            dudz = 0.0;
        }
        if (z_former >= 3)
        {
            dudz2 = 2.0 * (data[x][z + 1]->V - data[x][z]->V) / (dx * dx);
            dudz = 0.0;
        }
        if (z_latter >= 3 && z_former >= 3)
        {
            dudz2 = 0.0;
        }
    }

    // set up the diffusion

    dDiffusiondt = diffusion * (dudx2 + dudz2) + ddiffusion_x * dudx + ddiffusion_z * dudz;

    return dDiffusiondt;
}

int main()
{
    omp_set_num_threads(CPU_NUM);

    FILE *in_1, *in_2, *in_3, *out;    //, *activation, *X_direction, *Z_direction;
    char *str;
    int x, z, dd, **g;
    double ff, ** D, ** CM,  ** new_para, **old_para;
    double **t;
    int cnt = 0;
    int outputcounter = 0;

    double total_time, ddt, du;
    total_time = 1000.0;
    ddt = 0.005;

    clock_t time_begin, time_end;

    time_begin = clock();

    g = (int **)malloc((X) * sizeof(int *));
    D = (double **)malloc((X) * sizeof(double *));
    CM = (double **)malloc((X) * sizeof(double *));
    new_para = (double **)malloc((X) * sizeof(double *));
    old_para = (double **)malloc((X) * sizeof(double *));
    t = (double **)malloc((X) * sizeof(double *));

    for (x = 0; x < X; x++)
    {
        g[x]  = (int *)malloc((Z) * sizeof(int));
        D[x]  = (double *)malloc((Z) * sizeof(double));
        CM[x] = (double *)malloc((Z) * sizeof(double));
        new_para[x] = (double *)malloc((Z) * sizeof(double));
        old_para[x] = (double *)malloc((Z) * sizeof(double));
        t[x] = (double *)malloc((Z) * sizeof(double));

        for (z = 0; z < Z; z++)
        {
            g[x][z] = 0;
            D[x][z] = 0.0;
            CM[x][z] = 0.0;
            new_para[x][z] = 0.0;
            old_para[x][z] = 0.0;
            t[x][z] = 0.0;
        }
    }

    in_1 = fopen ("FILE/2D/geo_1.txt", "r");
    for (z = 0; z < Z; z++)
    {
        for (x = 0; x < X; x++)
        {
            fscanf(in_1, "%d ", &dd);
            g[x][z] = dd;
            /*            if (dd == 1)
                            g[x][z] = 2;*/
            if (dd >= 5)
                g[x][z] = 0;
        }
        fscanf(in_1, "\n");
    }
    fclose (in_1);

    in_2 = fopen ("FILE/2D/cm_geo.txt", "r");
    for (z = 0; z < Z; z++)
    {
        for (x = 0; x < X; x++)
        {
            fscanf(in_2, "%lf ", &ff);
            CM[x][z] = ff;
        }
        fscanf(in_2, "\n");
    }
    fclose (in_2);

    in_3 = fopen ("FILE/2D/diffusion_geo.txt", "r");
    for (z = 0; z < Z; z++)
    {
        for (x = 0; x < X; x++)
        {
            fscanf(in_3, "%lf ", &ff);
            D[x][z] = ff;
        }
        fscanf(in_3, "\n");
    }
    fclose (in_3);

    realtype tout;
    N_Vector y[X][Z];

    for (x = 0; x < X; x++)
        for (z = 0; z < Z; z++)
        {
            y[x][z] = NULL;
        }

    MY_CELL ***cell;
    cell = (MY_CELL ** *)malloc((X) * sizeof(MY_CELL **));
    for (x = 0; x < X; x++)
    {
        cell[x] = (MY_CELL **)malloc((Z) * sizeof(MY_CELL *));
        for (z = 0; z < Z; z++)
        {
            cell[x][z] = (MY_CELL *)malloc(sizeof(MY_CELL));
            cell[x][z]->cell_type = g[x][z];
            cell[x][z]->HCN = 0;
            if (cell[x][z]->cell_type >= 3)
            {
                cell[x][z]->capacitance = CM[x][z];
                y[x][z] = N_VNew_Serial(39);
                initial_for_cvode(y[x][z], cell[x][z]);
                prepare_for_cvode(ddt, y[x][z], cell[x][z]);
            }
            else if (cell[x][z]->cell_type == 2 || cell[x][z]->cell_type == 1)
            {
                y[x][z] = N_VNew_Serial(33);
                initial_for_cvode(y[x][z], cell[x][z]);
                prepare_for_cvode(ddt, y[x][z], cell[x][z]);
            }
            else
            {
                cell[x][z]->V = -79.48;
            }

            new_para[x][z] = cell[x][z]->V;
            old_para[x][z] = cell[x][z]->V;
        }
    }

    for (x = 0; x < X; x++)
    {
        free(g[x]);
        free(CM[x]);
    }
    free(g);
    free(CM);

    /*    FILE *diff_vtk;
        char *new_str = malloc (50 * sizeof(char));
        sprintf (new_str, "diffusion_Vtk.vtk");
        // diff_vtk = fopen("diffusion_Vtk.vtk", 'wt');
        write_to_vtk_2D(D, X,  Z, new_str);*/


#ifdef PARALLEL
    #pragma omp parallel default(none) \
    shared(new_para, old_para, cell, D, t, y, ddt, outputcounter, str, out, cnt, total_time) \
    private(x, z, du, tout)
#endif
    {
        for (tout = ddt; tout < total_time; tout += ddt)
        {
#ifdef PARALLEL
            #pragma omp for schedule(static)
#endif
            for (x = 1; x < X - 1; x++)
                for (z = 1; z < Z - 1; z++)
                {
                    if (cell[x][z]->cell_type > 0)
                    {
                        du = DIFFUSION_OP(x, z, D, cell, old_para);
                        new_para[x][z] = old_para[x][z] + (ddt / 2) * du;
                        cell[x][z]->V = new_para[x][z];
                    }
                } // end of for loop

#ifdef PARALLEL
            #pragma omp for schedule(static)
#endif
            for (x = 1; x < X - 1; x++)
                for (z = 1; z < Z - 1; z++)
                {
                    if (cell[x][z]->cell_type > 0)
                    {
                        CVode(cell[x][z]->cvode_mem, tout, y[x][z], &t[x][z], CV_NORMAL); // 1 time step solution.
                        new_para[x][z] = new_para[x][z] + ddt * cell[x][z]->dvdt;
                    }
                } // end of for loop

#ifdef PARALLEL
            #pragma omp for schedule(static)
#endif
            for (x = 1; x < X - 1; x++)
                for (z = 1; z < Z - 1; z++)
                {
                    if (cell[x][z]->cell_type > 0)
                    {
                        du = DIFFUSION_OP(x, z, D, cell, new_para);
                        old_para[x][z] = new_para[x][z] + (ddt / 2) * du;
                    }
                } // end of for loop
            // measurements here: beacuse you have old_state and new_state
            /*        #ifdef PARALLEL
                    #pragma omp parallel for default(none) shared(new_para, cell, act_time, sstime) private(x , z)
                    #endif
                    for(x=1; x<X; x++)
                        for(z=1; z<Z; z++)
                        {
                            if(cell[x][z]->cell_type>0 && cell[x][z]->V<=0.0 && new_para[x][z]>0)
                                act_time[x][z] = sstime;
                        }*/

            // renew the voltage
            /*#ifdef PARALLEL
                    #pragma omp parallel for default(none) shared(new_para, cell) private(x , z)
            #endif
                    for (x = 1; x < X - 1; x++)
                        for (z = 1; z < Z - 1; z++)
                        {
                            if (cell[x][z]->cell_type > 0)
                                cell[x][z]->V = new_para[x][z];
                        }*/
#ifdef PARALLEL
            #pragma omp single
#endif
            {
                outputcounter++;
                if (outputcounter % ((int)(1.0 / ddt)) == 1)
                {
                    // output the voltage.
                    str = (char *)malloc (12 * sizeof(char));
                    sprintf (str, "f%d.vtk", cnt++);
                    out = fopen (str, "wt");

                    fprintf (out, "# vtk DataFile Version 3.0\n");
                    fprintf (out, "vtk output\n");
                    fprintf (out, "ASCII\n");
                    fprintf (out, "DATASET STRUCTURED_POINTS\n");
                    fprintf (out, "DIMENSIONS %d %d %d\n", X, Z, 1);
                    fprintf (out, "SPACING 1 1 1\n");
                    fprintf (out, "ORIGIN 0 0 0\n");
                    fprintf (out, "POINT_DATA %d\n", X * 1 * Z);
                    fprintf (out, "SCALARS ImageFile float 1\n");
                    fprintf (out, "LOOKUP_TABLE default\n");

                    for (z = 0; z < Z; z++)
                    {
                        for (x = 0; x < X; x++)
                        {
                            if (cell[x][z]->cell_type > 0)
                                fprintf (out, "%5.2lf ", cell[x][z]->V);
                            else
                                fprintf (out, "%5.2lf", -100.0);
                            fprintf (out, "\n");
                        }
                    }
                    fclose (out);
                    free (str);

                    /*            X_direction = fopen("x_1d.dat","a+");
                                for ( x = 5; x < 192; x++)
                                {
                                    fprintf(X_direction,"%5.2f\t ", cell[x][115]->V);
                                }
                                fprintf(X_direction,"\n");
                                fclose(X_direction);

                                Z_direction = fopen("z_1d.dat","a+");
                                for ( z = 5; z < 125; z++)
                                {
                                    fprintf(Z_direction,"%5.2f\t ", cell[192][z]->V);
                                }
                                fprintf(Z_direction,"\n");
                                fclose(Z_direction);*/
                } // end of output.
            }
            #pragma omp barrier
        } // end of time loop
    }  // end of parallel

    // output activation time
    /*    activation = fopen("activation_time.dat", "a+");
        for (z = 1; z < Z; z++)
            for (x = 1; x < X; x++)
            {
                if (g[x][z] > 0)
                    fprintf(activation, "%d\t%d\t%5.2f\n", x, z, act_time[x][z]);
                else
                    fprintf(activation, "%d\t%d\t%5.2f\n", x, z, 0.0);
            }
        fclose(activation);*/

    for (x = 0; x < X; x++)
    {
        free(D[x]);
        free(new_para[x]);
        free(old_para[x]);
    }
    free(D);
    free(new_para);
    free(old_para);

    for (x = 0; x < Z; x++)
    {
        for (z = 0; z < Z; z++)
        {
            N_VDestroy_Serial(y[x][z]);
            CVodeFree(&cell[x][z]->cvode_mem);
            free(cell[x][z]);
        }
        free(cell[x]);
    }
    free(cell);

    time_end = clock();
    printf("CPU t1ime=%f seconds\n", ((double)(time_end - time_begin)) / CLOCKS_PER_SEC);

    return 0;
}