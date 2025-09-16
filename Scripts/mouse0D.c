#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "LIBRARY/for_cvode.h"

#define stim_offset     100
#define stim_period     200
#define stim_duration   4
#define stim_amplitude  12
#define number_of_apds  100

void measurement(double, double, double, double, double);
double min_potential[number_of_apds];
double tmin_potential[number_of_apds];
double max_potential[number_of_apds];
double dvdtmax[number_of_apds];
double vdvdtmax[number_of_apds];
double apd_start[number_of_apds];
double ddr[number_of_apds];
double top[number_of_apds];
double top_slope[number_of_apds];
double apd50[number_of_apds];
double apd_at_mins50[number_of_apds];
double apd90[number_of_apds];
double cycle_length[number_of_apds];
int param_counter = 0;
int start_output = 0;

int main(int argc, char *argv[]) {

    // for measurement
    for(int i=0;i<number_of_apds;i++)
    {
        min_potential[i] = 100000.0;
        max_potential[i] = -10000.0;
        dvdtmax[i]       = -10000.0;
        ddr[i]           = -10000.0;
        top[i]           =  10000.0;
        top_slope[i]     = -10000.0;
        apd50[i]         = -10000.0;
        apd_at_mins50[i]         = -10000.0;
        apd90[i]         = -10000.0;
        cycle_length[i]  = -10000.0;
    }

    FILE *outputcurrents, *outputmeasurements;
    // outputcurrents = fopen("Mouse_SAN.dat", "w+");
    outputmeasurements = fopen("MouseSAN_measurement.dat","w+");

    int current_trace_counter = 0;

    double Vnew, dvdtnew, total_time = 5000.0;
    double ddt = 0.1;
    double Istim, past;
    int NEQ;
    clock_t time_begin, time_end;
    double ach, iso;
    ach = (atof(argv[1])/100.0) * 0.012/*((0.04)*(atof(argv[1])/100.0))*/*1.0e-06; // (M)
    iso = pow(10, (-8.0+6.0*(atof(argv[2])/100.0)))/*pow(10, (-5.0+5.0*(atof(argv[2])/100.0)))*/*1.0e-06; // (M)


    time_begin = clock();

    MY_CELL *cell;
    cell = (MY_CELL *) malloc(sizeof(MY_CELL));
    cell->cell_type = 3;
    cell->capacitance = 0.025;
    cell->ACh = ach;
    cell->ISO = iso;
    // cell->state = (double *)malloc((40) * sizeof(double));
    cell->HCN = 1;
    cell->cvode_mem = NULL;
    cell->current = (double *) malloc(17 * sizeof(double));

    if (cell->cell_type >= 3) {
        NEQ = 41;
    } else {
        NEQ = 33;
    }

    realtype t, tout;
    N_Vector y;

    y = NULL;

    y = N_VNew_Serial(NEQ);
    initial_for_cvode(y, cell);
    prepare_for_cvode(ddt, y, cell);

    for (tout = ddt; tout < total_time; tout += ddt) {
        CVode(cell->cvode_mem, tout, y, &t, CV_NORMAL); // 1 time step solution.
        past =  floor(tout / stim_period) * stim_period;
        Istim = (cell->cell_type == 2 && tout - past >= stim_offset && tout - past <= stim_offset + stim_duration ? stim_amplitude : 0.00000);

        Vnew = cell->V + ddt * (cell->dvdt + Istim);
        measurement(Vnew, cell->V, cell->dvdt, dvdtnew, tout);
        cell->V = Vnew;
        dvdtnew = cell->dvdt;

/*        if (tout>(total_time-1000.0) && current_trace_counter % (int(1/ddt))  == 0) {
            fprintf(outputcurrents, "%10.2f\t", tout-(total_time-1000));
            fprintf(outputcurrents, "%10.2f\t", cell->V);
            // fprintf(outputcurrents, "%10.10f\t", (cell->current[8]/cell->capacitance), (cell->current[9]/cell->capacitance), cell->current[6]/cell->capacitance, cell->current[7]/cell->capacitance );
            // fprintf(outputcurrents, "%10.10f\t%10.10f\t%10.10f\t%10.10f\t", (cell->current[6]+cell->current[7])/cell->capacitance, cell->current[3]/cell->capacitance, (cell->current[10]+cell->current[11]+cell->current[12])/cell->capacitance, cell->current[4]/cell->capacitance );
            // fprintf(outputcurrents, "%10.10f\t", Ith(y, 28));

            // for(int i=0; i<NEQ; i++)
            // {
            //     fprintf(outputcurrents, "%10.10f\t", Ith(y,i));
            // }
            // for (int i = 0; i < 17; i++) {
            //     fprintf(outputcurrents, "%10.10f\t", cell->current[i]/cell->capacitance);
            // }
            // fprintf(outputcurrents, "%10.10f\t%10.10f\t%10.10f\t", Ith(y, 28), Ith(y, 31), Ith(y,32));
            fprintf(outputcurrents, "\n");
        }

        current_trace_counter++;*/
    } // end of time loop

    // fclose(outputcurrents);

    for(int i=1; i<param_counter; i++)
    {
        fprintf(outputmeasurements,"%10.5f\n",cycle_length[i]);
    }
    fclose(outputmeasurements);

    N_VDestroy_Serial(y);
    CVodeFree(&cell->cvode_mem);
    free(cell);

    time_end = clock();
    printf("CPU t1ime=%f seconds\n", ((double)(time_end - time_begin)) / CLOCKS_PER_SEC);

    return 0;
}

/***** parameters measurement *********************************************************/

void measurement(double Vnow, double Vold, double dvdt, double dvdtold, double sstime)
{
    if(dvdt>=0.0&&dvdtold<0.0)
    {
        min_potential[param_counter] = Vold;
        tmin_potential[param_counter] = sstime;
        start_output = 1;
    }

    if(dvdt>dvdtmax[param_counter]&&start_output>0)
    {
        dvdtmax[param_counter] = dvdt;
        apd_start[param_counter] = sstime;
        vdvdtmax[param_counter] = Vnow;
    }

    if(dvdtold>0.0&&dvdt<=0.0)
    {
        max_potential[param_counter] = Vold;
        top_slope[param_counter] = (max_potential[param_counter]-min_potential[param_counter])/(sstime - tmin_potential[param_counter]);
    }

    if((param_counter>0)&&(dvdtold<=top_slope[param_counter-1])&&(dvdt>top_slope[param_counter-1]))
    {
        top[param_counter] = Vold;
        ddr[param_counter] = (Vold - min_potential[param_counter])/(sstime - tmin_potential[param_counter]);
    }

    if(Vnow<=0.5*max_potential[param_counter]+0.5*min_potential[param_counter]&&Vold>0.5*max_potential[param_counter]+0.5*min_potential[param_counter])
    {
        if(apd_start[param_counter]>0.0)
            apd50[param_counter] = sstime - apd_start[param_counter];
    }

    if(Vnow<=-50.0 && Vold>-50.0)
    {
        if(apd_start[param_counter]>0.0)
            apd_at_mins50[param_counter] = sstime - apd_start[param_counter];
    }

    if(Vnow<=0.1*max_potential[param_counter]+0.9*min_potential[param_counter]&&Vold>0.1*max_potential[param_counter]+0.9*min_potential[param_counter])
    {
        if(apd_start[param_counter]>0.0)
        {
            apd90[param_counter] = sstime - apd_start[param_counter];
            cycle_length[param_counter] = apd_start[param_counter]-apd_start[param_counter-1];
            // printf("MDP=%10.2f\n",min_potential[param_counter]);
            // printf("OS is%10.2f\n",max_potential[param_counter]);
            // printf("dvdt_max=%10.2f\n",dvdtmax[param_counter]);
            // printf("APD50=%10.2f\n",apd50[param_counter]);
            // printf("APD90=%10.2f\n",apd90[param_counter]);
            // printf("APD(at -50 mV )=%10.2f\n",apd_at_mins50[param_counter]);
            printf("Cycle Length is%10.2f\n", cycle_length[param_counter]);
            printf("Heart Rate   is%10.2f\n", (60000.0/cycle_length[param_counter]));
            // printf("DDR is%10.2f\n",ddr[param_counter]);
            // printf("TOP is%10.2f\n",top[param_counter]);
            printf("\n");
            param_counter++;
        }
    }
}// end of measurement
