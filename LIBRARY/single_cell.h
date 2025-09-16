#ifndef SINGLE_CELL_H
#define SINGLE_CELL_H

#include "constants.h"

typedef struct {
    int cell_type, HCN;
    double ACh, ISO;
    double * state, * current;
    double V, dvdt;
    double capacitance;
    void * cvode_mem;
} MY_CELL;


double SAN(double * state, double * rate, MY_CELL * data);
double ATRIUM(double * state, double * rate, MY_CELL * data);

/***** mouse SAN model******************************************************/
double SAN(double * state, double * rate, MY_CELL * data) {
    int HCN = data->HCN;
    double capacitance = data->capacitance;
    double ach = data->ACh;
    double iso = data->ISO;

    double c_cent, c_peri;

    c_cent = 0.025;
    c_peri = 0.05;

    /****we can make the conductance gradients here****/

    double Fc, Fp;
    Fc = (capacitance - c_cent) / (c_peri - c_cent);
    Fp = (c_peri - capacitance) / (c_peri - c_cent);

    double gst_cent, gst_peri, gst;

    gst_cent = 0.00006;
    gst_peri = gst_cent * 2.0;
    gst = Fp * gst_cent + Fc * gst_peri;

    double gna_ttxs_cent, gna_ttxs_peri, gna_ttxs;

    gna_ttxs_cent = 0.1 * 5.925e-05;
    gna_ttxs_peri = gna_ttxs_cent; // *2.015;   //*10.0;
    gna_ttxs = Fp * gna_ttxs_cent + Fc * gna_ttxs_peri;

    double gna_ttxr_cent, gna_ttxr_peri, gna_ttxr;

    gna_ttxr_cent = 0.1 * 5.925e-05;
    gna_ttxr_peri = gna_ttxr_cent * 1.385 * 7.0; //10.0;
    gna_ttxr = Fp * gna_ttxr_cent + Fc * gna_ttxr_peri;

    double gcal12_cent, gcal12_peri, gcal12;
    double gcal13_cent, gcal13_peri, gcal13;

    gcal12_cent = 0.0010 * 4.0 * 1.5;
    gcal13_cent = 0.0030 * 4.0 * 1.5;

    gcal12_peri = gcal12_cent * 3.4; // 5.006;
    gcal13_peri = gcal13_cent * 5.4; // 6.878;

    gcal12 = Fp * gcal12_cent + Fc * gcal12_peri;
    gcal13 = Fp * gcal13_cent + Fc * gcal13_peri;

    double gcat_cent, gcat_peri, gcat;
    gcat_cent = 0.75 * 0.01862;
    gcat_peri = gcat_cent * 2.009;
    gcat = Fp * gcat_cent + Fc * gcat_peri;

    double gto_cent, gto_peri, gto;
    gto_cent = 0.00492;
    gto_peri = gto_cent * 2.3;
    gto = Fp * gto_cent + Fc * gto_peri;

    double gsus_cent, gsus_peri, gsus;
    gsus_cent = 0.00039060;
    gsus_peri = gsus_cent * 2.9;
    gsus = Fp * gsus_cent + Fc * gsus_peri;

    double gkr_cent, gkr_peri, gkr;
    gkr_cent = 0.8 * 0.002955;
    gkr_peri = gkr_cent * 2.0;
    gkr = Fp * gkr_cent + Fc * gkr_peri;

    double gks_cent, gks_peri, gks;
    gks_cent = 0.000299;
    gks_peri = gks_cent * 5.4; // 11.271*0.2;
    gks = Fp * gks_cent + Fc * gks_peri;

    double gk1_cent, gk1_peri, gk1;
    gk1_cent = 0.229 * 0.0039228 * 0.9;
    gk1_peri = gk1_cent * 5.558;
    gk1 = Fp * gk1_cent + Fc * gk1_peri;

    double ghk1_cent = 0.0015;
    double ghk2_cent = 0.000728;
    double ghk4_cent = 0.0025;
    double ghna1_cent = 0.001;
    double ghna2_cent = 0.000858;
    double ghna4_cent = 0.00363;

    if (HCN == 1) {
        ghk1_cent = 0.0;
        ghna1_cent = 0.0;
    } else if (HCN == 2) {
        ghk2_cent = 0.0;
        ghna2_cent = 0.0;
    } else if (HCN == 4) {
        ghk4_cent = 0.0;
        ghna4_cent = 0.0;
    }

    double ghk1_peri = ghk1_cent * 1.5;
    double ghk2_peri = ghk2_cent * 1.442;
    double ghk4_peri = ghk4_cent * 15.0; // 1.904;
    double ghna1_peri = ghna1_cent * 1.5;
    double ghna2_peri = ghna2_cent * 1.597;
    double ghna4_peri = ghna4_cent * 10.0; // 1.492;

    double ghk1;
    double ghk2;
    double ghk4;
    double ghna1;
    double ghna2;
    double ghna4;

    ghk1 = Fp * ghk1_cent + Fc * ghk1_peri;
    ghk2 = Fp * ghk2_cent + Fc * ghk2_peri;
    ghk4 = Fp * ghk4_cent + Fc * ghk4_peri;
    ghna1 = Fp * ghna1_cent + Fc * ghna1_peri;
    ghna2 = Fp * ghna2_cent + Fc * ghna2_peri;
    ghna4 = Fp * ghna4_cent + Fc * ghna4_peri;

    double gbna_cent, gbna_peri, gbna;
    gbna_cent = 0.0001215;
    gbna_peri = gbna_cent * 2.0;
    gbna = Fp * gbna_cent + Fc * gbna_peri;

    double gbca_cent, gbca_peri, gbca;
    gbca_cent = 0.000015;
    gbca_peri = gbca_cent * 2.0;
    gbca = Fp * gbca_cent + Fc * gbca_peri;

    double gbk_cent, gbk_peri, gbk;
    gbk_cent = 0.0000025;
    gbk_peri = gbk_cent * 2.0;
    gbk = Fp * gbk_cent + Fc * gbk_peri;

    double knaca_cent, knaca_peri, kNaCa;
    knaca_cent = 5.5;
    knaca_peri = knaca_cent * 2.182;
    kNaCa = Fp * knaca_cent + Fc * knaca_peri;

    double inakmax_cent, inakmax_peri, inakmax;
    inakmax_cent = 0.14245;
    inakmax_peri = inakmax_cent * 2.185 * 1.4; // *2.11;
    inakmax = Fp * inakmax_cent + Fc * inakmax_peri;

    double Pup_cent, Pup_peri, Pup;
    Pup_cent = 0.04;
    Pup_peri = Pup_cent * 2.182;
    Pup = Fp * Pup_cent + Fc * Pup_peri;

    double ks_cent, ks_peri, ks;
    ks_cent = 1300000.0;
    ks_peri = ks_cent; // *2.182*0.5;
    ks = Fp * ks_cent + Fc * ks_peri;

    double v                   = data->V;
    double dst                 = state[0];
    double fst                 = state[1];
    double dt                  = state[2];
    double ft                  = state[3];
    double ikr_act             = state[4];
    double ikr_inact           = state[5];
    double iks_act             = state[6];
    double fl12                = state[7];
    double dl12                = state[8];
    double fl13                = state[9];
    double dl13                = state[10];
    double r                   = state[11];
    double m_ttxr              = state[12];
    double h_ttxr              = state[13];
    double j_ttxr              = state[14];
    double m_ttxs              = state[15];
    double h_ttxs              = state[16];
    double j_ttxs              = state[17];
    double y_1                 = state[18];
    double carel               = state[19];
    double caup                = state[20];
    double casub               = state[21];
    double Ftc                 = state[22];
    double Ftmc                = state[23];
    double Ftmm                = state[24];
    double Fcms                = state[25];
    double Fcmi                = state[26];
    double Fcq                 = state[27];
    double cai                 = state[28];
    double q                   = state[29];
    double fca                 = state[30];
    double nai                 = state[31];
    double ki                  = state[32];
    double resting             = state[33];
    double open                = state[34];
    double inactivated         = state[35];
    double resting_inactivated = state[36];
    double y_2                 = state[37];
    double y_4                 = state[38];
    double Ach_j               = state[39];
    double Ach_k               = state[40];

    double ist;
    double ib;
    double ik1;
    double icat;
    double ikr;
    double iks;
    double ical12;
    double ical13;
    double ina_ttxs;
    double ina_ttxr;
    double ih;
    double ito;
    double isus;
    double inak;
    double inaca;
    double ikach;

    /*Reversal Potentials****************************************************/
    double ena, eca, ek, eks;

    ena = (R * T_SAN / F) * log(nao_SAN / nai);
    ek  = (R * T_SAN / F) * log(ko_SAN / ki);
    eks = ((R * T_SAN) / F) * log((ko_SAN + 0.12 * nao_SAN) / (ki + 0.12 * nai));
    eca = (R * T_SAN / (2 * F)) * log(cao_SAN / casub);

    /*Effect of ACh**********************************************************/

    /****ikach****/
    double alpha_j, beta_j, alpha_k, beta_k, Ach_j_inf, Ach_k_inf, tau_j, tau_k, kach;

    alpha_j = 0.763678 * 73.1e-03; //ms-1
    beta_j =  1.035639 * 120.0e-03 / (1.0 + exp(-(v + 1.295674 * 50.0) / (0.867184 * 15.0))); // ms-1
    tau_j = 1.0 / (alpha_j + beta_j);
    Ach_j_inf = alpha_j * tau_j;
    alpha_k = 6.362677 * 3.7e-03; // ms-1
    beta_k =  0.101428 * 5.820e-03 / (1.0 + exp(-(v + 1.295674 * 50.0) / (0.867184 * 15.0))); // ms-1
    tau_k = 1.0 / (alpha_k + beta_k);
    Ach_k_inf = alpha_k * tau_k;

    ikach = 0.0705 * Ach_j * Ach_j * (ach / (5.49e-08 + ach)) * (ko_SAN / (10.0 + ko_SAN)) * (v - ek) / (1 + exp((v - ek - 0.329262 * 140) * F / (2.5 * R * T_SAN)));

    /****fractional block of Ical****/
    double b_ical_ach;

    b_ical_ach =  0.56 * (ach / (0.2e-07 + ach)); // 0.56 is bmax, 0.9e-07 is kca.

    /****shift of ih activation curve*****/

    double s_ih_ach;

    s_ih_ach = -16.8 * (pow(ach, 0.69) / (pow(1.26e-08, 0.69) + pow(ach, 0.69))); // 1.26e-08 is K_0.5_f, 0.69 is the value of n_f

    /*Effect of ISO (isoprenaline) **********************************************/

    iso = iso * 10.0e6;

    /**** increase of gcal ****/
    double fcal_iso;               // 0.9 is fCaL_max, 15.0e-03 is K_0.5_Ca.

    fcal_iso = 0.9 * (iso / (15.0e-3 + iso));

    /**** shift of I-V of Ical ****/
    double scal_iso;               // 10.0 is sCaL_max, 15.0e-03 is K_0.5_Ca.

    scal_iso = 10.0 * (iso / (15.0e-3 + iso));

    /**** increase of gcat ****/
    double fcat_iso;               // 1.0 is fCaT_max, 15.0e-03 is K_0.5_Ca.

    fcat_iso = 1.0 * (iso / (15.0e-3 + iso));

    /**** shift of activation curve of if ****/
    double s_ih_iso;               // 13.6 is Sf_max, 13.5e-03 is K_0.5_f, 0.392 is the value of n_f

    s_ih_iso = 13.6 * (pow(iso, 0.392) / (pow(iso, 0.392) + pow(13.5e-03, 0.392)));

    /**** increase of gkr and gks ****/
    double fk_iso;                 // 1.87 is fk_max, 19.0e-03(uM) is K_0.5_gk

    fk_iso = 1.87 * (iso / (iso + 19.0e-03));

    /**** shift in activation of IKr ****/
    double skr_iso;                // 7.5e-03(uM) is K_0.5_Kacti

    skr_iso = 15.0 * (iso / (iso + 7.5e-03));

    /**** increase of gst ****/
    double fst_iso;

    fst_iso = 1.0 * (iso / (iso + 33.0e-03));

    /**** Phosphorylation level of RyR ****/
    double fryr_iso;

    fryr_iso = 1.0 / (1.0 + exp(-2.55 * (log10(iso) + 2.5)));

    /**** Phosphorylation level of PLB ****/
    double fplb_iso;

    fplb_iso = 1.0 / (1.0 + exp(-3.0 * (log10(iso) + 2.2)));

    /*Ist********************************************************************/

    double qa, qi, tauqa, tauqi, alphaqa, betaqa, alphaqi, betaqi;
    qa = 1.0 / (1.0 + exp(-(v + 67.0) / 5.0));
    alphaqa = 1.0 / (0.15 * exp(-(v) / 11.0) + 0.2 * exp(-(v) / 700.0));
    betaqa  =  1.0 / (16.0 * exp((v) / 8.0) + 15.0 * exp((v) / 50.0));
    tauqa = 1.0 / (alphaqa + betaqa);
    alphaqi = 0.15 * 1.0 / (3100.0 * exp((v + 10.0) / 13.0) + 700.3 * exp((v + 10.0) / 70.0));
    betaqi =  0.15 * 1.0 / (95.7 * exp(-(v + 10.0) / 10.0) + 50.0 * exp(-(v + 10.0) / 700.0)) + 0.000229 / (1 + exp(-(v + 10.0) / 5.0));
    qi = alphaqi / (alphaqi + betaqi);
    tauqi = 1.0 / (alphaqi + betaqi);

    ist = (1.0 + fst_iso) * gst * dst * fst * (v - eist_SAN);

    /* Ib ************************************************************************/

    double ibca, ibna, ibk;
    ibna = gbna * (v - ena);
    ibca = gbca * (v - eca);
    ibk  = gbk * (v - ek);

    ib = (ibna + ibca + ibk);

    /*IK1**********************************************************************/

    double xk1inf;
    xk1inf = 1.0 / (1.0 + exp(0.070727 * (v - ek)));

    ik1 = gk1 * xk1inf * (ko_SAN / (ko_SAN + 0.228880)) * (v - ek);

    /**ICaT Cav3.1**************************************************************/

    double dt_inf, tau_dt, ft_inf, tau_ft;
    tau_dt = 1.0 / (1.068 * exp((v + 26.3) / 30.0) + 1.068 * exp(-(v + 26.3) / 30.0));
    dt_inf = 1.0 / (1.0 + exp(-(v + 26.0) / 6.0));
    tau_ft = 1.0 / (0.0153 * exp(-(v + 61.7) / 83.3) + 0.015 * exp((v + 61.7) / 15.38));
    ft_inf = 1.0 / (1.0 + exp((v + 61.7) / 5.6));

    icat = (1.0 + fcat_iso) * gcat * ft * dt * (v - ecat_SAN);

    /*Ikr********************************************************************/

    double ikr_act_inf, tau_ikr_act, ikr_inact_inf, tau_ikr_inact;
    ikr_act_inf = 1.0 / (1.0 + exp(-(v + 21.173694 + skr_iso) / 9.757086));
    tau_ikr_act = 0.699821 / (0.003596 * exp((v) / 15.339290) + 0.000177 * exp(-(v) / 25.868423));
    ikr_inact_inf = 1.0 / (1.0 + exp((v + 20.758474 - 4.0) / (19.0)));
    tau_ikr_inact = 0.2 + 0.9 * 1.0 / (0.1 * exp(v / 54.645) + 0.656 * exp(v / 106.157));

    ikr = (1.0 + fk_iso) * gkr * ikr_act * ikr_inact * (v - ek);

    /**IKs********************************************************************/

    double iks_act_inf, tau_iks_act;
    iks_act_inf = 1.0 / (1.0 + exp(-(v - 20.876040) / 11.852723));
    tau_iks_act =  1000.0 / (13.097938 / (1.0 + exp(-(v - 48.910584) / 10.630272)) + exp(-(v) / 35.316539));

    iks = (1.0 + fk_iso) * gks * iks_act * iks_act * (v - eks);

    /*ICaL*******************************************************************/

    double alpha_dl, beta_dl, tau_dl, tau_fl12, tau_fl13;
    double dl12_inf, fl12_inf;
    double dl13_inf, fl13_inf;
    double fca_inf, taufca;
    if (fabs(v) <= 0.001) {
        alpha_dl  = -28.39 * (v + 35.0) / (exp(-(v + 35.0) / 2.5) - 1.0) + 408.173;
    } else if (fabs(v + 35.0) <= 0.001) {
        alpha_dl  = 70.975 - 84.9 * v / (exp(-0.208 * v) - 1.0);
    } else if (fabs(v) > 0.001 && fabs(v + 35.0) > 0.001) {
        alpha_dl  = -28.39 * (v + 35.0) / (exp(-(v + 35.0) / 2.5) - 1.0) - 84.9 * v / (exp(-0.208 * v) - 1.0);
    }

    if (fabs(v - 5.0) <= 0.001)
        beta_dl   = 28.575;
    else if (fabs(v - 5.0) > 0.001)
        beta_dl   = 11.43 * (v - 5.0) / (exp(0.4 * (v - 5.0)) - 1.0);

    tau_dl  = 2000.0 / (alpha_dl + beta_dl);
    dl13_inf = 1.0 / (1 + exp(-(v + 13.5 + scal_iso) / 6.0));
    fl13_inf = 1.0 / (1 + exp((v + 35.0 + scal_iso) / 7.3));
    tau_fl13 = (7.4 + 45.77 * exp(-0.5 * (v + 28.7) * (v + 28.7) / (11 * 11)));
    dl12_inf = 1.0 / (1 + exp(-(v + 3.0 + scal_iso) / 5.0));
    fl12_inf = 1.0 / (1 + exp((v + 36.0 + scal_iso) / 4.6));
    tau_fl12 = (7.4 + 45.77 * exp(-0.5 * (v + 24.7) * (v + 24.7) / (11 * 11)));
    fca_inf = kmfca_SAN / (kmfca_SAN + casub);
    taufca = fca_inf / alpha_fca_SAN;

    ical12 = (1.0 + fcal_iso) * (1.0 - b_ical_ach) * gcal12 * fl12 * dl12 * fca * (v - ecal_SAN);
    ical13 = (1.0 + fcal_iso) * (1.0 - b_ical_ach) * gcal13 * fl13 * dl13 * fca * (v - ecal_SAN);

    /**INa**********************************************************************/

    double m3_inf_ttxr, h_inf_ttxr;
    double m3_inf_ttxs, h_inf_ttxs;
    double m_inf_ttxr , m_inf_ttxs;
    double tau_mr, tau_hr, tau_jr;
    double fna, hs, hsr;
    fna = (9.52e-02 * exp(-6.3e-2 * (v + 34.4)) / (1.0 + 1.66 * exp(-0.225 * (v + 63.7)))) + 8.69e-2;
    m3_inf_ttxr = 1.0 / (1.0 + exp(-(v + 45.213705) / 7.219547));
    h_inf_ttxr = 1.0 / (1.0 + exp((v + 62.578120 ) / 6.084036));
    m3_inf_ttxs = 1.0 / (1.0 + exp(-(v + 31.097331) / 5.0));
    h_inf_ttxs = 1.0 / (1.0 + exp((v + 56.0) / 3.0));
    m_inf_ttxr = pow(m3_inf_ttxr, 0.333);
    m_inf_ttxs = pow(m3_inf_ttxs, 0.333);
    tau_mr = 1000.0 * ((0.6247e-03 / (0.832 * exp(-0.335 * (v + 56.7)) + 0.627 * exp(0.082 * (v + 65.01)))) + 0.0000492);
    tau_hr = 1000.0 * (((3.717e-06 * exp(-0.2815 * (v + 17.11))) / (1 + 0.003732 * exp(-0.3426 * (v + 37.76)))) + 0.0005977);
    tau_jr = 1000.0 * (((0.00000003186 * exp(-0.6219 * (v + 18.8))) / (1 + 0.00007189 * exp(-0.6683 * (v + 34.07)))) + 0.003556);
    hs = (1.0 - fna) * h_ttxs + fna * j_ttxs;
    hsr = (1.0 - fna) * h_ttxr + fna * j_ttxr;

    if (fabs(v) > 0.005)
        ina_ttxs = gna_ttxs * m_ttxs * m_ttxs * m_ttxs * hs * nao_SAN * (F * F / (R * T_SAN)) * ((exp((v - ena) * F / (R * T_SAN)) - 1.0) / (exp(v * F / (R * T_SAN)) - 1.0)) * v;
    else
        ina_ttxs = gna_ttxs * m_ttxs * m_ttxs * m_ttxs * hs * nao_SAN * F * ((exp((v - ena) * F / (R * T_SAN)) - 1.0));
    if (fabs(v) > 0.005)
        ina_ttxr = gna_ttxr * m_ttxr * m_ttxr * m_ttxr * hsr * nao_SAN * (F * F / (R * T_SAN)) * ((exp((v - enattxr_SAN) * F / (R * T_SAN)) - 1.0) / (exp(v * F / (R * T_SAN)) - 1.0)) * v;
    else
        ina_ttxr = gna_ttxr * m_ttxr * m_ttxr * m_ttxr * hsr * nao_SAN * F * ((exp((v - enattxr_SAN) * F / (R * T_SAN)) - 1.0));

    /**If**************************************************************************/

    double ihk, ihna, ihk_hcn1, ihk_hcn2, ihk_hcn4, ihna_hcn1, ihna_hcn2, ihna_hcn4;
    double y_inf_1, y_inf_2, y_inf_4, tau_y_1, tau_y_2, tau_y_4;

    y_inf_1 = 1.0 / (1.0 + exp((v + 75.2 - s_ih_ach - s_ih_iso) / 6.57));
    y_inf_2 = 1.0 / (1.0 + exp((v + 92.0 - s_ih_ach - s_ih_iso) / 9.5)); //y_inf_2 = 1.0/(1.0 + exp((v+102.0)/5.49));
    y_inf_4 = 1.0 / (1.0 + exp((v + 91.2 - s_ih_ach - s_ih_iso) / 10.5)); //y_inf_4 = 1.0/(1.0 + exp((v+101.0)/9.71));
    tau_y_1 = 0.0035017700 / (exp(-(v + 590.192) * 0.0208819) + exp((v - 85.4612) / 13.26));
    tau_y_2 = 437.116 / (0.000225 * exp(-v / 13.02) + 105.395 * exp(v / 13.02));
    tau_y_4 = 3121.55 / (2.758e-05 * exp(-v / 9.92) + 766.3 * exp(v / 9.92));

    ihk_hcn1 = ghk1 * y_1 * (v - ek);
    ihk_hcn2 = ghk2 * y_2 * (v - ek);
    ihk_hcn4 = ghk4 * y_4 * (v - ek);
    ihna_hcn1 = ghna1 * y_1 * (v - ena);
    ihna_hcn2 = ghna2 * y_2 * (v - ena);
    ihna_hcn4 = ghna4 * y_4 * (v - ena);

    ihk = ihk_hcn1 + ihk_hcn2 + ihk_hcn4;
    ihna = ihna_hcn1 + ihna_hcn2 + ihna_hcn4;

    ih = (ihk + ihna);

    /*Ito*************************************************************************/

    double r_inf, tau_r, q_inf, tau_q;
    q_inf = 1.0 / (1.0 + exp((v + 49.0) / 13.0));
    tau_q = (6.06 + 39.102 / (0.57 * exp(-0.08 * (v + 44.0)) + 0.065 * exp(0.1 * (v + 45.93)))) / 0.67;
    r_inf = 1.0 / (1.0 + exp(-(v - 19.3) / 15.0));
    tau_r = (2.75 + 14.40516 / (1.037 * exp(0.09 * (v + 30.61)) + 0.369 * exp(-0.12 * (v + 23.84)))) / 0.303;

    ito = gto * q * r * (v - ek);

    /*Isus***********************************************************************/

    isus = gsus * r * (v - ek);

    /*Inak***********************************************************************/

    inak = inakmax * (pow(ko_SAN, 1.2) / (pow(kmkp_SAN, 1.2) + pow(ko_SAN, 1.2))) * (pow(nai, 1.3) / (pow(kmnap_SAN, 1.3) + pow(nai, 1.3))) / (1.0 + exp(-(v - ena + 120.0) / 30.0));

    /****iNaCa*******************************************************************/

    double di, doo, k43, k12, k14, k41, k34, k21, k23, k32, x1, x2, x3, x4;
    di = 1 + (casub / Kci_SAN) * (1 + exp(-Qci_SAN * v * F / (R * T_SAN)) + nai / Kcni_SAN) + (nai / K1ni_SAN) * (1 + (nai / K2ni_SAN) * (1 + nai / K3ni_SAN));
    doo = 1 + (cao_SAN / Kco_SAN) * (1 + exp(Qco_SAN * v * F / (R * T_SAN))) + (nao_SAN / K1no_SAN) * (1 + (nao_SAN / K2no_SAN) * (1 + nao_SAN / K3no_SAN));
    k43 = nai / (K3ni_SAN + nai);
    k12 = (casub / Kci_SAN) * exp(-Qci_SAN * v * F / (R * T_SAN)) / di;
    k14 = (nai / K1ni_SAN) * (nai / K2ni_SAN) * (1 + nai / K3ni_SAN) * exp(Qn_SAN * v * F / (R * T_SAN * 2.0)) / di;
    k41 = exp(-Qn_SAN * v * F / (R * T_SAN * 2.0));
    k34 = nao_SAN / (K3no_SAN + nao_SAN);
    k21 = (cao_SAN / Kco_SAN) * exp(Qco_SAN * v * F / (R * T_SAN)) / doo;
    k23 = (nao_SAN / K1no_SAN) * (nao_SAN / K2no_SAN) * (1 + nao_SAN / K3no_SAN) * exp(-Qn_SAN * v * F / (R * T_SAN * 2.0)) / doo;
    k32 = exp(Qn_SAN * v * F / (R * T_SAN * 2));
    x1 = k34 * k41 * (k23 + k21) + k21 * k32 * (k43 + k41);
    x2 = k43 * k32 * (k14 + k12) + k41 * k12 * (k34 + k32);
    x3 = k43 * k14 * (k23 + k21) + k12 * k23 * (k43 + k41);
    x4 = k34 * k23 * (k14 + k12) + k21 * k14 * (k34 + k32);

    inaca = kNaCa * (k21 * x2 - k12 * x1) / (x1 + x2 + x3 + x4);

    /****Ca handling in the SR***************************************************/

    double Jrel, Jup, Jtr, kcasr, kosrca, kisrca, drdt, dodt, didt, dridt;
    double pumpkmf;

    pumpkmf = pumpkmf_SAN * (1.0 - 0.5 * fplb_iso);
    Jup = Pup * (pow(cai / pumpkmf, pumphill_SAN) - pow(caup / pumpkmr_SAN, pumphill_SAN)) / (1.0 + pow(cai / pumpkmf, pumphill_SAN) + pow(caup / pumpkmr_SAN, pumphill_SAN));
    Jtr  = (caup - carel) / Ttr_SAN;
    Jrel = ks * open * (carel - casub);
    kcasr = maxsr_SAN - (maxsr_SAN - minsr_SAN) / (1.0 + pow(eca50sr_SAN / carel, hsrr_SAN));
    kosrca = koca_SAN * (1.0 + fryr_iso) / kcasr;
    kisrca = kica_SAN * kcasr;

    drdt = kim_SAN * resting_inactivated - kisrca * casub * resting - kosrca * casub * casub * resting + kom_SAN * open;
    dodt = kosrca * casub * casub * resting - kom_SAN * open - kisrca * casub * open + kim_SAN * inactivated;
    didt = kisrca * casub * open - kim_SAN * inactivated - kom_SAN * inactivated + kosrca * casub * casub * resting_inactivated;
    dridt = kom_SAN * inactivated - kosrca * casub * casub * resting_inactivated - kim_SAN * resting_inactivated + kisrca * casub * resting;

    /****Ca diffusion ***********************************************************/

    double Jcadif;

    Jcadif = (casub - cai) / tdifca_SAN;

    /****Ca buffering ***********************************************************/

    double dFtcdt, dFtmcdt, dFtmmdt, dFcmsdt, dFcmidt, dFcqdt;

    dFtcdt = kfTC_SAN * cai * (1.0 - Ftc) - kbTC_SAN * Ftc;
    dFtmcdt = kfTMC_SAN * cai * (1.0 - Ftmc - Ftmm) - kbTMC_SAN * Ftmc;
    dFtmmdt = kfTMM_SAN * Mgi_SAN * (1.0 - Ftmc - Ftmm) - kbTMM_SAN * Ftmm;
    dFcmsdt = kfCM_SAN * casub * (1.0 - Fcms) - kbCM_SAN * Fcms;
    dFcmidt = kfCM_SAN * cai * (1.0 - Fcmi) - kbCM_SAN * Fcmi;
    dFcqdt = kfCQ_SAN * carel * (1.0 - Fcq) - kbCQ_SAN * Fcq;

    /****Intracellular ionic concentrations *************************************/

    double ca_flux, dcasubdt, dcaidt, dcareldt, dcaupdt, nai_tot, ki_tot;;

    ca_flux = (ical12 + ical13 + icat - 2.0 * inaca + ibca) / (2.0 * F);
    dcasubdt = (-ca_flux + Jrel * vrel_SAN) / vsub_SAN - Jcadif - ConcCM_SAN * dFcmsdt;
    dcaidt = (Jcadif * vsub_SAN - Jup * vup_SAN) / vi_SAN - (ConcCM_SAN * dFcmidt + ConcTC_SAN * dFtcdt + ConcTMC_SAN * dFtmcdt);
    dcareldt = Jtr - Jrel - ConcCQ_SAN * dFcqdt;
    dcaupdt = Jup - Jtr * vrel_SAN / vup_SAN;
    nai_tot = ihna + ina_ttxr + ina_ttxs + 3.0 * inak + 3.0 * inaca + ist + ibna;
    ki_tot = ihk + iks + ikr + ik1 + ibk - 2.0 * inak + isus + ito;

    /*****Membrane popential ****************************************************/

    data->dvdt = - (ih + ina_ttxr + ina_ttxs + ical12 + ical13 + iks + ikr + ik1 + ist + ib + icat + inak + isus + inaca + ito + ikach) / capacitance;

    /*****renew variables *******************************************************/

    rate[0]  = ((qa - dst) / tauqa);
    rate[1]  = ((qi - fst) / tauqi);
    rate[2]  = ((dt_inf - dt) / tau_dt);
    rate[3]  = ((ft_inf - ft) / tau_ft);
    rate[4]  = (ikr_act_inf - ikr_act) / tau_ikr_act;
    rate[5]  = (ikr_inact_inf - ikr_inact) / tau_ikr_inact;
    rate[6]  = (iks_act_inf - iks_act) / tau_iks_act;
    rate[7]  = (fl12_inf - fl12) / tau_fl12;
    rate[8]  = (dl12_inf - dl12) / tau_dl;
    rate[9]  = (fl13_inf - fl13) / tau_fl13;
    rate[10] = (dl13_inf - dl13) / tau_dl;
    rate[11] = ((r_inf - r) / tau_r);
    rate[12] = (m_inf_ttxr - m_ttxr) / tau_mr;
    rate[13] = (h_inf_ttxr - h_ttxr) / tau_hr;
    rate[14] = (h_inf_ttxr - j_ttxr) / tau_jr;
    rate[15] = (m_inf_ttxs - m_ttxs) / tau_mr;
    rate[16] = (h_inf_ttxs - h_ttxs) / tau_hr;
    rate[17] = (h_inf_ttxs - j_ttxs) / tau_jr;
    rate[18] = (y_inf_1 - y_1) / tau_y_1;
    rate[19] = dcareldt;
    rate[20] = dcaupdt;
    rate[21] = dcasubdt;
    rate[22] = dFtcdt;
    rate[23] = dFtmcdt;
    rate[24] = dFtmmdt;
    rate[25] = dFcmsdt;
    rate[26] = dFcmidt;
    rate[27] = dFcqdt;
    rate[28] = dcaidt;
    rate[29] = ((q_inf - q) / tau_q);
    rate[30] = (fca_inf - fca) / taufca;
    rate[31] = -(nai_tot) / (F * vi_SAN);
    rate[32] = -(ki_tot) / (F * vi_SAN);
    rate[33] = drdt;
    rate[34] = dodt;
    rate[35] = didt;
    rate[36] = dridt;
    rate[37] = (y_inf_2 - y_2) / tau_y_2;
    rate[38] = (y_inf_4 - y_4) / tau_y_4;
    rate[39] = (Ach_j_inf - Ach_j) / tau_j;
    rate[40] = (Ach_k_inf - Ach_k) / tau_k;

    double *current = data->current;
    current[0] = ist;
    current[1] = ib;
    current[2] = ik1;
    current[3] = icat;
    current[4] = ikr;
    current[5] = iks;
    current[6] = ical12;
    current[7] = ical13;
    current[8] = ina_ttxs;
    current[9] = ina_ttxr;
    current[10] = ihk_hcn1 + ihna_hcn1;
    current[11] = ihk_hcn2 + ihna_hcn2;
    current[12] = ihk_hcn4 + ihna_hcn4;
    current[13] = ito;
    current[14] = isus;
    current[15] = inak;
    current[16] = inaca;

    return 0;
}// end of SAN

/***** mouse atrium cell model**********************************************/
double ATRIUM(double * state, double * rate, MY_CELL * data) {
    double v       = data->V;
    double cai     = state[0];
    double ca_ss   = state[1];
    double ca_jsr  = state[2];
    double ca_nsr  = state[3];
    double LTRPNca = state[4];
    double HTRPNca = state[5];
    double O       = state[6];
    double C2      = state[7];
    double C3      = state[8];
    double C4      = state[9];
    double I1      = state[10];
    double I2      = state[11];
    double I3      = state[12];
    double P_c2    = state[13];
    double P_o1    = state[14];
    double P_o2    = state[15];
    double P_ryr   = state[16];
    double C_na1   = state[17];
    double C_na2   = state[18];
    double O_na    = state[19];
    double IF_na   = state[20];
    double I1_na   = state[21];
    double I2_na   = state[22];
    double IC_na2  = state[23];
    double IC_na3  = state[24];
    double nai     = state[25];
    double ki      = state[26];
    double a_tof   = state[27];
    double i_tof   = state[28];
    double n_ks    = state[29];
    double a_ur    = state[30];
    double i_ur    = state[31];
    double a_kss   = state[32];

    double Ical;
    double Ipca;
    double Inaca;
    double Icab;
    double Ina;
    double Inab;
    double Inak;
    double Iktof;
    double Ik1;
    double Iks;
    double Ikur;
    double Ikss;
    double Iclca;

    double Ena, Ek, Eca;
    Ena = (R * T_atrium / F) * log((0.9 * nao_atrium + 0.1 * ko_atrium) / (0.9 * nai + 0.1 * ki));
    Ek  = (R * T_atrium / F) * log(ko_atrium / ki);
    Eca = (R * T_atrium / (2 * F)) * log(cao_atrium / cai);

    /*----- Ical: Calcium Currents -----*/
    double C1, alpha, beta, gama, K_pcf, d_O, d_C2, d_C3, d_C4, d_I1, d_I2, d_I3;

    alpha = 0.4 * exp((v + 12.0) / 10.0) * (1 + 0.7 * exp(-pow((v + 40.0), 2) / 10.0) - 0.75 * exp(-pow((v + 20.0), 2) / 400.0)) / (1 + 0.12 * exp((v + 12.0) / 10.0));
    beta = 3 * 0.05 * exp(-(v + 30.0) / 13.0); //weijian
    gama = K_pc_max_atrium * ca_ss / (K_pc_half_atrium + ca_ss);
    K_pcf = 13.0 * (1 - exp(-pow((v + 22.5), 2) / 100.0));

    Ical = Gcal_atrium * O * (v - E_cal_atrium);

    C1 = 1 - (O + C2 + C3 + C4 + I1 + I2 + I3);
    d_O = alpha * C4 - 4 * beta * O + K_pcb_atrium * I1 - gama * O + 0.001 * (alpha * I2 - K_pcf * O);

    d_C2 = 4 * alpha * C1 - beta * C2 + 2 * beta * C3 - 3 * alpha * C2;
    d_C3 = 3 * alpha * C2 - 2 * beta * C3 + 3 * beta * C4 - 2 * alpha * C3;
    d_C4 = 2 * alpha * C3 - 3 * beta * C4 + 4 * beta * O - alpha * C4 + 0.01 * (4 * K_pcb_atrium * beta * I1 - alpha * gama * C4) + 0.002 * (4 * beta * I2 - K_pcf * C4) + 4 * beta * K_pcb_atrium * I3 - gama * K_pcf * C4;
    d_I1 = gama * O - K_pcb_atrium * I1 + 0.001 * (alpha * I3 - K_pcf * I1) + 0.01 * (alpha * gama * C4 - 4 * beta * K_pcb_atrium * I1);
    d_I2 = 0.001 * (K_pcf * O - alpha * I2) + K_pcb_atrium * I3 - gama * I2 + 0.002 * (K_pcf * C4 - 4 * beta * I2);
    d_I3 = 0.001 * (K_pcf * I1 - alpha * I3) + gama * I2 - K_pcb_atrium * I3 + gama * K_pcf * C4 - 4 * beta * K_pcb_atrium * I3;


    /*----- Ipca: Calcium pump current -----*/

    Ipca = Ipca_max_atrium * cai * cai / (K_m_pca_atrium * K_m_pca_atrium + cai * cai);

    /*----- Inaca: Na/Ca exchange current -----*/

    Inaca = k_naca_atrium / (pow(K_m_na_atrium, 3) + pow(nao_atrium, 3)) / (K_m_ca_atrium + cao_atrium) / (1 + K_sat_atrium * exp((eta_atrium - 1) * (v - 40.0) * F / R / T_atrium))
            * (exp(eta_atrium * (v - 40.0) * F / R / T_atrium) * pow(nai, 3) * cao_atrium - exp((eta_atrium - 1) * (v - 40.0) * F / R / T_atrium) * pow(nao_atrium, 3) * cai);

    /*----- Icab: Calcium background current -----*/

    Icab = Gcab_atrium * (v - Eca);

    /*----- Ina: Fast Na current -----*/
    double C_na3, d_C_na2, d_C_na1, d_O_na, d_IF_na, d_I1_na, d_I2_na, d_IC_na2, d_IC_na3;
    double a_na11, a_na12, a_na13, b_na11, b_na12, b_na13, a_na3, b_na3, a_na2, b_na2, a_na4, b_na4, a_na5, b_na5;

    Ina = Gna_atrium * O_na * (v - Ena);

    a_na11 = 3.802 / (0.1027 * exp(-(v + 2.5) / 17.0) + 0.2 * exp(-(v + 2.5) / 150.0));
    a_na12 = 3.802 / (0.1027 * exp(-(v + 2.5) / 15.0) + 0.23 * exp(-(v + 2.5) / 150.0));
    a_na13 = 3.802 / (0.1027 * exp(-(v + 2.5) / 12.0) + 0.25 * exp(-(v + 2.5) / 150.0));
    b_na11 = 0.1917 * exp(-(v + 2.5) / 20.3);
    b_na12 = 0.20 * exp(-(v - 2.5) / 20.3);
    b_na13 = 0.22 * exp(-(v - 7.5) / 20.3);
    a_na3 = 7.0e-7 * exp(-(v + 7.0) / 7.7);
    b_na3 = 0.0084 + 0.00002 * (v + 7.0);
    a_na2 = 1.0 / (0.188495 * exp(-(v + 7.0) / 16.6) + 0.393956);
    b_na2 = a_na13 * a_na2 * a_na3 / (b_na13 * b_na3);
    a_na4 = a_na2 / 1000.0;
    b_na4 = a_na3;
    a_na5 = a_na2 / 95000.0;
    b_na5 = a_na3 / 50.0;

    C_na3 = 1 - (O_na + C_na1 + C_na2 + IF_na + I1_na + I2_na + IC_na2 + IC_na3);
    d_C_na2 = a_na11 * C_na3 - b_na11 * C_na2 + b_na12 * C_na1 - a_na12 * C_na2 + a_na3 * IC_na2 - b_na3 * C_na2;
    d_C_na1 = a_na12 * C_na2 - b_na12 * C_na1 + b_na13 * O_na - a_na13 * C_na1 + a_na3 * IF_na - b_na3 * C_na1;
    d_O_na = a_na13 * C_na1 - b_na13 * O_na + b_na2 * IF_na - a_na2 * O_na;
    d_IF_na = a_na2 * O_na - b_na2 * IF_na + b_na3 * C_na1 - a_na3 * IF_na + b_na4 * I1_na - a_na4 * IF_na + a_na12 * IC_na2 - b_na12 * IF_na;
    d_I1_na = a_na4 * IF_na - b_na4 * I1_na + b_na5 * I2_na - a_na5 * I1_na;
    d_I2_na = a_na5 * I1_na - b_na5 * I2_na;
    d_IC_na2 = a_na11 * IC_na3 - b_na11 * IC_na2 + b_na12 * IF_na - a_na12 * IC_na2 + b_na3 * C_na2 - a_na3 * IC_na2;
    d_IC_na3 = b_na11 * IC_na2 - a_na11 * IC_na3 + b_na3 * C_na3 - a_na3 * IC_na3;

    /*----- Inab: Background Na current -----*/

    Inab = Gnab_atrium * (v - Ena);

    /*----- Iktof: Transient outward K current -----*/
    double d_a_tof, d_i_tof, a_a, b_a, a_i, b_i;

    a_a = 0.09564 * exp(0.028477 * (v + 40.0));
    b_a = 0.1156 * exp(-0.07237 * (v + 40.0));

    a_i = 0.000152 * exp(-(v + 13.5) / 7.0) / (0.067083 * exp(-(v + 33.5) / 7.0) + 1.0);
    b_i = 0.00095 * exp((v + 33.5) / 7.0) / (0.051335 * exp((v + 33.5) / 7.0) + 1.0);

    d_a_tof = a_a * (1 - a_tof) - b_a * a_tof;
    d_i_tof = a_i * (1 - i_tof) - b_i * i_tof;

    Iktof = Gktof_atrium * pow(a_tof, 3) * i_tof * (v - Ek);

    /*----- Ik1: Time-independent K current -----*/

    Ik1 = 1.35908 * (ko_atrium / (ko_atrium + 500.0)) * ((v - Ek) / (1 + exp(0.04886 * (v - Ek - 0.35) + 2.36)));

    /*----- Iks: Slow delayed rectifier K cureent -----*/
    double a_n, b_n, d_n_ks;

    Iks = Gks_atrium * pow(n_ks, 2) * (v - Ek);

    a_n = 0.00000481333 * (v + 26.5) / (1 - exp(-0.128 * (v + 26.5)));
    b_n = 0.0000953333 * exp(-0.038 * (v + 26.5));
    d_n_ks = a_n * (1 - n_ks) - b_n * n_ks;

    /*----- Ikur: Ultraraapidly activiting delayed rectifier K current -----*/
    double tau_iur, tau_aur, a_ss, i_ss, d_i_ur, d_a_ur;

    Ikur = Gkur_atrium * a_ur * i_ur * (v - Ek);

    tau_aur = 0.593 * exp(-0.00429 * v) + 2.058;
    tau_iur = 850.0 - 170.0 / (1 + exp((v + 45.2) / 5.7));

    a_ss = 1 / (1 + exp(-(v + 22.5) / 7.7));
    i_ss = 1 / (1 + exp((v + 45.2) / 5.7));

    d_a_ur = (a_ss - a_ur) / tau_aur;
    d_i_ur = (i_ss - i_ur) / tau_iur;

    /*----- Ikss: Non-inactivating steady-state K current -----*/
    double tau_kss, d_a_kss;

    Ikss = Gkss_atrium * a_kss * (v - Ek);

    tau_kss = 39.3 * exp(-0.00562 * v) + 13.17;

    d_a_kss = (a_ss - a_kss) / tau_kss;

    /*----- Inak: Na/K pump current -----*/
    double f_nak, rou;

    rou = 1 / 7 * (exp(nao_atrium / 67300.0) - 1.0);

    f_nak = 1 / (1 + 0.1245 * exp(-0.1 * v * F / R / T_atrium) + 0.0365 * rou * exp(-v * F / R / T_atrium));

    Inak = Inak_max_atrium * f_nak / (1.0 + pow(K_m_nai_atrium / nai, 3 / 2)) * ko_atrium / (ko_atrium + K_m_ko_atrium);

    /*----- Iclca: Ca-activated cl current -----*/
    double O_clca;

    O_clca = 0.2 / (1 + exp(-(v - 46.7) / 7.8));

    Iclca = Gclca_atrium * O_clca * cai / (cai + K_m_cl_atrium) * (v - E_cl_atrium);

    /*----- Calcium fluxes -----*/
    double d_P_ryr, J_rel, J_tr, J_xfer, J_leak, J_up, J_trpn;

    J_rel = v1_atrium_atrium * (P_o1 + P_o2) * (ca_jsr - ca_ss) * P_ryr;
    J_tr = (ca_nsr - ca_jsr) / tau_tr_atrium;
    J_xfer = (ca_ss - cai) / tau_xfer_atrium;
    J_leak = v2_atrium * (ca_nsr - cai);
    J_up = v3_atrium * cai * cai / (K_m_up_atrium * K_m_up_atrium + cai * cai);
    J_trpn = k_pos_htrpn_atrium * cai * (HTRPN_tot_atrium - HTRPNca) - k_neg_htrpn_atrium * HTRPNca + k_pos_ltrpn_atrium * cai * (LTRPN_tot_atrium - LTRPNca) - k_neg_ltrpn_atrium * LTRPNca;

    d_P_ryr = -0.04 * P_ryr - 0.1 * Ical / Ical_max_atrium * exp(-pow((v - 5.0), 2) / 648.0);

    /*----- Calcium concentration -----*/
    double B_jsr, B_ss, B_i, d_ca_nsr, d_ca_jsr, d_ca_ss, d_cai;

    B_jsr = pow((1 + (CSQN_tot_atrium * K_CSQN_m_atrium / pow((K_CSQN_m_atrium + ca_jsr), 2))), -1);
    B_ss = pow((1 + (CMDN_tot_atrium * K_CMDN_m_atrium / pow((K_CMDN_m_atrium + ca_ss), 2))), -1);
    B_i = pow((1 + (CMDN_tot_atrium * K_CMDN_m_atrium / pow((K_CMDN_m_atrium + cai), 2))), -1);

    d_cai = B_i * (J_leak + J_xfer - J_up - J_trpn - (Icab - 2 * Inaca + Ipca) * A_cap_atrium * Cm / (2 * V_myo_atrium * F));
    d_ca_ss = B_ss * (J_rel * V_jsr_atrium / V_ss_atrium - J_xfer * V_myo_atrium / V_ss_atrium - Ical * A_cap_atrium * Cm / (2 * V_ss_atrium * F));
    d_ca_jsr = B_jsr * (J_tr - J_rel);
    d_ca_nsr = (J_up - J_leak) * V_myo_atrium / V_nsr_atrium - J_tr * V_jsr_atrium / V_nsr_atrium;

    /*----- Calcium Buffering -----*/
    double d_LTRPNca, d_HTRPNca;

    d_LTRPNca = k_pos_ltrpn_atrium * cai * (LTRPN_tot_atrium - LTRPNca) - k_neg_ltrpn_atrium * LTRPNca;
    d_HTRPNca = k_pos_htrpn_atrium * cai * (HTRPN_tot_atrium - HTRPNca) - k_neg_htrpn_atrium * HTRPNca;

    /*----- Ryanodine receptors -----*/
    double P_c1, d_P_o1, d_P_o2, d_P_c2;

    P_c1 = 1 - (P_c2 + P_o1 + P_o2);

    d_P_o1 = k_pos_a_atrium * pow(ca_ss, N_atrium) * P_c1 - k_neg_a_atrium * P_o1 - k_pos_b_atrium * pow(ca_ss, M_atrium) * P_o1 + k_neg_b_atrium * P_o2 - k_pos_c_atrium * P_o1 + k_neg_c_atrium * P_c2;
    d_P_o2 = k_pos_b_atrium * pow(ca_ss, M_atrium) * P_o1 - k_neg_b_atrium * P_o2;
    d_P_c2 = k_pos_c_atrium * P_o1 - k_neg_c_atrium * P_c2;

    /*----- Na concentration -----*/
    double d_nai;

    d_nai = -(Ina + Inab + 3 * Inaca + 3 * Inak) * A_cap_atrium * Cm / (V_myo_atrium * F);

    /*----- K concentration -----*/
    double d_ki;

    d_ki = -(Iktof + Ik1 + Iks + Ikss + Ikur - 2 * Inak) * A_cap_atrium * Cm / (V_myo_atrium * F);

    /*----- Update values -----*/

    data->dvdt =  -(Ical + Ipca + Inaca + Icab + Ina + Inab + Inak + Iktof + Ik1 + Iks + Ikur + Ikss + Iclca) ; //(pA/pF)

    /*****renew variables ***********************************************/

    rate[0]  = d_cai;
    rate[1]  = d_ca_ss;
    rate[2]  = d_ca_jsr;
    rate[3]  = d_ca_nsr;
    rate[4]  = d_LTRPNca;
    rate[5]  = d_HTRPNca;
    rate[6]  = d_O;
    rate[7]  = d_C2;
    rate[8]  = d_C3;
    rate[9]  = d_C4;
    rate[10] = d_I1;
    rate[11] = d_I2;
    rate[12] = d_I3;
    rate[13] = d_P_c2;
    rate[14] = d_P_o1;
    rate[15] = d_P_o2;
    rate[16] = d_P_ryr;
    rate[17] = d_C_na1;
    rate[18] = d_C_na2;
    rate[19] = d_O_na;
    rate[20] = d_IF_na;
    rate[21] = d_I1_na;
    rate[22] = d_I2_na;
    rate[23] = d_IC_na2;
    rate[24] = d_IC_na3;
    rate[25] = d_nai;
    rate[26] = d_ki;
    rate[27] = d_a_tof;
    rate[28] = d_i_tof;
    rate[29] = d_n_ks;
    rate[30] = d_a_ur;
    rate[31] = d_i_ur;
    rate[32] = d_a_kss;

    return 0;
}// end of ATRIUM

#endif