#ifndef CONSTANTS_H
#define CONSTANTS_H

#define R 8.314472					// Universal gas constant
#define F 96.4845					// Faraday’s constant

/*SAN cell paramters. Central cell.*****************************************************/
#define T_SAN 310.5 				// Absolute temperature (in K)
#define vcell_SAN 3.0 				// Cell volume, pL
#define l_cell_SAN 66.3767257
#define r_cell_SAN 3.792956
#define vrel_SAN 0.0036				// Volume of the JSR, pL
#define vsub_SAN 0.03328117 		// Subspace volume, pL
#define vup_SAN 0.0348 				// Volume of the NSR, pL
#define vi_SAN 1.34671883			// Myoplasmic volume available for Ca2+ diffusion
#define Mgi_SAN 2.5 				// Intracellular Mg2ϩ concentration, mM
#define nao_SAN 140.0 				// extracellular Na+ concentration, mM
#define cao_SAN 1.8 				// extracellular Ca2+ concentration, mM
#define ko_SAN 5.4 					// extracellular K+ concentration, mM
#define eist_SAN 17.0 				// Reversal potential of Ist, mV
#define ecal_SAN 47.0 				// Reversal potential of ICaL,1.2 and ICaL,1.3, mV
#define kmfca_SAN 0.00035 			// Dissociation constant for Ca2ϩ-dependent inactivation of ICaL,1.2 and ICaL,1.3, mM
#define alpha_fca_SAN 0.021
#define ecat_SAN 45.0 				// Reversal potential of ICaT. mV	
#define enattxr_SAN 41.5761 		// Reversal potential of INa,1.5, mV
#define kmnap_SAN 14.0 				// Half-maximal [Na+]i for INak, mM
#define kmkp_SAN 1.4 				// Half-maximal [K+]o for INak, mM
#define K1ni_SAN 395.3 				// Dissociation constant for [Na+]i binding to the first site on the INaCa transporter, mM
#define K1no_SAN 1628 				// Dissociation constant for [Na+]o binding to the first site on the INaCa transporter, mM
#define K2ni_SAN 2.289 				// Dissociation constant for [Na+]i binding to the second site on the INaCa transporter, mM
#define K2no_SAN 561.4 				// Dissociation constant for [Na+]o binding to second site on the INaCa transporter, mM
#define K3ni_SAN 26.44 				// Dissociation constant for [Na+]i binding to the third site on the INaCa transporter, mM
#define K3no_SAN 4.663 				// Dissociation constant for [Na+]o binding to the third site on the INaCa transporter, mM
#define Kci_SAN 0.0207 				// Dissociation constant for [Ca2+]sub binding to the INaCa transporter, mM
#define Kco_SAN 3.663 				// Dissociation constant for [Ca2+]o binding to the INaCa transporter, mM
#define Kcni_SAN 26.44 				// Dissociation constant for [Na+]i and [Ca2+]sub simultaneous binding to the INaCa transporter, mM
#define Qci_SAN 0.1369 				// Fractional charge movement during the [Ca2+]sub occlusion reaction of the INaCa transporter
#define Qco_SAN 0.0 				// Fractional charge movement during the [Ca2ϩ]o occlusion reaction of the INaCa transporter
#define Qn_SAN 0.4315 				// Fractional charge movement during Naϩ occlusion reactions of the INaCa transporter
#define tdifca_SAN 0.04 			// ms
#define Ttr_SAN 40.0 				// ms
#define ConcTC_SAN 0.031 			// Total concentration of the troponin Ca2+ site, mM
#define ConcTMC_SAN 0.062 			// Total concentration of the troponin Mg2+ site, mM
#define kfTC_SAN 88.8 				// Ca2+ association constant of troponin, mM/ms
#define kfTMC_SAN 237.7 			// Ca2+ association constant of the troponin Mg2+ site, mM/ms
#define kbTC_SAN 0.446 				// Ca2+ dissociation constant of the troponin Ca2+ site, ms-1
#define kbTMC_SAN 0.00751 			// Ca2+ dissociation constant of the troponin Mg2+ site, ms-1
#define kfTMM_SAN 2.277 			// Mg2+ association constant of the troponin Mg2+ site, mM/ms
#define kbTMM_SAN 0.751 			// Mg2+ dissociation constant of the troponin Mg2+ site, mM/ms
#define ConcCM_SAN 0.045 			// Total calmodulin concentration, mM
#define kfCM_SAN 237.7 				// Ca2+ association constant of calmodulin, mM-1/ms
#define kbCM_SAN 0.542 				// Ca2+ dissociation constant of calmodulin, ms-1
#define ConcCQ_SAN 10.0 			// total calsequestrin concentration, mM 
#define kfCQ_SAN 0.534 				// Ca2+ association constant of calsequestrin, mM-1*ms-1
#define kbCQ_SAN 0.445 				// Ca2+ dissociation constant of calsequestrin, ms-1
#define koca_SAN 10.0 				// Baseline non-SR-dependent transition rate constant for the RyR, Mm-2
#define kom_SAN 0.06 				// Rate transition constant for RyR, ms-1
#define kica_SAN 0.5 				// Baseline non-SR-dependent transition rate constant for the RyR, mM-1*ms-1
#define kim_SAN 0.005 				// Rate transition constant for RyR, ms-1
#define eca50sr_SAN 0.45 			// EC50 for Ca2+ SR-dependent activation of SR Ca2+ release,
#define maxsr_SAN 15.0 				// Ca2+ modeling parameter, 
#define minsr_SAN 1.0 				// Ca2+ modeling parameter
#define hsrr_SAN 2.5 				// Parameter for Ca2+ dependent activation of SR Ca2+ release
#define pumphill_SAN 2.0 			// SR Ca2+ uptake and Hill coefficient
#define pumpkmf_SAN 0.00008
#define pumpkmr_SAN 4.5

/*Atrium Cell paramters************************************/
//Cell geometry parameters
#define A_cap_atrium 1.534e-4  //(cm2)
#define V_myo_atrium 25.84e-6
#define V_jsr_atrium 0.12e-6
#define V_nsr_atrium 2.098e-6
#define V_ss_atrium 1.485e-9   *0.7

//Extracellular ion concentrations
#define ko_atrium    5400.0
#define nao_atrium 140000.0    //(uM)
#define cao_atrium   1800.0

//SR parameter
#define v1_atrium_atrium 7.5
#define v2_atrium 1.74e-5
#define v3_atrium 0.65*1.7    //modulate cajsr and cansr
#define K_m_up_atrium 0.5
#define tau_tr_atrium 1.0  //20.0
#define tau_xfer_atrium 8.0
#define k_pos_a_atrium 0.006075
#define k_neg_a_atrium 0.07125
#define k_pos_b_atrium 0.00405
#define k_neg_b_atrium 0.965
#define k_pos_c_atrium 0.009
#define k_neg_c_atrium 0.0008
#define N_atrium 4.0
#define M_atrium 3.0

//L-type Ca2+ channel parameters
#define Gcal_atrium 0.056
#define E_cal_atrium 45.0
#define K_pc_max_atrium 0.23324
#define K_pc_half_atrium 20.0
#define K_pcb_atrium 0.0005
#define Ical_max_atrium 5.4310

//Buffering parameters
#define LTRPN_tot_atrium 70.0
#define HTRPN_tot_atrium 140.0
#define k_pos_htrpn_atrium 0.00237
#define k_neg_htrpn_atrium 3.2e-5
#define k_pos_ltrpn_atrium 0.0327
#define k_neg_ltrpn_atrium 0.0196
#define CMDN_tot_atrium 50.0
#define CSQN_tot_atrium 15000.0
#define K_CMDN_m_atrium 0.238
#define K_CSQN_m_atrium 800.0

#define Gktof_atrium 0.17556
#define Gkur_atrium 0.073125 //RA
//#define Gkur_atrium 0.14945    //LA
#define Gkss_atrium 0.0486

//Membrane current parameter
#define Cm 1.0 //(uF/cm2)
#define T_atrium 298.0
#define k_naca_atrium 878.4
#define K_m_na_atrium 87500.0
#define K_m_ca_atrium 1380.0
#define K_sat_atrium 0.1
#define eta_atrium 0.35
#define Inak_max_atrium 3.08
#define K_m_nai_atrium 21000.0
#define K_m_ko_atrium 1500.0
#define Ipca_max_atrium 1.0
#define K_m_pca_atrium 0.5
#define Gcab_atrium 0.0022754
#define Gna_atrium 13.0
#define Gnab_atrium 0.00572
#define Gks_atrium 0.00575 //(mS/uF)
#define Gclca_atrium 10.0
#define K_m_cl_atrium 10.0
#define E_cl_atrium -40.0

#endif
