// #include <Rdefines.h>
#include "rcpp_test.h"
#include <Rmath.h>

#define in_isNull   Rf_isNull
#undef in_isNull
#define in_isNull(s)    (TYPEOF(s) == NILSXP)

#include "BioEcoPar.h" // Class is defined in this file.

// using namespace Rcpp;

extern "C" {

void BioEcoPar::EcoDCF(SEXP list, int ind_t, int persCalc, double dr, int VERBOSE)
{

//ofstream fichier;
//if (ind_t ==4) fichier.open ("C:\\Users\\fbriton\\Dropbox\\These\\IAM_Dvt\\EcoDCF.txt", ios::out | ios::trunc);

//Rprintf("\nJ1\n");fichier << "J1" << endl;
    static SEXP att_dimCstF = Rf_install("DimCst");
    SEXP Flist;
    Rf_protect(Flist = getListElement(list, "Fleet"));

    Rf_protect(out_EcoDCF);

    SEXP dimCstF, DimF, dimnamesF, dimCstFM, dimCstFini, dimCstFMini, DimFM, DimFMini, dimnamesFM, dimnamesFMini; //formatage des objets resultats

    SEXP eFACTf, eFACTfm, elmt;

    SEXP    nbv_f, nbv_f_m, lc_f_m, lcd_f_m, tripLgth_f, tripLgth_f_m, nbTrip_f, nbTrip_f_m, nbds_f, nbds_f_m,
            /*effort1_f,*/ effort1_f_m, /*effort2_f,*/ effort2_f_m, Lref_f_m, cnb_f_m, ovcDCF_f_m, fc_f_m, vf_f_m, cshr_f_m, cshr_f, cnb_f, persc_f,
            eec_f, mwh_f, rep_f, gc_f, fixc_f, FTE_f, dep_f, ic_f, /*K_f, inv_f,*/ FTE_f_m, GVLref_f_m, /*ue_f,*/ ue_f_m;

    SEXP    dc_nbv_f, dc_nbv_f_m, dc_lc_f_m, dc_lcd_f_m, dc_tripLgth_f, dc_tripLgth_f_m, dc_nbTrip_f, dc_nbTrip_f_m, dc_nbds_f, dc_nbds_f_m,
            /*dc_effort1_f,*/ dc_effort1_f_m, /*dc_effort2_f,*/ dc_effort2_f_m, dc_Lref_f_m, dc_cnb_f_m, dc_ovcDCF_f_m, dc_fc_f_m, dc_vf_f_m, dc_cshr_f_m, dc_cshr_f, dc_cnb_f, dc_persc_f,
            dc_eec_f, dc_mwh_f, dc_rep_f, dc_gc_f, dc_fixc_f, dc_FTE_f, dc_dep_f, dc_ic_f, /*dc_K_f, dc_inv_f,*/ dc_FTE_f_m, dc_GVLref_f_m, /*dc_ue_f,*/ dc_ue_f_m;

    int     *dim_nbv_f, *dim_nbv_f_m, *dim_lc_f_m, *dim_tripLgth_f_m, *dim_nbTrip_f_m,
            *dim_Lref_f_m, *dim_cnb_f_m, *dim_ovcDCF_f_m, *dim_fc_f_m, *dim_vf_f_m, *dim_cshr_f_m, *dim_cshr_f, *dim_persc_f,
            *dim_eec_f, *dim_mwh_f, *dim_rep_f, *dim_gc_f, *dim_fixc_f, *dim_FTE_f, *dim_dep_f, *dim_ic_f, /* *dim_K_f, *dim_inv_f,*/ *dim_GVLref_f_m,
            *dim_ue_f_m;

    double  *r_nbv_f, *r_nbv_f_m, *r_lc_f_m, *r_lcd_f_m, *r_tripLgth_f_m, *r_nbTrip_f_m,
            *r_Lref_f_m, *r_cnb_f_m, *r_ovcDCF_f_m, *r_fc_f_m, *r_vf_f_m, *r_cshr_f_m, *r_cshr_f, *r_persc_f,
            *r_eec_f, *r_mwh_f, *r_rep_f, *r_gc_f, *r_fixc_f, *r_FTE_f, *r_dep_f, *r_ic_f, /* *r_K_f, *r_inv_f,*/ *r_GVLref_f_m;

    double  *r_ET_f_m_out,
            *r_GVLcom_f_m_e_out,*r_GVLst_f_m_e_out, *r_GVLtot_f_m_out, *r_GVLav_f_m_out, *r_GVLtot_f_out,
            *r_GVLav_f_out, *r_NGVLav_f_m_out, *r_NGVLav_f_out, *r_cnb_f_m_out, *r_cnb_f_out,
            *r_rtbs_f_m_out, *r_rtbs_f_out, *r_cshrT_f_m_out, *r_cshrT_f_out, *r_ncshr_f_out, *r_ocl_f_out, *r_cs_f_out, *r_csTot_f_out, *r_gva_f_out, *r_gvamargin_f_out,
            *r_gva_FTE_f_out, *r_ccw_f_out, *r_ccwCr_f_out, *r_wageg_f_out, *r_wagen_f_out, *r_wageg_FTE_f_out, *r_wageg_h_f_out, *r_gp_f_out, *r_gpmargin_f_out,
            *r_ncf_f_out, *r_np_f_out, *r_npmargin_f_out, *r_prof_f_out, *r_npmargin_trend_f_out, /**r_ssTot_f_out,*/ *r_ps_f_out, *r_sts_f_out, *r_BER_f_out, *r_CR_BER_f_out,
            *r_fuelEff_f_out, *r_ratio_fvol_gva_f_out, *r_ratio_gp_gva_f_out, /* *r_ratio_GVL_K_f_out, *r_ratio_gp_K_f_out, *r_RoFTA_f_out, *r_ROI_f_out, */
            /* *r_ratio_np_K_f_out, *r_ratio_GVL_cnb_ue_f_out, */
            *r_rtbsAct_f_out, *r_csAct_f_out, *r_gvaAct_f_out, *r_gpAct_f_out, *r_psAct_f_out, *r_stsAct_f_out, *r_QuotaExp_f_out;


    //Rprintf("Eco 2");fichier << "Eco2" << endl;
    //definition des dimensions


    Rf_protect(dimnamesF = Rf_allocVector(VECSXP,2));
    Rf_protect(dimnamesFM = Rf_allocVector(VECSXP,3));
    Rf_protect(dimnamesFMini = Rf_allocVector(VECSXP,2));

    SET_VECTOR_ELT(dimnamesF, 0, fleetList); SET_VECTOR_ELT(dimnamesF, 1, times);
    SET_VECTOR_ELT(dimnamesFM, 0, fleetList); SET_VECTOR_ELT(dimnamesFM, 1, metierListEco); SET_VECTOR_ELT(dimnamesFM, 2, times);
    SET_VECTOR_ELT(dimnamesFMini, 0, fleetList); SET_VECTOR_ELT(dimnamesFMini, 1, metierListEco);

    Rf_protect(dimCstF = _BioEcoPar_init_DimCst(nbF, 0, 0, nbT));
    Rf_protect(dimCstFini = _BioEcoPar_init_DimCst(nbF, 0, 0, 0));
    Rf_protect(dimCstFM = _BioEcoPar_init_DimCst(nbF, nbMe, 0, nbT));
    Rf_protect(dimCstFMini = _BioEcoPar_init_DimCst(nbF, nbMe, 0, 0));
    Rf_protect(DimF = _BioEcoPar_n_DimCst(dimCstF));
    Rf_protect(DimFM = _BioEcoPar_n_DimCst(dimCstFM));
    Rf_protect(DimFMini = _BioEcoPar_n_DimCst(dimCstFMini));

    // facteurs des indices generiques F/FM

    Rf_protect(eFACTf = iDim(INTEGER(dimCstF)));
    Rf_protect(eFACTfm = iDim(INTEGER(dimCstFM)));
    int *eF_f = INTEGER(eFACTf);
    int *eF_fm = INTEGER(eFACTfm);

    Rf_protect(nbv_f = getListElement(Flist, "nbv_f"));                Rf_protect(dc_nbv_f = iDim(INTEGER(Rf_getAttrib(nbv_f, att_dimCstF)))); //Rprintf("Eco 31");
    Rf_protect(nbv_f_m = getListElement(Flist, "nbv_f_m"));            Rf_protect(dc_nbv_f_m = iDim(INTEGER(Rf_getAttrib(nbv_f_m, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(lc_f_m = getListElement(Flist, "lc_f_m"));              Rf_protect(dc_lc_f_m = iDim(INTEGER(Rf_getAttrib(lc_f_m, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(lcd_f_m = getListElement(Flist, "lcd_f_m"));            Rf_protect(dc_lcd_f_m = iDim(INTEGER(Rf_getAttrib(lcd_f_m, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(tripLgth_f = getListElement(Flist, "tripLgth_f"));      Rf_protect(dc_tripLgth_f = iDim(INTEGER(Rf_getAttrib(tripLgth_f, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(tripLgth_f_m = getListElement(Flist, "tripLgth_f_m"));  Rf_protect(dc_tripLgth_f_m = iDim(INTEGER(Rf_getAttrib(tripLgth_f_m, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(nbTrip_f = getListElement(Flist, "nbTrip_f"));          Rf_protect(dc_nbTrip_f = iDim(INTEGER(Rf_getAttrib(nbTrip_f, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(nbTrip_f_m = getListElement(Flist, "nbTrip_f_m"));      Rf_protect(dc_nbTrip_f_m = iDim(INTEGER(Rf_getAttrib(nbTrip_f_m, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(nbds_f = getListElement(Flist, "nbds_f"));              Rf_protect(dc_nbds_f = iDim(INTEGER(Rf_getAttrib(nbds_f, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(nbds_f_m = getListElement(Flist, "nbds_f_m"));          Rf_protect(dc_nbds_f_m = iDim(INTEGER(Rf_getAttrib(nbds_f_m, att_dimCstF))));//Rprintf("Eco 31");
    // Rf_protect(effort1_f = getListElement(Flist, "effort1_f"));        Rf_protect(dc_effort1_f = iDim(INTEGER(Rf_getAttrib(effort1_f, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(effort1_f_m = getListElement(Flist, "effort1_f_m"));    Rf_protect(dc_effort1_f_m = iDim(INTEGER(Rf_getAttrib(effort1_f_m, att_dimCstF))));//Rprintf("Eco 31");
    // Rf_protect(effort2_f = getListElement(Flist, "effort2_f"));        Rf_protect(dc_effort2_f = iDim(INTEGER(Rf_getAttrib(effort2_f, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(effort2_f_m = getListElement(Flist, "effort2_f_m"));    Rf_protect(dc_effort2_f_m = iDim(INTEGER(Rf_getAttrib(effort2_f_m, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(Lref_f_m = getListElement(Flist, "Lref_f_m"));          Rf_protect(dc_Lref_f_m = iDim(INTEGER(Rf_getAttrib(Lref_f_m, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(cnb_f_m = getListElement(Flist, "cnb_f_m"));            Rf_protect(dc_cnb_f_m = iDim(INTEGER(Rf_getAttrib(cnb_f_m, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(ovcDCF_f_m = getListElement(Flist, "ovcDCF_f_m"));      Rf_protect(dc_ovcDCF_f_m = iDim(INTEGER(Rf_getAttrib(ovcDCF_f_m, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(fc_f_m = getListElement(Flist, "fc_f_m"));              Rf_protect(dc_fc_f_m = iDim(INTEGER(Rf_getAttrib(fc_f_m, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(vf_f_m = getListElement(Flist, "vf_f_m"));              Rf_protect(dc_vf_f_m = iDim(INTEGER(Rf_getAttrib(vf_f_m, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(cshr_f_m = getListElement(Flist, "cshr_f_m"));          Rf_protect(dc_cshr_f_m = iDim(INTEGER(Rf_getAttrib(cshr_f_m, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(cshr_f = getListElement(Flist, "cshr_f"));              Rf_protect(dc_cshr_f = iDim(INTEGER(Rf_getAttrib(cshr_f, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(cnb_f = getListElement(Flist, "cnb_f"));                Rf_protect(dc_cnb_f = iDim(INTEGER(Rf_getAttrib(cnb_f, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(persc_f = getListElement(Flist, "persc_f"));            Rf_protect(dc_persc_f = iDim(INTEGER(Rf_getAttrib(persc_f, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(eec_f = getListElement(Flist, "eec_f"));                Rf_protect(dc_eec_f = iDim(INTEGER(Rf_getAttrib(eec_f, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(mwh_f = getListElement(Flist, "mwh_f"));                Rf_protect(dc_mwh_f = iDim(INTEGER(Rf_getAttrib(mwh_f, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(rep_f = getListElement(Flist, "rep_f"));                Rf_protect(dc_rep_f = iDim(INTEGER(Rf_getAttrib(rep_f, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(gc_f = getListElement(Flist, "gc_f"));                  Rf_protect(dc_gc_f = iDim(INTEGER(Rf_getAttrib(gc_f, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(fixc_f = getListElement(Flist, "fixc_f"));              Rf_protect(dc_fixc_f = iDim(INTEGER(Rf_getAttrib(fixc_f, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(FTE_f = getListElement(Flist, "FTE_f"));                Rf_protect(dc_FTE_f = iDim(INTEGER(Rf_getAttrib(FTE_f, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(dep_f = getListElement(Flist, "dep_f"));                Rf_protect(dc_dep_f = iDim(INTEGER(Rf_getAttrib(dep_f, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(ic_f = getListElement(Flist, "ic_f"));                  Rf_protect(dc_ic_f = iDim(INTEGER(Rf_getAttrib(ic_f, att_dimCstF))));//Rprintf("Eco 31");
    // Rf_protect(K_f = getListElement(Flist, "K_f"));                    Rf_protect(dc_K_f = iDim(INTEGER(Rf_getAttrib(K_f, att_dimCstF))));//Rprintf("Eco 31");
    // Rf_protect(inv_f = getListElement(Flist, "inv_f"));                Rf_protect(dc_inv_f = iDim(INTEGER(Rf_getAttrib(inv_f, att_dimCstF))));//Rprintf("Eco 31");
    Rf_protect(FTE_f_m = getListElement(Flist, "FTE_f_m"));            Rf_protect(dc_FTE_f_m = iDim(INTEGER(Rf_getAttrib(FTE_f_m, att_dimCstF))));//Rprintf("Eco 31"); // TODO : not used
    Rf_protect(GVLref_f_m = getListElement(Flist, "GVLref_f_m"));      Rf_protect(dc_GVLref_f_m = iDim(INTEGER(Rf_getAttrib(GVLref_f_m, att_dimCstF))));//Rprintf("Eco 31");

    // ---> P = 14 + 35*2 = 84
    //Rprintf("Eco 4");fichier << "Eco4" << endl;

    // Rf_protect(ue_f = Rf_allocVector(REALSXP,nbF));
    // Rf_setAttrib(ue_f, R_DimSymbol, Rf_getAttrib(getListElement(Flist, "effort1_f"), R_DimSymbol));
    // Rf_setAttrib(ue_f, R_DimNamesSymbol, Rf_getAttrib(getListElement(Flist, "effort1_f"), R_DimNamesSymbol));
    // Rf_setAttrib(ue_f, att_dimCstF, Rf_getAttrib(getListElement(Flist, "effort1_f"), att_dimCstF));

    Rf_protect(ue_f_m = Rf_allocVector(REALSXP,nbF*nbMe));
    Rf_setAttrib(ue_f_m, R_DimSymbol, Rf_getAttrib(getListElement(Flist, "effort1_f_m"), R_DimSymbol));
    Rf_setAttrib(ue_f_m, R_DimNamesSymbol, Rf_getAttrib(getListElement(Flist, "effort1_f_m"), R_DimNamesSymbol));
    Rf_setAttrib(ue_f_m, att_dimCstF, Rf_getAttrib(getListElement(Flist, "effort1_f_m"), att_dimCstF));

    // double *r_ue_f = REAL(ue_f); double *reff1_f = REAL(getListElement(Flist, "effort1_f")) ; double *reff2_f = REAL(getListElement(Flist, "effort2_f"));
    double *r_ue_f_m = REAL(ue_f_m); double *reff1 = REAL(getListElement(Flist, "effort1_f_m")) ; double *reff2 = REAL(getListElement(Flist, "effort2_f_m"));

    for (int ind_f = 0 ; ind_f < nbF ; ind_f++) {
        // r_ue_f[ind_f] = reff1_f[ind_f]*reff2_f[ind_f];
        for (int ind_m = 0 ; ind_m < nbMe ; ind_m++) r_ue_f_m[ind_f + nbF*ind_m] = reff1[ind_f + nbF*ind_m]*reff2[ind_f + nbF*ind_m];
    }

    //Rprintf("Eco 5");fichier << "Eco5" << endl;
    // Rf_protect(dc_ue_f = iDim(INTEGER(Rf_getAttrib(ue_f, att_dimCstF))));
    Rf_protect(dc_ue_f_m = iDim(INTEGER(Rf_getAttrib(ue_f_m, att_dimCstF))));
    dim_ue_f_m = INTEGER(dc_ue_f_m);
    // ---> P = 84 + 4 = 88

    dim_nbv_f = INTEGER(dc_nbv_f);                          r_nbv_f = REAL(nbv_f);//Rprintf("Eco 51\n");
    dim_nbv_f_m = INTEGER(dc_nbv_f_m);                      r_nbv_f_m = REAL(nbv_f_m);//Rprintf("Eco 51\n");
    dim_lc_f_m = INTEGER(dc_lc_f_m);                        r_lc_f_m = REAL(lc_f_m);//Rprintf("Eco 51\n");
    /*dim_lcd_f_m = INTEGER(dc_lcd_f_m); */                 r_lcd_f_m = REAL(lcd_f_m);//Rprintf("Eco 51\n");
    //dim_tripLgth_f = INTEGER(dc_tripLgth_f);                r_tripLgth_f = REAL(tripLgth_f);//Rprintf("Eco 51\n");
    dim_tripLgth_f_m = INTEGER(dc_tripLgth_f_m);            r_tripLgth_f_m = REAL(tripLgth_f_m);//Rprintf("Eco 51\n");
    //dim_nbTrip_f = INTEGER(dc_nbTrip_f);                    r_nbTrip_f = REAL(nbTrip_f);//Rprintf("Eco 51\n");
    dim_nbTrip_f_m = INTEGER(dc_nbTrip_f_m);                r_nbTrip_f_m = REAL(nbTrip_f_m);//Rprintf("Eco 51\n");
    //dim_nbds_f = INTEGER(dc_nbds_f);                        r_nbds_f = REAL(nbds_f);//Rprintf("Eco 51\n");
    //dim_nbds_f_m = INTEGER(dc_nbds_f_m);                    r_nbds_f_m = REAL(nbds_f_m);//Rprintf("Eco 51\n");
    //dim_effort1_f = INTEGER(dc_effort1_f);                  r_effort1_f = REAL(effort1_f);//Rprintf("Eco 51\n");
    //dim_effort1_f_m = INTEGER(dc_effort1_f_m);              r_effort1_f_m = REAL(effort1_f_m);//Rprintf("Eco 51\n");
    //dim_effort2_f = INTEGER(dc_effort2_f);                  r_effort2_f = REAL(effort2_f);//Rprintf("Eco 51\n");
    //dim_effort2_f_m = INTEGER(dc_effort2_f_m);              r_effort2_f_m = REAL(effort2_f_m);//Rprintf("Eco 51\n");
    dim_Lref_f_m = INTEGER(dc_Lref_f_m);                    r_Lref_f_m = REAL(Lref_f_m);//Rprintf("Eco 51\n");
    dim_cnb_f_m = INTEGER(dc_cnb_f_m);                      r_cnb_f_m = REAL(cnb_f_m);//Rprintf("Eco 51\n");
    dim_ovcDCF_f_m = INTEGER(dc_ovcDCF_f_m);                r_ovcDCF_f_m = REAL(ovcDCF_f_m);//Rprintf("Eco 51\n");
    dim_fc_f_m = INTEGER(dc_fc_f_m);                        r_fc_f_m = REAL(fc_f_m);//Rprintf("Eco 51\n");
    dim_vf_f_m = INTEGER(dc_vf_f_m);                        r_vf_f_m = REAL(vf_f_m);//Rprintf("Eco 51\n");
    dim_cshr_f_m = INTEGER(dc_cshr_f_m);                    r_cshr_f_m = REAL(cshr_f_m);//Rprintf("Eco 51\n");
    dim_cshr_f = INTEGER(dc_cshr_f);                        r_cshr_f = REAL(cshr_f);//Rprintf("Eco 51\n");
    //dim_cnb_f = INTEGER(dc_cnb_f);                          r_cnb_f = REAL(cnb_f);//Rprintf("Eco 51\n");
    dim_persc_f = INTEGER(dc_persc_f);                      r_persc_f = REAL(persc_f);//Rprintf("Eco 51\n");
    dim_eec_f = INTEGER(dc_eec_f);                          r_eec_f = REAL(eec_f);//Rprintf("Eco 51\n");
    dim_mwh_f = INTEGER(dc_mwh_f);                          r_mwh_f = REAL(mwh_f);//Rprintf("Eco 51\n");
    dim_rep_f = INTEGER(dc_rep_f);                          r_rep_f = REAL(rep_f);//Rprintf("Eco 51\n");
    dim_gc_f = INTEGER(dc_gc_f);                            r_gc_f = REAL(gc_f);//Rprintf("Eco 51\n");
    dim_fixc_f = INTEGER(dc_fixc_f);                        r_fixc_f = REAL(fixc_f);//Rprintf("Eco 51\n");
    dim_FTE_f = INTEGER(dc_FTE_f);                          r_FTE_f = REAL(FTE_f);//Rprintf("Eco 51\n");
    dim_dep_f = INTEGER(dc_dep_f);                          r_dep_f = REAL(dep_f);//Rprintf("Eco 51\n");
    dim_ic_f = INTEGER(dc_ic_f);                            r_ic_f = REAL(ic_f);//Rprintf("Eco 51\n");
    // dim_K_f = INTEGER(dc_K_f);                              r_K_f = REAL(K_f);//Rprintf("Eco 51\n");
    // dim_inv_f = INTEGER(dc_inv_f);                          r_inv_f = REAL(inv_f);//Rprintf("Eco 51\n");
    //dim_FTE_f_m = INTEGER(dc_FTE_f_m);                    r_FTE_f_m = REAL(FTE_f_m);//Rprintf("Eco 51\n");
    dim_GVLref_f_m = INTEGER(dc_GVLref_f_m);                r_GVLref_f_m = REAL(GVLref_f_m);//Rprintf("Eco 51\n");


    int nbC=0;
    //int nbI=0;

    //Rprintf("Eco 6");fichier << "Eco6" << endl;

    if (ind_t==0) {

        SEXP ETini_f_m, fvolue_f_m, ovcDCFue_f_m, rtbsIni_f, ccwr_f, opersc_f, eco_names,
             GVLcom_f_m_e_out, GVLcom_f_m_eStat_out, GVLcom_f_m_e, GVLst_f_m_e_out, GVLst_f_m_eStat_out, GVLst_f_m_e, GVL_f_m_e_out, GVL_f_m_eStat_out, GVLtot_f_m_e,
             GVLtot_f_m_out, GVLav_f_m_out, /*GVLtot_f_out, GVLav_f_out,*/ NGVLav_f_m_out, /*NGVLav_f_out,*/ ET_f_m_out,
             cnb_f_m_out, /*cnb_f_out,*/ rtbs_f_m_out, /*rtbs_f_out, rtbsAct_f_out,*/ cshrT_f_m_out, 
            //  cshrT_f_out, ncshr_f_out, ocl_f_out, cs_f_out, csAct_f_out, csTot_f_out,
            //  gva_f_out, gvaAct_f_out, gvamargin_f_out, gva_FTE_f_out, ccw_f_out, ccwCr_f_out, wageg_f_out, wagen_f_out, wageg_FTE_f_out, wageg_h_f_out,
            //  gp_f_out, gpAct_f_out, gpmargin_f_out, ncf_f_out, np_f_out, npmargin_f_out, prof_f_out, npmargin_trend_f_out, 
            //  ssTot_f_out, ps_f_out, psAct_f_out, sts_f_out, stsAct_f_out, BER_f_out, CR_BER_f_out, fuelEff_f_out, 
            //  ratio_fvol_gva_f_out, ratio_gp_gva_f_out, ratio_GVL_K_f_out, ratio_gp_K_f_out, RoFTA_f_out, ROI_f_out, ratio_np_K_f_out, ratio_GVL_cnb_ue_f_out,
             QuotaExp_f_out;

        double  *r_ETini_f_m, *r_fvolue_f_m, *r_ovcDCFue_f_m, *r_rtbsIni_f, *r_ccwr_f, *r_opersc_f;

        //-------------------------
        // Stade preliminaire (temps initial)
        //-------------------------

            Rf_protect(ETini_f_m = Rf_allocVector(REALSXP,nbF*nbMe));                 r_ETini_f_m = REAL(ETini_f_m);  // TODO : numericMatrix here
            Rf_protect(fvolue_f_m = Rf_allocVector(REALSXP,nbF*nbMe));                r_fvolue_f_m = REAL(fvolue_f_m);
            Rf_protect(ovcDCFue_f_m = Rf_allocVector(REALSXP,nbF*nbMe));              r_ovcDCFue_f_m = REAL(ovcDCFue_f_m);
            Rf_protect(rtbsIni_f = Rf_allocVector(REALSXP,nbF));                      r_rtbsIni_f = REAL(rtbsIni_f);
            Rf_protect(ccwr_f = Rf_allocVector(REALSXP,nbF));                         r_ccwr_f = REAL(ccwr_f);
            Rf_protect(opersc_f = Rf_allocVector(REALSXP,nbF));                       r_opersc_f = REAL(opersc_f);
        // ---> P(t0) = 6
        //Rprintf("Eco 7");fichier << "Eco7" << endl;

        // on cree ETini
        //    double *rnbTrip = REAL(getListElement(Flist, "nbTrip_f_m"));
        //    double *rtripLgth = REAL(getListElement(Flist, "tripLgth_f_m"));
        //    double *rnbv = REAL(getListElement(Flist, "nbv_f_m"));
        //    double *rcnb = REAL(getListElement(Flist, "cnb_f_m"));

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++) {
            for (int ind_m = 0 ; ind_m < nbMe ; ind_m++) {
                r_ETini_f_m[ind_f + nbF*ind_m] =
                finite(
                r_Lref_f_m[ind_f*dim_Lref_f_m[0] + ind_m*dim_Lref_f_m[1] + 0*dim_Lref_f_m[2] + ind_t*dim_Lref_f_m[3]] /
                (r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]] *
                    r_nbTrip_f_m[ind_f*dim_nbTrip_f_m[0] + ind_m*dim_nbTrip_f_m[1] + 0*dim_nbTrip_f_m[2] + ind_t*dim_nbTrip_f_m[3]] *
                    r_tripLgth_f_m[ind_f*dim_tripLgth_f_m[0] + ind_m*dim_tripLgth_f_m[1] + 0*dim_tripLgth_f_m[2] + ind_t*dim_tripLgth_f_m[3]] *
                    r_cnb_f_m[ind_f*dim_cnb_f_m[0] + ind_m*dim_cnb_f_m[1] + 0*dim_cnb_f_m[2] + ind_t*dim_cnb_f_m[3]]));
                //if (ISNA(r_ETini_f_m[ind_f + nbF*ind_m])) r_ETini_f_m[ind_f + nbF*ind_m] = 0.0;
            }
        }


        //Rprintf("Eco 8");fichier << "Eco8" << endl;
        for (int e = 0 ; e < nbE+nbEstat ; e++) {

            if (e<nbE) {
             Rf_protect(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));
            } else {
             Rf_protect(elmt = getListElement(list, CHAR(STRING_ELT(sppListStat,e-nbE))));
            }

            //if (e<nbE) nbI = Rf_length(getListElement(elmt, "modI"));
            if (e<nbE) nbC = Rf_length(getListElement(elmt, "modC"));

            Rf_protect(GVLtot_f_m_e = Rf_allocVector(REALSXP,nbF*nbMe*nbT));
            Rf_protect(GVLcom_f_m_e = Rf_allocVector(REALSXP,nbF*nbMe*nbT));
            Rf_protect(GVLst_f_m_e = Rf_allocVector(REALSXP,nbF*nbMe*nbT));

            double *r_GVLtot_f_m_e = REAL(GVLtot_f_m_e);
            double *r_GVLcom_f_m_e = REAL(GVLcom_f_m_e);
            double *r_GVLst_f_m_e = REAL(GVLst_f_m_e);

            double *r_Lbio_f_m_e ,  *r_P_f_m_e, r_Pst_e=NA_REAL, *r_LD_efmc=&NA_REAL, *r_statLDor_efm=&NA_REAL, *r_statLDst_efm=&NA_REAL, r_theta_e;
            int *dim_Lbio_e, *dim_P_e;
            //Rprintf("Eco 9");fichier << "Eco9" << endl;
            if (e<nbE) {
                r_Lbio_f_m_e = REAL(VECTOR_ELT(out_L_efmct, e));
                //r_Lbio_f_sum_e = REAL(aggregObj(VECTOR_ELT(out_L_efmct, e),dimCstF));
                r_P_f_m_e = REAL(VECTOR_ELT(out_P_t, e));
                r_LD_efmc = REAL(VECTOR_ELT(out_LD_efmc, e));
                r_theta_e = REAL(getListElement(elmt, "theta_e"))[0];
                dim_Lbio_e = INTEGER(iDim(INTEGER(Rf_getAttrib(VECTOR_ELT(out_L_efmct, e), att_dimCstF))));
                dim_P_e = INTEGER(iDim(INTEGER(Rf_getAttrib(VECTOR_ELT(out_P_t, e), att_dimCstF))));
            } else {
                r_Lbio_f_m_e = REAL(VECTOR_ELT(out_Lstat, e-nbE));
                //r_Lbio_f_sum_e = REAL(aggregObj(VECTOR_ELT(out_Lstat, e-nbE),dimCstF));
                r_P_f_m_e = REAL(VECTOR_ELT(out_Pstat, e-nbE));
                r_statLDor_efm = REAL(VECTOR_ELT(out_statLDor_efm, e-nbE));
                r_statLDst_efm = REAL(VECTOR_ELT(out_statLDst_efm, e-nbE));
                r_theta_e = REAL(getListElement(elmt, "theta_e"))[0];
                dim_Lbio_e = INTEGER(iDim(INTEGER(Rf_getAttrib(VECTOR_ELT(out_Lstat, e-nbE), att_dimCstF))));
                dim_P_e = INTEGER(iDim(INTEGER(Rf_getAttrib(VECTOR_ELT(out_Pstat, e-nbE), att_dimCstF))));
                r_Pst_e = REAL(getListElement(elmt, "Pst_e"))[0];
            }

            //------------------------------
            //Equations de la table "p"
            //------------------------------
            //Rprintf("Eco 10");fichier << "Eco10" << endl;

            for (int ind_f = 0 ; ind_f < nbF ; ind_f++){   //on rappelle ici que ind_t est en fait egal a 0

            //double countGVLtotf = 0.0; //pour sommer GVLtot_f_m_e sur les metiers

                for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

                    //-- 3. GVLtot_f_m_e

                    double countCom = 0.0;

                    if (e<nbE) {

                        if (ISNA(r_theta_e)) r_theta_e = 1.0;

                        for (int ind_c = 0 ; ind_c < (nbC-1) ; ind_c++){ //sur les classes non sous-tailles

                            if (ISNA(r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + ind_c*dim_P_e[2] + ind_t*dim_P_e[3]])){
                                r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + ind_c*dim_P_e[2] + ind_t*dim_P_e[3]] = 0.0;
                            }

                            if (!ISNA(r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]])) {
                                countCom = countCom +
                                r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + ind_c*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 * //prix au kg
                                r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]] +
                                r_theta_e * r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + ind_c*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 * //prix au kg
                                r_LD_efmc[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]];
                            }

                        }

                        if (ISNA(r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + (nbC-1)*dim_P_e[2] + ind_t*dim_P_e[3]])){
                            r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + (nbC-1)*dim_P_e[2] + ind_t*dim_P_e[3]] = 0.0;
                        }

                        //if (ISNA(r_theta_e)) r_theta_e = 1.0;

                        if (!ISNA(r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + (nbC-1)*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]])) {

                            r_GVLst_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                            r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + (nbC-1)*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 * //prix au kg
                            r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + (nbC-1)*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]] +
                            r_theta_e * r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + (nbC-1)*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 * //prix au kg
                            r_LD_efmc[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + (nbC-1)*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]];

                        } else {
                            r_GVLst_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] = 0.0;
                        }

                        //if (e==1 & ind_f==0 & ind_m==4) PrintValue(ETini_f_m_out);

                    } else {

                        if (ISNA(r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + 0*dim_P_e[2] + ind_t*dim_P_e[3]])){
                            r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + 0*dim_P_e[2] + ind_t*dim_P_e[3]] = 0.0;
                        }

                        if (ISNA(r_theta_e)) r_theta_e = 1.0;

                        if (!ISNA(r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]])) {
                            countCom = r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + 0*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 * //prix au kg
                                r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]] +
                                r_theta_e * r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + 0*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 *
                                finite(r_statLDor_efm[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]]) ;
                        }

                        if (!ISNA(r_statLDst_efm[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]])) {

                            if (ISNA(r_Pst_e)) r_Pst_e = 0.0;

                            r_GVLst_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                                r_Pst_e * 1000 * r_statLDst_efm[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]];

                        }

                    }

                    r_GVLcom_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] = countCom;
                    r_GVLtot_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_GVLcom_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] +
                    r_GVLst_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

            //Rprintf("Eco 11");fichier << "Eco11" << endl;


                }
            }

            //on formatte le(s) resultat(s) et on les integre a 'eVar'
            //Rprintf("Eco 13");fichier << "Eco13" << endl;
            Rf_setAttrib(GVLtot_f_m_e, R_DimSymbol, DimFM);
            Rf_setAttrib(GVLtot_f_m_e, R_DimNamesSymbol, dimnamesFM);
            Rf_setAttrib(GVLtot_f_m_e, att_dimCstF, dimCstFM);
            if (e<nbE) SET_VECTOR_ELT(VECTOR_ELT(eVar, e), 41, GVLtot_f_m_e); else SET_VECTOR_ELT(VECTOR_ELT(eStatVar, e-nbE), 1, GVLtot_f_m_e);

            Rf_setAttrib(GVLcom_f_m_e, R_DimSymbol, DimFM);
            Rf_setAttrib(GVLcom_f_m_e, R_DimNamesSymbol, dimnamesFM);
            Rf_setAttrib(GVLcom_f_m_e, att_dimCstF, dimCstFM);
            if (e<nbE) SET_VECTOR_ELT(VECTOR_ELT(eVar, e), 246, GVLcom_f_m_e); else SET_VECTOR_ELT(VECTOR_ELT(eStatVar, e-nbE), 8, GVLcom_f_m_e);

            Rf_setAttrib(GVLst_f_m_e, R_DimSymbol, DimFM);
            Rf_setAttrib(GVLst_f_m_e, R_DimNamesSymbol, dimnamesFM);
            Rf_setAttrib(GVLst_f_m_e, att_dimCstF, dimCstFM);
            if (e<nbE) SET_VECTOR_ELT(VECTOR_ELT(eVar, e), 247, GVLst_f_m_e); else SET_VECTOR_ELT(VECTOR_ELT(eStatVar, e-nbE), 9, GVLst_f_m_e);

            //Rprintf("Eco 14");fichier << "Eco14" << endl;
            Rf_unprotect(4);

        }


        //Rprintf("Eco 8");fichier << "Eco8" << endl;

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

            double countRTBSnum = 0.0;

            for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

                countRTBSnum = countRTBSnum + finite(r_GVLref_f_m[ind_f*dim_GVLref_f_m[0] + ind_m*dim_GVLref_f_m[1] + 0*dim_GVLref_f_m[2] + ind_t*dim_GVLref_f_m[3]] *
                                            r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]] *
                                            (1 - finite(r_lc_f_m[ind_f*dim_lc_f_m[0] + ind_m*dim_lc_f_m[1] + 0*dim_lc_f_m[2] + ind_t*dim_lc_f_m[3]]))) -
                                            finite(r_ovcDCF_f_m[ind_f*dim_ovcDCF_f_m[0] + ind_m*dim_ovcDCF_f_m[1] + 0*dim_ovcDCF_f_m[2] + ind_t*dim_ovcDCF_f_m[3]]) -
                                            finite(r_fc_f_m[ind_f*dim_fc_f_m[0] + ind_m*dim_fc_f_m[1] + 0*dim_fc_f_m[2] + ind_t*dim_fc_f_m[3]]);

                //-- 4. fvolue_f_m

                r_fvolue_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                            finite( r_fc_f_m[ind_f*dim_fc_f_m[0] + ind_m*dim_fc_f_m[1] + 0*dim_fc_f_m[2] + ind_t*dim_fc_f_m[3]] /
                            (r_vf_f_m[ind_f*dim_vf_f_m[0] + ind_m*dim_vf_f_m[1] + 0*dim_vf_f_m[2] + ind_t*dim_vf_f_m[3]] *
                            r_ue_f_m[ind_f + ind_m*nbF]) );

                //-- 5. ovcDCFue_f_m

                r_ovcDCFue_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                            finite( r_ovcDCF_f_m[ind_f*dim_ovcDCF_f_m[0] + ind_m*dim_ovcDCF_f_m[1] + 0*dim_ovcDCF_f_m[2] + ind_t*dim_ovcDCF_f_m[3]] /
                            r_ue_f_m[ind_f + ind_m*nbF] );

            }

            //Rprintf("Eco 18");fichier << "Eco18" << endl;

            //-- 6. rtbsIni_f

            r_rtbsIni_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                countRTBSnum / r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];


            //-- 7. ccwr_f

            r_ccwr_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_persc_f[ind_f*dim_persc_f[0] + 0*dim_persc_f[1] + 0*dim_persc_f[2] + ind_t*dim_persc_f[3]] /
                r_rtbsIni_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

            //-- 8. opersc_f

            r_opersc_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_persc_f[ind_f*dim_persc_f[0] + 0*dim_persc_f[1] + 0*dim_persc_f[2] + ind_t*dim_persc_f[3]] -
                (0.01 * r_cshr_f[ind_f*dim_cshr_f[0] + 0*dim_cshr_f[1] + 0*dim_cshr_f[2] + ind_t*dim_cshr_f[3]] * //cshr_f en %
                r_rtbsIni_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]]);

        }

        //Rprintf("Eco 19\n");fichier << "Eco19" << endl;

        //on formatte le(s) resultat(s) et on integre a fVar

        Rf_setAttrib(fvolue_f_m, R_DimSymbol, DimFMini);
        Rf_setAttrib(fvolue_f_m, R_DimNamesSymbol, dimnamesFMini);
        Rf_setAttrib(fvolue_f_m, att_dimCstF, dimCstFMini);
        SET_VECTOR_ELT(fVar, 4, fvolue_f_m);

        Rf_setAttrib(ovcDCFue_f_m, R_DimSymbol, DimFMini);
        Rf_setAttrib(ovcDCFue_f_m, R_DimNamesSymbol, dimnamesFMini);
        Rf_setAttrib(ovcDCFue_f_m, att_dimCstF, dimCstFMini);
        SET_VECTOR_ELT(fVar, 10, ovcDCFue_f_m);

        Rf_setAttrib(ccwr_f, R_NamesSymbol, fleetList);
        Rf_setAttrib(ccwr_f, att_dimCstF, dimCstFini);
        SET_VECTOR_ELT(fVar, 27, ccwr_f);

        Rf_setAttrib(opersc_f, R_NamesSymbol,  fleetList);
        Rf_setAttrib(opersc_f, att_dimCstF, dimCstFini);
        SET_VECTOR_ELT(fVar, 28, opersc_f);

        SET_VECTOR_ELT(fVar, 31, rtbsIni_f);

        SET_VECTOR_ELT(fVar, 33, ETini_f_m);
        //Rprintf("Eco 20\n");fichier << "Eco20" << endl;
        //enfin, on initialise l'output


        Rf_protect(GVLcom_f_m_e_out = Rf_allocVector(VECSXP, nbE));
        Rf_setAttrib(GVLcom_f_m_e_out, R_NamesSymbol, sppList);
        SET_VECTOR_ELT(out_EcoDCF, 0, GVLcom_f_m_e_out);//Rprintf("Eco 20\n");

        Rf_protect(GVLcom_f_m_eStat_out = Rf_allocVector(VECSXP, nbEstat));
        Rf_setAttrib(GVLcom_f_m_eStat_out, R_NamesSymbol, sppListStat);
        SET_VECTOR_ELT(out_EcoDCF, 1, GVLcom_f_m_eStat_out);//Rprintf("Eco 20\n");

        Rf_protect(GVLst_f_m_e_out = Rf_allocVector(VECSXP, nbE));
        Rf_setAttrib(GVLst_f_m_e_out, R_NamesSymbol, sppList);
        SET_VECTOR_ELT(out_EcoDCF, 2, GVLst_f_m_e_out);//Rprintf("Eco 20\n");

        Rf_protect(GVLst_f_m_eStat_out = Rf_allocVector(VECSXP, nbEstat));
        Rf_setAttrib(GVLst_f_m_eStat_out, R_NamesSymbol, sppListStat);
        SET_VECTOR_ELT(out_EcoDCF, 3, GVLst_f_m_eStat_out);//Rprintf("Eco 20\n");

        Rf_protect(GVL_f_m_e_out = Rf_allocVector(VECSXP, nbE));
        Rf_setAttrib(GVL_f_m_e_out, R_NamesSymbol, sppList);
        SET_VECTOR_ELT(out_EcoDCF, 4, GVL_f_m_e_out);//Rprintf("Eco 20\n");

        Rf_protect(GVL_f_m_eStat_out = Rf_allocVector(VECSXP, nbEstat));
        Rf_setAttrib(GVL_f_m_eStat_out, R_NamesSymbol, sppListStat);
        SET_VECTOR_ELT(out_EcoDCF, 5, GVL_f_m_eStat_out);//Rprintf("Eco 20\n");

        Rf_protect(GVLtot_f_m_out = Rf_allocVector(REALSXP,nbF*nbMe*nbT));
        Rf_setAttrib(GVLtot_f_m_out, R_DimSymbol, DimFM);
        Rf_setAttrib(GVLtot_f_m_out, R_DimNamesSymbol, dimnamesFM);
        Rf_setAttrib(GVLtot_f_m_out, att_dimCstF, dimCstFM);
        SET_VECTOR_ELT(out_EcoDCF, 6, GVLtot_f_m_out);//Rprintf("Eco 20\n");

        Rf_protect(GVLav_f_m_out = Rf_allocVector(REALSXP,nbF*nbMe*nbT));
        Rf_setAttrib(GVLav_f_m_out, R_DimSymbol, DimFM);
        Rf_setAttrib(GVLav_f_m_out, R_DimNamesSymbol, dimnamesFM);
        Rf_setAttrib(GVLav_f_m_out, att_dimCstF, dimCstFM);
        SET_VECTOR_ELT(out_EcoDCF, 7, GVLav_f_m_out);//Rprintf("Eco 20\n");

        // Rf_protect(GVLtot_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(GVLtot_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(GVLtot_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(GVLtot_f_out, att_dimCstF, dimCstF);
        SET_VECTOR_ELT(out_EcoDCF, 8, _BioEcoPar_init_matrix(nbF, nbT,  dimnamesF, dimCstF));//Rprintf("Eco 20\n");

        // Rf_protect(GVLav_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(GVLav_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(GVLav_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(GVLav_f_out, att_dimCstF, dimCstF);
        SET_VECTOR_ELT(out_EcoDCF, 9, _BioEcoPar_init_matrix(nbF, nbT,  dimnamesF, dimCstF));//Rprintf("Eco 20\n");

        Rf_protect(NGVLav_f_m_out = Rf_allocVector(REALSXP,nbF*nbMe*nbT));
        Rf_setAttrib(NGVLav_f_m_out, R_DimSymbol, DimFM);
        Rf_setAttrib(NGVLav_f_m_out, R_DimNamesSymbol, dimnamesFM);
        Rf_setAttrib(NGVLav_f_m_out, att_dimCstF, dimCstFM);
        SET_VECTOR_ELT(out_EcoDCF, 10, NGVLav_f_m_out);//Rprintf("Eco 20\n");

        // Rf_protect(NGVLav_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(NGVLav_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(NGVLav_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(NGVLav_f_out, att_dimCstF, dimCstF);
        SET_VECTOR_ELT(out_EcoDCF, 11, _BioEcoPar_init_matrix(nbF, nbT,  dimnamesF, dimCstF));//Rprintf("Eco 20\n");

        Rf_protect(ET_f_m_out = Rf_allocVector(REALSXP,nbF*nbMe*nbT));
        Rf_setAttrib(ET_f_m_out, R_DimSymbol, DimFM);
        Rf_setAttrib(ET_f_m_out, R_DimNamesSymbol, dimnamesFM);
        Rf_setAttrib(ET_f_m_out, att_dimCstF, dimCstFM);
        SET_VECTOR_ELT(out_EcoDCF, 12, ET_f_m_out);//Rprintf("Eco 20\n");

        Rf_protect(cnb_f_m_out = Rf_allocVector(REALSXP,nbF*nbMe*nbT));
        Rf_setAttrib(cnb_f_m_out, R_DimSymbol, DimFM);
        Rf_setAttrib(cnb_f_m_out, R_DimNamesSymbol, dimnamesFM);
        Rf_setAttrib(cnb_f_m_out, att_dimCstF, dimCstFM);
        SET_VECTOR_ELT(out_EcoDCF, 13, cnb_f_m_out);//Rprintf("Eco 20\n");

        // Rf_protect(cnb_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(cnb_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(cnb_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(cnb_f_out, att_dimCstF, dimCstF);
        SET_VECTOR_ELT(out_EcoDCF, 14, _BioEcoPar_init_matrix(nbF, nbT,  dimnamesF, dimCstF));//Rprintf("Eco 20\n");

        Rf_protect(rtbs_f_m_out = Rf_allocVector(REALSXP,nbF*nbMe*nbT));
        Rf_setAttrib(rtbs_f_m_out, R_DimSymbol, DimFM);
        Rf_setAttrib(rtbs_f_m_out, R_DimNamesSymbol, dimnamesFM);
        Rf_setAttrib(rtbs_f_m_out, att_dimCstF, dimCstFM);
        SET_VECTOR_ELT(out_EcoDCF, 15, rtbs_f_m_out);

        // Rf_protect(rtbs_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(rtbs_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(rtbs_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(rtbs_f_out, att_dimCstF, dimCstF);
        SET_VECTOR_ELT(out_EcoDCF, 16, _BioEcoPar_init_matrix(nbF, nbT,  dimnamesF, dimCstF));//Rprintf("Eco 20\n");

        // Rf_protect(rtbsAct_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(rtbsAct_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(rtbsAct_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(rtbsAct_f_out, att_dimCstF, dimCstF);
        SET_VECTOR_ELT(out_EcoDCF, 17, _BioEcoPar_init_matrix(nbF, nbT,  dimnamesF, dimCstF));

        Rf_protect(cshrT_f_m_out = Rf_allocVector(REALSXP,nbF*nbMe*nbT));
        Rf_setAttrib(cshrT_f_m_out, R_DimSymbol, DimFM);
        Rf_setAttrib(cshrT_f_m_out, R_DimNamesSymbol, dimnamesFM);
        Rf_setAttrib(cshrT_f_m_out, att_dimCstF, dimCstFM);
        SET_VECTOR_ELT(out_EcoDCF, 18, cshrT_f_m_out);

        // TODO : rcpp replace loop from here !
        {
        // Rf_protect(cshrT_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(cshrT_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(cshrT_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(cshrT_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 19, cshrT_f_out);//Rprintf("Eco 20\n");

        // Rf_protect(ncshr_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(ncshr_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(ncshr_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(ncshr_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 20, ncshr_f_out);//Rprintf("Eco 20\n");

        // Rf_protect(ocl_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(ocl_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(ocl_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(ocl_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 21, ocl_f_out);//Rprintf("Eco 20\n");

        // Rf_protect(cs_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(cs_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(cs_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(cs_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 22, cs_f_out);//Rprintf("Eco 20\n");

        // Rf_protect(csAct_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(csAct_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(csAct_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(csAct_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 23, csAct_f_out);

        // Rf_protect(csTot_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(csTot_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(csTot_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(csTot_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 24, csTot_f_out);

        // Rf_protect(gva_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(gva_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(gva_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(gva_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 25, gva_f_out);//Rprintf("Eco 20\n");

        // Rf_protect(gvaAct_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(gvaAct_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(gvaAct_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(gvaAct_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 26, gvaAct_f_out);//Rprintf("Eco 20\n");

        // Rf_protect(gvamargin_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(gvamargin_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(gvamargin_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(gvamargin_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 27, gvamargin_f_out);//Rprintf("Eco 20\n");

        // Rf_protect(gva_FTE_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(gva_FTE_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(gva_FTE_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(gva_FTE_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 28, gva_FTE_f_out);//Rprintf("Eco 20\n");

        // Rf_protect(ccw_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(ccw_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(ccw_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(ccw_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 29, ccw_f_out);//Rprintf("Eco 20\n");

        // Rf_protect(ccwCr_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(ccwCr_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(ccwCr_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(ccwCr_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 30, ccwCr_f_out);//Rprintf("Eco 20\n");

        // Rf_protect(wageg_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(wageg_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(wageg_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(wageg_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 31, wageg_f_out);//Rprintf("Eco 20\n");

        // Rf_protect(wagen_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(wagen_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(wagen_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(wagen_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 32, wagen_f_out);//Rprintf("Eco 20\n");

        // Rf_protect(wageg_FTE_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(wageg_FTE_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(wageg_FTE_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(wageg_FTE_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 33, wageg_FTE_f_out);//Rprintf("Eco 20\n");

        // Rf_protect(wageg_h_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(wageg_h_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(wageg_h_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(wageg_h_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 34, wageg_h_f_out);//Rprintf("Eco 20\n");

        // Rf_protect(gp_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(gp_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(gp_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(gp_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 35, gp_f_out);//Rprintf("Eco 20\n");

        // Rf_protect(gpAct_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(gpAct_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(gpAct_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(gpAct_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 36, gpAct_f_out);//Rprintf("Eco 20\n");

        // Rf_protect(gpmargin_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(gpmargin_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(gpmargin_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(gpmargin_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 37, _BioEcoPar_init_matrix(nbF, nbT,  dimnamesF, dimCstF));//Rprintf("Eco 20\n");

        // Rf_protect(ncf_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(ncf_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(ncf_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(ncf_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 38, ncf_f_out);//Rprintf("Eco 20\n");

        // Rf_protect(np_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(np_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(np_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(np_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 39, np_f_out);//Rprintf("Eco 20\n");

        // Rf_protect(npmargin_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(npmargin_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(npmargin_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(npmargin_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 40, npmargin_f_out);//Rprintf("Eco 20\n");

        // Rf_protect(prof_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(prof_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(prof_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(prof_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 41, prof_f_out);//Rprintf("Eco 20\n");

        // Rf_protect(npmargin_trend_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(npmargin_trend_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(npmargin_trend_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(npmargin_trend_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 42, npmargin_trend_f_out);//Rprintf("Eco 20\n");
    
        // Rf_protect(ssTot_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(ssTot_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(ssTot_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(ssTot_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 43, ssTot_f_out);
        // SET_VECTOR_ELT(out_EcoDCF, 43, _BioEcoPar_init_matrix(nbF, nbT,  dimnamesF, dimCstF));

        // Rf_protect(ps_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(ps_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(ps_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(ps_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 44, ps_f_out);

        // Rf_protect(psAct_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(psAct_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(psAct_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(psAct_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 45, _BioEcoPar_init_matrix(nbF, nbT,  dimnamesF, dimCstF));

        // Rf_protect(sts_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(sts_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(sts_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(sts_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 46, sts_f_out);

        // Rf_protect(stsAct_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(stsAct_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(stsAct_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(stsAct_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 47, stsAct_f_out);

        // Rf_protect(BER_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(BER_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(BER_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(BER_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 48, BER_f_out);

        // Rf_protect(CR_BER_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(CR_BER_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(CR_BER_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(CR_BER_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 49, CR_BER_f_out);

        // Rf_protect(fuelEff_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(fuelEff_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(fuelEff_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(fuelEff_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 50, fuelEff_f_out);

        // Rf_protect(ratio_fvol_gva_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(ratio_fvol_gva_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(ratio_fvol_gva_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(ratio_fvol_gva_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 51, ratio_fvol_gva_f_out);

        // Rf_protect(ratio_gp_gva_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(ratio_gp_gva_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(ratio_gp_gva_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(ratio_gp_gva_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 52, ratio_gp_gva_f_out);

        // Rf_protect(ratio_GVL_K_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(ratio_GVL_K_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(ratio_GVL_K_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(ratio_GVL_K_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 53, ratio_GVL_K_f_out);

        // Rf_protect(ratio_gp_K_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(ratio_gp_K_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(ratio_gp_K_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(ratio_gp_K_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 54, ratio_gp_K_f_out);

        // Rf_protect(RoFTA_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(RoFTA_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(RoFTA_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(RoFTA_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 55, RoFTA_f_out);

        // Rf_protect(ROI_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(ROI_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(ROI_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(ROI_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 56, ROI_f_out);

        // Rf_protect(ratio_np_K_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(ratio_np_K_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(ratio_np_K_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(ratio_np_K_f_out, att_dimCstF, dimCstF);

        // Rf_protect(ratio_GVL_cnb_ue_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        // Rf_setAttrib(ratio_GVL_cnb_ue_f_out, R_DimSymbol, DimF);
        // Rf_setAttrib(ratio_GVL_cnb_ue_f_out, R_DimNamesSymbol, dimnamesF);
        // Rf_setAttrib(ratio_GVL_cnb_ue_f_out, att_dimCstF, dimCstF);
        // SET_VECTOR_ELT(out_EcoDCF, 58, ratio_GVL_cnb_ue_f_out);
    }

        for (int rcpp = 19; rcpp < 59; rcpp++){
            SET_VECTOR_ELT(out_EcoDCF, rcpp, _BioEcoPar_init_matrix(nbF, nbT,  dimnamesF, dimCstF));
        }

        Rf_protect(QuotaExp_f_out = Rf_allocVector(REALSXP,nbF*nbT));
        Rf_setAttrib(QuotaExp_f_out, R_DimSymbol, DimF);
        Rf_setAttrib(QuotaExp_f_out, R_DimNamesSymbol, dimnamesF);
        Rf_setAttrib(QuotaExp_f_out, att_dimCstF, dimCstF);
        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
            for (int ind_tt = 0 ; ind_tt < nbT ; ind_tt++) REAL(QuotaExp_f_out)[ind_f + nbF*ind_tt] = 0.0;
        SET_VECTOR_ELT(out_EcoDCF, 59, QuotaExp_f_out);

        //Rprintf("Eco 20.8\n");fichier << "Eco 20.8" << endl;
        //on nomme les elements de out_EcoDCF




        const char *namesEco[60] = {"GVLcom_f_m_e_out","GVLcom_f_m_eStat_out","GVLst_f_m_e_out","GVLst_f_m_eStat_out","GVL_f_m_e_out","GVL_f_m_eStat_out","GVLtot_f_m_out",
                                    "GVLav_f_m_out","GVLtot_f_out","GVLav_f_out","NGVLav_f_m_out","NGVLav_f_out","ET_f_m_out","cnb_f_m_out","cnb_f_out","rtbs_f_m_out","rtbs_f_out",
                                    "rtbsAct_f_out","cshrT_f_m_out","cshrT_f_out","ncshr_f_out","ocl_f_out","cs_f_out","csAct_f_out","csTot_f_out","gva_f_out","gvaAct_f_out",
                                    "gvamargin_f_out","gva_FTE_f_out","ccw_f_out","ccwCr_f_out","wageg_f_out","wagen_f_out","wageg_FTE_f_out","wageg_h_f_out","gp_f_out",
                                    "gpAct_f_out","gpmargin_f_out","ncf_f_out","np_f_out","npmargin_f_out","prof_f_out","npmargin_trend_f_out","ssTot_f_out","ps_f_out",
                                    "psAct_f_out","sts_f_out","stsAct_f_out","BER_f_out","CR_BER_f_out","fuelEff_f_out","ratio_fvol_gva_f_out","ratio_gp_gva_f_out",
                                    "ratio_GVL_K_f_out","ratio_gp_K_f_out","RoFTA_f_out","ROI_f_out","ratio_np_K_f_out","ratio_GVL_cnb_ue_f_out","QuotaExp_f_out"};

        //Rprintf("Eco 20.9\n");fichier << "Eco 20.9" << endl;
            Rf_protect(eco_names = Rf_allocVector(STRSXP, 60));

            for(int ct = 0; ct < 60; ct++) SET_STRING_ELT(eco_names, ct, Rf_mkChar(namesEco[ct])); //PrintValue(out_EcoDCF);
        //Rprintf("Eco 20.95\n");fichier << "Eco20.95" << endl;
            Rf_setAttrib(out_EcoDCF, R_NamesSymbol, eco_names);
        //Rprintf("Eco 21\n");fichier << "Eco21" << endl;

        // ---> P(t0) = 6 + 60 = 66
    }


    // on importe les outputs afin de les mettre a jour a l'instant ind_t

    // r_GVLcom_f_m_e_out = REAL(VECTOR_ELT(out_EcoDCF, 0));//Rprintf("Eco 20\n");
    // r_GVLcom_f_m_eStat_out = REAL(VECTOR_ELT(out_EcoDCF, 1));//Rprintf("Eco 20\n");
    // r_GVLst_f_m_e_out = REAL(VECTOR_ELT(out_EcoDCF, 2));//Rprintf("Eco 20\n");
    // r_GVLst_f_m_eStat_out = REAL(VECTOR_ELT(out_EcoDCF, 3));//Rprintf("Eco 20\n");
    // r_GVL_f_m_e_out = REAL(VECTOR_ELT(out_EcoDCF, 4));//Rprintf("Eco 20\n");
    // r_GVL_f_m_eStat_out = REAL(VECTOR_ELT(out_EcoDCF, 5));//Rprintf("Eco 20\n");
    r_GVLtot_f_m_out = REAL(VECTOR_ELT(out_EcoDCF, 6));//Rprintf("Eco 20\n");
    r_GVLav_f_m_out = REAL(VECTOR_ELT(out_EcoDCF, 7));//Rprintf("Eco 20\n");
    r_GVLtot_f_out = REAL(VECTOR_ELT(out_EcoDCF, 8));//Rprintf("Eco 20\n");
    r_GVLav_f_out = REAL(VECTOR_ELT(out_EcoDCF, 9));//Rprintf("Eco 20\n");
    r_NGVLav_f_m_out = REAL(VECTOR_ELT(out_EcoDCF, 10));//Rprintf("Eco 20\n");
    r_NGVLav_f_out = REAL(VECTOR_ELT(out_EcoDCF, 11));//Rprintf("Eco 20\n");
    r_ET_f_m_out = REAL(VECTOR_ELT(out_EcoDCF, 12));//Rprintf("Eco 20\n");
    r_cnb_f_m_out = REAL(VECTOR_ELT(out_EcoDCF, 13));//Rprintf("Eco 20\n");
    r_cnb_f_out = REAL(VECTOR_ELT(out_EcoDCF, 14));//Rprintf("Eco 20\n");
    r_rtbs_f_m_out = REAL(VECTOR_ELT(out_EcoDCF, 15));
    r_rtbs_f_out = REAL(VECTOR_ELT(out_EcoDCF, 16));//Rprintf("Eco 20\n");
    r_rtbsAct_f_out = REAL(VECTOR_ELT(out_EcoDCF, 17));
    r_cshrT_f_m_out = REAL(VECTOR_ELT(out_EcoDCF, 18));
    r_cshrT_f_out = REAL(VECTOR_ELT(out_EcoDCF, 19));//Rprintf("Eco 20\n");
    r_ncshr_f_out = REAL(VECTOR_ELT(out_EcoDCF, 20));//Rprintf("Eco 20\n");
    r_ocl_f_out = REAL(VECTOR_ELT(out_EcoDCF, 21));//Rprintf("Eco 20\n");
    r_cs_f_out = REAL(VECTOR_ELT(out_EcoDCF, 22));//Rprintf("Eco 20\n");
    r_csAct_f_out = REAL(VECTOR_ELT(out_EcoDCF, 23));
    r_csTot_f_out = REAL(VECTOR_ELT(out_EcoDCF, 24));
    r_gva_f_out = REAL(VECTOR_ELT(out_EcoDCF, 25));//Rprintf("Eco 20\n");
    r_gvaAct_f_out = REAL(VECTOR_ELT(out_EcoDCF, 26));//Rprintf("Eco 20\n");
    r_gvamargin_f_out = REAL(VECTOR_ELT(out_EcoDCF, 27));//Rprintf("Eco 20\n");
    r_gva_FTE_f_out = REAL(VECTOR_ELT(out_EcoDCF, 28));//Rprintf("Eco 20\n");
    r_ccw_f_out = REAL(VECTOR_ELT(out_EcoDCF, 29));//Rprintf("Eco 20\n");
    r_ccwCr_f_out = REAL(VECTOR_ELT(out_EcoDCF, 30));//Rprintf("Eco 20\n");
    r_wageg_f_out = REAL(VECTOR_ELT(out_EcoDCF, 31));//Rprintf("Eco 20\n");
    r_wagen_f_out = REAL(VECTOR_ELT(out_EcoDCF, 32));//Rprintf("Eco 20\n");
    r_wageg_FTE_f_out = REAL(VECTOR_ELT(out_EcoDCF, 33));//Rprintf("Eco 20\n");
    r_wageg_h_f_out = REAL(VECTOR_ELT(out_EcoDCF, 34));//Rprintf("Eco 20\n");
    r_gp_f_out = REAL(VECTOR_ELT(out_EcoDCF, 35));//Rprintf("Eco 20\n");
    r_gpAct_f_out = REAL(VECTOR_ELT(out_EcoDCF, 36));//Rprintf("Eco 20\n");
    r_gpmargin_f_out = REAL(VECTOR_ELT(out_EcoDCF, 37));//Rprintf("Eco 20\n");
    r_ncf_f_out = REAL(VECTOR_ELT(out_EcoDCF, 38));//Rprintf("Eco 20\n");
    r_np_f_out = REAL(VECTOR_ELT(out_EcoDCF, 39));//Rprintf("Eco 20\n");
    r_npmargin_f_out = REAL(VECTOR_ELT(out_EcoDCF, 40));//Rprintf("Eco 20\n");
    r_prof_f_out = REAL(VECTOR_ELT(out_EcoDCF, 41));//Rprintf("Eco 20\n");
    r_npmargin_trend_f_out = REAL(VECTOR_ELT(out_EcoDCF, 42));//Rprintf("Eco 20\n");
    // r_ssTot_f_out = REAL(VECTOR_ELT(out_EcoDCF, 43));
    r_ps_f_out = REAL(VECTOR_ELT(out_EcoDCF, 44));
    r_psAct_f_out = REAL(VECTOR_ELT(out_EcoDCF, 45));
    r_sts_f_out = REAL(VECTOR_ELT(out_EcoDCF, 46));
    r_stsAct_f_out = REAL(VECTOR_ELT(out_EcoDCF, 47));
    r_BER_f_out = REAL(VECTOR_ELT(out_EcoDCF, 48));
    r_CR_BER_f_out = REAL(VECTOR_ELT(out_EcoDCF, 49));
    r_fuelEff_f_out = REAL(VECTOR_ELT(out_EcoDCF, 50));
    r_ratio_fvol_gva_f_out = REAL(VECTOR_ELT(out_EcoDCF, 51));
    r_ratio_gp_gva_f_out = REAL(VECTOR_ELT(out_EcoDCF, 52));
    // r_ratio_GVL_K_f_out = REAL(VECTOR_ELT(out_EcoDCF, 53));
    // r_ratio_gp_K_f_out = REAL(VECTOR_ELT(out_EcoDCF, 54));
    // r_RoFTA_f_out = REAL(VECTOR_ELT(out_EcoDCF, 55));
    // r_ROI_f_out = REAL(VECTOR_ELT(out_EcoDCF, 56));
    // r_ratio_np_K_f_out = REAL(VECTOR_ELT(out_EcoDCF, 57));
    // r_ratio_GVL_cnb_ue_f_out = REAL(VECTOR_ELT(out_EcoDCF, 58));
    r_QuotaExp_f_out = REAL(VECTOR_ELT(out_EcoDCF, 59));

    // Rprintf("Eco 22\n");fichier << "Eco22" << endl;

    double *r_fvolue_f_m2 = REAL(VECTOR_ELT(fVar,4));
    double *r_ovcDCFue_f_m2 = REAL(VECTOR_ELT(fVar,10));
    double *r_ccwr_f2 = REAL(VECTOR_ELT(fVar,27));
    double *r_opersc_f2 = REAL(VECTOR_ELT(fVar,28));

    SEXP countLf;
    Rf_protect(countLf = Rf_allocVector(REALSXP,nbF)); // --> 67
    double *r_countLf = REAL(countLf);
    for (int INd_f = 0 ; INd_f < nbF ; INd_f++) r_countLf[INd_f] = 0.0; // pour le calcul de 'fuelEff'

    // ---> P = 88 + 1 = 89



    double *rnbv = REAL(getListElement(Flist, "nbv_f_m"));
    double *rnbv_f = REAL(getListElement(Flist, "nbv_f"));
    double *rnbTrip = REAL(getListElement(Flist, "nbTrip_f_m"));
    double *rtripLgth = REAL(getListElement(Flist, "tripLgth_f_m"));
    double *rnbTrip_f = REAL(getListElement(Flist, "nbTrip_f"));
    double *rtripLgth_f = REAL(getListElement(Flist, "tripLgth_f"));


    //Rprintf("Eco 23\n");fichier << "Eco23" << endl;


    //initialisation de cnb, GVLtot et NGVLav, et remplissage de ET


    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        r_cnb_f_out[ind_f + nbF*ind_t] = 0.0;

        for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

            if (ind_t==0) {
              r_ET_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = REAL(VECTOR_ELT(fVar,33))[ind_f + nbF*ind_m];
            } else {
              r_ET_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = r_ET_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*(ind_t-1)];
            }

            r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = 0.001;

            r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] = 0.0;

            r_NGVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] = 0.0;
        }
    }



    // indicateurs especes ---------------------------------------------------------

    for (int e = 0 ; e < nbE+nbEstat ; e++) {//on assume qu'il y a au moins une espece modelisee, qu'elle soit dynamique ou non --> pas de if

        if (e<nbE) {
         Rf_protect(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e)))); //espece dynamique
        } else {
         Rf_protect(elmt = getListElement(list, CHAR(STRING_ELT(sppListStat,e-nbE)))); //espece statique
        }
        //Rprintf("Eco X10\n");fichier << "EcoX10" << endl;
        //if (e<nbE) nbI = Rf_length(getListElement(elmt, "modI"));
        if (e<nbE) nbC = Rf_length(getListElement(elmt, "modC"));

        double *r_Lbio_f_sum_t_e=&NA_REAL, *r_GVLtot_f_m_e2=&NA_REAL, *r_Lbio_f_m_e=&NA_REAL, *r_P_f_m_e=&NA_REAL, *r_LD_efmc=&NA_REAL, r_theta_e=NA_REAL, *r_statLDor_efm=&NA_REAL,
               *r_statLDst_efm=&NA_REAL, r_Pst_e=NA_REAL, *r_LD_f_sum_t_e=&NA_REAL, *r_statLDor_f_sum_t_e=&NA_REAL, *r_statLDst_f_sum_t_e=&NA_REAL;
        int *dim_Lbio_e, *dim_P_e;

        //Rprintf("Eco X11\n");fichier << "EcoX11" << endl;
        if ((nbE>0) & (e<nbE)) {
            SEXP gg1=R_NilValue, gg2=R_NilValue, Pgg1=R_NilValue, Pgg2=R_NilValue;
            Rf_protect(Pgg1=VECTOR_ELT(out_L_efmct, e));
            Rf_protect(Pgg2=VECTOR_ELT(out_LD_efmc, e));
            Rf_protect(gg1=aggregObj(Pgg1,dimCstF));
            Rf_protect(gg2=aggregObj(Pgg2,dimCstF));

            r_GVLtot_f_m_e2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, e),41));
            r_GVLcom_f_m_e_out = REAL(VECTOR_ELT(VECTOR_ELT(eVar, e),246));
            r_GVLst_f_m_e_out = REAL(VECTOR_ELT(VECTOR_ELT(eVar, e),247));
            r_Lbio_f_m_e = REAL(VECTOR_ELT(out_L_efmct, e));
            r_Lbio_f_sum_t_e = REAL(gg1);
            r_LD_efmc = REAL(VECTOR_ELT(out_LD_efmc, e));
            r_LD_f_sum_t_e = REAL(gg2);
            r_theta_e = REAL(getListElement(elmt, "theta_e"))[0];
            //r_Lref_f_e = REAL(getListElement(elmt, "Lref_f_e"));
            r_P_f_m_e = REAL(VECTOR_ELT(out_P_t, e));
            dim_Lbio_e = INTEGER(iDim(INTEGER(Rf_getAttrib(VECTOR_ELT(out_L_efmct, e), att_dimCstF))));
            dim_P_e = INTEGER(iDim(INTEGER(Rf_getAttrib(VECTOR_ELT(out_P_t, e), att_dimCstF))));
            //Rprintf("Eco X12\n");fichier << "EcoX12" << endl;
            Rf_unprotect(4);
        }
        if ((nbEstat>0) & (e>=nbE)) {

            SEXP gg1=R_NilValue, gg2=R_NilValue, gg3=R_NilValue, Pgg1=R_NilValue, Pgg2=R_NilValue, Pgg3=R_NilValue;
            Rf_protect(Pgg1=VECTOR_ELT(out_Lstat, e-nbE));
            Rf_protect(Pgg2=VECTOR_ELT(out_statLDor_efm, e-nbE));
            Rf_protect(Pgg3=VECTOR_ELT(out_statLDst_efm, e-nbE));
            Rf_protect(gg1=aggregObj(Pgg1,dimCstF));
            Rf_protect(gg2=aggregObj(Pgg2,dimCstF));
            Rf_protect(gg3=aggregObj(Pgg3,dimCstF));

            r_GVLtot_f_m_e2 = REAL(VECTOR_ELT(VECTOR_ELT(eStatVar, e-nbE),1));
            r_GVLcom_f_m_e_out = REAL(VECTOR_ELT(VECTOR_ELT(eStatVar, e-nbE),8));
            r_GVLst_f_m_e_out = REAL(VECTOR_ELT(VECTOR_ELT(eStatVar, e-nbE),9));
            r_Lbio_f_m_e = REAL(VECTOR_ELT(out_Lstat, e-nbE));
            r_Lbio_f_sum_t_e = REAL(gg1);
            r_statLDor_efm = REAL(VECTOR_ELT(out_statLDor_efm, e-nbE));
            r_statLDor_f_sum_t_e = REAL(gg2);
            r_statLDst_efm = REAL(VECTOR_ELT(out_statLDst_efm, e-nbE));
            r_statLDst_f_sum_t_e = REAL(gg3);
            r_theta_e = REAL(getListElement(elmt, "theta_e"))[0];
            //r_Lref_f_e = REAL(getListElement(elmt, "Lref_f_e"));
            r_P_f_m_e = REAL(VECTOR_ELT(out_Pstat, e-nbE));
            dim_Lbio_e = INTEGER(iDim(INTEGER(Rf_getAttrib(VECTOR_ELT(out_Lstat, e-nbE), att_dimCstF))));
            dim_P_e = INTEGER(iDim(INTEGER(Rf_getAttrib(VECTOR_ELT(out_Pstat, e-nbE), att_dimCstF))));
            r_Pst_e = REAL(getListElement(elmt, "Pst_e"))[0];
            //Rprintf("Eco X13\n");fichier << "EcoX13" << endl;
            Rf_unprotect(6);
        }



        //---------------------
        //Equations de la table "t"
        //---------------------

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

            if (e<nbE) {

            r_countLf[ind_f] = r_countLf[ind_f] + finite(r_Lbio_f_sum_t_e[ind_f]) + finite(r_LD_f_sum_t_e[ind_f]);

            } else {

            r_countLf[ind_f] = r_countLf[ind_f] + finite(r_Lbio_f_sum_t_e[ind_f]) + finite(r_statLDor_f_sum_t_e[ind_f]) + finite(r_statLDst_f_sum_t_e[ind_f]);

            }
            ////Rprintf("Eco X131\n");fichier << "EcoX131" << endl;
            for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

                //-- 1. GVL_f_m_e

                double countCom = 0.0;
                ////Rprintf("Eco X132\n");fichier << "EcoX132" << endl;

                if (e<nbE) {

                    if (ISNA(r_theta_e)) r_theta_e = 1.0;

                    for (int ind_c = 0 ; ind_c < (nbC-1) ; ind_c++){ //sur les classes non sous-tailles

                        if (ISNA(r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + ind_c*dim_P_e[2] + ind_t*dim_P_e[3]])){
                            r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + ind_c*dim_P_e[2] + ind_t*dim_P_e[3]] = 0.0;
                        }

                        if (!ISNA(r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]])) {

                            countCom = countCom +
                                r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + ind_c*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 * //prix au kg
                                r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]] +
                                r_theta_e * r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + ind_c*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 * //prix au kg
                                finite(r_LD_efmc[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]]);

                            r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] +
                                r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]] +
                                finite(r_LD_efmc[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]]);

                        }
                    }

                    ////Rprintf("Eco X133\n");fichier << "EcoX133" << endl;
                    if (ISNA(r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + (nbC-1)*dim_P_e[2] + ind_t*dim_P_e[3]])){
                        r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + (nbC-1)*dim_P_e[2] + ind_t*dim_P_e[3]] = 0.0;
                    }

                    ////Rprintf("Eco X1331\n");fichier << "EcoX1331" << endl;
                    if (!ISNA(r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + (nbC-1)*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]])) {
                        ////Rprintf("Eco X1332\n");fichier << "EcoX1332" << endl;
                        r_GVLst_f_m_e_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                            r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + (nbC-1)*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 * //prix au kg
                            r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + (nbC-1)*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]] +
                            r_theta_e * r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + (nbC-1)*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 * //prix au kg
                            finite(r_LD_efmc[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + (nbC-1)*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]]);

                        r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] +
                            r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + (nbC-1)*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]] +
                            finite(r_LD_efmc[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + (nbC-1)*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]]);
                        ////Rprintf("Eco X1333\n");fichier << "EcoX1333" << endl;
                    } else {
                        ////Rprintf("Eco X1334\n");fichier << "EcoX1334" << endl;
                        r_GVLst_f_m_e_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] = 0.0;
                        ////Rprintf("Eco X1335\n");fichier << "EcoX1335" << endl;
                    }

                } else {

                    ////Rprintf("Eco X134\n");fichier << "EcoX134" << endl;

                    if (ISNA(r_theta_e)) r_theta_e = 1.0;

                    if (ISNA(r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + 0*dim_P_e[2] + ind_t*dim_P_e[3]])){
                        r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + 0*dim_P_e[2] + ind_t*dim_P_e[3]] = 0.0;
                    }

                    if (!ISNA(r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]])){

                        countCom = r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + 0*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 * //prix au kg
                            r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]] +
                            r_theta_e * r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + 0*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 *
                            finite(r_statLDor_efm[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]]) ;

                        r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] +
                            r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]] +
                            finite(r_statLDor_efm[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]]);
                    }

                    if (!ISNA(r_statLDst_efm[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]])){

                        if (ISNA(r_Pst_e)) r_Pst_e = 0.0;

                        r_GVLst_f_m_e_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                            r_Pst_e * 1000 * finite(r_statLDst_efm[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]]);

                        r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] +
                            finite(r_statLDst_efm[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]]);

                    } else {
                        r_GVLst_f_m_e_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] = 0.0;
                    }
                    ////Rprintf("Eco X135\n");fichier << "EcoX135" << endl;
                }

                ////Rprintf("Eco X1338\n");fichier << "EcoX1338" << endl;
                r_GVLcom_f_m_e_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] = countCom;
                ////Rprintf("Eco X1339\n");fichier << "EcoX1339" << endl;
                r_GVLtot_f_m_e2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                r_GVLcom_f_m_e_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] + r_GVLst_f_m_e_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];
                ////Rprintf("Eco X136\n");fichier << "EcoX136" << endl;

                //-- 2. GVLtot_f_m


                r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] +
                    finite(r_GVLtot_f_m_e2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]]);

                double LC = 0.0, LCD = 0.0;
                if (!ISNA(r_lc_f_m[ind_f*dim_lc_f_m[0] + ind_m*dim_lc_f_m[1] + 0*dim_lc_f_m[2] + ind_t*dim_lc_f_m[3]])){
                    LC = r_lc_f_m[ind_f*dim_lc_f_m[0] + ind_m*dim_lc_f_m[1] + 0*dim_lc_f_m[2] + ind_t*dim_lc_f_m[3]];
                }
                if (!ISNA(r_lcd_f_m[ind_f*dim_lc_f_m[0] + ind_m*dim_lc_f_m[1] + 0*dim_lc_f_m[2] + ind_t*dim_lc_f_m[3]])){
                    LCD = r_lcd_f_m[ind_f*dim_lc_f_m[0] + ind_m*dim_lc_f_m[1] + 0*dim_lc_f_m[2] + ind_t*dim_lc_f_m[3]];
                }
                        ////Rprintf("Eco X137\n");fichier << "EcoX137" << endl;
                r_NGVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                r_NGVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] +
                finite(r_GVLcom_f_m_e_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]]) * (1 - 0.01*LC) +
                finite(r_GVLst_f_m_e_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]]) * (1 - 0.01*LCD);


            }
        }
        ////Rprintf("Eco X14\n");fichier << "EcoX14" << endl;

        ////Rprintf("aa1"); fichier << "aa1" << endl;
        //if (ind_t==4 & e==26) PrintValue(VECTOR_ELT(out_EcoDCF, 13));
        ////Rprintf("aa1");


        if (e<nbE) {
            SET_VECTOR_ELT(VECTOR_ELT(out_EcoDCF,4), e, VECTOR_ELT(VECTOR_ELT(eVar, e),41));//Rprintf("Eco X1\n");
            SET_VECTOR_ELT(VECTOR_ELT(out_EcoDCF,0), e, VECTOR_ELT(VECTOR_ELT(eVar, e),246));//Rprintf("Eco X2\n");
            SET_VECTOR_ELT(VECTOR_ELT(out_EcoDCF,2), e, VECTOR_ELT(VECTOR_ELT(eVar, e),247));//Rprintf("Eco X3\n");
        } else {
            SET_VECTOR_ELT(VECTOR_ELT(out_EcoDCF,5), e-nbE, VECTOR_ELT(VECTOR_ELT(eStatVar, e-nbE),1));//Rprintf("Eco X4\n");
            SET_VECTOR_ELT(VECTOR_ELT(out_EcoDCF,1), e-nbE, VECTOR_ELT(VECTOR_ELT(eStatVar, e-nbE),8));//Rprintf("Eco X5\n");
            SET_VECTOR_ELT(VECTOR_ELT(out_EcoDCF,3), e-nbE, VECTOR_ELT(VECTOR_ELT(eStatVar, e-nbE),9));//Rprintf("Eco X6\n");
        }

        Rf_unprotect(1);

        //Rprintf("\nJ3\n");fichier << "J3" << endl;

        ////Rprintf("aa1"); fichier << "aa1" << endl;
        //if (ind_t==4 & e==26) PrintValue(VECTOR_ELT(out_EcoDCF, 13));
        ////Rprintf("aa1");
        //if (ind_t==4) PrintValue(VECTOR_ELT(out_EcoDCF, 13));
    }

    // Calcul quota costs
    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){ //Reinitialisation si plusieurs appels au module
        r_QuotaExp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = 0.0;
    }

    SEXP nam_eQuota, PQuot_et,QuotaTrade_fe;
    double *r_PQuot_et, *r_QuotaTrade_fe;
    for (int eQuota = 0 ; eQuota < nbEQuotaMarket ; eQuota++) {

        Rf_protect(nam_eQuota = STRING_ELT(sppListQM,eQuota));
        Rf_protect(PQuot_et = getListElement(out_PQuot_et, CHAR(nam_eQuota)));

        // if(!in_isNull(getListElement(out_QuotaTrade_fe, CHAR(nam_eQuota)))){ //explicitely traded: amount traded = landings - holdings
        //     Rf_protect(QuotaTrade_fe = getListElement(out_QuotaTrade_fe, CHAR(nam_eQuota)));
        //     r_QuotaTrade_fe = REAL(QuotaTrade_fe);
        // } else { //otherwise amount traded = landings
        if (!in_isNull (getListElement(out_L_efmit, CHAR(nam_eQuota)))){ // espece dyn
            Rf_protect(QuotaTrade_fe = aggregObj(getListElement(out_L_efmit, CHAR(nam_eQuota)),dimCstF));
        } else{ // espece stat
            Rf_protect(QuotaTrade_fe = aggregObj(getListElement(out_Lstat, CHAR(nam_eQuota)),dimCstF));
        }
        r_QuotaTrade_fe = REAL(QuotaTrade_fe);
        // }

        r_PQuot_et = REAL(PQuot_et);

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
            r_QuotaExp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = r_QuotaExp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] +
                                                                                        r_PQuot_et[ind_t] * r_QuotaTrade_fe[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] * 1000;
            // if (ind_f==6){
            //     fichier << "e = " << CHAR(nam_eQuota) <<
            //                 "; Pquot = " << r_PQuot_et[ind_t] <<
            //                 "; Traded amount = " << r_QuotaTrade_fe[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]]*1000 <<
            //                 "; QuotaExp = " << r_QuotaExp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] << endl;
            // }

        }
        Rf_unprotect(3);
    }


    //Rprintf("aa");
    //if (ind_t==4) PrintValue(VECTOR_ELT(out_EcoDCF, 13));
    //Rprintf("bb");//Rprintf("%f\n",r_cnb_f_m_out[0 + nbF*7 + nbF*nbMe*4]);//Rprintf("%f\n",r_ET_f_m_out[0 + nbF*7 + nbF*nbMe*4]);
    ////Rprintf("%f\n",rnbv[0 + nbF*7]);//Rprintf("%f\n",rnbTrip[0 + nbF*7]);//Rprintf("%f\n",rtripLgth[0 + nbF*7]);
    //if (ind_t==1) PrintValue(VECTOR_ELT(out_EcoDCF, 12));
    //Rprintf("cc");

    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){
            ////Rprintf("Eco X1336\n");fichier << "EcoX1336" << endl; //a ce moment, cnb contient les debarquements totaux par flottille et metier
            r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = finite(r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] /
                (r_ET_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] * rnbv[ind_f + nbF*ind_m] *
                rnbTrip[ind_f + nbF*ind_m] * rtripLgth[ind_f + nbF*ind_m]));

            ////Rprintf("Eco X1337\n");fichier << "EcoX1337" << endl;//calcul du numerateur de cnb_f

            r_cnb_f_out[ind_f + nbF*ind_t] = r_cnb_f_out[ind_f + nbF*ind_t] +
                finite(r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] * rnbv[ind_f + nbF*ind_m] *
                rnbTrip[ind_f + nbF*ind_m] * rtripLgth[ind_f + nbF*ind_m]);
        }
    }

    //if (ind_t==1) PrintValue(VECTOR_ELT(out_EcoDCF, 13));

    // --------------------------------------------------------------------------------------




    // A ce stade, plus de consideration d'espece pour les indicateurs

    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
        double NGVLtot_f = 0.0, RTBStot_f = 0.0;

        for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

            //-- 5. GVLav_f_m

            r_GVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] /
                r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]];

            //-- 6. GVLtot_f & NGVLav_f_m

            if (ind_m==0) {

                if (!ISNA(r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]])) {

                    r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                } else {

                    r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = 0.0;
                }

            } else {

                if (!ISNA(r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]])) {

                    r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] +
                        r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];
                }
            }


            if (!ISNA(r_NGVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]])) {

                NGVLtot_f = NGVLtot_f +
                    r_NGVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                    r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]];

            }

            r_NGVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                r_NGVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] /
                r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]];



            //-- 11. rtbs_f_m

            r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                r_NGVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] -
                ((finite(r_ovcDCFue_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]]) +
                finite(r_fvolue_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]]) *
                r_vf_f_m[ind_f*dim_vf_f_m[0] + ind_m*dim_vf_f_m[1] + 0*dim_vf_f_m[2] + ind_t*dim_vf_f_m[3]]) *
                r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] / pow(1+0.0,ind_t));

            if (!ISNA(r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]])) {

                RTBStot_f = RTBStot_f +
                    r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                    r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]];

            }

            if (persCalc<2) {

                r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    0.01 * r_cshr_f_m[ind_f*dim_cshr_f_m[0] + ind_m*dim_cshr_f_m[1] + 0*dim_cshr_f_m[2] + ind_t*dim_cshr_f_m[3]] *
                    r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

            } else if (persCalc==5){

                r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    0.01 * r_cshr_f_m[ind_f*dim_cshr_f_m[0] + ind_m*dim_cshr_f_m[1] + 0*dim_cshr_f_m[2] + ind_t*dim_cshr_f_m[3]] *
                    r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

            }else {

                r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] = NA_REAL;

            }

        } //on sort de la boucle sur les niveaux metiers


        //-- 7. GVLav_f

        r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
            r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];


        //-- 8.5. NGVLav_f

        r_NGVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            NGVLtot_f / r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];

        //-- 10. cnb_f

        r_cnb_f_out[ind_f + nbF*ind_t] =
            finite(r_cnb_f_out[ind_f + nbF*ind_t] / (rnbv_f[ind_f] * rtripLgth_f[ind_f] * rnbTrip_f[ind_f]));


        //-- 12. rtbs_f

        r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            RTBStot_f / r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];



        //version actualisee
        r_rtbsAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);


            //-- 14. cshrT_f


        if (persCalc==0) {  //salaires par marin fixes

            r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                0.01*r_cshr_f[ind_f*dim_cshr_f[0] + 0*dim_cshr_f[1] + 0*dim_cshr_f[2] + ind_t*dim_cshr_f[3]] *
                r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] *
                r_cnb_f_out[ind_f + nbF*ind_t] / r_cnb_f_out[ind_f + nbF*0];

        }

        if (persCalc==1) {  //part equipage constante (RAP)

            r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                0.01*r_cshr_f[ind_f*dim_cshr_f[0] + 0*dim_cshr_f[1] + 0*dim_cshr_f[2] + ind_t*dim_cshr_f[3]] *
                r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

        }

        if (persCalc==2) {  //part equipage constante calculee (RAP) - ccwr

            r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                0.01*r_ccwr_f2[ind_f] *
                r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

        }


        if (persCalc==3) {  //part equipage constante (RAP) + salaire marin supplementaire fixe

            r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                0.01*r_cshr_f[ind_f*dim_cshr_f[0] + 0*dim_cshr_f[1] + 0*dim_cshr_f[2] + ind_t*dim_cshr_f[3]] *
                (r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] +
                r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] *
                (r_cnb_f_out[ind_f + nbF*0 + nbF*ind_t] - r_cnb_f_out[ind_f + nbF*0 + nbF*0]) /
                r_cnb_f_out[ind_f + nbF*0 + nbF*0]);

        }


        if (persCalc==4) {  //part equipage constante calculee (RAP)- salaires marin supplementaire fixe

            r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_ccwr_f2[ind_f] *
                (r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] +
                r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] *
                (r_cnb_f_out[ind_f + nbF*0 + nbF*ind_t] - r_cnb_f_out[ind_f + nbF*0 + nbF*0]) /
                r_cnb_f_out[ind_f + nbF*0 + nbF*0]);

        }

        if (persCalc==5) {  //part equipage constante (GVL)

            r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                0.01*r_cshr_f[ind_f*dim_cshr_f[0] + 0*dim_cshr_f[1] + 0*dim_cshr_f[2] + ind_t*dim_cshr_f[3]] *
                r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] ;

        }
        //-- 15. ncshr_f

        r_ncshr_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
            r_eec_f[ind_f*dim_eec_f[0] + 0*dim_eec_f[1] + 0*dim_eec_f[2] + ind_t*dim_eec_f[3]];

        //-- 16. ocl_f

        r_ocl_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_mwh_f[ind_f*dim_mwh_f[0] + 0*dim_mwh_f[1] + 0*dim_mwh_f[2] + 0*dim_mwh_f[3]] *
            r_cnb_f_out[ind_f + nbF*0 + nbF*ind_t] * rtripLgth_f[ind_f] * rnbTrip_f[ind_f];


        //-- 17. cs_f

        r_cs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_ncshr_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
            r_ocl_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

        //version actualisee
        r_csAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_cs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);


        //-- 18. csTot_f

        r_csTot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_cs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] * rnbv_f[ind_f];


        //-- 19. gva_f

        r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
            (r_rep_f[ind_f*dim_rep_f[0] + 0*dim_rep_f[1] + 0*dim_rep_f[2] + ind_t*dim_rep_f[3]] +
            r_fixc_f[ind_f*dim_fixc_f[0] + 0*dim_fixc_f[1] + 0*dim_fixc_f[2] + ind_t*dim_fixc_f[3]] +
            r_gc_f[ind_f*dim_gc_f[0] + 0*dim_gc_f[1] + 0*dim_gc_f[2] + ind_t*dim_gc_f[3]])/ pow(1+0.0,ind_t) ;//+
            // r_QuotaExp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+0.0,ind_t) ;


        //version actualisee
        r_gvaAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
           r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);


        //-- 20. gvamargin_f

        r_gvamargin_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
           r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
            r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];


        //-- 21. gva_FTE_f

        r_gva_FTE_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
           r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
           r_FTE_f[ind_f*dim_FTE_f[0] + 0*dim_FTE_f[1] + 0*dim_FTE_f[2] + ind_t*dim_FTE_f[3]];


        //-- 22. ccw_f

        r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
           r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

        if ( (persCalc==0) | (persCalc==1) | (persCalc==3) ) {

            r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] + r_opersc_f2[ind_f];
        }


        //-- 23. ccwCr_f

        r_ccwCr_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / r_cnb_f_out[ind_f + nbF*ind_t];


        //-- 24. wageg_f

        r_wageg_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / r_cnb_f_out[ind_f + nbF*ind_t];


        //-- 25. wagen_f

        r_wagen_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_ncshr_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / r_cnb_f_out[ind_f + nbF*ind_t];

        //-- 26. wageg_FTE_f

        r_wageg_FTE_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / r_FTE_f[ind_f*dim_FTE_f[0] + 0*dim_FTE_f[1] + 0*dim_FTE_f[2] + ind_t*dim_FTE_f[3]];
           // r_wageg_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / r_FTE_f[ind_f*dim_FTE_f[0] + 0*dim_FTE_f[1] + 0*dim_FTE_f[2] + ind_t*dim_FTE_f[3]];


        //-- 27. wageg_h_f

        r_wageg_h_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_wageg_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / (rtripLgth_f[ind_f] * rnbTrip_f[ind_f]);
           // r_wageg_FTE_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / (rtripLgth_f[ind_f] * rnbTrip_f[ind_f]);

        //-- 28. gp_f

        r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] - r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];


        //version actualisee
        r_gpAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
           r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);


        //-- 29. gpmargin_f

        r_gpmargin_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

        //-- 30. ncf_f

        r_ncf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] - r_dep_f[ind_f*dim_dep_f[0] + 0*dim_dep_f[1] + 0*dim_dep_f[2] + ind_t*dim_dep_f[3]];


        //-- 31. np_f

        r_np_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_ncf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] - r_ic_f[ind_f*dim_ic_f[0] + 0*dim_ic_f[1] + 0*dim_ic_f[2] + ind_t*dim_ic_f[3]];

        //-- 32. npmargin_f

        r_npmargin_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_np_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];


        //-- 33. prof_f

        r_prof_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = -1.0;
        if (r_npmargin_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]]>=0) r_prof_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = 0.0;
        if (r_npmargin_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]]>0.1) r_prof_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = 1.0;


        //-- 34. npmargin_trend_f

        r_npmargin_trend_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = -1.0;
        if (ind_t>=5) {
            double devTrend;
            devTrend = r_npmargin_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                (0.2 * (r_npmargin_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + (ind_t-5)*eF_f[3]] +
                r_npmargin_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + (ind_t-4)*eF_f[3]] +
                r_npmargin_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + (ind_t-3)*eF_f[3]] +
                r_npmargin_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + (ind_t-2)*eF_f[3]] +
                r_npmargin_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + (ind_t-1)*eF_f[3]]));
            if (devTrend>(-0.05)) r_npmargin_trend_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = 0.0;
            if (devTrend>0.05) r_npmargin_trend_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = 1.0;
        }

        //-- 35. ssTot_f

        // r_ssTot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
        //     r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] * rnbv_f[ind_f];

        //-- 36. ps_f

        r_ps_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            rnbv_f[ind_f] * (r_cs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] + r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]]);

        //version actualisee
        r_psAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_ps_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);


        //-- 37. sts_f

        r_sts_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = 0.0;

        for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

            r_sts_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = r_sts_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] +
                (finite(r_lc_f_m[ind_f*dim_lc_f_m[0] + ind_m*dim_lc_f_m[1] + 0*dim_lc_f_m[2] + ind_t*dim_lc_f_m[3]]) *
                r_GVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] * rnbv[ind_f + ind_m*nbF]);

        }

        //version actualisee
        r_stsAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_sts_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);



        //-- 38. ber_f

        r_BER_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] *
            (r_fixc_f[ind_f*dim_fixc_f[0] + 0*dim_fixc_f[1] + 0*dim_fixc_f[2] + ind_t*dim_fixc_f[3]] +
            r_dep_f[ind_f*dim_dep_f[0] + 0*dim_dep_f[1] + 0*dim_dep_f[2] + ind_t*dim_dep_f[3]] +
            r_ic_f[ind_f*dim_ic_f[0] + 0*dim_ic_f[1] + 0*dim_ic_f[2] + ind_t*dim_ic_f[3]]) /
            r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] ;

        //-- 39. CR_BER_f

        r_CR_BER_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
            r_BER_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];


        //-- 40. fuelEff_f
        double numFuelEff = 0.0;
        for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){
            numFuelEff = numFuelEff + finite(r_fvolue_f_m2[ind_f + ind_m*nbF] * r_ue_f_m[ind_f+ ind_m*nbF] * r_nbv_f_m[ind_f+ ind_m*nbF]);
        }
        r_fuelEff_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = numFuelEff / r_countLf[ind_f];


        //-- 41. ratio_fvol_GVA_f
        double numFvolGVA = 0.0;
        for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){
            numFvolGVA = numFvolGVA + finite(r_fvolue_f_m2[ind_f + ind_m*nbF] * r_ue_f_m[ind_f+ ind_m*nbF]);
        }
        r_ratio_fvol_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = numFvolGVA / r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

        //-- 42. ratio_gp_GVA_f

        r_ratio_gp_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
            r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
            r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

        //-- 43. ratio_GVL_K_f

        // r_ratio_GVL_K_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
        //     r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
        //     (r_K_f[ind_f*dim_K_f[0] + 0*dim_K_f[1] + 0*dim_K_f[2] + ind_t*dim_K_f[3]] / pow(1+0.0,ind_t) );

        //-- 44. ratio_gp_K_f

        // r_ratio_gp_K_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
        //     r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
        //     (r_K_f[ind_f*dim_K_f[0] + 0*dim_K_f[1] + 0*dim_K_f[2] + ind_t*dim_K_f[3]] / pow(1+0.0,ind_t) );

        //-- 45. RoFTA_f

        // r_RoFTA_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
        //     r_ncf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
        //     (r_K_f[ind_f*dim_K_f[0] + 0*dim_K_f[1] + 0*dim_K_f[2] + ind_t*dim_K_f[3]] / pow(1+0.0,ind_t) );

        //-- 46. ROI_f
        // r_ROI_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = NA_REAL;
        // if (finite(r_inv_f[ind_f*dim_inv_f[0] + 0*dim_inv_f[1] + 0*dim_inv_f[2] + ind_t*dim_inv_f[3]])>0) {
        //     r_ROI_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
        //         (r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] - r_inv_f[ind_f*dim_inv_f[0] + 0*dim_inv_f[1] + 0*dim_inv_f[2] + ind_t*dim_inv_f[3]]) /
        //         finite(r_inv_f[ind_f*dim_inv_f[0] + 0*dim_inv_f[1] + 0*dim_inv_f[2] + ind_t*dim_inv_f[3]]);
        // }

        //-- 47. ratio_np_K_f
        // r_ratio_np_K_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
        // // r_np_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] * 0;
        //     r_np_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
        //     r_K_f[ind_f*dim_K_f[0] + 0*dim_K_f[1] + 0*dim_K_f[2] + ind_t*dim_K_f[3]];



        //-- 48. ratio_GVL_cnb_ue_f

        // r_ratio_GVL_cnb_ue_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
        //     r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
        //     (r_cnb_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] * r_ue_f[ind_f] );


    }

    SEXP out_stripe;
    Rf_protect(out_stripe = _BioEcoPar_stripe_ecoCDF(out_EcoDCF, Flist, ind_t));
    SEXP ssTot = getListElement(out_stripe, "ssTot");                          SET_VECTOR_ELT(out_EcoDCF, 43, ssTot);
    SEXP ratio_GVL_K = getListElement(out_stripe, "ratio_GVL_K");              SET_VECTOR_ELT(out_EcoDCF, 53, ratio_GVL_K);
    SEXP ratio_gp_K = getListElement(out_stripe, "ratio_gp_K");                SET_VECTOR_ELT(out_EcoDCF, 54, ratio_gp_K);
    SEXP RoFTA = getListElement(out_stripe, "RoFTA");                          SET_VECTOR_ELT(out_EcoDCF, 55, RoFTA);
    SEXP ROI = getListElement(out_stripe, "ROI");                              SET_VECTOR_ELT(out_EcoDCF, 56, ROI);
    SEXP ratio_np_K = getListElement(out_stripe, "ratio_np_K");                SET_VECTOR_ELT(out_EcoDCF, 57, ratio_np_K);
    SEXP ratio_GVL_cnb_ue = getListElement(out_stripe, "ratio_GVL_cnb_ue");    SET_VECTOR_ELT(out_EcoDCF, 58, ratio_GVL_cnb_ue);
    
    // static SEXP att_dimCstF = att_dimCstF;
    Rf_unprotect(1);


    if (ind_t==0) Rf_unprotect(67-46);
    Rf_unprotect(85-6);

    //Rprintf("\nJ2\n");fichier << "J2" << endl;

    //if (ind_t>10) //PrintValue(VECTOR_ELT(out_EcoDCF, 44));
    //if (ind_t>10) //PrintValue(VECTOR_ELT(out_EcoDCF, 2));
    //if (ind_t>10) //PrintValue(VECTOR_ELT(fVar, 10));
    //if (ind_t>10) //PrintValue(VECTOR_ELT(fVar, 4));
    //if (ind_t>10) //PrintValue(vf_f_m);
    //if (ind_t>10) //PrintValue(ue_f_m);
    //if (ind_t==(nbT-1)) {//PrintValue(VECTOR_ELT(fVar,4));//PrintValue(VECTOR_ELT(fVar,10));}

    //fichier.close();

}}

