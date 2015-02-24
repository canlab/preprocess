%function f99_make_desMtx


load KUL_stat_gen
load KUL_stat_stoch
load KUL_stat_nonstoch
load KUL_stat_modul
load KUL_stat_basisfunct

f99_spm_fmri_design_ma(RT,nscan,kulRES,STOC,rep,n_cond,cond_name,nulev,soa_stoc,oc_prob,prob_type,Sstr,soa_n_stoc,onset,Ptype,Pstr,Etype,h_exp,h_pol,sel,p_mod,Rov,Cov,pst_ev,h,m,HRF,TD,W,VOLT,regr_us,regr_us_na,durat,hrffunct)