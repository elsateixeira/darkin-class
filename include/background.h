/** @file background.h Documented includes for background module */

#ifndef __BACKGROUND__
#define __BACKGROUND__

#include "common.h"
#include "quadrature.h"
#include "growTable.h"
#include "arrays.h"
#include "dei_rkck.h"
#include "parser.h"

/** list of possible parametrisations of the DE equation of state */

enum equation_of_state {CLP,EDE};

/** ET: scalar field potential types */
enum scf_potential_type {scf_potential_exp, scf_potential_double_exp};

/** ET: scalar field conformal coupling function types */
enum scf_conformal_form_type {scf_conformal_exp, scf_conformal_matter_exp};

/** ET: scalar field coupling types */
enum scf_coupling_type {
  scf_coupling_none,
  scf_coupling_conformal,
  scf_coupling_disformal,
  scf_coupling_mixed,
  scf_coupling_entropy,
  scf_coupling_momentum
};

/** ET: scalar field shooting target (for explicit scf_* inputs) */
enum scf_shooting_target_type {
  scf_shoot_none,
  scf_shoot_V0,
  scf_shoot_lambda,
  scf_shoot_V0_2,
  scf_shoot_lambda_2,
  scf_shoot_C0,
  scf_shoot_beta,
  scf_shoot_alpha,
  scf_shoot_D0,
  scf_shoot_phi_ini,
  scf_shoot_phi_prime_ini
};


/** list of possible parametrizations of the varying fundamental constants */

enum varconst_dependence {varconst_none,varconst_instant};

/** list of formats for the vector of background quantities */

enum vecback_format {short_info, normal_info, long_info};

/** list of interpolation methods: search location in table either
    by bisection (inter_normal), or step by step starting from given
    index (inter_closeby) */

enum interpolation_method {inter_normal, inter_closeby};

/**
 * background structure containing all the background information that
 * other modules need to know.
 *
 * Once initialized by the backgound_init(), contains all necessary
 * information on the background evolution (except thermodynamics),
 * and in particular, a table of all background quantities as a
 * function of time and scale factor, used for interpolation in other
 * modules.
 */

struct background
{
  /** @name - input parameters initialized by user in input module
   *  (all other quantities are computed in this module, given these parameters
   *   and the content of the 'precision' structure)
   *
   * The background cosmological parameters listed here form a parameter
   * basis which is directly usable by the background module. Nothing
   * prevents from defining the input cosmological parameters
   * differently, and to pre-process them into this format, using the input
   * module (this might require iterative calls of background_init()
   * e.g. for dark energy or decaying dark matter). */

  //@{

  double H0; /**< \f$ H_0 \f$: Hubble parameter (in fact, [\f$H_0/c\f$]) in \f$ Mpc^{-1} \f$ */
  double h;  /**< reduced Hubble parameter */

  double Omega0_g; /**< \f$ \Omega_{0 \gamma} \f$: photons */
  double T_cmb;    /**< \f$ T_{cmb} \f$: current CMB temperature in Kelvins */

  double Omega0_b; /**< \f$ \Omega_{0 b} \f$: baryons */

  double Omega0_ur; /**< \f$ \Omega_{0 \nu r} \f$: ultra-relativistic neutrinos */

  double Omega0_cdm;      /**< \f$ \Omega_{0 cdm} \f$: cold dark matter */

  double Omega0_idm; /**< \f$ \Omega_{0 idm} \f$: interacting dark matter with photons, baryons, and idr */

  double Omega0_qcdm; /**< \f$ \Omega_{0 qcdm} \f$: dark matter coupled to the scalar field */

  double Omega_ini_qcdm; /**< \f$ ET: \Omega_{ini,qcdm} \f$: coupled cold dark matter initial value */

  double omega_ini_qcdm; /**< \f$ ET: \omega_{ini,qcdm} \f$: coupled cold dark matter initial value */

  double Omega_ini_idm; /**< \f$ ET: \Omega_{ini,idm} \f$: coupled cold dark matter initial value */

  double omega_ini_idm; /**< \f$ ET: \Omega_{ini,idm} \f$: coupled cold dark matter initial value */

  double Omega0_idr; /**< \f$ \Omega_{0 idr} \f$: interacting dark radiation */
  double T_idr;      /**< \f$ T_{idr} \f$: current temperature of interacting dark radiation in Kelvins */

  double Omega0_dcdmdr;   /**< \f$ \Omega_{0 dcdm}+\Omega_{0 dr} \f$: decaying cold dark matter (dcdm) decaying to dark radiation (dr) */
  double Omega_ini_dcdm;  /**< \f$ \Omega_{ini,dcdm} \f$: rescaled initial value for dcdm density (see 1407.2418 for definitions) */
  double Gamma_dcdm;      /**< \f$ \Gamma_{dcdm} \f$: decay constant for decaying cold dark matter */
  double tau_dcdm;

  int N_ncdm;                            /**< Number of distinguishable ncdm species */
  /* the following parameters help to define tabulated ncdm p-s-d passed in file */
  char * ncdm_psd_files;                 /**< list of filenames for tabulated p-s-d */
  int * got_files;                       /**< list of flags for each species, set to true if p-s-d is passed through file */
  /* the following parameters help to define the analytical ncdm phase space distributions (p-s-d) */
  double * ncdm_psd_parameters;          /**< list of parameters for specifying/modifying ncdm p.s.d.'s, to be customized for given model
                                            (could be e.g. mixing angles) */
  double * M_ncdm;                       /**< vector of masses of non-cold relic: dimensionless ratios m_ncdm/T_ncdm */
  double * m_ncdm_in_eV;                 /**< list of ncdm masses in eV (inferred from M_ncdm and other parameters above) */
  double * Omega0_ncdm, Omega0_ncdm_tot; /**< Omega0_ncdm for each species and for the total Omega0_ncdm */
  double * T_ncdm,T_ncdm_default;        /**< list of 1st parameters in p-s-d of non-cold relics: relative temperature
                                            T_ncdm1/T_gamma; and its default value */
  double * ksi_ncdm, ksi_ncdm_default;   /**< list of 2nd parameters in p-s-d of non-cold relics: relative chemical potential
                                            ksi_ncdm1/T_ncdm1; and its default value */
  double * deg_ncdm, deg_ncdm_default;    /**< vector of degeneracy parameters in factor of p-s-d: 1 for one family of neutrinos
                                             (= one neutrino plus its anti-neutrino, total g*=1+1=2, so deg = 0.5 g*); and its
                                             default value */
  int * ncdm_input_q_size; /**< Vector of numbers of q bins */
  double * ncdm_qmax;      /**< Vector of maximum value of q */

  double Omega0_k;         /**< \f$ \Omega_{0_k} \f$: curvature contribution */

  double Omega0_lambda;    /**< \f$ \Omega_{0_\Lambda} \f$: cosmological constant */
  double Omega0_fld;       /**< \f$ \Omega_{0 de} \f$: fluid */
  double Omega0_scf;       /**< \f$ \Omega_{0 scf} \f$: scalar field */
  short use_ppf; /**< flag switching on PPF perturbation equations instead of true fluid equations for perturbations. It could have been defined inside
                    perturbation structure, but we leave it here in such way to have all fld parameters grouped. */
  double c_gamma_over_c_fld; /**< ppf parameter defined in eq. (16) of 0808.3125 [astro-ph] */
  enum equation_of_state fluid_equation_of_state; /**< parametrisation scheme for fluid equation of state */
  double w0_fld;   /**< \f$ w0_{DE} \f$: current fluid equation of state parameter */
  double wa_fld;   /**< \f$ wa_{DE} \f$: fluid equation of state parameter derivative */
  double cs2_fld;  /**< \f$ c^2_{s~DE} \f$: sound speed of the fluid in the frame comoving with the fluid (so, this is
                      not [delta p/delta rho] in the synchronous or newtonian gauge!) */
  double Omega_EDE;        /**< \f$ wa_{DE} \f$: Early Dark Energy density parameter */
  double * scf_parameters; /**< list of parameters describing the scalar field potential */
  // ET: Added extra parameters for the scalar field coupling and potential functions. The exact meaning of these parameters depends on the specific model considered, and they are defined in background.c.
  int scf_parameters_size; /**< size of scf_parameters */
  short use_scf_parameters; /**< ET: whether to use scf_parameters rather than explicit scf_* inputs */
  enum scf_potential_type scf_potential; /**< ET: scalar field potential type */
  enum scf_conformal_form_type scf_conformal_form; /**< ET: conformal coupling function type */
  enum scf_coupling_type scf_coupling; /**< ET: scalar field coupling type */
  enum scf_shooting_target_type scf_shooting_target; /**< ET: shooting target when using explicit scf_* inputs */
  short scf_ic_from_today; /**< ET: reconstruct scalar/qcdm initial conditions from today boundary values */
  short scf_allow_multiple_couplings; /**< ET: allow more than one scalar-field coupling sector */
  short scf_use_conformal; /**< ET: requested conformal scalar-field coupling */
  short scf_use_disformal; /**< ET: requested disformal scalar-field coupling */
  short scf_use_entropy; /**< ET: requested entropy scalar-field coupling */
  short scf_use_momentum; /**< ET: requested momentum scalar-field coupling */
  short has_idm_de; /**< ET: legacy IDM/SCF coupling flag; kept off when using qcdm */
  short has_idm_de_q; /**< ET: legacy Q_scf-sector coupling flag; kept off when using qcdm */
  short has_qcdm_de; /**< ET: coupling between qcdm and scalar field */
  short has_qcdm_de_q; /**< ET: qcdm Q_scf-sector coupling active (conformal/disformal/mixed) */
  short has_scf_conformal; /**< ET: conformal coupling active */
  short has_scf_disformal; /**< ET: disformal coupling active */
  short has_scf_entropy; /**< ET: entropy coupling active */
  short has_scf_momentum; /**< ET: momentum coupling active */
  double V0_scf;        /**< ET: scalar field potential amplitude */
  double lambda_scf;    /**< ET: scalar field potential slope */
  double V0_scf_2;      /**< ET: second exponential amplitude */
  double lambda_scf_2;  /**< ET: second exponential slope */
  double beta_scf;      /**< ET: conformal coupling strength */
  double C0_scf;        /**< ET: conformal coupling amplitude */
  double alpha_scf;     /**< ET: disformal coupling strength */
  double D0_scf;        /**< ET: disformal coupling scale (in meV^-1) */
  double phi_today_scf; /**< ET: scalar field value today for reconstructed ICs */
  double w_today_scf; /**< ET: scalar field equation of state today for reconstructed ICs */
  double phi_prime_today_sign_scf; /**< ET: sign of phi' today for reconstructed ICs */
  double scf_gamma0;      /**< ET: momentum-coupling strength (Type-3 gamma) */
  /* ET: entropy-coupling/source parameters (coupled through QCDM) */
  double g0_scf;        /**< entropy force coefficient entering g_scf(phi) */
  double h0_scf;        /**< entropy mixing coefficient entering h_scf(phi) */
  double As_scf;        /**< amplitude of entropy source mode delta_s */
  double ns_scf;        /**< tilt of entropy source mode delta_s */
  double kp_scf;        /**< pivot scale for entropy source mode delta_s */
  double kc_scf;        /**< cutoff scale for entropy source mode delta_s */
  double pc_scf;        /**< cutoff power for entropy source mode delta_s */
  short attractor_ic_scf;  /**< whether the scalar field has attractor initial conditions */
  int scf_tuning_index;    /**< index in scf_parameters used for tuning */
  double phi_ini_scf;      /**< \f$ \phi(t_0) \f$: scalar field initial value */
  double phi_prime_ini_scf;/**< \f$ d\phi(t_0)/d\tau \f$: scalar field initial derivative wrt conformal time */
  double varconst_alpha; /**< finestructure constant for varying fundamental constants */
  double varconst_me; /**< electron mass for varying fundamental constants */
  enum varconst_dependence varconst_dep; /**< dependence of the varying fundamental constants as a function of time */
  double varconst_transition_redshift; /**< redshift of transition between varied fundamental constants and normal fundamental constants in the 'varconst_instant' case*/

  //@}


  /** @name - related parameters */

  //@{

  double age; /**< age in Gyears */
  double conformal_age; /**< conformal age in Mpc */
  double K; /**< \f$ K \f$: Curvature parameter \f$ K=-\Omega0_k*a_{today}^2*H_0^2\f$; */
  int sgnK; /**< K/|K|: -1, 0 or 1 */
  double Neff; /**< so-called "effective neutrino number", computed at earliest time in interpolation table */
  double Omega0_dcdm; /**< \f$ \Omega_{0 dcdm} \f$: decaying cold dark matter */
  double Omega0_dr; /**< \f$ \Omega_{0 dr} \f$: decay radiation */
  double Omega0_m;  /**< total non-relativistic matter today */
  double Omega0_r;  /**< total ultra-relativistic radiation today */
  double Omega0_de; /**< total dark energy density today, currently defined as 1 - Omega0_m - Omega0_r - Omega0_k */
  double Omega0_nfsm; /**< total non-free-streaming matter, that is, cdm, baryons and wdm */
  double a_eq;      /**< scale factor at radiation/matter equality */
  double H_eq;      /**< Hubble rate at radiation/matter equality [Mpc^-1] */
  double z_eq;      /**< redshift at radiation/matter equality */
  double tau_eq;    /**< conformal time at radiation/matter equality [Mpc] */

  //@}


  /** @name - all indices for the vector of background (=bg) quantities stored in table */

  //@{

  int index_bg_a;             /**< scale factor (in fact (a/a_0), see
                                 normalisation conventions explained
                                 at beginning of background.c) */
  int index_bg_H;             /**< Hubble parameter in \f$Mpc^{-1}\f$ */
  int index_bg_H_prime;       /**< its derivative w.r.t. conformal time */

  /* end of vector in short format, now quantities in normal format */

  int index_bg_rho_g;         /**< photon density */
  int index_bg_rho_b;         /**< baryon density */
  int index_bg_rho_cdm;       /**< cdm density */
  int index_bg_rho_idm;       /**< idm density */
  int index_bg_rho_qcdm;      /**< qcdm density */
  int index_bg_rho_lambda;    /**< cosmological constant density */
  int index_bg_rho_fld;       /**< fluid density */
  int index_bg_w_fld;         /**< fluid equation of state */
  int index_bg_rho_idr;       /**< density of interacting dark radiation */
  int index_bg_rho_ur;        /**< relativistic neutrinos/relics density */
  int index_bg_rho_dcdm;      /**< dcdm density */
  int index_bg_rho_dr;        /**< dr density */

  int index_bg_phi_scf;       /**< scalar field value */
  int index_bg_phi_prime_scf; /**< scalar field derivative wrt conformal time */
  int index_bg_V_scf;         /**< scalar field potential V */
  int index_bg_dV_scf;        /**< scalar field potential derivative V' */
  int index_bg_ddV_scf;       /**< scalar field potential second derivative V'' */
  /* ET: momentum indices */
  int index_bg_mom_scf;       /**< momentum variable \f$ Z=-\phi'/a \f$ used in momentum coupling */
  int index_bg_gamma_scf;      /**< momentum-coupling function gamma_scf(Z) */
  int index_bg_dgamma_scf;     /**< first derivative dgamma_scf/dZ */
  int index_bg_ddgamma_scf;    /**< second derivative ddgamma_scf/dZ2 */
  /* ET: entropy indices */
  int index_bg_g_scf;         /**< entropy function g_scf(phi) */
  int index_bg_dg_scf;        /**< first derivative of g_scf(phi) */
  int index_bg_ddg_scf;       /**< second derivative of g_scf(phi) */
  int index_bg_h_scf;         /**< entropy function h_scf(phi) */
  int index_bg_dh_scf;        /**< first derivative of h_scf(phi) */
  int index_bg_ddh_scf;       /**< second derivative of h_scf(phi) */
  int index_bg_As_scf;        /**< entropy source amplitude */
  int index_bg_ns_scf;        /**< entropy source tilt */
  int index_bg_kp_scf;        /**< entropy source pivot scale */
  int index_bg_kc_scf;        /**< entropy source cutoff scale */
  int index_bg_pc_scf;        /**< entropy source cutoff power */
  int index_bg_rho_scf;       /**< scalar field energy density */
  int index_bg_p_scf;         /**< scalar field pressure */
  int index_bg_p_prime_scf;         /**< scalar field pressure */

  /* ET: conformal and disformal indices */
  int index_bg_C_scf;         /**< scalar field conformal factor C */
  int index_bg_dC_scf;        /**< scalar field conformal factor derivative C' */
  int index_bg_ddC_scf;       /**< scalar field conformal factor second derivative C'' */
  int index_bg_D_scf;         /**< scalar field disformal factor D */
  int index_bg_dD_scf;        /**< scalar field disformal factor derivative D' */
  int index_bg_ddD_scf;       /**< scalar field disformal factor second derivative D'' */
  int index_bg_Q_scf;         /**< scalar field coupling function */
  int index_bg_B_cff_scf;     /**< scf scalar field perturbed coupling denominator */
  int index_bg_B1_scf;        /**< scf scalar field perturbed coupling delta_c factor */
  int index_bg_B2_scf;        /**< scf scalar field perturbed coupling Phi' factor */
  int index_bg_B3_scf;        /**< scf scalar field perturbed coupling Psi factor */
  int index_bg_B4_scf;        /**< scf scalar field perturbed coupling delta phi' factor */
  int index_bg_B5_scf;        /**< scf scalar field perturbed coupling delta phi factor*/

  int index_bg_rho_ncdm1;     /**< density of first ncdm species (others contiguous) */
  int index_bg_p_ncdm1;       /**< pressure of first ncdm species (others contiguous) */
  int index_bg_pseudo_p_ncdm1;/**< another statistical momentum useful in ncdma approximation */

  int index_bg_rho_tot;       /**< Total density */
  int index_bg_p_tot;         /**< Total pressure */
  int index_bg_p_tot_prime;   /**< Conf. time derivative of total pressure */

  int index_bg_Omega_r;       /**< relativistic density fraction (\f$ \Omega_{\gamma} + \Omega_{\nu r} \f$) */

  /* end of vector in normal format, now quantities in long format */

  int index_bg_rho_crit;      /**< critical density */
  int index_bg_Omega_m;       /**< non-relativistic density fraction (\f$ \Omega_b + \Omega_cdm + \Omega_{\nu nr} \f$) */
  int index_bg_conf_distance; /**< conformal distance (from us) in Mpc */
  int index_bg_ang_distance;  /**< angular diameter distance in Mpc */
  int index_bg_lum_distance;  /**< luminosity distance in Mpc */
  int index_bg_time;          /**< proper (cosmological) time in Mpc */
  int index_bg_rs;            /**< comoving sound horizon in Mpc */

  int index_bg_D;             /**< scale independent growth factor D(a) for CDM perturbations */
  int index_bg_f;             /**< corresponding velocity growth factor [dlnD]/[dln a] */

  int index_bg_varc_alpha;    /**< value of fine structure constant in varying fundamental constants */
  int index_bg_varc_me;      /**< value of effective electron mass in varying fundamental constants */

  int bg_size_short;  /**< size of background vector in the "short format" */
  int bg_size_normal; /**< size of background vector in the "normal format" */
  int bg_size;        /**< size of background vector in the "long format" */

  //@}


  /** @name - background interpolation tables */

  //@{

  int bt_size;               /**< number of lines (i.e. time-steps) in the four following array */
  double * loga_table;       /**< vector loga_table[index_loga] with values of log(a) (in fact \f$ log(a/a0) \f$, logarithm of relative scale factor compared to today) */
  double * tau_table;        /**< vector tau_table[index_loga] with values of conformal time \f$ \tau \f$ (in fact \f$ a_0 c tau \f$, see normalisation conventions explained at beginning of background.c) */
  double * z_table;          /**< vector z_table[index_loga] with values of \f$ z \f$ (redshift) */
  double * background_table; /**< table background_table[index_tau*pba->bg_size+pba->index_bg] with all other quantities (array of size bg_size*bt_size) **/

  //@}


  /** @name - table of their second derivatives, used for spline interpolation */

  //@{

  double * d2tau_dz2_table; /**< vector d2tau_dz2_table[index_loga] with values of \f$ d^2 \tau / dz^2 \f$ (conformal time) */
  double * d2z_dtau2_table; /**< vector d2z_dtau2_table[index_loga] with values of \f$ d^2 z / d\tau^2 \f$ (conformal time) */
  double * d2background_dloga2_table; /**< table d2background_dtau2_table[index_loga*pba->bg_size+pba->index_bg] with values of \f$ d^2 b_i / d\log(a)^2 \f$ */

  //@}


  /** @name - all indices for the vector of background quantities to be integrated (=bi)
   *
   * Most background quantities can be immediately inferred from the
   * scale factor. Only few of them require an integration with
   * respect to conformal time (in the minimal case, only one quantity needs to
   * be integrated with time: the scale factor, using the Friedmann
   * equation). These indices refer to the vector of
   * quantities to be integrated with time.
   * {B} quantities are needed by background_functions() while {C} quantities are not.
   */

  //@{

  int index_bi_rho_dcdm;/**< {B} dcdm density */
  int index_bi_rho_idm;  /**< ET: {B} extra index for idm density evolution */
  int index_bi_rho_qcdm; /**< ET: {B} extra index for qcdm density evolution */
  int index_bi_rho_dr;  /**< {B} dr density */
  int index_bi_rho_fld; /**< {B} fluid density */
  int index_bi_phi_scf;       /**< {B} scalar field value */
  int index_bi_phi_prime_scf; /**< {B} scalar field derivative wrt conformal time */

  int index_bi_time;    /**< {C} proper (cosmological) time in Mpc */
  int index_bi_rs;      /**< {C} sound horizon */
  int index_bi_tau;     /**< {C} conformal time in Mpc */
  int index_bi_D;       /**< {C} scale independent growth factor D(a) for CDM perturbations. */
  int index_bi_D_prime; /**< {C} D satisfies \f$ [D''(\tau)=-aHD'(\tau)+3/2 a^2 \rho_M D(\tau) \f$ */

  int bi_B_size;        /**< Number of {B} parameters */
  int bi_size;          /**< Number of {B}+{C} parameters */

  //@}

  /** @name - flags describing the absence or presence of cosmological
      ingredients
      *
      * having one of these flag set to zero allows to skip the
      * corresponding contributions, instead of adding null contributions.
      */


  //@{

  short has_cdm;       /**< presence of cold dark matter? */
  short has_idm;       /**< ET: presence of interacting dark matter with photons, baryons, and idr */
  short has_qcdm;      /**< ET: presence of dark matter interacting with the scalar field */
  short has_dcdm;      /**< presence of decaying cold dark matter? */
  short has_dr;        /**< presence of relativistic decay radiation? */
  short has_scf;       /**< presence of a scalar field? */
  short has_ncdm;      /**< presence of non-cold dark matter? */
  short has_lambda;    /**< presence of cosmological constant? */
  short has_fld;       /**< presence of fluid with constant w and cs2? */
  short has_ur;        /**< presence of ultra-relativistic neutrinos/relics? */
  short has_idr;       /**< presence of interacting dark radiation? */
  short has_curvature; /**< presence of global spatial curvature? */
  short has_varconst;  /**< presence of varying fundamental constants? */

  //@}


  /**
   *@name - arrays related to sampling and integration of ncdm phase space distributions
   */

  //@{

  int * ncdm_quadrature_strategy; /**< Vector of integers according to quadrature strategy. */
  double ** q_ncdm_bg;  /**< Pointers to vectors of background sampling in q */
  double ** w_ncdm_bg;  /**< Pointers to vectors of corresponding quadrature weights w */
  double ** q_ncdm;     /**< Pointers to vectors of perturbation sampling in q */
  double ** w_ncdm;     /**< Pointers to vectors of corresponding quadrature weights w */
  double ** dlnf0_dlnq_ncdm; /**< Pointers to vectors of logarithmic derivatives of p-s-d */
  int * q_size_ncdm_bg; /**< Size of the q_ncdm_bg arrays */
  int * q_size_ncdm;    /**< Size of the q_ncdm arrays */
  double * factor_ncdm; /**< List of normalization factors for calculating energy density etc.*/

  //@}

  /** @name - technical parameters */

  //@{

  short shooting_failed;  /**< flag is set to true if shooting failed. */
  ErrorMsg shooting_error; /**< Error message from shooting failed. */

  short background_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  short is_allocated; /**< flag is set to true if allocated */
  //@}
};


/**
 * temporary parameters and workspace passed to the background_derivs function
 */

struct background_parameters_and_workspace {

  /* structures containing fixed input parameters (indices, ...) */
  struct background * pba;

  /* workspace */
  double * pvecback;

};

/**
 * temporary parameters and workspace passed to phase space distribution function
 */

struct background_parameters_for_distributions {

  /* structures containing fixed input parameters (indices, ...) */
  struct background * pba;

  /* Additional parameters */

  /* Index of current distribution function */
  int n_ncdm;

  /* Used for interpolating in file of tabulated p-s-d: */
  int tablesize;
  double *q;
  double *f0;
  double *d2f0;
  int last_index;

};

/**************************************************************/
/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int background_at_z(
                      struct background *pba,
                      double a_rel,
                      enum vecback_format return_format,
                      enum interpolation_method inter_mode,
                      int * last_index,
                      double * pvecback
                      );

  int background_at_tau(
                        struct background *pba,
                        double tau,
                        enum vecback_format return_format,
                        enum interpolation_method inter_mode,
                        int * last_index,
                        double * pvecback
                        );

  int background_tau_of_z(
                          struct background *pba,
                          double z,
                          double * tau
                          );

  int background_z_of_tau(
                          struct background *pba,
                          double tau,
                          double * z
                          );

  int background_functions(
                           struct background *pba,
                           double a_rel,
                           double * pvecback_B,
                           enum vecback_format return_format,
                           double * pvecback
                           );

  int background_w_fld(
                       struct background * pba,
                       double a,
                       double * w_fld,
                       double * dw_over_da_fld,
                       double * integral_fld);

  int background_varconst_of_z(
                               struct background* pba,
                               double z,
                               double* alpha,
                               double* me
                               );

  int background_init(
                      struct precision *ppr,
                      struct background *pba
                      );

  int background_free(
                      struct background *pba
                      );

  int background_free_noinput(
                              struct background *pba
                              );

  int background_free_input(
                            struct background *pba
                            );

  int background_indices(
                         struct background *pba
                         );

  int background_ncdm_distribution(
                                   void *pba,
                                   double q,
                                   double * f0
                                   );

  int background_ncdm_test_function(
                                    void *pba,
                                    double q,
                                    double * test
                                    );

  int background_ncdm_init(
                           struct precision *ppr,
                           struct background *pba
                           );

  int background_ncdm_momenta(
                              double * qvec,
                              double * wvec,
                              int qsize,
                              double M,
                              double factor,
                              double z,
                              double * n,
                              double * rho,
                              double * p,
                              double * drho_dM,
                              double * pseudo_p
                              );

  int background_ncdm_M_from_Omega(
                                   struct precision *ppr,
                                   struct background *pba,
                                   int species
                                   );

  int background_checks(
                        struct precision * ppr,
                        struct background *pba
                        );

  int background_solve(
                       struct precision *ppr,
                       struct background *pba
                       );

  int background_initial_conditions(
                                    struct precision *ppr,
                                    struct background *pba,
                                    double * pvecback,
                                    double * pvecback_integration,
                                    double * loga_ini
                                    );

  int background_find_equality(
                               struct precision *ppr,
                               struct background *pba
                               );


  int background_output_titles(struct background * pba,
                               char titles[_MAXTITLESTRINGLENGTH_]
                               );

  int background_output_data(
                             struct background *pba,
                             int number_of_titles,
                             double *data);

  int background_derivs(
                        double loga,
                        double * y,
                        double * dy,
                        void * parameters_and_workspace,
                        ErrorMsg error_message
                        );

  int background_sources(
                         double loga,
                         double * y,
                         double * dy,
                         int index_loga,
                         void * parameters_and_workspace,
                         ErrorMsg error_message
                         );

  int background_timescale(
                           double loga,
                           void * parameters_and_workspace,
                           double * timescale,
                           ErrorMsg error_message
                           );

  int background_output_budget(
                               struct background* pba
                               );

  /* ET: Scalar field potential and its derivatives */
  double V_scf(
               struct background *pba,
               double phi
               );

  double dV_scf(
		struct background *pba,
		double phi
		);

  double ddV_scf(
                  struct background *pba,
                  double phi
                );

  /* ET: momentum coupling helper function gamma_scf(scf_mom) and derivatives */
  double gamma_scf(
               struct background *pba,
               double scf_mom
               );

  double dgamma_scf(
                struct background *pba,
                double scf_mom
                );

  double ddgamma_scf(
                 struct background *pba,
                 double scf_mom
                 );

  double dddgamma_scf(
                   struct background *pba,
                   double scf_mom
                   );

  /* ET: entropy coupling/source helper functions */
  double g_scf(
               struct background *pba,
               double phi
               );

  double dg_scf(
                struct background *pba,
                double phi
                );

  double ddg_scf(
                 struct background *pba,
                 double phi
                 );

  double h_scf(
               struct background *pba,
               double phi
               );

  double dh_scf(
                struct background *pba,
                double phi
                );

  double ddh_scf(
                 struct background *pba,
                 double phi
                 );

  double As_scf(
                struct background *pba
                );

  double ns_scf(
                struct background *pba
                );

  double kp_scf(
                struct background *pba
                );

  double kc_scf(
                struct background *pba
                );

  double pc_scf(
                struct background *pba
                );

  /* ET: conformal factor and its derivatives */

  double C_scf(
                struct background *pba,
                double phi
                );

  double dC_scf(
      struct background *pba,
      double phi
    );

  double ddC_scf(
                  struct background *pba,
                  double phi
                  );
  
  /* ET: disformal factor and its derivatives */

  double D_scf(
                struct background *pba,
                double phi
                );

  double dD_scf(
      struct background *pba,
      double phi
    );

  double ddD_scf(
                  struct background *pba,
                  double phi
                  );

  /* ET: Coupling funtion and perturbed coupling coefficients */

  double Q_scf(
              struct background *pba,
              double phi,
              double phi_prime,
              double rho_idm,
              double a,
              double *pvecback);

  double B_cff_scf(
                struct background *pba,
                double phi,
                double phi_prime,
                double rho_idm,
                double a,
                double * pvecback);

  double B1_scf(
                struct background *pba,
                double phi,
                double phi_prime,
                double rho_idm,
                double a,
                double * pvecback);

  double B2_scf(
                struct background *pba,
                double phi,
                double phi_prime,
                double rho_idm,
                double a,
                double * pvecback);
  
  double B3_scf(
                struct background *pba,
                double phi,
                double phi_prime,
                double rho_idm,
                double a,
                double * pvecback);
  
  double B4_scf(
                struct background *pba,
                double phi,
                double phi_prime,
                double rho_idm,
                double a,
                double * pvecback);

  double B5_scf(
                struct background *pba,
                double phi,
                double phi_prime,
                double rho_idm,
                double a,
                double * pvecback);

#ifdef __cplusplus
}
#endif

/**************************************************************/

/**
 * @name Some conversion factors and fundamental constants needed by background module:
 */

//@{

#define _Mpc_over_m_ 3.085677581282e22  /**< conversion factor from meters to megaparsecs */
/* remark: CAMB uses 3.085678e22: good to know if you want to compare  with high accuracy */

#define _Gyr_over_Mpc_ 3.06601394e2 /**< conversion factor from megaparsecs to gigayears
                                       (c=1 units, Julian years of 365.25 days) */
#define _c_ 2.99792458e8            /**< c in m/s */
#define _G_ 6.67428e-11             /**< Newton constant in m^3/Kg/s^2 */
#define _eV_ 1.602176487e-19        /**< 1 eV expressed in J */

/* parameters entering in Stefan-Boltzmann constant sigma_B */
#define _k_B_ 1.3806504e-23
#define _h_P_ 6.62606896e-34
/* remark: sigma_B = 2 pi^5 k_B^4 / (15h^3c^2) = 5.670400e-8
   = Stefan-Boltzmann constant in W/m^2/K^4 = Kg/K^4/s^3 */

//@}

/**
 * @name Some limits on possible background parameters
 */

//@{

#define _h_BIG_ 1.5            /**< maximal \f$ h \f$ */
#define _h_SMALL_ 0.3         /**< minimal \f$ h \f$ */
#define _omegab_BIG_ 0.039    /**< maximal \f$ omega_b \f$ */
#define _omegab_SMALL_ 0.005  /**< minimal \f$ omega_b \f$ */

//@}

/**
 * @name Some limits imposed in other parts of the module:
 */

//@{

#define _SCALE_BACK_ 0.1  /**< logarithmic step used when searching
                             for an initial scale factor at which ncdm
                             are still relativistic */

#define _PSD_DERIVATIVE_EXP_MIN_ -30 /**< for ncdm, for accurate computation of dlnf0/dlnq, q step is varied in range specified by these parameters */
#define _PSD_DERIVATIVE_EXP_MAX_ 2  /**< for ncdm, for accurate computation of dlnf0/dlnq, q step is varied in range specified by these parameters */

#define _zeta3_ 1.2020569031595942853997381615114499907649862923404988817922 /**< for quandrature test function */
#define _zeta5_ 1.0369277551433699263313654864570341680570809195019128119741 /**< for quandrature test function */

//@}


#endif
/* @endcond */
