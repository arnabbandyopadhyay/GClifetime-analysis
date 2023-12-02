/*
 * cellthis.h
 * --------------------------------------------------
 * Deklaration von spezifischen Zelltypen unter Verwendung der Routinen  in cell.h
 * --------------------------------------------------------------------------------
 *
 * Struktur der Deklaration:
 *
 * 1) Eintragen des neuen Typs (XX) in den states-enumerator in point.h
 *
 * 2) Erstelle neue Klasse cellXX in dieser Datei:
 *
 * class cellXX : public cell {...} // Zelle aus einem Fragment
 *
 * oder
 *
 * class cellXX : public frag_cell { // Zelle mit mehr Fragmenten
 * public:
 *  cellXX() { }               // notwendiger Konstruktur
 *  cellXX(const cellXX& x);   // notwendiger Selbst-Konstruktur mit Zuordnung von Werten
 *  cellXX(cellYY& x);         // optional fuer Zuordnung zwischen verwandten Zelltypen
 *  ~cellXX();                 // notwendiger Destruktur
 *  void ini(const long& i, const long& li, lattice& l, AffinitySpace& shape);
 *                             // schreibt eine neue Zelle auf lattice und AffinitySpace ein.
 *                 // Dies koennte in den Konstruktur wenn jede erzeugte
 *                 // Zelle auch auf lattice und AffinitySpace erscheinen soll.
 * }
 *
 * 3) Optional Erstellung eines zellspezifischen Zell-enumerators zur Unterscheidung
 * verschiedener Zellzustaende:
 *
 * enum XXstates{state1, state2};
 *
 * und in die cellXX-Klasse Zustands variable einfuehren:
 *
 * XXstates state;
 *
 * 4) Einfuehrung der statischen Variablen, die die Prozesse
 * quantitativ steuern (Wahrscheinlichkeiten p_):
 *
 * private:
 * static double p_pro,             // Proliferation
 *              max_pro_distance,  // fuer 1-Fragment: pushing distance
 *              p_mut,             // Mutation (verwendet AffinitySpace)
 *      p_apo,             // Apoptose
 *      p_difu,            // Diffusion/Motility
 *      p_difu_width,
 *      distance_tolerace, // Formparameter bei Diffusion
 *      p_grow,            // Wachstum
 *      p_shrink,          // Schrumpfung
 *      p_dif;             // Differenzierung
 * public:
 * static int target_volume; // Zielzellvolumen
 * static void set_statics(const Parameter& par, lattice& l, ofstream& ana);
 * static void set_statics(const double& time, const Parameter& par);
 *
 * Es ist zu beachten, dass nur die fuer den Zelltyp relevanten
 * p_-Variablen zu verwenden sind.
 *
 * target_volume ist nicht unbedingt von aussen sichtbar zu definieren.
 *
 * set_statics() muss von aussen zu steuern sein. Die Argumente haengen
 * vom neuen Typ ab. Die erste Variante dient der Initialisierung, wobei
 * die Resultate in ana rausgeschrieben werden. Mit der zweiten Variante
 * koennen die Parameter waehrend der Laufzeit geaendert werden.
 *
 * Will man die p_ fuer einige Prozesse dynamisch gestalten, sind die
 * entsprechenden p_ nicht static zu deklarieren.
 *
 * 4a) Die Werte p_ und andere werden aus Parameter eingelesen. Die Parameter
 * muessen dort deklariert und definiert sein. Dies ist in setparam.x zu tun.
 *
 * 5) Einfuegung der erwuenschten Prozesse (Bsp nur aus cell und frag_cell).
 * Bei einigen Prozessen ist die Abfrage, ob diese durchgefuehrt werden
 * getrennt zu programmieren, etwa wie durch ask_mitosis().
 *
 * public:
 * // Diffusion fuer mehrfragmentige Zellen
 * double diffuse(const long& li, lattice& l)
 *  { return fragdiffuse(XX,distance_tolerance,li,p_difu,l); }
 * // Diffusion fuer einfragmentige Zellen
 * short diffuse(const long& li, lattice& l)
 *  { return do_diffuse(XX,li,l); }
 *
 * // Proliferation
 * long ask_mitosis(long* pp, lattice& l) {
 *  if (volume==1 && target_volume==1) {
 *    return find_mitosis_place(p_pro,max_pro_distance,pp,l);
 *  } // im Fall der Klasse cell genuegt find_mitosis_place(...)
 *  else { // case of more fragment object
 *    // probabilistic decision if proliferation is done
 *    if ( volume>0.9*target_volume && drandom(1.) < p_pro ) return 0;
 *    // Proliferation is allowed if the total volume is near the target-volume of a cell!
 *    return 1;
 *  }
 * }
 * short mitosis(const long& i, const long& li,
 *      const long& newli, frag_cell& newCell,
 *      lattice& l)
 *  { return do_mitosis(XX,i,li,newli,newCell,l); }
 *
 * // Mutation
 * short mutate(AffinitySpace& shape)
 *  { return do_mutate(shape,p_mut); }
 *
 * // Wachstum und Schrumpfung (### letzteres nicht programmiert)
 * short grow(const long& li, lattice& l)
 *  { if (target_volume==1) return 1; // d.h. kein Wachstum!
 *    else return do_grow(XX,li,target_volume,p_grow,p_shrink,l); }
 *
 * // Kontaktfragen
 * short find_contact(states celltype, const long& i, lattice& l); //cell
 * short check_for_border(const long& i, lattice& l); //cell
 * short contact(states celltype, lattice& l); //frag_cell
 *
 * // Signal Produktion
 * void signal_production(const long& i, lattice& l) {
 *  // Produce signal for differentiation of CB to CC
 *  signal_secretion(i,sig_differ2CC,p_mksignal,l);
 *  // other signals may be included here
 * }
 *
 * 5a) Initialisierung einer neuen Zellliste fuer XX.
 *
 * 5b) Ansteuerung der Prozesse in cellman.x durch eine neue Routine
 * calc_XX() und eventuelle weitere Routinen, die Zelltypumwandlung
 * oder Zellerzeugung -- also Zell-Listen-relevante Prozesse betreffen.
 *
 * 6) Das sollte es sein.
 */

#ifndef i_cellthis
#define i_cellthis
#include "setparam.h"
#include "antibody.h"
#include "cell.h"
#include "ss.h"
#include "sequencespace.h"
#include "space.h"
#include "ode.h"
#include "trackfate.h"
#include "histo_monitor.h"
#include "hhstat.h"

// ===============================================================
// ===============================================================
// ===============================================================
// =================== Specific cell types =======================
// ===============================================================
// ===============================================================

// ===============================================================
// ======================= cellCB ================================
// ===============================================================

//static const int MAX_GCS=1000;
class cellCC;
class cellFDC;
class cellTC;

class immunoglobulin_class {
  private:
   static const int matrix_dimension = nIg_classes * nIg_classes;
   static int get_matrix_index(const Ig_classes &i, const Ig_classes &j);
   static double switch_matrix[matrix_dimension];
  public:
   immunoglobulin_class();
   ~immunoglobulin_class();
   static short do_switch_classes;
   static void load_matrix(const Parameter &par, ofstream &ana);
   Ig_classes Ig_class;
   void class_switch();
   void set_class(const immunoglobulin_class &c);
};

// ### dies in cellCB integrieren und dann ueberall cellCB:: davor schreiben
// +++ muss enum static gemacht werden? Bekommt jedes Objekt ein neues enum?
// +++ nein, laut blog im Netz ist das nur eine Typendeklaration und jede Instanz
// +++ der Klasse belegt dafuer keinen Speicher.
enum centroblasts {
   cb_normal,cb_differentiate,
   cb_G1,cb_G0,cb_S,cb_G2,cb_M,        // keep these together!
   cb_divide,cb_stop_dividing,
   cb_statenumber
};

class cellCB: public frag_cell {
  private:
   static double
     limit_volume,
     p_mut,
     p_mut_after_selection,
     p_mut_after_dec_selection,
     p_mut_affinity_exponent,
     p_difu,
     p_difu_width,
     p_tension,
     p_grow,
     p_shrink,
     //p_dif,
     p_dif_target,
     start_differentiate,
     max_pro_distance,
     diffusion_tolerance_min,
     diffusion_tolerance_steepness,
   // for receptors:
     receptors,
     receptor_activation,
     receptor_binding,
     receptor_dissociation,
   // movement and shape
     tension,
     elongation,
     K_elongation,
     smoothmove,
     persistence,
     v_slow_factor,
     p_switch_v,
     max_adhesion;
   static double dtphase[cb_statenumber];
   //static const short phase_number = 5;
   //static const int dtphase_resolution = 21;
   //static int dtphase_frequency[phase_number][dtphase_resolution][MAX_GCS];
   static double fraction_of_phase_width;
   static short v_modi,n_v_states;
   // static long cell_cycle_delay,DEC205_cell_cycle_delay;
   static bool p_mut_affinity_dependent;
   static double total_n_of_DEC205_divisions;
   static double total_n_of_divisions;
   static bool transmit_CCdelay2cellcycle;
   centroblasts progress_cycle_phase();

   double cycle_state_time;
   double receptor_ligand;               // receptors with ligand bound

 public:
   //static const int max_n_of_divisions = 12;
  // static long cummulative_attributed_n_of_divisions[max_n_of_divisions + 1][MAX_GCS];
 private:
   //static long attributed_n_of_divisions[max_n_of_divisions + 1][MAX_GCS];

   //static const int mutation_bins = 20;
   //static long attributed_mutation_prob[mutation_bins][MAX_GCS];
   //static long cummulative_attributed_mutation_prob[mutation_bins][MAX_GCS];

  public:
   static int target_volume;
   static short receptor_use;
   static bool ag_loaded_CB_stop_mutation;
   centroblasts state;                   // Differentiation state of CB
   void transmit_CCdelay2cycle(double waited_time, hstat &hs);
   bool DEC205,DEC205_ova,pMHCdeficient,diff2output;
   // +++++++++++++++ OPTION pMHCdeficient ++++++++++++++++++++++++++++++++
   static const double p_pMHCdeficient;
   // +++++++++++++++ OPTION pMHCdeficient ++++++++++++++++++++++++++++++++
   void attribute_DEC205(double fracofpos);
   double retained_ag;
   vector<int> collected_ag_portions;
   bool iamhighag;
   static bool ag_loaded_CB_diff2output;
   static double asymmetric_polarity_index,smooth_PI;
   static double p_pro,delta_p_pro,p_CXCR4down;
   //static double average_seeder_affinity[MAX_GCS];

   cellCB();
   cellCB(const cellCB &x);
   cellCB(const cellCC &x);
   ~cellCB();
   void destruct() { delete[] fragments; } // Destruktor des Fragment-Feldes
   void ini(const long &i, const long &li, const double &t, space &l, AffinitySpace &shape);
   void make_CB_new(hstat &hs);

   static void set_statics(const Parameter &par, space &l, ofstream &ana, hstat &hs);
   static void set_statics(const double &time, const Parameter &par, hstat &hs);
   static void set_differentiation(const double &time, hstat &hs);
   static bool SMOOTH_DIFFERENTIATION;
   static double smooth_differentiation_time;

   static double total_cell_cycle_duration();
   static bool fixed_number_of_divisions();
   double time_of_cycle_state_switch;
   int n_divisions2do;
   static void show_number_of_divisions(double&, ofstream&, hstat &hs);
   static void show_cummulative_number_of_divisions(hstat &hs);
   static void show_mutation_prob(double&, double&, double&, ofstream&, hstat &hs);
   static void show_cummulative_mutation_prob(hstat &hs);
   static void show_cell_cycle_phase_duration(hstat &hs);

   immunoglobulin_class IgX;
   static double IgE_factor_cellcycle, IgE_factor_divisions;

   // Aktionen:
   // =========
   void set_p_move();
   void set_CBfromCC(cellCB &cb, cellCC &cc, hstat &hs);
   double move(const long &li, space &l, sigs &s, TRACK &td, double &time, hstat &hs);

   void set_CXCR4expression();
   void resensitise4CXCL12(sigs &s);
   void adapt_specific_ag(double factor);

   static double ag_preloaded;
   void preload_with_ag();

   void set_remaining_divisions();
   /*  double inverse_erf(double x);
    * double get_sample_from_normal(const double& mean, const double& width);
    * double get_positive_sample_from_normal(const double& mean, const double& width);
    */
   double set_cycle_state_duration(centroblasts &s, hstat &hs);
   static centroblasts get_virtual_cell_cycle_phase(double waited);
   static bool shiftCCdelay2CBcycle();
   long ask_mitosis(long * pp, space &l);
   short mitosis(const long &i, const long &li,
                 const long &newli, frag_cell &newCell,
                 space &l)
   { return do_mitosis(CB,i,li,newli,newCell,l); }

   void set_mutation_after_TC(AffinitySpace &shape, hstat &hs);
   short mutate(AffinitySpace &shape)
   { return do_mutate(shape); }

   short grow(const long &li, space &l) {
      if (target_volume == 1) {
         return 1;                    // d.h. kein Wachstum!
      } else { return do_grow(CB,li,target_volume,p_grow,p_shrink,l); }
   }
   void get_new_state(const long &i, double &dt, space &l, sigs &s, hstat &hs);
   centroblasts set2differentiation();
   void set_adhesion() { get_adhesion(max_adhesion); }
   short ask_differentiate(hstat &hs);

   cellCB&operator =(const cellCB &x);
   cellCB&operator =(const cellCC &x);
};

// ====================================================
// ================== cellCC ==============================
// ====================================================

enum centrocytes {
   unselected,contact,FDCselected,TCcontact,selected,apoptosis
};

class cellCC: public cell {
  private:
   static double
     p_apo,
     p_apo4FDCselected,
     p_mph,
     p_sel,
     TCselect_prob,
     p_FDCsignalling,
     p_dif,
     // p_dif_DEC,
    // p_dif2out,
     p_dif2out_target,
    // p_dif2out_DEC,
     p_dif2out_DEC_target,
     p_final_differentiation,
     start_differentiate,
     p_difu,
     p_difu_width,
     persistence,
     v_slow_factor,
     p_switch_v;
   static short v_modi,n_v_states;
   static short apoptotic_motility_mode;
   static double p_apo_randomwalk;
   static short TC_CC_selection;
   static bool force_selection_at_TFHthreshold, negativeTCselection, BCstaysonTCbyTCtime;
   static short TFHsignal_delivery_mode;
   static double tc_time, tc_dec205ova_binding_time, tc_time_width,
     BTtime_K, BTtime_min, BTtime_max, BTtime_n,
     TFHsignal_delivery_min, TFHsignal_delivery_max, 
     TFHsignal_delivery_KpMHC, TFHsignal_delivery_n;
   static double TFHsignal_decay;
   static double get_TFHsignal_delivery_factor(double, bool);
   static const int TFHsignal_delivery_max_pMHC = 200;
   static double TFHsignal_delivery_factor[TFHsignal_delivery_max_pMHC];
   static short mode_of_setting_tc_time;
   static double max_tc_signal_portion, rescue_signal_threshold, SST_tc_signal;
   static bool ag_deleted_in_fresh_CC;
   static bool inhibitFDCapoptosis;
   static double prob2kill_noFDCcontactBCs;
   short get_tc_selected(AffinitySpace &shape, hstat &hs);
   void progress_selection_state(AffinitySpace &shape, hstat &hs);
   void make_apoptosis(const double& time, AffinitySpace &shape, hstat &hs);
   void return2unselected(AffinitySpace &shape, hstat &hs);
   int bound_ag_index;
   void add_collected_ag();
   void add_tc_signal();
   double FDCselected_clock;
   double tc_clock, tc_signal, TFHsignalAUC, tc_interaction_time;
   bool SSTactive;
   static double AgThreshold4Selection;
   static short tc_search_duration_mode;
   static double tc_search_duration_per_FDCcontact, tc_search_duration_fixed;
   double tc_search_duration;
   double get_tc_search_duration();
   double get_pMHC_presentation();
   double get_tc_interaction_time();
   double set_selected_CC_delay();
   double individual_dif_delay;

   static bool pMHC_dependent_division, SIND;
   static short mode_of_DND_of_Tfh;
   static double pMHC_dependent_P_standard;
   static double pMHC_dependent_P_min;
   static double pMHC_dependent_P_max;
   static double pMHC_dependent_K, TFHsignal_dependent_K;
   static double pMHC_dependent_nHill;
   static double pMHC_dependent_pMHC_of_2divisions, TFHsignal_of_P0divisions;
   static double TFHgradient_dependent_K, TFHgradient_of_P0divisions;
   static double pMHC_of_DEC205_ova;
   double get_pMHC_dependent_division();
   double get_signal_induced_number_of_divisions();
  // static const int max_n_of_ag_portions = 100;
   static short present_specific_ag2TC;
   int get_max_collected_ag(bool returnindex);
   static short outputfiles;
   static double write_trackfate;

   /* Signal molecules */
   static short ICOSL_upregulation_mode; 
   static bool ICOSL_dependent_Tfh_signals, ICOSL_memory;
   static double ICOSL_upregulation_time;
   static short FoxO_mode, dT_FoxO_start, dT_FoxO_reg;
   static bool stopFoxOonTFH;
   static double FoxO_ini, FoxO_production, 
     FoxOup_min, FoxOup_max, FoxOup_K, FoxOup_n, KFoxO, nFoxO;
   static double mTORC1_production;
   void inhibitFoxO(double& pMHClevel);
   double get_FoxO_rate(double& pMHClevel);
   double ICOSL;
   double get_ICOSL_expression();
   double mTORC1, FoxO, FoxO_upregulation;
   bool hadcontact2Tfh;
   fateTRACK fatetracker;
   static int fatetracker_ndt;
   int fatetracker_n;
   bool came_from_FDCselected;
   double TimeOfLastTbinding;

   static void antigen_collection_statistics(int, bool, hstat &hs);
   // count the frequency of the duration of the period of pure search for Ag on FDC
  /* static histo_monitor monitor_pureFDCsearch[MAX_GCS];
   // count the frequency of different numbers of CC-TC interactions
   static histo_monitor monitor_nTCcontacts_selected[];
   static histo_monitor monitor_nTCcontacts_deleted[MAX_GCS];
   // count the frequency of acquired BC-Tfh search times
   static histo_monitor monitor_BsearchTtime_selected[MAX_GCS];
   static histo_monitor monitor_BsearchTtime_deleted[MAX_GCS];
   // count the frequency of acquired BC-Tfh interaction times
   static histo_monitor monitor_BinteractTtime[MAX_GCS];
   // count the frequency of times between two Tfh bindings
   static histo_monitor monitor_TimeBetweenTbindings[MAX_GCS];
   // count the frequency of the amount of integrated Tfh signal at time of selection:
   static histo_monitor monitor_TfhSignalAtSelection[MAX_GCS];

   static histo_monitor monitor_TfhSignalSpeed_selected[MAX_GCS];
   static histo_monitor monitor_TfhSignalSpeed_deleted[MAX_GCS];
   // count the frequency of the AUC of Tfh signal at time of selection:
   static histo_monitor monitor_TFHsignalAUC[MAX_GCS];*/
   // count the frequency of Tfh-signal-speed in BCs
   double calc_TfhSignalSpeed();
  public:
   cellCC();
   cellCC(const cellCC &x);
   cellCC(cellCB &x);
   ~cellCC();
   static void set_statics(const Parameter &par, space &l, ofstream &ana, hstat &hs);
   static void set_statics(const double &time, const Parameter &par, hstat &hs);
   void trackfate_initialize(hstat &hs);
   void trackfate(double, bool);
   void trackfate_show(hstat &hs);
   void writetrackfate(hstat &hs);
   //static const int fdc_max_encounters = 50;
   //static int fdc_encounters[fdc_max_encounters][MAX_GCS];
   static unsigned short CC_FDC_selection;
   static double p_CXCR5down;
   static void set_differentiation(const double &time, hstat &hs);
   static bool SMOOTH_DIFFERENTIATION;
   static double smooth_differentiation_time;
   static bool collectFDCsignals;
   static bool multipleTFHcontacts;
   static int reset_antigen_after_collection;
   static double ignore_affinity;
   static double dif_delay,dif_delay_DEC;
   static bool ag_loaded_CC_directly2TFH;
   static long test_delay, ICAM_delay;
   static int AgThreshold4Diff;
   static double collectFDCperiod;
   static double IgE_BCRlevel, IgE_prob_CXCR5down;
   static void show_cummulative_ag_collection(hstat &hs);
   //static long cummulative_ag_collection_all[max_n_of_ag_portions + 1][MAX_GCS];
   //static long cummulative_ag_collection_selected[max_n_of_ag_portions + 1][MAX_GCS];
   static bool simultaneousTfhFDC;

   short selectable,mobile;
   centrocytes state;
   double affinity;
   double selected_clock;
   bool selected4output;
   bool DEC205,DEC205_ova,pMHCdeficient;
   long tc_index, last_tc_index;
   double fdc_clock;
   int nFDCcontacts,nTCcontacts;
   vector<int> collected_ag_portions;
   immunoglobulin_class IgX;
   double BCRexpression;
   short CXCR5failure;
   double pMHC_dependent_number_of_divisions;
   double get_FoxO();

   void set_p_move();
   short move(const long &li, space &l, sigs &s, TRACK &td, double &time);
   void set_selectable();
   void set_CXCR5expression();
   void resensitise4CXCL13(sigs &s);
   void set_CXCR4expression();
   void resensitise4CXCL12(sigs &s);
   bool set_apoptotic_motility(sigs &s);
   void go2TCselection(AffinitySpace &shape);
   void delete_antigen();
   void process_antigen();
   long contact2FDC(space &l);
   short bind_antigen(cellFDC &fdc, int frag_index, int &ag_index,
                      AffinitySpace &shape, double &threshold);
   short select(AffinitySpace &shape, hstat &hs);
   short stop_collecting_FDCsignals(AffinitySpace &shape, const double& time, double &dt, hstat &hs);
   short dif2OUTorCB(double&, hstat &hs);
   bool final_differentiation();
   void attribute_tc_search_duration();
   bool same_TC_as_before(long new_tc_index);
   bool try2findTFH();
   void bind_TC(const double& time, cellTC &tcell, space &l, AffinitySpace &shape, hstat &hs);
   void progress_signalling(); 
   void decay_TFHsignal();
   void addTFHsignalAUC(double& dt);
   short got_tc_signals(const double &time, const double &dt, cellTC &tcell, 
            space &l, AffinitySpace &shape, hstat &hs);
   bool do_selection();
   short try_selection(const double &time, const double &dt, AffinitySpace &shape, hstat &hs);
   static short mode_of_apoptosis;
   short apoptose(const double &time, AffinitySpace &shape, hstat &hs);
   short macrophagocyte(AffinitySpace &shape, hstat &hs);
   static void update_histograms(hstat &hs);
   static void write_histograms(hstat &hs);

   cellCC&operator =(const cellCB &x);
   cellCC&operator =(const cellCC &x);
};

// ====================================================
// ================== cellTC ==========================
// ====================================================

enum tcells {
   TCnormal,TC_CCcontact,TCdivide
};

// class cellTC: public frag_cell {
/* This can be activated together with a change in cellTC::cellTC(..) in cellthis.cpp.
 * It generates a normal run (tested with bcinflow09rand0.par).
 * But I want to check what exactly happens to all the additional parameters
 * introduced by frag_cell. Probably they are just there and do nothing. Better to check.
 * ### no frag_cell yet for TC!!!
 */
class cellTC: public cell {
  private:
   static double p_difu, p_difu_width, persistence;
   static double proliferation, meancycle, cyclewidth, max_distance2divide;
   static int Ndivisions;

   long * CC_nn;
   short n_CC_nn;
   double * CC_affinity;
   short changed_nn;
   double cell_cycle_clock,cell_cycle_duration;

   static double dTrepolarise, dTsignal;
   double TsinceBCpolarisation;
   bool able2repolarise;
   long CurrentSignalTarget;
   void ini_polarisation_state();

  public:
   static bool do_division;
   // state of TC
   tcells state;
   bool able2signal;

   cellTC();
   cellTC(const cellTC &x);
   ~cellTC();
   void ini(const long &i, const long &li, const double &t, space &l, AffinitySpace &shape);
   static void set_statics(const Parameter &par, space &l, ofstream &ana);
   //  static void set_statics(const double& time, const Parameter& par);
   void evolve_polarisation_state(double&);
   void set_changed_nn(short x);

   // Aktionen:
   // =========
   void set_p_move();
   short move(const long &li, space &l, sigs &s, TRACK &td, double &time);
   void make_tc_cc_link(const long &index, const double &nFDCcontacts);
   void make_tc_cc_link(const long &index,
                        const long &CCpos,
                        int ag_index,
                        AffinitySpace &shape,
                        const bool &highag);
   void liberateCC(const long &index);
   void set_polarity(space &l);

   // Division:
   double set_cell_cycle_duration();
   void reset_cycle_times();
   void ask_enter_cell_cycle();
   bool progress_cell_cycle(double &dt);
   long ask_mitosis(long * pp, space &l);

   cellTC&operator =(const cellTC &x);
};

// ===============================================================
// ======================= cellFDC ===============================
// ===============================================================

enum FDCstates {
   none,soma,dendrite
};

class cellFDC: public frag_cell {
  private:
   //static double p_mksignal;
   static double p_mkCXCL13,p_mkSEMA4D;
   static short vesicle;
   static double ic_k_on,ic_k_off;
   static int FDCmaxFrags,n_Antigen_max,n_Antigen_dim_factor;
   static vector<double> ag_fraction;
   static short ag_distribution_mode, ag_detection_mode;
   int get_highest_amount_ag(const int &frag_index, AffinitySpace &AS);
   int get_highest_affinity_ag(const int &frag_index, const long &BCRposAS, AffinitySpace &AS);
   double get_total_antigen_at_site(int frag);
  public:
   cellFDC();
   cellFDC(const cellFDC&);
   ~cellFDC();
   static vector<double> ini_ag_fraction();
   static void set_statics(const Parameter &par, space &l, ofstream &ana, hstat &hs);
   static void set_statics(const double &time, const Parameter &par, hstat &hs);
   void signal_production(const long &i, sigs &l, hstat &hs);
   void add_antigen(int n_Antigen);
   double get_total_antigen();
   double get_total_antigen(int ag_index);
   int get_voidsites();
   int get_fragment_index(const long &fragpos);
   short consume_ag(const long &frag_index, const int &ag_index);
   int local_interaction_with_ag(const int &frag_index, const long &BCRpos_ss,
                                 AffinitySpace &shape);
   double get_total_immune_complex();
   double get_total_immune_complex(const int &length);
   double get_total_immune_complex(const int &ag_index, const int &length);

   FDCstates state;
   /*multiAg: evtl. do double** antigen_amount, where the 2nd dim runs over the different antigens,
    * i.e. antigen_amount[f][a] with f running over FDC fragments and a running over antigens.
    * As a vector one could add antigen later on. But not necessary, as it will be known how many
    * antigens
    * will be presented in the time of the GC, so dimension of array can be set accordingly.
    */
   double * * antigen_amount;
   /*multiAg: Now: ic_amount[f][b] with f running over FDC fragments and b running over affinity
    * bins.
    * For each FDC fragment we have a fixed Ag and a set of Ab affinity bins associated with it.
    * This recollects all Abs in the GC.
    * For multiple Ags, each Ag requires its own Ab affinity bin set, in order to calculate and save
    * the ic_amount,
    * thus, we need ic_amount[f][b][a] with a running over the antigens.
    */
   double * * * ic_amount;
   static unsigned short use_antigen;
   static double antigen_amount_per_FDC,antigen_saturation;
   short antigen_presence(const long &fragpos, int &ag_index);
   void mk_immune_complex(sigs &l);
   void mk_immune_complex(const double &d_t, sigs &l, hstat &hs);
   void mk_immune_complex(const double &d_t, AntibodyDyn &ABS, AffinitySpace &AS, hstat &hs);
   static double ag_threshold;
   //static long ab_sign_errors,ag_sign_errors,ic_sign_errors,ic_calculations;
   void set_antigen_amount(double time, double dt);

   static int DendriteLength;

   cellFDC&operator =(const cellFDC &x);
};

// ===============================================================
// ======================= cellOUT ===============================
// ===============================================================

class cellOUT: public frag_cell {
  private:
   static double 
     p_difu,
     p_difu_width,
     persistence,
     v_slow_factor,
     p_switch_v;
   static short v_modi,n_v_states;
   static short vesicle;
   static bool exit2tz;

  public:
   cellOUT();
   cellOUT(const cellOUT &x);
   ~cellOUT();
   static void set_statics(const Parameter &par, space &l, ofstream &ana, hstat &hs);
  /* static double average_affinity[MAX_GCS],max_affinity[MAX_GCS];

  */
   static double initial_ab_affinity;
   static short use_threshold;
   static double p_mk_ab;

   immunoglobulin_class IgX;
   bool DEC205;
   short try_eject(long i, space &l, AffinitySpace &shape);
   void set_p_move();
   short move(const long &li, space &l, sigs &s, TRACK &td, double &time);
   void signal_production(const long &i, sigs &l);

   cellOUT&operator =(const cellOUT &x);
   cellOUT&operator =(const cellCC &x);
   cellOUT&operator =(const cellCB &x);
};

// ===============================================================
// ======================= cellbeta ==============================
// ===============================================================

// enum betacells{};

class cellbeta: public frag_cell {
  private:
   // constants:
   static const double Avogadro;// = 6.02205e+23; // mol^-1
   static const double MFaraday;// = 9.6485309e-02; // Faraday constant in C/(micromol)
   // 9.6e+04 C/mol = 9.6e+04*1e-06 C / 1e-06mol = 9.6e-02 C/micromol
   static const double Faraday;// = 9.6485309e+04; // Faraday constant in C/(mol)
   static const double Rydberg;// = 8.315; // in J/(K*mol)
   static const double pi;// = 3.141592654;

   enum beta_currents {
      NaK,K_ATP,K_V,K_Ca,sK_Ca,Na_V,fNa_V,NCX,PMCA,Ca_L,Ca_T,
      SERCA,cIP3,N_currents
   };
   static double I[N_currents];
   static void get_currents(double * y,
                            double &Vbar_K, double &Vbar_Na, double &Vbar_Ca, double &Vbar_ER,
                            double &C_IP3_inh);
   static void get_current_factors(double t);

   static double
      max_pro_distance,
      p_tension,
      p_grow,
      p_shrink,
      p_difu,
      diffusion_tolerance_min,
      diffusion_tolerance_steepness,
   // movement and shape
      tension,
      elongation,
      K_elongation,
      smoothmove,
      persistence,
      v_slow_factor,
      p_switch_v,
      max_adhesion;
   static short v_modi,n_v_states;

   // needed for electrophysiology calculation:
   static double xi,xi_ER,xi_ERC,xi_S_ERC;
   static double J_K,J_Na,J_Ca,J_ions;
   static double V_ER_0;

   static double get_2sigmoidal(double &t, double t_a, double t_b,
                                double rest, double factor,
                                double kappa_a, double kappa_b);
   double get_glucose(const long &i, sigs &s);
   static double get_IP3(double &t);
   static double get_tau_IP3(double &ip3);
   static double get_sigmoidal(double &half, double &x, double &kappa);
   static double get_inactivation(double &half, double &x, double &kappa);
   static double get_Hill(double &x, double &half, double &coefficient);
   static double get_ca_buffer(double &b_0, double &c, double &dissociation);
   static double get_Ca_free_fraction(double &c);
   static double get_Ca_ER_free_fraction(double &c);
   static double get_dynamic_half(double &x, double a, double b);
   static double get_tau_V(double &x, double a, double b, double c, double Vx);
   static double get_tau_K_V_sherman88(double &x, double &VbarK);
   static double get_K_ext(double &t, double &k_0);
   static void get_all_Nernst(double &t, double &K, double &Na, double &Ca, double &Ca_ER,
                              double &Vbar_K, double &Vbar_Na, double &Vbar_Ca, double &Vbar_ER);
   static void insert(double &ion, double &voltage, double I_load, double valence_sign);

   void show(double t, double * y_n, ofstream &file, ofstream &file2);
   void show_rev(double t, double * v, ofstream &file);
   void show_c(double t, ofstream &file);
   void show_cr(double t, ofstream &file);
   void show_gap(double t, double * y, ofstream &file);

   // Declare cell specific variables
   // betacell specific variables for each cell
   double * y_n1; // try these to be static! ###
   double * y_n_old;
   static double betadt;
   static long betandt;
   void set_initial_values();

   // solver related functions:
   // static const ode_method method=RungeKutta_4th;     // available methods are listed in ode.h
   ode_method method;
   static ode solver;
   void (*prhs)(double, double*, double*); // Pointer on rhs(...)
   void get_gap_junction(dynarray<cellbeta> &bl, space &l);

   // each cell gets its own files
   ofstream result,result2,currents,currents_rho,reversal,gapfile;
   void open_files();
   static long n_write;
   long step_count;

  public:
   // ++++++++++++++++++++++++++++ OPTION +++++++++++++++++++++++++++++
   // if the number of betacells is larger than 200 this flag shall be true
   static const bool LOCAL_FILES = false;
   // choose true here, if gap-junctions shall also be treated with Runge-Kutta
   static const bool FULL_RUNGE = true;
   // set whether the sequence of cells in cellman::calc_CB shall be randomised
   static const bool RANDOMISE_SEQUENCE = false;
   // ++++++++++++++++++++++++ end OPTION +++++++++++++++++++++++++++++
   static const double z_K;// = 1.0; // valence of potassium ions
   static const double z_Na;// = 1.0; // valence of sodium ions
   static const double z_Ca;// = 2.0; // valence of calcium ions
   static betaWerte p;
   static double p_pro;
   static double glucose_rest;
   static int target_volume;
   enum beta_quantities {
      V,K,Na,Ca,                                                  // 0-3
      g_K_ATP,g_K_V,g_Na_V,g_Ca_L,g_Ca_T,                 // 4-8
      h_K_V,h_Na_V,h_Ca_L,h_Ca_T,                         // 9-12
      g_K_Ca,C_K_Ca,                                      // 13-14
      Ca_ER,V_ER,IP3,g_IP3,h_IP3,                         // 15-19
      g_sK_Ca,glu,                                        // 20-21
      gap_K,gap_K_0,gap_K_1,gap_K_2,gap_K_3,gap_K_4,gap_K_5,                 // 22-28
      gap_Na,gap_Na_0,gap_Na_1,gap_Na_2,gap_Na_3,gap_Na_4,gap_Na_5,          // 29-35
      gap_Ca,gap_Ca_0,gap_Ca_1,gap_Ca_2,gap_Ca_3,gap_Ca_4,gap_Ca_5,          // 36-42
      g_fNa_V,                                                               // 43
      N_equations
   };                                                                        // 44

   void randomise_protein_expression();

   // The actual quantities for each cell
   double * y_n;
   // double rho_gap,rho_K_ATP;
   double * rho;

   cellbeta();
   cellbeta(const cellbeta &x);
   ~cellbeta();
   void destruct() { delete[] fragments; } // Destruktor des Fragment-Feldes
   void ini(const long &i, const long &li, const double &t, space &l);
   static void set_statics(const Parameter &par, space &l, ofstream &ana);
   void synchronise();

   // Aktionen:
   // =========

   static double get_Nernst(double &x, double &x_ext, double valence);
   static void rhs(double t, double * y, double * derivative);
   void electrophysiology(double thr, double dthr, sigs &s,
                          dynarray<cellbeta> &bl, space &l);

   void set_p_move();
   double move(const long &li, space &l, sigs &s, TRACK &td, double &time, hstat &hs);
   long ask_mitosis(long * pp, space &l);
   short mitosis(const long &i, const long &li,
                 const long &newli, frag_cell &newCell,
                 space &l)
   { return do_mitosis(CB,i,li,newli,newCell,l); }
   short grow(const long &li, space &l) {
      if (target_volume == 1) {
         return 1;                    // d.h. kein Wachstum!
      } else { return do_grow(CB,li,target_volume,p_grow,p_shrink,l); }
   }
   void get_new_state(const long &i, sigs &s);
   void set_adhesion() { get_adhesion(max_adhesion); }
   void show_all(double t);
   void show_all(double t,
                 ofstream &f_a,ofstream &f_b,ofstream &f_i,
                 ofstream &f_r,ofstream &f_n,ofstream &f_g);
   static suffix beta_index;

   cellbeta&operator =(const cellbeta &x);
};

// ===========================================
// (Vergleichs)-Operatoren
char operator ==(const cellCB &a, const cellCB &b);
char operator !=(const cellCB &a, const cellCB &b);

#endif
