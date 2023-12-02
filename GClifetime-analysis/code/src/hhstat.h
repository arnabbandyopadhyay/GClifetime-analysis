#ifndef HHSTAT_H
#define HHSTAT_H

#include "histo_monitor.h"

class hstat{
public:

    static const int mutation_bins = 20;
    static const int max_n_of_divisions = 12;
    static const int max_mutation_bin = 30;
    static const short phase_number = 5;
    static const int dtphase_resolution = 21;
    static const int max_gcs=1000;
    static const int fdc_max_encounters = 50;
    static const int max_n_of_ag_portions = 100;

    struct cb {

        double p_dif;
        double average_seeder_affinity;
        long attributed_mutation_prob[mutation_bins];
        long cummulative_attributed_mutation_prob[mutation_bins];
        int dtphase_frequency[phase_number][dtphase_resolution]; //cellCB
        long cummulative_attributed_n_of_divisions[max_n_of_divisions + 1];//cellCB
        long attributed_n_of_divisions[max_n_of_divisions + 1];//cellCB
    } cb_obj;

     struct fdc {

         double p_mksignal;
         long ab_sign_errors,ag_sign_errors,ic_sign_errors,ic_calculations;

     } fdc_obj;

     struct cc {

         long cummulative_ag_collection_all[max_n_of_ag_portions + 1];
         long cummulative_ag_collection_selected[max_n_of_ag_portions + 1];

         double p_dif2out,p_dif2out_DEC;

         int fdc_encounters[fdc_max_encounters];

     } cc_obj;

     struct out {

         double average_affinity,max_affinity; //cellOUT

     } out_obj;

    // alpha_mean averaged over all cells
    double frac_average;

    // OPTIONAL: look for key MOVE_ANALYSIS for activation (in files cell.C and hyphasma.C)
    // counts the total number of tries and successes of movement (summed over all cells)
     long n_try_move,n_move_done,n_move_removed,n_move_forbidden,n_move_self,
      n_move_forbidden_back,n_move_self_back; // from cell class : This was not checked.


    double antibody_production_factor;// AntibodyDyn class
    int mutation_frequency[max_mutation_bin]; //cellman

    long abs_cell_index; //trackfate
    long abs_track_index; //trackfate

    histo_monitor monitor_pureFDCsearch;
    // count the frequency of different numbers of CC-TC interactions
    histo_monitor monitor_nTCcontacts_selected;
    histo_monitor monitor_nTCcontacts_deleted;
    // count the frequency of acquired BC-Tfh search times
    histo_monitor monitor_BsearchTtime_selected;
    histo_monitor monitor_BsearchTtime_deleted;
    // count the frequency of acquired BC-Tfh interaction times
    histo_monitor monitor_BinteractTtime;
    // count the frequency of times between two Tfh bindings
    histo_monitor monitor_TimeBetweenTbindings;
    // count the frequency of the amount of integrated Tfh signal at time of selection:
    histo_monitor monitor_TfhSignalAtSelection;
    histo_monitor monitor_TfhSignalSpeed_selected;
    histo_monitor monitor_TfhSignalSpeed_deleted;
    // count the frequency of the AUC of Tfh signal at time of selection:
    histo_monitor monitor_TFHsignalAUC;

    hstat();
    };

#endif // HHSTAT_H
