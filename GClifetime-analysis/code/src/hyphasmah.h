#ifndef HYPHASMAH_H
#define HYPHASMAH_H
#include "cellman.h"
#include<sstream>
#include<string.h>


class hyphasma{
public:

    int var;
    ofstream ana;
    hstat hs;
    space l;
    sigs s;

    hyphasma(){}
    hyphasma(Parameter &par);

};

class hyphasma2: public hyphasma
{
public:
//hstat hs; // hyphasma class
SS newshape;
AntibodyDyn Ab;
cellman c;
static int GCref;
hyphasma2(){}

hyphasma2(Parameter &par);

static void antibody_sum(hyphasma2 **hyp); //updates antibodies
static double add_absum(int ag_no); // antibody total from all bins of all GCs
static void initialize_abs(); // clear ab_sum variables
static string prefix_files(string file);
static const int max_ags=100;
static const int max_resolution=11;
static int newvar;

// Apoptosis of plasma cells
static bool PCsUndergoApoptosis;

//extrafollicular pcs
static bool includeEFPCs;
static double nEFPC,newEFPCstoptime; // number of EFPC and stop time to introduce new EFPCs
static double EFPC_affinity; // constant affinity for all EFPCs
static double EFPC_antibody; // antibody conc. from all EFPCs
static double EFPC_abDrate; // EFPC antibody deg.rate
static double EFPC_abPrate; // EFPC antibody prod.rate
static double EFPC_Prate,EFPC_Apoprate; // EFPC formation rate and apoptosis time

// variables to store antibodies from all GCs
static double ab_sum[max_resolution][max_ags];

static int ngcs_tt; //total number of gc objects
static int ngcs; //number of gcs with vol>0, calculated using N
static int small,medium,large; // number of small, medium &large gcs
static int num_ags; // number of antigens
static int DZ_count,LZ_count; // number of DZ and LZ cells excluding apoptotic

// overall readouts (from all GCs)
static double IP_tot[max_ags]; // total IP
static double gcvol_tot;   // total gcvolume
//static double abs_tot;
static double pcs_tot,out_tot,cbs_tot,ccs_tot,outext_tot;   // total number of plasma cells, out cells,cbs, ccs
static double PCaff_tot[max_ags];  // mean affinities of all plasma cells
static double OUTaff_tot[max_ags]; // includes outext and out
static double CBaff_tot[max_ags];
static double CCaff_tot[max_ags];

static double calculate_IP(int agno);
double calculate_gcvol(hyphasma2 **hyp,int nogc);
double calculate_freeag(hyphasma2 **hyp,int agno); // only for simulated GCs
double calculate_abs();
double calculate_nCells(hyphasma2 **hyp,int nogc,cells celltype);
double count_output(hyphasma2 **hyp,int nogc,space &l);
int calc_dzcount(hyphasma2 **hyp,int nogc,space &l);
int calc_lzcount(hyphasma2 **hyp,int nogc,space &l);
void calculate_PCaff_allAG(hyphasma2 **hyp,int Hno);
void calculate_cellaff_allAG(hyphasma2 **hyp,int Hno,double time,space &l,SS &newshape);
double calculate_sumaff(hyphasma2 **hyp,int agno,int nogc,cells celltype);
double calculate_outaff(hyphasma2 **hyp,int agno,int nogc,space &l,SS &newshape);
double antibodysum_eachgc(hyphasma2 **hyp,int nogc);

//other readouts from previous version
void write_sel_analyze(Parameter &par);
static void write_sel_readme();
void data_analysis(); //for analysis after simulation

//EFPC functions
// Extrafollicular plasmacells
static void update_newEFPC(double dt); // update no of EFPCs
static void update_EFPCapoptosis(double dt);
static void update_EFPCantibody(double dt); // update the antibody concentrations from EFPC
static void add_EFPCab_tobin(double dt);

};

#endif // HYPHASMAH_H
