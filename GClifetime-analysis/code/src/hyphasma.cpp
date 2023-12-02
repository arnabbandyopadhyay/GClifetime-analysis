// #include <math.h>
#include <stdlib.h>
#include <ctime>
#include <stdio.h>
#include <math.h>
#include "hyphasmah.h"
#include "hhstat.h"


hstat::hstat(){

    cerr<<endl<<"hstat constructor";
    for(int i=0;i<max_n_of_divisions + 1;i++){
        cb_obj.cummulative_attributed_n_of_divisions[i]=0;
        cb_obj.attributed_n_of_divisions[i]=0;

    }
    for(int i=0;i<mutation_bins;i++){
        cb_obj.attributed_mutation_prob[i]=0;
        cb_obj.cummulative_attributed_mutation_prob[i]=0;
    }

    cb_obj.p_dif=0.0;
    cb_obj.average_seeder_affinity=0.0;
    fdc_obj.p_mksignal=0.0;
    fdc_obj.ab_sign_errors=0;
    fdc_obj.ag_sign_errors=0;
    fdc_obj.ic_sign_errors=0;
    fdc_obj.ic_calculations=0;

    frac_average=0;
    antibody_production_factor=1e-17;
    abs_cell_index=0;
    abs_track_index=0;

    out_obj.average_affinity=0;
    out_obj.max_affinity=0;

    cc_obj.p_dif2out=0;
    cc_obj.p_dif2out_DEC=0;

    for(int i=0;i<hstat::phase_number;i++){
        for(int j=0;j<hstat::dtphase_resolution;j++){
            cb_obj.dtphase_frequency[i][j]=0;
        }
    }

    for(int i=0;i<hstat::max_n_of_ag_portions+1;i++){
        cc_obj.cummulative_ag_collection_all[i]=0;
        cc_obj.cummulative_ag_collection_selected[i]=0;

    }
    for(int i=0;i<hstat::fdc_max_encounters;i++){
       cc_obj.fdc_encounters[i]=0;
    }

    for(int i=0;i<hstat::max_mutation_bin;i++){
        mutation_frequency[i]=0;
    }

}

hyphasma::hyphasma(Parameter &par)
     :ana(hyphasma2::prefix_files("ana_ini.out")),
     l(par, ana),
     s(l,par,ana)

{ 
    cerr<<"Hyphasma constructor";

    immunoglobulin_class::load_matrix(par, ana);
    cell::set_statics(par, ana);
    cellCB::set_statics(par, l, ana, hs);
    cellCC::set_statics(par, l, ana, hs);
    cellTC::set_statics(par, l, ana);
    cellOUT::set_statics(par, l, ana, hs);
    cellFDC::set_statics(par, l, ana, hs);
    cellbeta::set_statics(par, l, ana);
    AntibodyDyn::set_statics(par, ana, hs);
    arupProtein::set_statics(par, ana);
}


double hyphasma2::ab_sum[max_resolution][max_ags]={0};
double hyphasma2::gcvol_tot=0.0;
double hyphasma2::IP_tot[max_ags]={0.0};
double hyphasma2::PCaff_tot[max_ags]={0.0};
double hyphasma2::CBaff_tot[max_ags]={0.0};
double hyphasma2::CCaff_tot[max_ags]={0.0};
double hyphasma2::OUTaff_tot[max_ags]={0.0};

double hyphasma2::pcs_tot=0.0;
double hyphasma2::cbs_tot=0.0;
double hyphasma2::ccs_tot=0.0;
double hyphasma2::out_tot=0.0;
double hyphasma2::outext_tot=0.0;

//double hyphasma2::abs_tot=0.0;

int hyphasma2::newvar=0;//not used
int hyphasma2::GCref=0;//reference number for each hyphasma objects
int hyphasma2::ngcs_tt=0; //total number of hyphasma objects
int hyphasma2::ngcs=0; //total number of GCs calculated from N of all instances
int hyphasma2::small=0;
int hyphasma2::medium=0;
int hyphasma2::large=0;

int hyphasma2::DZ_count=0;
int hyphasma2::LZ_count=0;

int hyphasma2::num_ags=1; //Number of antigens
bool hyphasma2::PCsUndergoApoptosis=false;

//EFPCS
bool hyphasma2::includeEFPCs=false;
double hyphasma2::EFPC_abDrate=log(2.0)/(24.0*14);
double hyphasma2::EFPC_abPrate=(5.0e-18)/((1.0e-2)*(1.0e+15)); // distributing abs in a vol 10ml
double hyphasma2::EFPC_affinity=0.0; // add to bin0
double hyphasma2::EFPC_antibody=0.0;
double hyphasma2::EFPC_Apoprate=log(2.0)/24.0;
double hyphasma2::EFPC_Prate=log(2.0)/0.005;
double hyphasma2::nEFPC=0.0; //start with no EFPCs
double hyphasma2::newEFPCstoptime=5.0*24.0; //days to hours

hyphasma2::hyphasma2(Parameter &par)
    :hyphasma(par),
    newshape(par,ana),
    Ab(par,newshape,ana),
    c(par,l,s,newshape,ana, hs)
{


    cerr<<"hyphasma2 constructor";
    c.show_ag_BrdU(true);

    if (par.Value.show_mode == islet) {
       c.show_BETA();
    }


    if (sigs::TEST_MODE) {
       s.write_TEST(c.time);
    }

}

// adds prefix to all files depending on the GC ref no.
string hyphasma2::prefix_files(string file){
    std::stringstream filename;
    filename << "H" << hyphasma2::GCref << "/" << file;
    file = filename.str();
    return file;

}

// sum of abs from all the bins of all the GCs
double hyphasma2::add_absum(int ag_no){
    double total=0;
    for(int i=0;i<AntibodyDyn::antibodies_resolution+1;i++)//for all ab bins
    {
     total+=ab_sum[i][ag_no];
    }
    return total;
}


// update the antibodies
void hyphasma2::antibody_sum(hyphasma2 **hyp){

    for(int j=0;j<num_ags;j++){// for all antigens

        for(int i=0;i<AntibodyDyn::antibodies_resolution+1;i++)//all ab bins
        {
         ab_sum[i][j]=ab_sum[i][j]+(*(*hyp)).Ab.ab_bins[j].antibodies[i];
        }
    }
}


void hyphasma2::initialize_abs(){
    newvar=0;

    for(int j=0;j<num_ags;j++){// for all antigens

        for(int i=0;i<AntibodyDyn::antibodies_resolution+1;i++)
        {
            ab_sum[i][j]=0;
        }
    }
}

double hyphasma2::calculate_gcvol(hyphasma2 **hyp, int nogc){

    double gcvol;
    gcvol=(*(*hyp)).c.gcvol();
    gcvol*=nogc;

    return gcvol;
}

double hyphasma2::calculate_freeag(hyphasma2 **hyp, int agno){

    double freeag;
    freeag=(*(*hyp)).c.return_freeag(agno);

    return freeag;
}


double hyphasma2::antibodysum_eachgc(hyphasma2 **hyp,int nogc){

    double abtotal=0.0;

    for(int j=0;j<1;j++){// for only 1 antigen

        for(int i=0;i<AntibodyDyn::antibodies_resolution+1;i++)//all ab bins
        {
         abtotal=abtotal+(*(*hyp)).Ab.ab_bins[j].antibodies[i];
        }
    }

    abtotal=abtotal/nogc; // antibodies will be from N GCs

    return abtotal;
}


double hyphasma2::calculate_nCells(hyphasma2 **hyp,int nogc,cells celltype){
    double nCells;
    nCells=(*(*hyp)).newshape.nCells(celltype);
    nCells*=nogc;
    return nCells;

}

double hyphasma2::count_output(hyphasma2 **hyp,int nogc, space &l){
  double nCells;
   nCells=(*(*hyp)).c.count_out(l);
    nCells*=nogc;
    return nCells;
  
}

int hyphasma2::calc_dzcount(hyphasma2 **hyp,int nogc, space &l){
  int nCells;
   nCells=(*(*hyp)).c.count_dz(l);
    nCells*=nogc;
    return nCells;

}

int hyphasma2::calc_lzcount(hyphasma2 **hyp,int nogc, space &l){
  int nCells;
   nCells=(*(*hyp)).c.count_lz(l);
    nCells*=nogc;
    return nCells;

}


void hyphasma2::calculate_PCaff_allAG(hyphasma2 **hyp,int Hno){
    // this function calls a function from shape space which writes the PCaff for all ags of this particular GC
     (*(*hyp)).newshape.write_PCaff_allAGS(Hno);
}

void hyphasma2::calculate_cellaff_allAG(hyphasma2 **hyp,int Hno,double time,space &l,SS &newshape){
    // (*(*hyp)).newshape.write_cellaff_allAGS(Hno,time);
    (*(*hyp)).c.write_cellaff_allAGS(Hno,time,l,newshape);
}

double hyphasma2::calculate_sumaff(hyphasma2 **hyp, int agno,int nogc,cells celltype){
    double sum_aff;
    sum_aff=(*(*hyp)).newshape.Cellaff_sum(agno,nogc,celltype);

    return sum_aff;
}

double hyphasma2::calculate_outaff(hyphasma2 **hyp, int agno, int nogc, space &l, SS &newshape){
  double sum_outaff;

  sum_outaff=(*(*hyp)).c.calc_outaff_sum(agno, nogc, l, newshape);

  return sum_outaff;
}

double hyphasma2::calculate_IP(int agno){
    double tt_antigen=1e-5; // antigen concentration in M
     double ic_allbins=0.0;

        for(int i=0;i<AntibodyDyn::antibodies_resolution+1;i++){
            double ki=0;
            ki=pow(10,-5.5-0.4*i);

            double ic=0.0;

            ic=(ab_sum[i][agno]*1.e+15*tt_antigen)/(ki+tt_antigen);
            if(ic<0){
                ic=0;
            }

            ic_allbins+=ic;
        }

       return ic_allbins/tt_antigen;

    }

 void hyphasma2::write_sel_readme(){

     ofstream sel_ana_text("sel_readme.out");
     sel_ana_text
        << "Definition of the values stored in the file sel_analyze.out:\n"
        << "------------------------------------------------------------\n"
        << "\n"
        << "The values shown in a single line are the following\n"
        << "from left to right:\n\n"
        << "column 01: Number of FDC sites\n"
        << "column 02: Number of FDC\n"
        << "column 03: Length of FDC arms in microns\n"
        << "column 04: Antigen presentation per FDC (in Binding-Quanta)\n"
        <<
        "column 05: Antigen amount per FDC fragment below which binding probability is linearly reduced\n"
        << "column 06: Number of TC\n"
        << "column 07: Duration of CC-TC contact (hours)\n"
        << "column 08: Minimum polarisation duration of TC towards CC for selection\n"
        << "column 09: Antibody production of OUTPUT in Mol per hour and cell\n"
        << "column 10: IC=Ab-Ag dynamics: k_on\n"
        << "column 11: IC=Ab-Ag dynamics: k_off\n"
        << "column 12: Threshold of Ag concentration for CC binding in Mol\n"
        << "column 13: CB2CC differentiation time in hours\n"
        << "column 14: Differentiation signal production in Quanta per FDC and hour\n"
        // and now results:
        << "column 15: First time of population below 100 cells (hours)\n"
        << "column 16: Number of cells at the end\n"
        << "column 17: Fraction of high affinity cells at the end (average CB+CC)\n"
        << "column 18: Number of produced OUTPUT cells\n"
        << "column 19: Average binding probability of produced OUTPUT cells\n"
        << "column 20: End of dark zone (hours)\n"
        << "column 21: Ratio of CC2CB at t=288 hours\n"
        << "column 22: Ratio of output cells at day 12 to day 6\n"
        << "column 23: Standard deviation from measured GC kinetics\n";
     sel_ana_text.close();
 }

 void hyphasma2::write_sel_analyze(Parameter &par){
     // +++ OPTION: Use an additional analysis-file with special results
     // The file is designed for several runs stored in the same file
        ofstream sel_analyze(hyphasma2::prefix_files("sel_analyze.out"), std::ios_base::app);

          if (par.Value.FDClength >= par.Value.dx) {
              sel_analyze << 12 * int (par.Value.FDClength / par.Value.dx + 0.5)
                 * par.Value.FDCnumber;
           } else {
              sel_analyze << 4 * par.Value.FDCnumber;
           }
           sel_analyze << "   " << par.Value.FDCnumber << "   " << par.Value.FDClength << "   ";
           if (par.Value.ag_per_FDC < 0) {
              sel_analyze << "-1   ";
           } else {
              sel_analyze << par.Value.ag_per_FDC << "   ";
           }
           sel_analyze << par.Value.ag_saturation_FDC << "   " << par.Value.totalTC << "   ";
           if (par.Value.TC_CC_selection == 1) {
              sel_analyze << par.Value.TC_time << "   " << par.Value.TC_rescue_time << "   ";
           } else {
              sel_analyze << "0   -1   ";
           }
           if (par.Value.mk_ab < 0.) {
              sel_analyze << "0   ";
           } else {
              sel_analyze << par.Value.mk_ab << "   ";
           }
           sel_analyze << par.Value.ic_k_on << "   " << par.Value.ic_k_off << "   "
                       << par.Value.ag_threshold << "   "
                       << log(2.) / par.Value.tolight << "   " << par.Value.mksignal << "   ";


 }

 void hyphasma2::data_analysis(){
     if (c.outputfiles == 0) {
        // if at least one cell is left do this analysis
        if (c.CB_list.benutzt() > 0) {
           // deformation parameter:
           long nmoves = 0;
           hs.frac_average = 0.;
           for (long n = 0; n < c.CB_list.benutzt(); ++n) {
              nmoves += c.CB_list[n].n_moves;
              hs.frac_average += (c.CB_list[n].n_moves * c.CB_list[n].alpha_mean);
           }
           hs.frac_average /= double (nmoves);
           ana << "Average fraction of considered free neighbors for movement = "
                   << hs.frac_average
                   << "\n";
        }

        /*// If no CB is left add one artificially to allow access to the data
        if (c.CB_list.benutzt() == 0) {
           cellCB tmpcell;
           c.CB_list.add(tmpcell);
        }*/

        // OPTION: MOVE_ANALYSIS
        /*
          ana <<"Fraction of successful movement tries = "
          <<double(hs.n_move_done)/double(hs.n_try_move)
          <<"\n";
          ana <<"Fraction of movements on itself = "
          <<double(hs.n_move_self)/double(hs.n_try_move)
          <<"\n";
          ana <<"Fraction of movements on itself of back fragments = "
          <<double(hs.n_move_self_back)/double(hs.n_try_move)
          <<"\n";
          ana <<"Fraction of movement tries removed by previous movement = "
          <<double(hs.n_move_removed)/double(hs.n_try_move)
          <<"\n";
          ana <<"Fraction of movement tries that are forbidden = "
          <<double(hs.n_move_forbidden)/double(hs.n_try_move)
          <<"\n";
          ana <<"Fraction of movement tries that are forbidden for back fragments= "
          <<double(hs.n_move_forbidden_back)/double(hs.n_try_move)
          <<"\n";
          double factor=-1.0;
          if (hs.n_try_move!=0
          && hs.n_move_self!=0
          && hs.n_move_forbidden!=0) {
          double denom=double(hs.n_move_done)/double(hs.n_try_move)
          +(1.0-double(hs.n_move_self_back)/double(hs.n_move_self))
          double(hs.n_move_self)/double(hs.n_try_move)
         +(1.0-double(hs.n_move_forbidden_back)/double(hs.n_move_forbidden))
          double(hs.n_move_forbidden)/double(hs.n_try_move);
          if (denom>0.)
          factor=(1.0-double(hs.n_move_removed)/double(hs.n_try_move))/denom;
          }
          ana<<"Correction factor for diffusion at the end = ";
          if (factor==-1.0) ana<<"error\n";
          else ana<<factor<<"\n";*/
         // end of OPTION: MOVE_ANALYSIS

        ana << "Ratio of performed to aimed move per time step (average over "
                << c.CB_list[0].n_directed_moves
                << " moves) = " << c.CB_list[0].performed2aimed_move << "\n";
     }

     ana << "Number of CBs at the end = " << c.CB_end << "\n";
     c.CB_average /= 48.;
     ana << "Average number of CBs during the last 2 days = " << c.CB_average << " +- "
             << sqrt(c.CB_variance / 48. - c.CB_average * c.CB_average) << "\n";
     double _OUT_haffinity, _OUT_steepness, _CB_haffinity, _CC_haffinity;
     newshape.Call_getStatistics(_OUT_haffinity, _OUT_steepness, _CB_haffinity, _CC_haffinity);
     //getStatistics is a private function
     //newshape.getStatistics(_OUT_haffinity, _OUT_steepness, _CB_haffinity, _CC_haffinity);
     ana << "High affinity of CBs at the end = " << 100. * _CB_haffinity << "%\n";
     ana << "High affinity of CCs at the end = " << 100. * _CC_haffinity << "%\n";
     ana << "High affinity of OUTs at the end = " << 100. * _OUT_haffinity << "%\n";
     ana << "Total number of OUTs at the end = " << newshape.get_sum_cell(sout) << "\n";
     ana << "Steepness of OUT at times 288h/144h = " << _OUT_steepness << "\n";
     ana << "Effectivity index out/mid_CB = " << double (newshape.get_sum_cell(sout))
        / c.CB_average << "\n";
     ana << "End of dark zone at " << c.t_dark_end << " h\n";
     if (hs.fdc_obj.ic_calculations > 0) {
        double ic_errs
           = double (hs.fdc_obj.ab_sign_errors + hs.fdc_obj.ag_sign_errors
                     + hs.fdc_obj.ic_sign_errors)
             / double (hs.fdc_obj.ic_calculations);
        ana << "Fraction of errors in mk_immune_complex(..): " << ic_errs << "\n";
        ana << "  ab: " << double (hs.fdc_obj.ab_sign_errors)
           / double (hs.fdc_obj.ic_calculations)
                << ", ag: " << double (hs.fdc_obj.ag_sign_errors)
           / double (hs.fdc_obj.ic_calculations)
                << ", ic: " << double (hs.fdc_obj.ic_sign_errors)
           / double (hs.fdc_obj.ic_calculations) << ".\n";
        if (ic_errs > 0.05) {
           cout << "Note large number of errors in cellFDC::mk_immune_complex(...)!\n";
        }
     }
     ana << "... end of data analysis.\n";
     ana.close();

     // +++ OPTION: results of optional analysis
     ofstream sel_analyze(hyphasma2::prefix_files("sel_analyze.out"), std::ios_base::app);

     sel_analyze << c.t_1st_under100 << "   " << newshape.get_sum_cell(sCB)
        + newshape.get_sum_cell(sCC) << "   ";
     sel_analyze << (double (newshape.get_sum_cell(sCB)) * _CB_haffinity
                     + double (newshape.get_sum_cell(sCC)) * _CC_haffinity)
        / double (newshape.get_sum_cell(sCB) + newshape.get_sum_cell(sCC)) << "   "
                 << newshape.get_sum_cell(sout) << "   " << _OUT_haffinity << "   "
                 << c.t_dark_end << "   "
                 << c.CC2CBratio << "   " << _OUT_steepness << "   ";
     // if (c.n_chi_values>1) sel_analyze<<c.chi_sum/(double(c.n_chi_values));
     double sigma2 = c.kinetics.get_sigma2_GCvolume();
     if (sigma2 >= 0.) {
        sel_analyze << sigma2;
     } else {
        sel_analyze << "*";
     }
     sel_analyze << "\n";
     sel_analyze.close();
     // end OPTION.
 }


 void hyphasma2::update_newEFPC(double dt){
     // this function is called only until the stoptime for adding newEFPCs
     nEFPC=nEFPC+(EFPC_Prate*dt);
 }

 void hyphasma2::update_EFPCapoptosis(double dt){
     double apop_pc;
     apop_pc=nEFPC*EFPC_Apoprate*dt;

     if(apop_pc>0.0){
         nEFPC=nEFPC-apop_pc;
     }

 }

 void hyphasma2::update_EFPCantibody(double dt){
     // antibody production by EFPC
     double efpc_ab;

     efpc_ab=(EFPC_abPrate*nEFPC*dt)-(EFPC_abDrate*EFPC_antibody*dt);
     EFPC_antibody=EFPC_antibody+efpc_ab;

 }

 void hyphasma2::add_EFPCab_tobin(double dt){
    // add EFPC antibodies to bin 0
     // only for the 1st antigen
     ab_sum[0][0]=ab_sum[0][0]+EFPC_antibody;

 }

    int main(int argn, char * * argument) {

        // to run the simulation
        // hyphasma Parameters/x1 Par2files/x2

        hyphasma2::PCsUndergoApoptosis=false;
        hyphasma2::includeEFPCs=false;

         // read the new parameter file
          string parfile2=argument[2];
          ifstream par2=ifstream(parfile2.c_str());

          string Pnames[3]={"Number of GCs","Starting time for each GC(hours)","N for each GC"};

          int tt_gcs=0; 

          while(!par2.eof()){
             string name, fromfile;
             name=Pnames[0];

             getline(par2,fromfile);
            // cerr<<endl<<fromfile;

             if(!fromfile.compare(name)){

                   // read the total number of GCs to be simulated
                   par2>>tt_gcs;
                   hyphasma2::ngcs_tt=tt_gcs;
                 }
          }

          cerr<<endl<<"Total number of simulated gcs="<<tt_gcs;

          //initiation time in hours
          double gc_time_av[tt_gcs],gc_time[tt_gcs],gc_time_sd;
          gc_time_sd=24.0;

          // value of N = Number of GCs represented by each simulated gc
          int Ngcs[tt_gcs];
          ofstream nGCs,ini_time;
          nGCs.open("nGCs.out");
          ini_time.open("ini_time.out");

          for(int i=1;i<3;i++){
            ifstream p2=ifstream(parfile2.c_str());

              while(!p2.eof()){

                 string name, fromfile;
                 name=Pnames[i];

                 getline(p2,fromfile);
                 if(!fromfile.compare(name)){


                     if(i==1){
                         cerr<<endl<<"starting times";
                         for(int j=0;j<tt_gcs;j++){
                             p2>>gc_time_av[j];
                             cerr<<endl<<gc_time_av[j];
                         }
                     }

                     if(i==2){
                         cerr<<endl<<"value of N";
                         for(int j=0;j<tt_gcs;j++){
                             p2>>Ngcs[j];
                             cerr<<endl<<Ngcs[j];
                             nGCs<<Ngcs[j]<<endl;
                         }
                     }

                 }

              }
          }
          // end of reading the new parameter file
         nGCs.close();

          Parameter par[tt_gcs];
          for(int i=0;i<tt_gcs;i++)
          {
              // read parameter file: same parameter file read for each GC
              par[i].wahl(argument[1],true,true);

          }

          /*cerr<<"<GClifetime paper> settings for Figure S3"<<endl;

         // founder distance with respect to specific antigen
          par[0].Value.founder_closeto=0;
          par[0].Value.min_seeder_dist=1;
          par[0].Value.max_seeder_dist=3;

          par[1].Value.founder_closeto=0;
          par[1].Value.min_seeder_dist=1;
          par[1].Value.max_seeder_dist=3;

          par[2].Value.founder_closeto=0;
          par[2].Value.min_seeder_dist=1;
          par[2].Value.max_seeder_dist=3;

          par[3].Value.founder_closeto=1;
          par[3].Value.min_seeder_dist=1;
          par[3].Value.max_seeder_dist=3;

          par[4].Value.founder_closeto=2;
          par[4].Value.min_seeder_dist=1;
          par[4].Value.max_seeder_dist=3;

          par[5].Value.founder_closeto=2;
          par[5].Value.min_seeder_dist=1;
          par[5].Value.max_seeder_dist=3;

          par[6].Value.founder_closeto=2;
          par[6].Value.min_seeder_dist=-1;
          par[6].Value.max_seeder_dist=-1;

          par[7].Value.founder_closeto=2;
          par[7].Value.min_seeder_dist=-1;
          par[7].Value.max_seeder_dist=-1;*/


          /*cerr<<"<GClifetime paper> settings for Figure 5"<<endl;

          // founder distance with respect to specific antigen

          par[0].Value.founder_closeto=2;
          par[0].Value.min_seeder_dist=-1;
          par[0].Value.max_seeder_dist=-1;

          par[1].Value.founder_closeto=2;
          par[1].Value.min_seeder_dist=-1;
          par[1].Value.max_seeder_dist=-1;

          par[2].Value.founder_closeto=2;
          par[2].Value.min_seeder_dist=1;
          par[2].Value.max_seeder_dist=3;

          par[3].Value.founder_closeto=1;
          par[3].Value.min_seeder_dist=1;
          par[3].Value.max_seeder_dist=3;

          par[4].Value.founder_closeto=0;
          par[4].Value.min_seeder_dist=1;
          par[4].Value.max_seeder_dist=3;

          par[5].Value.founder_closeto=0;
          par[5].Value.min_seeder_dist=1;
          par[5].Value.max_seeder_dist=3;

          par[6].Value.founder_closeto=0;
          par[6].Value.min_seeder_dist=1;
          par[6].Value.max_seeder_dist=3;

          par[7].Value.founder_closeto=0;
          par[7].Value.min_seeder_dist=1;
          par[7].Value.max_seeder_dist=3;*/


          hyphasma2::num_ags=par[0].Value.APeakNumber;//number of antigens from parfile

          if (par[0].Value.ini_random == -1) {
             srand(time(NULL));
          }


          for(int i=0;i<tt_gcs;i++){
              gc_time[i]=get_positive_sample_from_normal(gc_time_av[i],gc_time_sd);
               cerr<<gc_time[i]<<endl;
               ini_time<<gc_time[i]<<endl;

               // initiation time dependent or randomly sampled ag concentrations

               /* cerr << " <GClifetime paper> settings for Figure 6 " << endl;
               par[i].Value.ag_per_FDC=get_positive_sample_from_normal(1300,1100);
               */

               /* cerr << " <GClifetime paper> settings for Figure 4 " << endl;
               par[i].Value.ag_per_FDC=20000*exp(-0.026*gc_time[i]);
               */

               /* cerr << " <GClifetime paper> settings for Figure 7 " << endl;
               par[i].Value.ag_per_FDC=2000*exp(-0.01*gc_time[i]);
               */

          }

          ini_time.close();

         // exit(1);

          double time_i;
          double time_f;
          double time_p;

          time_i=par[0].Value.tmin; // Day 0
          time_p=0;// present time (between time_i and time_f)
          time_f=par[0].Value.tmax; //final time in hours
          double deltat=par[0].Value.deltat;

          // number of timesteps from time_i to time_f
          long int nmax2 = long ((time_f - time_i) / deltat + 0.5);

          double t_gap = 24.;

          short int write_x = 0;
          short write_t = 0, write_day = 0, write_5days = 0;

          // write zone-files each day
          short int write_z = 0;

          long write_int = long (nmax2 / (time_f - time_i));
          if (write_int == 0) {
             write_int = 1;
          }

          long write12_int = long (12 * nmax2 / (time_f - time_i));
          long write24_int = long (24 * nmax2 / (time_f - time_i));
          long write120_int =long (120 * nmax2 / (time_f - time_i));
          short int inject = 0;

          // Initialisierung des Dateinamens fuer den Output:
          suffix tnr = "0000";
          if (time_i > 0) {
             for (int push = 0; push < int (double (time_i) / t_gap + 0.5); push++) {
                addchar(tnr);
             }
          }

          //open files

          hyphasma2** h=new hyphasma2*[tt_gcs];

          for(int i=0;i<tt_gcs;i++)
          {
              hyphasma2::GCref=i; // set GC ref number

              ofstream gcn;
              gcn.open(hyphasma2::prefix_files("GC_new.out"));
              gcn.close();
              ofstream pb;
              pb.open(hyphasma2::prefix_files("plasma.out"));
              pb.close();
              ofstream pb_aff;
              pb_aff.open(hyphasma2::prefix_files("npc.out"));
              pb_aff.close();

              ofstream oall_aff;
              oall_aff.open(hyphasma2::prefix_files("outputall_aff.out"));
              oall_aff.close();

              ofstream self_ab;

              self_ab.open(hyphasma2::prefix_files("self_ab.out"));
              self_ab.close();
              ofstream self_ic_pergc;
              self_ic_pergc.open(hyphasma2::prefix_files("self_ic_pergc.out"));
              self_ic_pergc.close();
              ofstream self_ab_pergc;
              self_ab_pergc.open(hyphasma2::prefix_files("self_ab_pergc.out"));
              self_ab_pergc.close();

              ofstream IP,relIP,check,ab_tc1,self_ic;
              IP.open(hyphasma2::prefix_files("IP.out"));
              IP.close();

              relIP.open(hyphasma2::prefix_files("relIP.out"));
              relIP.close();

              check.open(hyphasma2::prefix_files("check.out"));
              check.close();

              ab_tc1.open(hyphasma2::prefix_files("ab_tc1.out"));
              ab_tc1.close();

              self_ic.open(hyphasma2::prefix_files("self_ic.out"));
              self_ic.close();

              par[i].Value.N_GC=Ngcs[i]; // overwrite N values in parameter class

              h[i]=new hyphasma2(par[i]);

               // write sel_analyze.out file
              h[i]->write_sel_analyze(par[i]);

              ofstream CBaff,CCaff,OUTaff,PCaff;
              stringstream fileCB,fileCC,fileOUT,filePC;
              fileCB << "CBaff_GC" << i << ".out";
              fileCC << "CCaff_GC" << i << ".out";
              fileOUT << "OUTaff_GC" << i << ".out";
              filePC << "PCaff_GC" << i << ".out";

              CBaff.open(fileCB.str());
              CCaff.open(fileCC.str());
              OUTaff.open(fileOUT.str());
              PCaff.open(filePC.str());

              CBaff.close();
              CCaff.close();
              OUTaff.close();
              PCaff.close();
          }

              // The stored values are explained in a separate text-file:
              hyphasma2::write_sel_readme();

              // overall readout for all GCs
              ofstream tt_ab,tt_gcvol,tt_ip,tt_PCaff,tt_CBaff,tt_CCaff,tt_OUTaff,tt_pcs,tt_cbs,tt_ccs,tt_out,tt_ngcs,
              GCvolume,EFPC_count,EFPC_antibody,Plasmacells,Antibodies,CBs,CCs,OUT,DZcount,LZcount;

              string filenames[]={"PCaffinity","CBaffinity","CCaffinity","OUTaffinity","Freeag"};


              for(int i=0;i<hyphasma2::num_ags;i++){

                  for(int j=0;j<5;j++){
                      string fname;
                      stringstream filename;

                      filename << filenames[j] << "_Ag" << i <<".out";
                      fname=filename.str();
                      ofstream openfile;
                      openfile.open(fname);
                      openfile.close();
                  }

              }

              DZcount.open("DZcount.out");
              LZcount.open("LZcount.out");
              tt_ab.open("tt_ab.out");
              tt_gcvol.open("tt_gcvol.out");
              tt_ip.open("tt_ip.out");
              tt_PCaff.open("tt_PCaff.out");
              tt_CBaff.open("tt_CBaff.out");
              tt_CCaff.open("tt_CCaff.out");
              tt_OUTaff.open("tt_OUTaff.out");
              tt_pcs.open("tt_pcs.out");
              tt_cbs.open("tt_cbs.out");
              tt_ccs.open("tt_ccs.out");
              tt_out.open("tt_out.out");

              tt_ngcs.open("tt_ngcs.out");
              Plasmacells.open("Plasmacells.out");
              CBs.open("CBs.out");
              CCs.open("CCs.out");
              OUT.open("OUT.out");

              Antibodies.open("Antibodies.out");
              GCvolume.open("GCvolume.out");
              EFPC_count.open("EFPC_count.out");
              EFPC_antibody.open("EFPC_antibody.out");

              DZcount.close();
              LZcount.close();
              tt_ab.close();
              tt_gcvol.close();
              tt_ip.close();
              tt_PCaff.close();
              tt_CBaff.close();
              tt_CCaff.close();
              tt_OUTaff.close();
              tt_pcs.close();
              tt_cbs.close();
              tt_ccs.close();
              tt_out.close();
              tt_ngcs.close();
              Plasmacells.close();
              CBs.close();
              CCs.close();
              OUT.close();

              Antibodies.close();
              GCvolume.close();
              EFPC_count.close();
              EFPC_antibody.close();


                long write_beta = long (par[0].betaValue.dt_output / (3600. * par[0].Value.deltat) + 0.5);
                long beta_count = 0;
                if ((par[0].Value.show_mode == islet) && (cellbeta::LOCAL_FILES == false) && (write_beta == 0)) {
                   write_beta = 1;
                   cout << "WARNING: beta-output-dt is too small for the external timestep!\n";
                }

             // calculated by using 1 par file
                long switch_toGCphase = -1;

                long switch_toDifferentiation;
                if (cellCB::SMOOTH_DIFFERENTIATION) {
                   switch_toDifferentiation = 1;
                } else {
                   switch_toDifferentiation
                      = long ((par[0].Value.Start_Differentiation - time_i) / deltat + 0.5) + 1;
                }
                long switch_toOutput;
                if (cellCC::SMOOTH_DIFFERENTIATION) {
                   switch_toOutput = 1;
                } else {
                   switch_toOutput = long ((par[0].Value.StartOutput - time_i) / deltat + 0.5) + 1;
                }

                suffix tnz = "0000";
                if (par[0].Value.tmin > 0) {
                   for (int push = 0; push < int (double (par[0].Value.tmin) / t_gap + 0.5); push++) {
                      addchar(tnz);
                   }
                }


                //array of n., each index for one gc
                long int *n;

                // length of array n=total number of simulated GC
                n= new long int[tt_gcs];

                // Initialize all to 0. 0 means the GC has not been initialized yet
                for (int i=0;i<tt_gcs;i++){
                    n[i]=0;
                }

                 //Number of gcs (w/o considering N) active at a timepoint
                int n_gcs=0;


                for (int nloop = 1; (nloop <= nmax2) ; nloop++) {
                  
                  long write_tt_1h = ((nloop == (nloop / write_int) * write_int));
                  long write_tt_day = ((nloop == (nloop / write24_int) * write24_int));
                  long write_tt_5days = (nloop == (nloop / write120_int) * write120_int);
                  
                   // every timestep
                 n_gcs=0;
                 long *m;
                 time_p+=deltat;


                    for(int i1=0;i1<tt_gcs;i1++){
                        // in case it is not initialised
                        if(n[i1]==0){
                            if(gc_time[i1]<=time_p){
                                //then it needs to be initialised
                                n[i1]=1;
                                // does not need GCref number?
                                h[i1]->c.set_pars(par[0], h[i1]->l, h[i1]->hs);

                            }
                        }
                    }

                    for(int i2=0;i2<tt_gcs;i2++){
                        if(n[i2]>0){
                            n_gcs=n_gcs+1;
                        }
                    }

                 m=new long[tt_gcs];
                 for(int i=0;i<tt_gcs;i++){
                 m[i]=i;
                 }

                 random2_sequence(m,tt_gcs);// not necessary?

                    for(int j=0;j<tt_gcs;j++)
                    {

                        hyphasma2::GCref=m[j];// set GC ref number

                        if(n[m[j]]!=0){

                            write_x = (n[m[j]] == (n[m[j]] / par[m[j]].Value.ToFileStep) * par[m[j]].Value.ToFileStep);
                            write_t = ((n[m[j]] == (n[m[j]] / write_int) * write_int) || (n[m[j]] == nmax2));
                            inject = (n[m[j]] == (n[m[j]] / write12_int) * write12_int);
                            write_day = (n[m[j]] == (n[m[j]] / write24_int) * write24_int);
                            write_5days = (n[m[j]] == (n[m[j]] / write120_int) * write120_int);
                            write_z = ((n[m[j]] == (n[m[j]] / (int (t_gap) * write_int)) * int (t_gap) * write_int) || (n[m[j]] == nmax2));

                            h[m[j]]->c.time += deltat;


                            if ((n[m[j]] == switch_toGCphase)
                                ||   //	  n==switch_toMutation ||
                                (n[m[j]] == switch_toDifferentiation)
                                || (n[m[j]] == switch_toOutput))
                           {

                               h[m[j]]->c.set_pars(par[m[j]], h[m[j]]->l, h[m[j]]->hs);

                           }

                          hyphasma2::num_ags=h[m[j]]->Ab.ab_bins.size();//number of antigens from number of Ab bins

                        h[m[j]]->c.time_step(write_t, write_day, h[m[j]]->l, h[m[j]]->s, h[m[j]]->Ab, h[m[j]]->newshape, h[m[j]]->ana, par[m[j]], h[m[j]]->hs);

                        if (write_x == 1) {
                           addchar(tnr);
                           h[m[j]]->c.movie << "$1\"" << tnr << ".gif\" ";
                           // c.mkgifs<<"ppmtogif "<<"xy"<<tnr<<".ppm > xy"<<tnr<<".gif\n";
                           h[m[j]]->c.xfiles(tnr, h[m[j]]->l);
                           h[m[j]]->c.write_mutations(h[m[j]]->c.time);
                           h[m[j]]->newshape.write_gcbc_hamming(h[m[j]]->c.time);
                           h[m[j]]->newshape.write_gcbc_affinity(h[m[j]]->c.time);
                           h[m[j]]->s.write_files(tnr, false);
                           if (sigs::TEST_MODE) {
                              h[m[j]]->s.write_TEST(h[m[j]]->c.time);
                           }

                           if (par[m[j]].Value.safety_checks == 1) {
                              h[m[j]]->c.check_all(h[m[j]]->l, h[m[j]]->newshape);
                              h[m[j]]->Ab.check4consistency(h[m[j]]->newshape, h[m[j]]->ana);
                           }
                        }
                        if ((h[m[j]]->c.outputfiles < 2) && (write_z == 1) && (h[m[j]]->c.show_mode != islet)) {
                           addchar(tnz);
                           h[m[j]]->c.zone_files(tnz, h[m[j]]->l);
                        }
                        if (inject == 1) {
                           // every 12 hours
                           if (n[m[j]] >= switch_toDifferentiation - 1) {
                              h[m[j]]->c.inject_Ki67();
                           }
                           h[m[j]]->c.write_log_bc_aff(h[m[j]]->newshape);
                        }
                        if ((cellbeta::LOCAL_FILES == false) && (beta_count == write_beta)) {
                           h[m[j]]->c.show_BETA();
                           beta_count = 0;
                        }


                        }



                    }

                    hyphasma2::initialize_abs(); // clear the new antibody bins

                    // add antibodies from all GCs
                    for(int i=0;i<tt_gcs;i++)
                    {
                       hyphasma2::antibody_sum(&h[i]);
                    }

                    // add EFPC abs
                    if(hyphasma2::includeEFPCs==true) {
                        if(time_p<=hyphasma2::newEFPCstoptime){
                            hyphasma2::update_newEFPC(deltat);
                        }
                        hyphasma2::update_EFPCapoptosis(deltat);
                        hyphasma2::update_EFPCantibody(deltat);

                        hyphasma2::add_EFPCab_tobin(deltat);
                    }


                    if(write_tt_1h) //write_t
                    {
                        hyphasma2::gcvol_tot=0;
                        hyphasma2::pcs_tot=0;
                        hyphasma2::cbs_tot=0;
                        hyphasma2::ccs_tot=0;
                        hyphasma2::out_tot=0;
                        hyphasma2::outext_tot=0;
                        hyphasma2::ngcs=0;
                        hyphasma2::small=0;
                        hyphasma2::medium=0;
                        hyphasma2::large=0;
                        hyphasma2::DZ_count=0;
                        hyphasma2::LZ_count=0;

                        for(int i=0;i<tt_gcs;i++){

                            // to classify gcs
                            hyphasma2::gcvol_tot +=(h[i]->calculate_gcvol(&h[i],Ngcs[i]));
                            // cutoff for small, medium and large GCs
                            if((h[i]->calculate_gcvol(&h[i],Ngcs[i])/Ngcs[i])>100){
                                hyphasma2::ngcs +=Ngcs[i];

                                if(((h[i]->calculate_gcvol(&h[i],Ngcs[i]))/Ngcs[i])<=2000){
                                    hyphasma2::small +=Ngcs[i];
                                }
                                else if(((h[i]->calculate_gcvol(&h[i],Ngcs[i]))/Ngcs[i])>2000 && ((h[i]->calculate_gcvol(&h[i],Ngcs[i]))/Ngcs[i])< 5000){
                                    hyphasma2::medium +=Ngcs[i];
                                }
                                else{
                                    hyphasma2::large +=Ngcs[i];
                                }
                            }

                            hyphasma2::pcs_tot +=(h[i]->calculate_nCells(&h[i],Ngcs[i],soutextproduce));
                            hyphasma2::cbs_tot += (h[i]->calculate_nCells(&h[i],Ngcs[i],sCB));
                            hyphasma2::ccs_tot += (h[i]->calculate_nCells(&h[i],Ngcs[i],sCC));
                            //hyphasma2::out_tot += (h[i]->calculate_nCells(&h[i],Ngcs[i],sout));
                            hyphasma2::out_tot += (h[i]->count_output(&h[i],Ngcs[i],h[i]->l));
                           // hyphasma2::outext_tot += (h[i]->calculate_nCells(&h[i],Ngcs[i],soutext));

                            hyphasma2::DZ_count += (h[i]->calc_dzcount(&h[i],Ngcs[i],h[i]->l));
                            hyphasma2::LZ_count += (h[i]->calc_lzcount(&h[i],Ngcs[i],h[i]->l));

                        }

                        for(int g=0;g<hyphasma2::num_ags;g++){
                            hyphasma2::PCaff_tot[g]=0.0;
                            hyphasma2::CBaff_tot[g]=0.0;
                            hyphasma2::CCaff_tot[g]=0.0;
                            hyphasma2::OUTaff_tot[g]=0.0;

                            for(int i=0;i<tt_gcs;i++){
                                //  calculate_sumaff only returns the sum of affinities of all cells ofcelltype
                                hyphasma2::PCaff_tot[g] +=(h[i]->calculate_sumaff(&h[i],g,Ngcs[i],soutextproduce));
                                hyphasma2::CBaff_tot[g] +=(h[i]->calculate_sumaff(&h[i],g,Ngcs[i],sCB));
                                hyphasma2::CCaff_tot[g] +=(h[i]->calculate_sumaff(&h[i],g,Ngcs[i],sCC));
                                //hyphasma2::OUTaff_tot[g] +=(h[i]->calculate_sumaff(&h[i],g,Ngcs[i],sout));
                                hyphasma2::OUTaff_tot[g] +=(h[i]->calculate_outaff(&h[i],g,Ngcs[i],h[i]->l,h[i]->newshape)); 
                               // hyphasma2::OUTaff_tot[g] +=(h[i]->calculate_sumaff(&h[i],g,Ngcs[i],soutext));

                            }

                            if(hyphasma2::pcs_tot>0) hyphasma2::PCaff_tot[g]/=hyphasma2::pcs_tot;
                            else hyphasma2::PCaff_tot[g]=0;


                            if(hyphasma2::cbs_tot>0) hyphasma2::CBaff_tot[g]/=hyphasma2::cbs_tot;
                            else hyphasma2::CBaff_tot[g]=0;


                            if(hyphasma2::ccs_tot>0) hyphasma2::CCaff_tot[g]/=hyphasma2::ccs_tot;
                            else hyphasma2::CCaff_tot[g]=0;


                            if((hyphasma2::out_tot)>0) hyphasma2::OUTaff_tot[g]/=(hyphasma2::out_tot);
                            else hyphasma2::OUTaff_tot[g]=0;
                        }

                    }

                    // write readouts of every gc, efpcs
                    if(write_tt_1h) //write_t
                    {

                    GCvolume.open("GCvolume.out",ofstream::app);
                    Plasmacells.open("Plasmacells.out",ofstream::app);
                    CBs.open("CBs.out",ofstream::app);
                    CCs.open("CCs.out",ofstream::app);
                    OUT.open("OUT.out",ofstream::app);
                    Antibodies.open("Antibodies.out",ofstream::app);
                    EFPC_count.open("EFPC_count.out",ofstream::app);
                    EFPC_antibody.open("EFPC_antibody.out",ofstream::app);

                    GCvolume << time_p/24 << " ";
                    Plasmacells << time_p/24 << " ";
                    CBs << time_p/24 << " ";
                    CCs << time_p/24 << " ";
                    OUT << time_p/24 << " ";

                    Antibodies << time_p/24 << " ";
                    EFPC_count<< time_p/24 << " ";
                    EFPC_antibody<< time_p/24 <<" ";

                    for(int i=0;i<tt_gcs;i++){

                         GCvolume<< h[i]->calculate_gcvol(&h[i],1) << " ";
                         Plasmacells<< h[i]->calculate_nCells(&h[i],1,soutextproduce) << " ";
                         CBs<< h[i]->calculate_nCells(&h[i],1,sCB) << " ";
                         CCs<< h[i]->calculate_nCells(&h[i],1,sCC) << " ";
                         //OUT<< (h[i]->calculate_nCells(&h[i],1,sout) << " ";
                         OUT<< (h[i]->count_output(&h[i],1,h[i]->l)) << " ";

                         double pc=h[i]->calculate_nCells(&h[i],1,soutextproduce);
                         Antibodies<< h[i]->antibodysum_eachgc(&h[i],Ngcs[i])*1e+15 << " ";

                    }

                    EFPC_count << hyphasma2::nEFPC;
                    EFPC_antibody << hyphasma2::EFPC_antibody*1e+15;

                    GCvolume << endl;
                    Plasmacells << endl;
                    CBs << endl;
                    CCs << endl;
                    OUT << endl;

                    Antibodies << endl;

                    EFPC_count << endl;
                    EFPC_antibody << endl;

                    //write PCaffinity_AgX, CBaffinity_AgX, CCaffinity_AgX and OUTaffinity_AgX
                    for(int i1=0;i1<hyphasma2::num_ags;i1++ ){
                        ofstream PCaffinity,CCaffinity,CBaffinity,OUTaffinity,Freeag;
                        string fnamePC,fnameCC,fnameCB,fnameOUT,fnameAg;
                        stringstream filenamePC,filenameCC,filenameCB,filenameOUT,filenameAg;

                        filenamePC<< "PCaffinity_Ag" << i1 << ".out";
                        filenameCB<< "CBaffinity_Ag" << i1 << ".out";
                        filenameCC<< "CCaffinity_Ag" << i1 << ".out";
                        filenameOUT<< "OUTaffinity_Ag" << i1 << ".out";
                        filenameAg<< "Freeag_Ag" << i1 << ".out";

                        fnamePC=filenamePC.str();
                        fnameCB=filenameCB.str();
                        fnameCC=filenameCC.str();
                        fnameOUT=filenameOUT.str();
                        fnameAg=filenameAg.str();

                        PCaffinity.open(fnamePC,ofstream::app);
                        CBaffinity.open(fnameCB,ofstream::app);
                        CCaffinity.open(fnameCC,ofstream::app);
                        OUTaffinity.open(fnameOUT,ofstream::app);
                        Freeag.open(fnameAg,ofstream::app);

                        PCaffinity << time_p/24 << " ";
                        CBaffinity << time_p/24 << " ";
                        CCaffinity << time_p/24 << " ";
                        OUTaffinity << time_p/24 << " ";
                        Freeag << time_p/24 << " ";


                        for(int i=0;i<tt_gcs;i++){

                            double pcaff,pc,cb,cbaff,cc,ccaff,out,outaff;

                            pc=h[i]->calculate_nCells(&h[i],1,soutextproduce);
                            cb=h[i]->calculate_nCells(&h[i],1,sCB);
                            cc=h[i]->calculate_nCells(&h[i],1,sCC);
                            //out=(h[i]->calculate_nCells(&h[i],1,sout);
                            out=(h[i]->count_output(&h[i],1,h[i]->l)); 


                            if(pc>0)
                            {
                                pcaff=(h[i]->calculate_sumaff(&h[i],i1,1,soutextproduce))/(pc);
                            }else{
                                pcaff=0;
                            }

                            if(cb>0)
                            {
                                cbaff=(h[i]->calculate_sumaff(&h[i],i1,1,sCB))/(cb);
                            }else{
                                cbaff=0;
                            }

                            if(cc>0)
                            {
                                ccaff=(h[i]->calculate_sumaff(&h[i],i1,1,sCC))/(cc);
                            }else{
                                ccaff=0;
                            }

                            if(out>0)
                            {
                                outaff=(h[i]->calculate_outaff(&h[i],i1,1,h[i]->l,h[i]->newshape))/(out);
                            }else{
                                outaff=0;
                            }


                            PCaffinity << pcaff << " ";
                            CBaffinity << cbaff << " ";
                            CCaffinity << ccaff << " ";
                            OUTaffinity << outaff << " ";
                            Freeag << h[i]->calculate_freeag(&h[i],i1) << " " ;

                        }
                        PCaffinity << endl;
                        CBaffinity << endl;
                        CCaffinity << endl;
                        OUTaffinity << endl;
                        Freeag << endl;

                        PCaffinity.close();
                        CBaffinity.close();
                        CCaffinity.close();
                        OUTaffinity.close();
                        Freeag.close();

                    }

                    }

                    GCvolume.close();
                    Plasmacells.close();
                    CBs.close();
                    CCs.close();
                    OUT.close();

                 //   PCaffinity.close();
                    Antibodies.close();
                    EFPC_count.close();
                    EFPC_antibody.close();

                    if(write_tt_1h) //write_t
                    {  //writing readouts in file

                        ofstream tt_gcvol,tt_ip,tt_PCaff,tt_CBaff,tt_CCaff,tt_OUTaff,tt_pcs,tt_ab, tt_ngcs,
                                tt_cbs,tt_ccs,tt_out,DZcount,LZcount;

                        DZcount.open("DZcount.out",ofstream::app);
                        LZcount.open("LZcount.out",ofstream::app);

                        tt_ab.open("tt_ab.out",ofstream::app);
                        tt_ip.open("tt_ip.out",ofstream::app);
                        tt_pcs.open("tt_pcs.out",ofstream::app);
                        tt_cbs.open("tt_cbs.out",ofstream::app);
                        tt_ccs.open("tt_ccs.out",ofstream::app);
                        tt_out.open("tt_out.out",ofstream::app);

                        tt_PCaff.open("tt_PCaff.out",ofstream::app);
                        tt_CBaff.open("tt_CBaff.out",ofstream::app);
                        tt_CCaff.open("tt_CCaff.out",ofstream::app);
                        tt_OUTaff.open("tt_OUTaff.out",ofstream::app);
                        tt_ngcs.open("tt_ngcs.out",ofstream::app);

                        tt_ngcs<<time_p/24<<" "<<hyphasma2::ngcs<<" "<<hyphasma2::small<<" "<<hyphasma2::medium<<" "<<hyphasma2::large<<endl;

                        DZcount << time_p/24 << " " << hyphasma2::DZ_count << endl;
                        LZcount << time_p/24 << " " << hyphasma2::LZ_count << endl;

                        tt_ab<<time_p/24<<" ";

                        for(int c1=0;c1<hyphasma2::num_ags;c1++){
                            tt_ab<<hyphasma2::add_absum(c1)*1e+15<<" ";
                        }

                        tt_ab<<endl;


                        tt_gcvol.open("tt_gcvol.out",ofstream::app);

                        tt_gcvol<<time_p/24<<" "<<hyphasma2::gcvol_tot<<endl;
                        tt_ip<<time_p/24<<" ";
                        for(int g=0;g<hyphasma2::num_ags;g++){
                            tt_ip<<hyphasma2::calculate_IP(g)<<" ";
                        }
                        tt_ip<<endl;

                        tt_PCaff<<time_p/24<<" ";
                        tt_CBaff<<time_p/24<<" ";
                        tt_CCaff<<time_p/24<<" ";
                        tt_OUTaff<<time_p/24<<" ";

                        for(int i=0;i<hyphasma2::num_ags;i++){
                            tt_PCaff<<hyphasma2::PCaff_tot[i]<<" ";
                            tt_CBaff<<hyphasma2::CBaff_tot[i]<<" ";
                            tt_CCaff<<hyphasma2::CCaff_tot[i]<<" ";
                            tt_OUTaff<<hyphasma2::OUTaff_tot[i]<<" ";
                        }

                        tt_PCaff<<endl;
                        tt_CBaff<<endl;
                        tt_CCaff<<endl;
                        tt_OUTaff<<endl;

                        tt_pcs<<time_p/24<<" "<<hyphasma2::pcs_tot<<endl;
                        tt_cbs<<time_p/24<<" "<<hyphasma2::cbs_tot<<endl;
                        tt_ccs<<time_p/24<<" "<<hyphasma2::ccs_tot<<endl;
                        tt_out<<time_p/24<<" "<<hyphasma2::out_tot<<endl;

                        tt_gcvol.close();
                        tt_ip.close();
                        tt_PCaff.close();
                        tt_pcs.close();
                        tt_cbs.close();
                        tt_ccs.close();
                        tt_out.close();
                        tt_ab.close();
                        DZcount.close();
                        LZcount.close();

                    }

                    // increment n for already initialized gcs
                    for(int i3=0;i3<tt_gcs;i3++){
                        if(n[i3]!=0){
                            n[i3]+=1;
                        }
                    }

                    if(write_tt_5days==1){ //write_5days
                        for(int i=0;i<tt_gcs;i++){
                            if(n[i]!=0){
                                h[i]->calculate_cellaff_allAG(&h[i],i,time_p/24,h[i]->l,h[i]->newshape);
                            }
                        }
                    }

                }
                // end of all timesteps


                for(int i=0;i<tt_gcs;i++){
                    hyphasma2::GCref=i; // again set GCref number

                    // write the Plasma cell affinity file for all GCs

                    h[i]->calculate_PCaff_allAG(&h[i],i);

                    h[i]->c.write_final(h[i]->l.dx, h[i]->hs);
                    h[i]->c.write_cell_specific(h[i]->hs);
                    cout << "Close movie-files ... \n";
                    h[i]->c.movie << "> $1\".gif\"\n";
                    h[i]->c.movie.close();
                    // c.mkgifs.close();
                    // cout<<"done.\n";
                    h[i]->c.check_all(h[i]->l, h[i]->newshape);
                    h[i]->Ab.check4consistency(h[i]->newshape, h[i]->ana);

                    h[i]->ana << " ... end of calculation.\n\n";

                    h[i]->ana << "Data analysis:\n";

                    h[i]->data_analysis();
                    // +++++++++++++++ OPTION ++++++++++++++++++++++++++++++++++++++++++++++
                    // Write all signals to a file
                    bool writeallsignals = true;
                    if (writeallsignals) {
                       cout
                          <<
                          "Write signals to signal.out that can be used as input if copied to signal.sig.\n";
                       h[i]->s.write_all_signals(cell::chemo_max, cell::chemo_steep, cell::chemo_half);
                    } else {
                       cout << "Do not write signals to signal.out that could be used as input signal.sig.\n";
                    }
                    // +++++++++++ end OPTION ++++++++++++++++++++++++++++++++++++++++++++++

                    h[i]->c.close_files(); // close cellman files
                    h[i]->newshape.close_files(); // close affinityspace files

                    h[i]->c.trackdata.Write_files();
                    h[i]->c.trackmutations.save_to_file();
                    // activate this for brainbow analysis
                    //c.trackmutations.analysis();
                    // The following was used for the analysis in Cell Reports 2012 Fig. S1
                    h[i]->c.show_BrdU();
                    h[i]->c.show_ag_loaded();
                    h[i]->c.show_cell_in_BrdU();
                    h[i]->c.show_ag_BrdU(false);  // FACS: ag versus BrdU

                    if ((h[i]->c.show_mode == islet) && (cellbeta::LOCAL_FILES == false)) {
                       h[i]->c.mk_single_beta_files(h[i]->l);
                    }
                }

                return 0;


    }




