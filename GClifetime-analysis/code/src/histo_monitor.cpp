#include "histo_monitor.h"
#include "hyphasmah.h"



histo_monitor::histo_monitor() { };
histo_monitor::~histo_monitor() { }
void histo_monitor::initialize(double max_value, int set_bin_max, 
			       double tmin, double tmax, 
			       std::string set_filename,
			       std::string set_quantity_name) {
  //  const unsigned int bin_max_tmp = set_bin_max;
  ndays = int((tmax - tmin) / (24.0)) + 1;
  thisday = int((tmin / 24.0)) + 0;

  bin_max = set_bin_max;
  dvalue = max_value / double(bin_max);
  /*std::vector<int> vtottmp(bin_max + 1);
  histo_total = vtottmp;
  vtottmp.clear();
  std::vector <std::vector<int> > vtmp(ndays, std::vector<int>(bin_max + 1));
  histo_day = vtmp;
  vtmp.clear();*/
 //-----------------
  histo_total=new int[bin_max+1];

  histo_day = new int*[ndays];
  for (unsigned int f = 0; f < ndays; f++) {
     histo_day[f] = new int[bin_max+1];
  }
//------------------

  for (unsigned int day = 0; day < ndays; day++) {
    for (unsigned int nc = 0; nc < bin_max+1; nc++) {
      histo_day[day][nc] = 0;
    }
  }
  for (unsigned int nc = 0; nc < bin_max+1; nc++) {
    histo_total[nc] = 0;
  }
  filename = set_filename;
  quantity_name = set_quantity_name;
}
int histo_monitor::get_bin(double value) {
  int bin = int(value / dvalue + 0.5);
  if (bin > int(bin_max)) {
    std::cerr << "\n                               "
	      << "In int histo_monitor::get_bin(value="
	      << value << ", " << quantity_name << "): bin = " 
	      << bin << " larger than max = "
	      << bin_max << " was generated.\n";
    bin = bin_max;
  }
  return bin;
}
double histo_monitor::get_value(int bin) {
  double value = double(bin) * dvalue;
  return value;
}
void histo_monitor::add(double value) {

  ++histo_day[thisday][get_bin(value)];
}
void histo_monitor::rm(double value) {
  --histo_day[thisday][get_bin(value)];
}
void histo_monitor::show_all() {
  std::cout << "bin : value : total_frey : day 0 : ... : day last\n";
  for (unsigned ni = 0; ni <= bin_max; ++ni) {
    std::cout << ni << "  " << get_value(ni) << "    " << histo_total[ni] << "    ";
    for (unsigned ti = 0; ti < ndays; ++ti)
      { std::cout << histo_day[ti][ni] << "  "; }
    std::cout << "\n";
  }
}  
void histo_monitor::update_day() {
  // save this day in histo_total
  for (unsigned nc = 0; nc < bin_max+1; nc++)
    { histo_total[nc] += histo_day[thisday][nc]; }
  // switch to next day
  ++thisday;
}
void histo_monitor::write2file() {
  std::ofstream outfile(hyphasma2::prefix_files(filename.c_str()));
  outfile << "! bin : " << quantity_name << " : total freq : day 0 : ... : day last\n";
  for (unsigned ni = 0; ni <= bin_max; ++ni) {
    outfile << ni << "  " << get_value(ni) << "    " << histo_total[ni] << "    ";
    for (unsigned ti = 0; ti < ndays; ++ti) 
      { outfile << histo_day[ti][ni] << "  "; }
    outfile << "\n";
  }
  outfile.close();
}
