#ifndef i_histo_monitor
#define i_histo_monitor
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

class histo_monitor {
 public:
  histo_monitor();
  ~histo_monitor();
  void initialize(double, int, double, double, std::string, std::string);
  void add(double value);
  void rm(double value);
  void update_day();
  void write2file();
  void show_all();
 private:
  //std::vector <std::vector<int> > histo_day;
  //std::vector<int> histo_total;
  int* histo_total;
  int** histo_day;

  unsigned int bin_max;
  double dvalue;
  unsigned int thisday, ndays;
  int get_bin(double value);
  double get_value(int bin);
  std::string filename, quantity_name;
};

#endif
