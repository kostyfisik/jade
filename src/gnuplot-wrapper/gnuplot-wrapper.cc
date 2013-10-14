///
/// @file   gnuplot-wrapper.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @date   Fri Sep 27 12:44:12 2013
/// @copyright 2013 Ladutenko Konstantin
///
/// gnuplot-wrapper is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// gnuplot-wrapper is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with gnuplot-wrapper.  If not, see <http://www.gnu.org/licenses/>.
///
///
/// @brief Wrapper class around gnuplot for ease of use. It produces
/// data file, plot description for gnuplot and shell script to call
/// gnuplot.
///
//        gnuplot-wrapper.h
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <vector>
#include <string>
#include "gnuplot-wrapper.h"
namespace gnuplot{  
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void GnuplotWrapper::PrintPlotFile() {
    FILE *fp;
    std::string fname = plot_name_ + ".plt";
    fp = fopen(fname.c_str(), "w");
    std::string plt_format = "# " + plot_name_ + "\n" +
      //"set terminal png nocrop size " +
      "set terminal svg size " +
      std::to_string(plot_size_x_) + ","
      + std::to_string(plot_size_y_) + " dynamic enhanced\n" +
      + "set object 1 rect from screen 0, 0, 0 to screen 1, 1, 0 behind\n" +
      + "set object 1 rect fc  rgb \"white\"  fillstyle solid 1.0\n"
      //"set output \"" + plot_name_ + ".png\"\n"+
      "set output \"" + plot_name_ + ".svg\"\n"+
      "set xlabel \"" + x_label_name_ + "\"\n"+
      "set ylabel \"" + y_label_name_ + "\"\n" +
      "set key below\n";
    if (x_range_[0] != x_range_[1])
      plt_format += "set xrange [" + std::to_string(x_range_[0])
        + ":" + std::to_string(x_range_[1]) + "]\n";
      // #set yrange [-0.1:0.1]
      // #set yrange [-1:1]
      // #set linestyle 1 lt 2 lw 3
    plt_format += "set title \"" + plot_name_ + "\"\n"
      +"plot \\\n";
    for (int i = 1; i < column_names_.size(); ++i)
      plt_format += "\"" + plot_name_ + ".dat\" using 1:"
        + std::to_string(i+1) + " title '" +
        column_names_[i] + "' " +
        plot_draw_style_ +",\\\n";
    plt_format.erase (plt_format.end()-3, plt_format.end());
    plt_format += "\n";
    fprintf(fp, "%s", plt_format.c_str());
    fclose(fp);
  }  // end of void GnuplotWrapper::PrintPlotFile()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void GnuplotWrapper::PrintShellFile() {
    FILE *fp;
    std::string fname = "run-gnuplot-" + plot_name_ + ".sh";
    fp = fopen(fname.c_str(), "w");
    std::string shell = std::string("#!/bin/bash\n")
      + "gnuplot " + plot_name_ + ".plt\n";
    fprintf(fp, "%s", shell.c_str());
    fclose(fp);
  }  // end of void GnuplotWrapper::PrintShellFile()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void GnuplotWrapper::PrintDataFile() {
    FILE *fp;
    std::string fname = plot_name_ + ".dat";
    fp = fopen(fname.c_str(), "w");
    std::string header_line = "#";
    for (auto name : column_names_) header_line += name + "\t";
    header_line.push_back('\n');
    fprintf(fp, "%s", header_line.c_str());
    header_line = "# "+ plot_name_ + "\n";
    fprintf(fp, "%s", header_line.c_str());    
    for (auto row : data_) {
      for (auto cell : row) fprintf(fp, "%.19g\t", cell);
      fprintf(fp, "\n");
    }  // end of for each row
    fclose(fp);
  }  // end of void GnuplotWrapper::PrintDataFile()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void GnuplotWrapper::MakeOutput() {
    int number_or_columns = column_names_.size();
    for (auto row : data_)
      if (number_or_columns != row.size())
        throw std::invalid_argument("All rows should be the same size!");
    PrintDataFile();
    PrintPlotFile();
    PrintShellFile();
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void GnuplotWrapper::AddMultiPoint(std::vector<double> point) {
    data_.push_back(point);
  }  // end of void GnuplotWrapper::AddMultiPoint(std::vector<double> point)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  ///GnuplotWrapper::
}  // end of namespace gnuplot
