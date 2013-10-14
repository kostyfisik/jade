#ifndef SRC_GNUPLOT_WRAPPER_GNUPLOT_WRAPPER_H_
#define SRC_GNUPLOT_WRAPPER_GNUPLOT_WRAPPER_H_
///
/// @file   gnuplot-wrapper.h
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @date   Fri Sep 27 12:40:59 2013
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
/// @brief Wrapper class around gnuplot for ease of use. It produces
/// data file, plot description for gnuplot and shell script to call
/// gnuplot.
/// 
///
#include <vector>
#include <string>
namespace gnuplot {
  class GnuplotWrapper {
   public:
    void SetPlotName(std::string name) {plot_name_ = name;};
    void SetDrawStyle(std::string style) {plot_draw_style_ = style;};
    void SetXLabelName(std::string name) {x_label_name_ = name;};
    void SetYLabelName(std::string name) {y_label_name_ = name;};
    void SetXRange(std::vector<double> range) {x_range_ = range;};
    void AddMultiPoint(std::vector<double> point);
    void AddColumnName(std::string name) {column_names_.push_back(name);};
    void MakeOutput();
    
   private:
    void PrintDataFile();
    void PrintPlotFile();
    void PrintShellFile();
    int plot_size_x_ = 800;
    int plot_size_y_ = 600;
    /* int plot_size_x_ = 1366; */
    /* int plot_size_y_ = 768; */
    std::vector<double> x_range_ = {0, 0};
    std::string plot_name_ = "gnuplot";
    std::string x_label_name_ = "x";
    std::string y_label_name_ = "y";
    std::string plot_draw_style_ = "w l lw 2";
    std::vector< std::string> column_names_;
    std::vector< std::vector<double> > data_;
  };  // end of class GnuplotWRapper
}  // end of namespace gnuplot
#endif  // SRC_GNUPLOT_WRAPPER_GNUPLOT_WRAPPER_H_
