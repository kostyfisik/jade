/**
 * @file   read-spectra.cc
 * @author Konstantin Ladutenko <kostyfisik at gmail (.) com>
 * @date   Wed Mar 11 11:51:26 2015
 * 
 * @copyright 2015 Konstantin Ladutenko
 *
 * @brief  Read complex spectra from file in format 'WL real imag'
 * 
 * read-spectra is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * read-spectra is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with read-spectra.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#include <algorithm>
#include <complex>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <iostream>
#include <vector>
#include "read-spectra.h"
namespace read_spectra {
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void ReadSpectra::ReadFromFile(std::string fname) {
    //std::cout<<"Reading file: "<< fname << std::endl;
    std::ifstream infile(fname.c_str());
    data_.clear();
    std::string line;
    while (std::getline(infile, line))
      {
	if (line.front() == '#') continue; //do not read comments
	if (line.find('#') != std::string::npos) 
	  throw std::invalid_argument("Error! Comments should be marked with # in the begining of the line!\n");
	std::istringstream iss(line);	
	double wl, re, im;
	if (!(iss >> wl >> re >> im)) throw std::invalid_argument("Error! Unexpected format of the line!\n");
	data_.push_back(std::vector<double>({wl,re,im}));
	//std::cout<<wl<<' '<<re<<' '<<im<<std::endl;
      }  // end of wile reading file 
    std::sort(data_.begin(), data_.end(),
	      [](const std::vector<double>& a, const std::vector<double>& b) {
		return a.front() < b.front();
	      });

  }  // end of void ReadSpectra::ReadFromFile(std::string fname)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void ReadSpectra::PrintData() {
    for (auto row : data_) {
      for (auto cell : row) printf("%g\t", cell);
      printf("\n");
    }  // end of for each row
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //

}  // end of namespace read_spectra

