#ifndef SRC_COMMON_H_
#define SRC_COMMON_H_
///
/// @file   common.h
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @copyright 2013 Ladutenko Konstantin
/// @section LICENSE
/// This file is part of JADE++.
///
/// JADE++ is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// JADE++ is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with JADE++.  If not, see <http://www.gnu.org/licenses/>.
/// @date   Wed Aug 14 13:42:51 2013
///
/// @brief  Global variables, enumerations, macros, types, etc..
#include <stdint.h>  // Including stdint.h to use int64_t.
#include <sstream>
#include <string>
namespace jade {
  // ********************************************************************** //
  // **********************     Constants           *********************** //
  // ********************************************************************** //
  const int kOutput = 0;
  /// @brief Error codes
  ///
  /// Error codes used with onza.
  enum Errors {
    /// no error
    kDone = 0,
    /// Unspecified (pending to be described).
    kError,
    /// Dublicated key during reading config file.
    kErrorConfigFileDublicatedKey,
    /// Was not able to open config file during init of
    /// HaloExchangeProcess
    kErrorConfigFileWasNotAbleToOpen,
    /// After passing from main() start parameters to
    /// HaloExchangeProcess::Init(int argc, char *argv[])
    /// config file name was not found
    kErrorConfigFileNameWasNotDefined,
    /// Error while converting string to positive decimal.
    kErrorConfgiFileWrongDecimalInput,
    /// exchange buffer was initiated to be not contiguous with
    /// BasicSimulationCore::Init() or was detected to be not
    /// contiguous with HaloToExchange::Init().
    kErrorExchangeBufferIsNotContiguous,
    /// After splitting model to Cartesian grid current procces was
    /// not found in it.
    kErrorProcessNotInGrid,
    /// Trying to start running timer.
    kErrorProfilingTimerSecondStart,
    /// Some error during preseting the timer.
    kErrorProfilingTimerWrongPreset,
    /// Some error during reseting the timer.
    kErrorProfilingTimerWrongReset,
    /// Some error during stoping profiling timer.
    kErrorProfilingTimerWrongStop,
    /// Base must be >= 0!
    kErrorProfilingTimerZeroBase,
    /// Printing a printer that's not specified
    kErrorProfilingTimerPrintUnknownTimer,
    /// Buffers to send and recieve  halo has different sizes.
    /// Checked with HaloToExchange::Init().
    kErrorSendAndReceiveBuffersHasDifferentSizes,
    /// Grid size in any dimension should be more or equal to 1.
    /// @see GridInputConfig::set_total_grid_length()
    kErrorSettingWrongGridSize,
    /// Number of elements in subdomain should be more than 1.
    kErrorSingleElementSubdomain,
    /// case of too many processes during decomposition
    kErrorSubdomainSizeLessThanHaloWidth,
    /// onza::SimulationInputConfig forced to set total pml width in
    /// some dimension bigger than size of this dimension.
    kErrorTooWidePml,
    /// found uninitiated simulation core when
    /// HaloExchangeProcess::RunSimulation()
    kErrorUninitiatedSimulationCore,
    /// Using InputConfig too early
    kErrorUsingInputConfigTooEarly,
    /// Should use some FDTD algorithm.
    kErrorWrongAlgorithm,
    /// FDTD algorithm`s should have some subdomain halo.
    kErrorWrongHaloWidth,
    /// wrong index while checking with
    /// HaloExchangeProcess::CheckSubdomainIndexes()
    kErrorWrongIndexDifference,
    /// FDTD algorithm`s time depth should be >= 1!
    kErrorWrongTimeDepth
  };
  /// @brief Simulation status
  enum SimulationStatus {kSimulationStatusFinished = 0,
                         kSimulationStatusRunning,
                         kSimulationStatusInitiated,
                         kSimulationStatusUndefined
                         };
  /// @brief Algorithm selection
  enum Algorithm { kAlgorithmSimpleX1D = 0, kAlgorithmSimpleY1D,
                   kAlgorithmSimpleZ1D,
                   kAlgorithmSimpleTMz2D,
                   kAlgorithmSimple3D};
  /// @brief Simulation core and halo border predefined names for data
  /// components
  enum DataComponents {kEz = 0, kHy, kCezh, kSrcEz,  // 1D axis x
                       kCeze, kChyh, kChye,
                       kHx, kChxh, kChxe,          // 2D TM axis z
                       kEx, kCexe, kCexh,
                       kEy, kCeye, kCeyh,
                       kHz, kChzh, kChze};  // 3D FDTD
  /// @brief Status of reading input config.
  ///
  /// For use with onza::SimulationInputConfig::status_
  enum InputConfig {kInputConfigUndefined = 0,
                    kInputConfigErrorWidePml,
                    kInputConfigAllDone = 10000};
  /// @brief Type of boundary conditions for simulation model.
  ///
  /// Currently only PML condition is valid
  enum BoundaryCondition {
    /// Uninitialized boundary condition.
    kBoundaryConditionUndefined = 0,
    /// Boundary is periodical.
    kBoundaryConditionPeriodical,
    /// Boundary is reduced.
    kBoundaryConditionReduced,
    /// There is PML region near this boundary inside current domain.
    kBoundaryConditionPML
  };  // end of enum BoundaryCondition
  /// @brief Total number of Cartesian axes.
  ///
  /// Should allways be 3!!! To define 1D or 2D model set
  /// corresponding grid length(s) to 1.
  const int kDimensions = 3;
  /// @brief Names for Cartesian axes.
  ///
  /// Use them to access dimesion relative arrays.
  enum Axis {kAxisX = 0, kAxisY, kAxisZ};
  /// @brief Position of halo borders to exchange relative to current process.
  ///
  /// Order in enum is used by GetOppositeBorder()
  /// Simple math (due to order in enums): border_axis = border % kDimensions;
  enum BorderPosition {kBorderLeft = 0, kBorderBottom, kBorderBack,
                       kBorderRight, kBorderTop, kBorderFront};
  /// @brief Some values for tags for MPI messages.
  enum MpiTag {kMpiTagCheckIndex = 1};
  // ********************************************************************** //
  // **********************         Global functions       **************** //
  // ********************************************************************** //
  /// @brief Returns border opposite to input border
  ///
  /// Calculation of opposite border is based on definition of #BorderPosition.
  /// @param[in] input_border
  /// @return Border opposite to input_border
  inline BorderPosition GetOppositeBorder(BorderPosition input_border) {
    return static_cast<BorderPosition>((input_border+3)%6);
  }
  /// @brief Mutiply components of vector
  ///
  /// Used to get product of any vector type, e.g.
  /// @code
  /// int array_of_int[] = {1, 2, 3};
  /// // 6
  /// printf("Product = %i\n", onza::MultiplyComponents(array_of_int));
  /// double array_of_double[] = { 0.2, 0.3};
  /// // 0.06
  /// printf("Product = %g\n", onza::MultiplyComponents(array_of_double));
  /// @endcode
  /// @param[in] Some vector
  /// @return Product of vector components.
  /// @todo3 check speed of template realization of MultiplyComponents function.
  /// as far as it is belived it has very small performance overhead
  /// (with compiler optimizing flag -O2 or -O3).
  template <class VectorType, int dimensions> inline
      VectorType MultiplyComponents(VectorType(& vector)[dimensions]) {
    if (dimensions == 3) return  vector[0]*vector[1]*vector[2];
    if (dimensions == 2) return  vector[0]*vector[1];
    if (dimensions == 1) return  vector[0];
    if (dimensions > 3) {
      VectorType product = vector[0];
      for (int i = 1; i < dimensions; ++i) {
        product *= vector[i];
      }  // end of for
      return product;
    }  // end of if dimensions
  }  // end of template MultiplyComponents
  /// @brief Sum up components of vector
  ///
  /// Used to get sum of any vector type components.
  /// @see MultiplyComponents()
  /// @param[in] Some vector
  /// @return Sum of vector components.
  /// @todo check speed of template realization of SumUpComponents function.
  /// as far as it is belived it has very small performance overhead
  /// (with compiler optimizing flag -O2 or -O3).
  template <class VectorType, int dimensions> inline
      VectorType SumUpComponents(VectorType(& vector)[dimensions]) {
    if (dimensions == 3) return  vector[0] + vector[1] + vector[2];
    if (dimensions == 2) return  vector[0] + vector[1];
    if (dimensions == 1) return  vector[0];
    if (dimensions > 3) {
      VectorType product = vector[0];
      for (int i = 1; i < dimensions; ++i) {
        product += vector[i];
      }  // end of for
      return product;
    }  // end of if dimensions
  }  // end of template SumUpComponents
  /// @brief Square some value
  template<class T> inline T pow2(const T value) {return value*value;}
  /// @brief Convert string to positive integer.
  ///
  /// Return zero for any error.
  /// Modified anser of Nawaz from
  /// http://stackoverflow.com/questions/7370887/
  ///            convert-an-ascii-string-to-long-long
  template<class DecimalType> DecimalType
      convert_to_positive(const std::string &s, DecimalType &output) {
    DecimalType Error = 0;
    if (s.size() == 0) return Error;
    DecimalType v = 0;
    size_t i = 0;
    char sign = (s[0] == '-' || s[0] == '+') ? (++i, s[0]) : '+';
    for (; i < s.size(); ++i) {
      if (s[i] < '0' || s[i] > '9') return Error;
      v = v * 10 + s[i] - '0';
    }
    if (sign == '-') return Error;
    output = v;
    return output;
  }  // end of DecimalType convert(const std::string &s);
  /// @bief Convert anything to std::string
  ///
  /// Modified  answer from Konrad Rudolph from
  /// http://stackoverflow.com/questions/332111/
  ///           how-do-i-convert-a-double-into-a-string-in-c
  template <typename T>
      std::string to_string(T const& value) {
    std::stringstream sstr;
    // @todo3 Optimized for double (fraction in %) to string conversion
    // Need work to check template is still usable for other conversions.
    sstr.width(14);
    sstr.precision(11);
    sstr.fill('0');
    sstr.setf(std::ios::fixed);
    sstr << value;
    return sstr.str();
  }
}  // end of namespace onza
#endif  // SRC_COMMON_H_
