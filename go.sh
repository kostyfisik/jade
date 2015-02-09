#!/bin/bash
# @file   go.sh
# @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
# @copyright 2013 Ladutenko Konstantin
# @section LICENSE
# This file is part of JADE++.
#
# JADE++ is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# JADE++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with JADE++.  If not, see <http://www.gnu.org/licenses/>.
#
# @brief  Compile JADE++ and run tests.
usage_msg="\nUsage: JADE++ can be (re-)compiled and run with call \n
$./go.sh [mode]\n
Possible modes:\n
new     \t - Default. Make new build and run.
old     \t - Re-build and run.
test    \t - Re-build and run JADE++ on a set of test functions.
single  \t - Re-build and run JADE++ on a single test function.
Possible enviroment parameters:\n
\n
JADE_MPI_size   \t- total number of MPI processes.\n
JADE_MPI_nodes  \t- total number of MPI nodes for cluster enviroment.\n
\n
In the case of enviroment parameters are not (or set with value 'unset')\n
set some default vale (depending on executing host) will be used.\n"
#############################################################################
if [[ ! $JADE_MPI_size ]]; then JADE_MPI_size="unset"; fi
if [[ ! $JADE_MPI_nodes ]]; then JADE_MPI_nodes="unset"; fi
MPI_options=
#MPI_options="--bind-to-core"
#############################################################################
#   Parse input parameters   
#############################################################################
mode=$1; config_file=$2; wrong_params=$3
#HOST=`cat /etc/hostname`
HOST=`hostname`
if [[ $config_file ]]; then
    echo ================ !ERROR! =================
    echo Should be no more than single  input parameter
    echo ================ !ERROR! =================
    echo -e  $usage_msg
    exit 1
fi
# # Mode names
mode_new1="new"
# ? clang mode 
mode_new2="new2"  # gcc
# mode_new3="new3"  # gcc with -O3
mode_old1="old";     mode_old2="old2";         
mode_test="test"; # mode_prof="prof";     mode_old2prof="old2prof"; mode_pgo="pgo"
mode_single="single";
# mode_custom="custom"; mode_debug="debug";       mode_build="build"
# Default values
yes="yes";        no="no"
compiler_gcc="gcc";
#compiler_clang="clang"
usedCompiler=$compiler_gcc # or clang
useGCC48=$yes  # use gcc 4.8 if it is available in build area of scripts folder
isNew=$yes;
isTest=$no ;  isProfile=$no ; isPGO=$no
isBuildOnly=$no;
# if [[ $mode = $mode_build || $config_file = $mode_build ]]; then
#     isBuildOnly=$yes    
# fi
path_jade=$PWD;   path_bin=$path_jade/bin;   path_build=$path_jade/build
path_src=$path_jade/src
# Should be same as in cmake config in $path_src

if [[ ! $mode ]]; then  mode=$mode_new1; fi
if [[ $mode = "help" || $mode = "--help" || $mode = "-h" ]]; then
 echo -e $usage_msg
 exit 0
fi 
#jade_bin="run-size-sweep"
#jade_bin="run-optimize-cloak"
jade_bin="run-coating-w-sweep"
#jade_bin="run-optimize-meander-cloak"
#jade_bin="run-quasi-pec-spectra"
if [[ $mode = $mode_test ]]; then
    echo Run JADE++ test!
    jade_bin="run-jade-test"
fi 
if [[ $mode = $mode_single ]]; then
    echo Run JADE++ test for single function!
    jade_bin="run-jade-single"
fi 
# # Check mode
# if [[ $mode != $mode_new1 \
#  # && $mode != $mode_new2 && $mode != $mode_new3 && \
#  #    $mode != $mode_old1 && $mode != $mode_old2 && \
#  #    $mode != $mode_test && \
#  #    $mode != $mode_prof && $mode != $mode_old2prof && \
#  #    $mode != $mode_pgo && \
#  #    $mode != $mode_custom && \
#  #    $mode != $mode_build && \
#  #    $mode != $mode_debug
#  ]]; then
#     # So mode may be miss spelled or contains config file path.
#     if [[ $config_file ]]; then
#         echo ================ !ERROR! =================
#         echo Undefined mode
#         echo ================ !ERROR! =================
#         echo -e  $usage_msg
#         exit 1
#     fi  
#     echo Using default mode: Full build with clang compiler.
#     config_file=$mode 
#     mode=$mode_new1
# fi 
# if [[ ( $mode = $mode_test || $mode = $mode_build ) && $config_file ]]; then
#     echo ================ !ERROR! =================
#     echo Test and build modes do not support external config file
#     echo ================ !ERROR! =================
#     exit 1
# fi
# # Check config file(s)
# path_default_config=$path_jade/data/default-jade.config
# # Test containing "X1D-zero" in its name is a special optin in TuneJADEOptionsMPI
# tests="self-test-X1D-zero self-test-TMz2D-speedup self-test-3D-simple"
# if [[ $mode = $mode_test ]]; then
#     echo Check for tests config files...
#     for test in $tests; do
#         test_name=${test//-/_}
#         path_test_config="path_${test_name}_config"
#         eval path_${test_name}_config=$path_jade/data/$test.config  
#         echo ${!path_test_config}
#         if [[ ! -r ${!path_test_config} ]];
#         then
#             echo ================ !ERROR! =================
#             echo Test mode was not able to access some of config files for tests
#             echo ${!path_test_config}
#             echo ================ !ERROR! =================
#             exit 1
#         fi
#     done
# fi
# # Should be tested after checking $test_mode 
# if [[ ! $config_file && $mode != $mode_build ]]; then
#     echo Setting default config file $path_default_config
#     config_file=$path_default_config
# fi
# if  [[ ! -a $config_file && $isBuildOnly != $yes ]]; then
#     echo ================ !ERROR! =================
#     echo Was not able to found config file using path: $config_file
#     echo ================ !ERROR! =================
#     exit 1
# fi
# # Convert relative path for custom config file to absolute.
# firstChar=${config_file:0:1}
# if [[ $firstChar != "/" && $isBuildOnly != $yes ]]; then
#     config_file=$path_jade/$config_file
#     echo Change config file path to absolute path: $config_file
# fi
echo ============ Current script settings =============
echo mode: $mode
# if [[ $mode != $mode_test && $isBuildOnly != $yes ]]; then
#     echo config file: $config_file
# fi
echo base dir path: $path_jade
# Check directory structure
if [[ ! -d $path_src ]]; then
    echo ================ !ERROR! =================
    echo No source folder $path_src
    echo ================ !ERROR! =================
    exit
fi
force_new_build=$no
if [[ -a $path_build ]]; then
    echo Found build folder.
else
    echo Creating build folder...
    mkdir build
    force_new_build=$yes
fi
if [[ -a $path_bin ]]; then
    echo Found bin folder.
else
    echo Creating build folder...
    mkdir bin
fi
#############################################################################
#   Compile settings
#############################################################################
echo ============ Compile settings =============
if [[ ( $mode = $mode_old1 || $mode = $mode_old2 || $mode = $mode_test || $mode=$mode_single ) && $force_new_build = $no ]]; then
    isNew=$no
    echo Recompile mode is on.
else
    echo Cleaning build path.
    rm -r $path_build/* >/dev/null 2>&1
fi 
echo Cleaning bin
rm -r $path_bin/* >/dev/null 2>&1
# Profiling mode should use gcc compiler.
if [[ $mode = $mode_prof || $mode = $mode_old2prof ]]; then
    echo Using gprof.
    isProfile=$yes
    flag_cmake_profile="-DCMAKE_CXX_FLAGS=-pg -DCMAKE_EXE_LINKER_FLAGS=-pg"
fi 
# PGO mode should use gcc compiler.
if [[ $mode = $mode_pgo ]]; then
    echo Using pgo.
    isPGO=$yes
fi 
if  [[ $HOST == "rh-lum.metalab.ifmo.ru" ]]; then
    echo Setting MPI path on rh-lum.metalab.ifmo.ru !
    ompi_path_bin=/usr/lib64/openmpi/bin/
    ompi_path_lib=/usr/lib64/openmpi/lib/
    if [ -d "$ompi_path_bin" ] && [[ ":$PATH:" != *":$ompi_path_bin:"* ]]; then
        PATH="${PATH:+"$PATH:"}$ompi_path_bin"
    fi
    export LD_LIBRARY_PATH=$ompi_path_lib
fi
# Set compiler

if [[ $mode = $mode_new1 || $mode = $mode_old1 || $mode = $mode_test || $mode=$mode_single ]]; then
    echo Using \'clang\' compiler.
    usedCompiler=$compiler_clang
    #path_clang33=/home/mmedia/soft/clang/clang+llvm-3.3-amd64-debian6/bin/
    if [[ $HOST == "head.phoif.ifmo.ru" ]]; then
        path_clang33=/home/nfs-shared/tig/clang
    elif  [[ $HOST == "dmmrkovich-birzha" ]]; then
        path_clang33=
    elif  [[ $HOST == "deb00" || $HOST == "tig-laptop2" ]]; then
        path_clang33=/home/mmedia/soft/clang/clang33        
    elif  [[ $HOST == "rh-lum.metalab.ifmo.ru" ]]; then
        path_clang33=/home/tig/clang/clang33        
	ompi_path_bin=/usr/lib64/openmpi/bin/
	ompi_path_lib=/usr/lib64/openmpi/lib/
	if [ -d "$ompi_path_bin" ] && [[ ":$PATH:" != *":$ompi_path_bin:"* ]]; then
            PATH="${PATH:+"$PATH:"}$ompi_path_bin"
	fi
	export LD_LIBRARY_PATH=$ompi_path_lib
    else
        path_clang33=
    fi
    echo path_clang: $path_clang33
    export OMPI_CC=$path_clang33/bin/clang
    export OMPI_CXX="$path_clang33/bin/clang++ -I$path_clang33/include -stdlib=libc++"
    if  [[ $HOST == "tig-laptop3" ]]; then
	path_clang33=/usr 
	export OMPI_CC=$path_clang33/bin/clang
	export OMPI_CXX="$path_clang33/bin/clang++ -I$path_clang33/inclde/c++/v1 -stdlib=libc++"
    fi

else
    path_gcc48=/home/mmedia/soft/gcc/gcc48
    if [[ -a $path_gcc48 && $useGCC48 = $yes ]]; then
        echo Using \'gcc-4.8\' compiler.
        export OMPI_CC=$path_gcc48/bin/gcc
        #export OMPI_CXX="$path_gcc48/bin/g++"
        export OMPI_CXX="$path_gcc48/bin/g++ -I$path_gcc48/include/c++/4.8.1 -Wl,-rpath,$path_gcc48/lib64"
    else
        echo Using gcc compiler.
        path_gcc48=
        export OMPI_CC=gcc
        export OMPI_CXX=g++
    fi
fi 

# Select OMPI_CXXFLAGS
#debug
#flags_O2="-std=c++11"
flags_O2="-O2 -Wall -std=c++11"
#flags_O2="-O2 -ftemplate-depth-30 -Wall -std=c++11"
flags_debug="-ftemplate-depth-30 -Wall -std=c++11  -stdlib=libc++"
flags_O3="-O3 -ffast-math -ftemplate-depth-30 -march=native -mtune=native -mfpmath=both -malign-double -mstackrealign -ftree-vectorize -msse2 -ftree-vectorizer-verbose=5  -Wall  -std=c++11 -stdlib=libc++" 
# TODO option -flto   -- Do we need it?
export OMPI_CXXFLAGS=$flags_O2
if [[ $mode = $mode_debug ]]; then
    echo Using debug mode.
    export OMPI_CXXFLAGS=$flags_debug
fi
if [[ $mode = $mode_new3 ]]; then
    export OMPI_CXXFLAGS=$flags_O3
fi
#############################################################################
#   Build
#############################################################################
function BuildJADE {
    echo "Building JADE++..."
    local path_tmp=`pwd`
    cd $path_build
    # if [[ $flag_cmake_profile && $flag_cmake_pgo ]]; then
    #     echo ================ !ERROR! =================
    #     echo Config with cmake flags $flag_cmake_profile and $flag_cmake_pgo
    #     echo ================ !ERROR! =================
    #     exit        
    # fi
    flag_cmake=
    # if [[ $flag_cmake_profile ]]; then
    #     flag_cmake=$flag_cmake_profile
    # fi
    # if [[ $flag_cmake_pgo ]]; then
    #     flag_cmake=$flag_cmake_pgo
    # fi
    if [[ $isNew = $yes ]]; then
        echo $PWD
        echo "Doing new build.."
        CC=mpicc CXX=mpic++ VERBOSE=1 cmake $path_jade \
            -DCMAKE_INSTALL_PREFIX="$path_bin" \
            $flag_cmake
    fi
    make -j4 
    make install
    cd $path_tmp
}  # end of function BuildJADE
BuildJADE
#############################################################################
#   Run
#############################################################################
echo "Executing $jade_bin on host -> $HOST <-"
if [[ $isBuildOnly = $yes ]]; then
    cd $path_bin
    mv run-jade-test jade-test.bin
    exit 0; 
fi
if [[ $mode = $mode_test ]]; then isTest=$yes; fi
if [[ $isProfile = $yes ]]; then 
    # Grpof output for each process in separate file
    export GMON_OUT_PREFIX='gmon.out'
fi
# Task and host dependant MPI parameters tuning.
function TuneJADEOptionsMPI {
    # if [[ $JADE_MPI_size = "unset" && \
    #      "$config_file" = *"X1D-zero"* ]]; then
    #     JADE_MPI_size=2
    #     JADE_MPI_nodes=1
    # fi
    if [[ $HOST == "head.phoif.ifmo.ru" ]]; then
        if [[ $JADE_MPI_size = "unset" ]]; then JADE_MPI_size=16; fi
        if [[ $JADE_MPI_nodes = "unset" ]]; then JADE_MPI_nodes=8; fi
    elif  [[ $HOST == "deb00" || $HOST == "dmmrkovich-birzha" ]]; then
        if [[ $JADE_MPI_size = "unset" ]]; then JADE_MPI_size=4; fi            
    elif  [[ $HOST == "rh-lum.metalab.ifmo.ru" ]]; then
        if [[ $JADE_MPI_size = "unset" ]]; then JADE_MPI_size=12; fi            
    else
        if [[ $JADE_MPI_size = "unset" ]]; then JADE_MPI_size=2; fi
    fi
}  # end of TuneJADEOptionsMPI
# Host dependant MPI run of JADE
function RunJADE {
    local path_tmp=`pwd`
    cd $path_bin
    # cp $config_file $path_bin/jade.config
    if [[ $HOST == "head.phoif.ifmo.ru" ]]; then
        echo "Waiting for shared file system to distibute files..."
        sleep 2
        echo "(1) Nodes $JADE_MPI_nodes procs $JADE_MPI_size"
        salloc -N $JADE_MPI_nodes -n $JADE_MPI_size -p max1day \
            mpirun $MPI_options ./$jade_bin #jade.config
    else
        echo "(1) Nodes 1  procs $JADE_MPI_size"
        #time mpirun -np $JADE_MPI_size $MPI_options ./$jade_bin #jade.config
        mpirun -np $JADE_MPI_size $MPI_options ./$jade_bin #jade.config
    fi    
    cd $path_tmp
}  # end of RunJADE
if [[ $isTest = $no ]]; then
    TuneJADEOptionsMPI
    time RunJADE
    if [[ $isPGO = $yes ]]; then
        # export OMPI_CXXFLAGS="$flags_O2 -fprofile-generate"
        # export OMPI_LDFLAGS=$OMPI_CXXFLAGS
        flag_cmake_pgo="-DCMAKE_CXX_FLAGS=-fprofile-generate -DCMAKE_EXE_LINKER_FLAGS=-fprofile-generate"
        rm -r $path_bin/*
        BuildJADE
        echo ============================================
        echo Run PGO profiling
        sleep 2
        RunJADE

        cp $path_build/src/mpi-decomposition/CMakeFiles/mpi-decomposition.dir/__/simulation-core/basic-fdtd.cc.gcda $path_build/src/simulation-core/CMakeFiles/simulation-core.dir/basic-fdtd.cc.gcda
        cp $path_build/src/mpi-decomposition/CMakeFiles/mpi-decomposition.dir/__/profiling/timer.cc.gcda $path_build/src/profiling/CMakeFiles/profiling.dir/timer.cc.gcda
        # export OMPI_CXXFLAGS="$flags_O2 -fprofile-use"
        # export OMPI_LDFLAGS=$OMPI_CXXFLAGS
        flag_cmake_pgo="-DCMAKE_CXX_FLAGS=-fprofile-use -DCMAKE_EXE_LINKER_FLAGS=-fprofile-use"
        rm -r $path_bin/*
        BuildJADE
        echo ============================================
        echo Run after PGO compiled
        sleep 2
        RunJADE 
    fi
fi  # end of if [[ $isTest = $no ]]
# if [[ $isProfile = $yes ]]; then 
#     gprof  $jade_bin gmon.out* | grep '^index'
#     for file in gmon.out*;
#     do
#         echo $file
#         gprof  $jade_bin $file | grep 'DoBorderStep' | grep '^\['
#         gprof  $jade_bin $file | grep 'DoStep' | grep '^\['
#     done
#     echo Average
#     gprof  $jade_bin gmon.out* | grep 'DoBorderStep' | grep '^\['
#     gprof  $jade_bin gmon.out* | grep 'DoStep' | grep '^\['
#     gprof --no-flat-profile $jade_bin gmon.out* > average-profile
#     gprof $jade_bin gmon.out* > flat-average-profile
#     rm gmon.out*
# fi  # end of if isProfile
if [[ $isTest = $yes ]]; then
    TuneJADEOptionsMPI
    RunJADE
fi
# if [[ $isTest = $yes ]]; then
#     echo "Prepare files for tests ..."
#     cd $path_bin
#     for test in $tests; do
#         # For each test make dir, copy config and binary to it, run.
#         cd $path_bin
#         mkdir $test
#         test_name=${test//-/_}
#         #path_test_config="path_${test_name}_config"
#         path_test="$path_bin/$test"  
#         #path_test_config="path_${test}_config"
#         #cp ${!path_test_config} ${path_test}/jade.config
#         cp $jade_bin $path_test
#         cd $path_test
#         backup_JADE_MPI_size=$JADE_MPI_size
#         backup_JADE_MPI_nodes=$JADE_MPI_nodes
#         if [[ $test = "self-test-X1D-zero" ]]; then
#             JADE_MPI_size=2
#             JADE_MPI_nodes=1
#         fi
#         echo "============ Running test $test ============="
#         if [[ $HOST == "head.phoif.ifmo.ru" ]]; then
#             echo "Waiting for shared file system to distibute files..."
#             sleep 2
#             if [[ $JADE_MPI_size = "unset" ]]; then JADE_MPI_size=16; fi
#             if [[ $JADE_MPI_nodes = "unset" ]]; then JADE_MPI_nodes=8; fi
#             echo "(1) Nodes $JADE_MPI_nodes procs $JADE_MPI_size"
#             salloc -N $JADE_MPI_nodes -n $JADE_MPI_size -p max1hour \
#                 mpirun $MPI_options ./$jade_bin #jade.config
#             if [[ $test = "self-test-TMz2D-speedup" ]]; then
#                 echo "(*******) Nodes 16 procs 128 (1024 x 1024, step 1000, ~15.4s)"
#                 salloc -N 16 -n 128 -p max1hour \
#                     mpirun $MPI_options ./$jade_bin #jade.config
#                 echo "(*******) Nodes 16 procs 16 (1024 x 1024, step 1000, ~7s)"
#                 salloc -N 16 -n 16 -p max1hour \
#                     mpirun $MPI_options ./$jade_bin #jade.config
#                 echo "(*******) Nodes 1 procs 8 (1024 x 1024, step 1000, ~30.2s)"
#                 salloc -N 1 -n 8 -p max1hour \
#                     mpirun $MPI_options ./$jade_bin #jade.config
#                 echo "(*******) Nodes 1 procs 1 (1024 x 1024, step 1000, ~47.6s)"
#                 salloc -N 1 -n 1 -p max1hour \
#                     mpirun $MPI_options ./$jade_bin #jade.config
#             fi

#         elif  [[ $HOST == "deb00" || $HOST == "dmmrkovich-birzha" ]]; then
#             echo "(*******) Procs 1"
#             mpirun -np 1 $MPI_options ./$jade_bin #jade.config
#             echo "(*******) Procs 2"
#             mpirun -np 2 $MPI_options ./$jade_bin #jade.config
#             echo "(*******) Procs 4"
#             mpirun -np 4 $MPI_options ./$jade_bin #jade.config
#         else
#             if [[ $JADE_MPI_size = "unset" ]]; then JADE_MPI_size=2; fi
#             echo "(1) Nodes 1  procs $JADE_MPI_size"
#             mpirun -np $JADE_MPI_size $MPI_options ./$jade_bin #jade.config
#         fi
#         JADE_MPI_size=$backup_JADE_MPI_size
#         JADE_MPI_nodes=$backup_JADE_MPI_nodes
#         rm *.jade  >/dev/null  2>&1
#     done  # end of for test in $tests; do
# fi  # end of if [[ $isTest = $yes ]]
# if [[ $config_file = $path_self_test_X1D_zero_config ]]; then
#     echo "Prepare *.png from gnuplot ..."
#     cp $path_jade/data/gnuplot/* ./
#     ./gnuplot-all.sh >/dev/null  2>&1
#     # mkdir tmpdir
#     # cp *0241-* $path_bin/tmpdir
#     # rm $path_bin/*    
#     # cp $path_bin/tmpdir/* ./
#     # rm -r $path_bin/tmpdir
# fi
#rm *.jade  >/dev/null  2>&1
cd bin
chmod +x run-gnuplot*.sh
for file in `ls run-gnuplot*`; do 
 ./$file
done   

cp $path_jade/scripts/prepare-overview.py $path_bin
cp $path_jade/scripts/filter.py $path_bin
cd path_bin
./filter.py
./prepare-overview.py
