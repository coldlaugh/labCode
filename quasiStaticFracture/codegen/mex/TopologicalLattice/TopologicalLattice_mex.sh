MATLAB="/Applications/MATLAB_R2016a.app"
Arch=maci64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/Users/leyou/.matlab/R2016a"
OPTSFILE_NAME="./setEnv.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for TopologicalLattice" > TopologicalLattice_mex.mki
echo "CC=$CC" >> TopologicalLattice_mex.mki
echo "CFLAGS=$CFLAGS" >> TopologicalLattice_mex.mki
echo "CLIBS=$CLIBS" >> TopologicalLattice_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> TopologicalLattice_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> TopologicalLattice_mex.mki
echo "CXX=$CXX" >> TopologicalLattice_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> TopologicalLattice_mex.mki
echo "CXXLIBS=$CXXLIBS" >> TopologicalLattice_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> TopologicalLattice_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> TopologicalLattice_mex.mki
echo "LD=$LD" >> TopologicalLattice_mex.mki
echo "LDFLAGS=$LDFLAGS" >> TopologicalLattice_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> TopologicalLattice_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> TopologicalLattice_mex.mki
echo "Arch=$Arch" >> TopologicalLattice_mex.mki
echo OMPFLAGS= >> TopologicalLattice_mex.mki
echo OMPLINKFLAGS= >> TopologicalLattice_mex.mki
echo "EMC_COMPILER=Xcode with Clang" >> TopologicalLattice_mex.mki
echo "EMC_CONFIG=optim" >> TopologicalLattice_mex.mki
"/Applications/MATLAB_R2016a.app/bin/maci64/gmake" -B -f TopologicalLattice_mex.mk
