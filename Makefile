# to use:  make -f make_beda  step1

#PETSC_DIR = /home/gelmut/petsc-3.0.0-p12
PETSC_DIR = /home/thanhnguyen/soft/petsc-3.0.0-p12
#PETSC_DIR = 

#PETSC_ARCH = linux-gnu-c-debug

#PETSC_ARCH = /home/thanhnguyen/soft/petsc-3.0.0-p12
PETSC_ARCH = 

include ${PETSC_DIR}/conf/base
     
     
CC = gcc
CXX = g++

#CFLAGS   = -O3 -xHOST -g

CXXFLAGS = -O3 -xHOST 
#-Wno-deprecated 
#-Wno-write-strings 
#-Wnon-template-friend
 
CPPFLAGS +=  -I${MV_INCLUDE} -I. -Wno-non-template-friend
#-Wno-deprecated
#-Wno-write-strings 
#-Wno-format

FFLAGS = -I../

BOPT = g++

#MV_DIR = /home/gelmut/mv
#MV_INCLUDE = /home/gelmut/mv/include
#MPI_INCLUDE=/home/gelmut/petsc-3.0.0-p12/linux-gnu-cxx-debug/include/mpiuni
#PETSC_ADDITIONAL_INCLUDE=/home/gelmut/petsc-3.0.0-p12/linux-gnu-c-debug/include
#MPI_INCLUDE = /home/gelmut/petsc-3.0.0-p12/include/mpiuni

MPI_INCLUDE = ${PETSC_DIR}/include/mpiuni

CPPFLAGS += -I$(MV_INCLUDE)
CPPFLAGS += -I$(MPI_INCLUDE)

CPPFLAGS += -I.
CPPFLAGS += -fpermissive

# Allow unsupported header files from old version of gcc
# (like <strstream> which is not at all included in gcc 3.4,
# for some good reasons
#CFLAGS += -idirafter /usr/include/c++/3.2.3/backward/

OBJGRID = WavesGridB.o \
	WavesNeighborFE.o \
	WavesSparseDS.o \
	WavesFEcomp.o 


OBJGRID_impl = WavesGridB.o \
	       WavesNeighborFE.o \
               WavesSparseDS.o \
	       WavesFEcomp_new.o
 
test_optional = WavesOptional.o \
 WavesSDGeometry.o \
 WavesSDOperator.o \
 WavesSDMaskIndexes.o \
 Wavesutil2.o \
 WavesSDIndexes.o 
 

OBJ6 = WavesSDGeometry.o \
       WavesSDOperator.o \
       WavesSDIndexes.o \
       WavesSDMaskIndexes.o \
       WavesExchangeOperator.o \
       Wavesutil2.o \
       Wavesutil3.o \
       Wavesinitgeom.o \
       Wavesoverlap.o  
  

#Mesh generation:

mesh_gen:  mesh_gen3D.o   $(test_optional) WavesOutputs.o    WavesfindBoundaryNodes.o  $(OBJGRID) chkopts
	g++ -o  mesh_gen mesh_gen3D.o  $(test_optional)  WavesOutputs.o  WavesfindBoundaryNodes.o  $(OBJGRID) $(PETSC_LIB) $(CCFLAGS)


inp2inp:  inp2inp.o   $(test_optional) WavesOutputs.o    WavesfindBoundaryNodes.o  $(OBJGRID) chkopts
	g++ -o  inp2inp inp2inp.o  $(test_optional)  WavesOutputs.o  WavesfindBoundaryNodes.o  $(OBJGRID) $(PETSC_LIB) $(CCFLAGS)

vec_inp2inp:  vec_inp2inp.o   $(test_optional) WavesOutputs.o    WavesfindBoundaryNodes.o  $(OBJGRID) chkopts
	g++ -o  vec_inp2inp vec_inp2inp.o  $(test_optional)  WavesOutputs.o  WavesfindBoundaryNodes.o  $(OBJGRID) $(PETSC_LIB) $(CCFLAGS)

# forward solver and related functions:

forward_solver3D: forward_solver3D.o  $(test_optional) WavesReflectedWave.o  WavesScalarEqOpOpt.o   WavesPlaneWaveOpexpl.o WavesBCOperators.o WavesFuncDefs.o WavesOutputs.o  WavesWaveEquationOp.o WavesWaveEqOpSeismic.o  WavesWaveEqOpPar.o  $(OBJGRID) chkopts
		${CXX} -o  forward_solver3D  forward_solver3D.o  $(test_optional)  WavesReflectedWave.o WavesScalarEqOpOpt.o WavesPlaneWaveOpexpl.o   WavesBCOperators.o WavesFuncDefs.o WavesOutputs.o WavesWaveEquationOp.o WavesWaveEqOpSeismic.o  WavesWaveEqOpPar.o  $(OBJGRID) $(PETSC_LIB) $(CCFLAGS)

# SAR_3D: temporary not available.
sar_3d: SAR_3D.o  $(test_optional) WavesBCOperators.o WavesFuncDefs.o WavesOutputs.o $(OBJGRID) chkopts
		${CXX} -o  sar_3d  SAR_3D.o  $(test_optional)  WavesBCOperators.o WavesFuncDefs.o WavesOutputs.o $(OBJGRID) $(PETSC_LIB) $(CCFLAGS)
		
laplace :  laplace_transform.o  $(test_optional) WavesBCOperators.o WavesFuncDefs.o WavesOutputs.o $(OBJGRID) chkopts
	g++ -o  laplace  laplace_transform.o  $(test_optional)  WavesBCOperators.o WavesFuncDefs.o WavesOutputs.o $(OBJGRID) $(PETSC_LIB) $(CCFLAGS)

extract_scatteredwave :  extract_scatteredwave.o  $(test_optional) WavesBCOperators.o WavesFuncDefs.o WavesOutputs.o  $(OBJGRID) chkopts
	g++ -o  extract_scatteredwave  extract_scatteredwave.o  $(test_optional)  WavesBCOperators.o WavesFuncDefs.o WavesOutputs.o  $(OBJGRID) $(PETSC_LIB) $(CCFLAGS)

create_exchangeMask :  create_exchangeMask.o  $(test_optional) WavesReflectedWave.o  WavesScalarEqOpOpt.o   WavesPlaneWaveOpexpl.o WavesBCOperators.o WavesFuncDefs.o WavesOutputs.o  WavesWaveEquationOp.o  WavesWaveEqOpPar.o WavesWaveEqOp123D.o WavesWaveEqOpSeismic.o $(OBJGRID) chkopts
	g++ -o  create_exchangeMask  create_exchangeMask.o  $(test_optional)  WavesReflectedWave.o WavesScalarEqOpOpt.o WavesPlaneWaveOpexpl.o   WavesBCOperators.o WavesFuncDefs.o WavesOutputs.o WavesWaveEquationOp.o  WavesWaveEqOpPar.o WavesWaveEqOp123D.o WavesWaveEqOpSeismic.o  $(OBJGRID) $(PETSC_LIB) $(CCFLAGS)


# globally convergent algorithm: 

dat4gca_scatwave :  dat4gca_scatwave.o  $(test_optional)  WavesGridFEHier.o   WavesfindBoundaryNodes.o  WavesBCOperators.o WavesFuncDefs.o WavesOutputs.o $(OBJGRID) chkopts
	g++ -o  dat4gca_scatwave  dat4gca_scatwave.o $(test_optional)  WavesGridFEHier.o   WavesfindBoundaryNodes.o  WavesBCOperators.o WavesFuncDefs.o WavesOutputs.o $(OBJGRID) $(PETSC_LIB) $(CCFLAGS)

dat4gca_totalwave :  dat4gca_totalwave.o  $(test_optional)  WavesGridFEHier.o   WavesfindBoundaryNodes.o  WavesBCOperators.o WavesFuncDefs.o WavesOutputs.o $(OBJGRID) chkopts
	g++ -o  dat4gca_totalwave  dat4gca_totalwave.o $(test_optional)  WavesGridFEHier.o   WavesfindBoundaryNodes.o  WavesBCOperators.o WavesFuncDefs.o WavesOutputs.o $(OBJGRID) $(PETSC_LIB) $(CCFLAGS)

gca_totalwave :  gca_totalwave.o  $(test_optional)  WavesWaveEqOpElliptic_log.o  WavesElliptic_GlobConvAlg.o   WavesPoissonEqOp.o	  WavesGridFEHier.o   WavesfindBoundaryNodes.o  WavesReflectedWave.o  WavesScalarEqOpOpt.o   WavesPlaneWaveOpexpl.o WavesBCOperators.o WavesFuncDefs.o WavesOutputs.o  WavesWaveEquationOp.o  WavesWaveEqOpPar.o WavesWaveEqOp123D.o WavesWaveEqOpSeismic.o $(OBJGRID) chkopts
	g++ -o  gca_totalwave  gca_totalwave.o $(test_optional)  WavesWaveEqOpElliptic_log.o  WavesElliptic_GlobConvAlg.o   WavesPoissonEqOp.o  WavesGridFEHier.o   WavesfindBoundaryNodes.o  WavesReflectedWave.o WavesScalarEqOpOpt.o WavesPlaneWaveOpexpl.o   WavesBCOperators.o WavesFuncDefs.o WavesOutputs.o WavesWaveEquationOp.o  WavesWaveEqOpPar.o WavesWaveEqOp123D.o WavesWaveEqOpSeismic.o  $(OBJGRID) $(PETSC_LIB) $(CCFLAGS)

gca_new_totalwave :  gca_new_totalwave.o  $(test_optional)  WavesWaveEqOpElliptic_log.o  WavesElliptic_GlobConvAlg.o   WavesPoissonEqOp.o	  WavesGridFEHier.o   WavesfindBoundaryNodes.o  WavesReflectedWave.o  WavesScalarEqOpOpt.o   WavesPlaneWaveOpexpl.o WavesBCOperators.o WavesFuncDefs.o WavesOutputs.o  WavesWaveEquationOp.o  WavesWaveEqOpPar.o WavesWaveEqOp123D.o WavesWaveEqOpSeismic.o $(OBJGRID) chkopts
	g++ -o  gca_new_totalwave  gca_new_totalwave.o $(test_optional)  WavesWaveEqOpElliptic_log.o  WavesElliptic_GlobConvAlg.o   WavesPoissonEqOp.o  WavesGridFEHier.o   WavesfindBoundaryNodes.o  WavesReflectedWave.o WavesScalarEqOpOpt.o WavesPlaneWaveOpexpl.o   WavesBCOperators.o WavesFuncDefs.o WavesOutputs.o WavesWaveEquationOp.o  WavesWaveEqOpPar.o WavesWaveEqOp123D.o WavesWaveEqOpSeismic.o  $(OBJGRID) $(PETSC_LIB) $(CCFLAGS)

gca_scatwave :  gca_scatwave.o  $(test_optional)  WavesWaveEqOpElliptic_log.o  WavesElliptic_GlobConvAlg.o   WavesPoissonEqOp.o	  WavesGridFEHier.o   WavesfindBoundaryNodes.o  WavesReflectedWave.o  WavesScalarEqOpOpt.o   WavesPlaneWaveOpexpl.o WavesBCOperators.o WavesFuncDefs.o WavesOutputs.o  WavesWaveEquationOp.o  WavesWaveEqOpPar.o WavesWaveEqOp123D.o WavesWaveEqOpSeismic.o $(OBJGRID) chkopts
	g++ -o  gca_scatwave  gca_scatwave.o $(test_optional)  WavesWaveEqOpElliptic_log.o  WavesElliptic_GlobConvAlg.o   WavesPoissonEqOp.o  WavesGridFEHier.o   WavesfindBoundaryNodes.o  WavesReflectedWave.o WavesScalarEqOpOpt.o WavesPlaneWaveOpexpl.o   WavesBCOperators.o WavesFuncDefs.o WavesOutputs.o WavesWaveEquationOp.o  WavesWaveEqOpPar.o WavesWaveEqOp123D.o WavesWaveEqOpSeismic.o  $(OBJGRID) $(PETSC_LIB) $(CCFLAGS)

# gca2_scatwave: the simplified algorithm
gca2_scatwave :  gca2_scatwave.o  $(test_optional)  WavesWaveEqOpElliptic_log.o  WavesElliptic_GlobConvAlg.o   WavesPoissonEqOp.o	  WavesGridFEHier.o   WavesfindBoundaryNodes.o  WavesReflectedWave.o  WavesScalarEqOpOpt.o   WavesPlaneWaveOpexpl.o WavesBCOperators.o WavesFuncDefs.o WavesOutputs.o  WavesWaveEquationOp.o  WavesWaveEqOpPar.o WavesWaveEqOp123D.o WavesWaveEqOpSeismic.o $(OBJGRID) chkopts
	g++ -o  gca2_scatwave  gca2_scatwave.o $(test_optional)  WavesWaveEqOpElliptic_log.o  WavesElliptic_GlobConvAlg.o   WavesPoissonEqOp.o  WavesGridFEHier.o   WavesfindBoundaryNodes.o  WavesReflectedWave.o WavesScalarEqOpOpt.o WavesPlaneWaveOpexpl.o   WavesBCOperators.o WavesFuncDefs.o WavesOutputs.o WavesWaveEquationOp.o  WavesWaveEqOpPar.o WavesWaveEqOp123D.o WavesWaveEqOpSeismic.o  $(OBJGRID) $(PETSC_LIB) $(CCFLAGS)


all:	$(OBJ) $(OBJGRID) chkopts

сlean:
	rm $(OBJ)
