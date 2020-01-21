LIB_SRC=/opt/openfoam7/src
FOAM_USER_LIBBIN=/home/octavio/OpenFOAM/octavio-7/platforms/linux64GccDPInt32Opt/lib
gcc -fdirectives-only -E porousScalarTransportFoam.C -I$LIB_SRC/finiteVolume/lnInclude
    -I$LIB_SRC/meshTools/lnInclude \
    -I$LIB_SRC/thermophysicalModels/specie/lnInclude \
    -I$LIB_SRC/thermophysicalModels/basic/lnInclude \
    -I$LIB_SRC/thermophysicalModels/reactionThermo/lnInclude \
    -I../../libraries/porousBoundaryConditions/lnInclude \
    -I../../libraries/porousModels/lnInclude \
    -I../../libraries/reactionModels/lnInclude \
    -I../../libraries/toolsGIS/lnInclude \
    -I../../libraries/numericalMethods/lnInclude
    -lfiniteVolume \
    -lspecie \
    -lfluidThermophysicalModels \
    -lreactionThermophysicalModels \
    -L$FOAM_USER_LIBBIN \
    -lPMFnumericalMethods \
    -lporousBoundaryConditions \
    -lporousModels \
    -lreactionModels \
    -ltoolsGIS
