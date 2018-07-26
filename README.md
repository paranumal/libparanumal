# holmes
It's elementary.

# 0. Code block diagram
<img src="http://www.math.vt.edu/people/tcew/libParanumalDiagramLocal-crop-V2.png" width="600" >

## 1. Clone: Holmes
git clone https://github.com/tcew/holmes

## 2. OCCA dependency (currently OCCA 1.0 forked by Noel Chalmers)
git clone https://github.com/noelchalmers/occa

## 3. Build OCCA
cd occa

export OCCA_DIR=``pwd``

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OCCA_DIR/lib

make -j

cd ../

## 4. Build holmes elliptic example
cd holmes/solvers/elliptic

make -j

## 5-1. Run elliptic example with provided quadrilateral set up file on a single device:
./ellipticMain setups/setupQuad2D.rc

## 5-2. Run the same example with two devices:
mpiexec -n 2 ./ellipticMain setups/setupQuad2D.rc
