source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup cmake v3_9_0
setup gcc v6_4_0
setup eigen v3_3_5
setup root v6_28_10b -q e20:p3915:prof
setup python v2_7_3

export PYTHONPATH=${PYTHONPATH}:${PWD}/lib/

which gcc
which root
which python
# Optional
#setup root v6_12_06a -q e15:prof
#setup root v6_28_10b -q e20:p3915:prof
