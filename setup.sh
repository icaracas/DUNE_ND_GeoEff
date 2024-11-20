source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup cmake v3_9_0
setup gcc v6_4_0
setup eigen v3_3_5
setup python v2_7_3

export PYTHONPATH=${PYTHONPATH}:${PWD}/lib/

# Optional
setup root v6_12_06a -q e15:prof

