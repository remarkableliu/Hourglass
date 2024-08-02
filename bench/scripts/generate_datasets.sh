#! /bin/bash

#
# This file is part of Grafite <https://github.com/marcocosta97/grafite>.
# Copyright (C) 2023 Marco Costa.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters, usage: generate_datasets.sh <build_path> <real_datasets_path>"
fi

BUILD_PATH=$(realpath $1)
if [ ! -d "$BUILD_PATH" ]; then
  echo "Grafite build path does not exist"
  exit 1
fi
REAL_DATASETS_PATH=$(realpath $2)
if [ ! -d "$REAL_DATASETS_PATH" ]; then
  echo "Real datasets path does not exist"
  exit 1
fi

WORKLOAD_GEN_PATH=$(realpath $BUILD_PATH/bench/workload_gen)
if [ ! -f "$WORKLOAD_GEN_PATH" ]; then
  echo "Workload generator does not exist"
  exit 1
fi

generate_corr_test() {
  i=0
  x=0.0

  while [ $i -le 10 ]
  do
    $WORKLOAD_GEN_PATH --kdist kuniform --qdist qcorrelated --corr-degree ${x}
    mv kuniform/ kuniform_${i}/
    x=$(echo $x + 0.1 | bc)
    i=$(($i + 1))
  done
}

generate_adapt_test() {
  i=0  
  x=0

  while [ $i -le 4 ]
  do
    $WORKLOAD_GEN_PATH --kdist kuniform --qdist qadapt-uniform --uq-ratio ${x}
    mv kuniform/ kuniform_${x}/
    $WORKLOAD_GEN_PATH --kdist kuniform --qdist qadapt-corr --uq-ratio ${x}
    mv kuniform/ kuniform_qcorr${x}/
    x=$(echo $x + 5 | bc)
    i=$(($i + 1))
  done
}

generate_constr_time_test() {
  x=15938355

  $WORKLOAD_GEN_PATH --kdist kuniform --qdist quniform --range-size 10 -n ${x} -q $(echo "($x * 0.0001)/1" | bc)
  mv kuniform/ kuniform_quniform/

  $WORKLOAD_GEN_PATH --kdist knormal --qdist quniform --range-size 10 -n ${x} -q $(echo "($x * 0.0001)/1" | bc)
  mv knormal/ knormal_quniform/

  $WORKLOAD_GEN_PATH --kdist kuniform --qdist qcorrelated --range-size 10 -n ${x} -q $(echo "($x * 0.0001)/1" | bc)
  mv kuniform/ kuniform_qcorrelated/

  $WORKLOAD_GEN_PATH --kdist knormal --qdist qcorrelated --range-size 10 -n ${x} -q $(echo "($x * 0.0001)/1" | bc)
  mv knormal/ knormal_qcorrelated/

  done
}

generate_query_time_test() {
  x=15938355

  $WORKLOAD_GEN_PATH --kdist kuniform --qdist quniform --range-size 10 -n ${x} -q $(echo "($x * 0.05)/1" | bc)
  mv kuniform/ kuniform_quniform/

  $WORKLOAD_GEN_PATH --kdist knormal --qdist qcorrelated --range-size 10 -n ${x} -q $(echo "($x * 0.05)/1" | bc)
  mv kuniform/ kuniform_qcorrelated/

  done
}

mkdir -p ../adapt_test && cd ../adapt_test || exit 1
if ! generate_adapt_test ; then
  echo "[!!] adapt_test generation failed"
  exit 1
fi
echo "[!!] adapt_test dataset generated"

mkdir -p ../corr_test && cd ../corr_test || exit 1
if ! generate_corr_test ; then
  echo "[!!] corr_test generation failed"
  exit 1
fi
echo "[!!] corr_test dataset generated"

mkdir -p ../fpr_test && cd ../fpr_test || exit 1
if ! $WORKLOAD_GEN_PATH --kdist knormal; then
  echo "[!!] fpr_test generation failed"
  exit 1
fi
if ! $WORKLOAD_GEN_PATH --kdist kuniform; then
  echo "[!!] fpr_test generation failed"
  exit 1
fi
echo "[!!] fpr_test (figure 3) dataset generated"

mkdir -p ../fpr_real_test && cd ../fpr_real_test || exit 1
if ! $WORKLOAD_GEN_PATH --binary-keys $REAL_DATASETS_PATH/books_200M_uint64  $REAL_DATASETS_PATH/osm_cellids_200M_uint64 ; then
  echo "[!!] fpr_real_test generation failed"
  exit 1
fi
echo "[!!] fpr_real_test dataset generated"

mkdir -p ../constr_time_test && cd ../constr_time_test || exit 1
if ! generate_constr_time_test ; then
  echo "[!!] constr_time_test generation failed"
  exit 1
fi
echo "[!!] constr_time_test dataset generated"

mkdir -p ../query_time_test && cd ../query_time_test || exit 1
if ! generate_query_time_test ; then
  echo "[!!] query_time_test generation failed"
  exit 1
fi
echo "[!!] query_time_test dataset generated"

mkdir -p ../true_test && cd ../true_test || exit 1
if ! $WORKLOAD_GEN_PATH --qdist qtrue ; then
  echo "[!!] true_test generation failed"
  exit 1
fi
echo "[!!] true_test dataset generated"

echo "[!!] success, all datasets generated"