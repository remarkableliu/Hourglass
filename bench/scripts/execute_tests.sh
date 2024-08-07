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
    echo "Illegal number of parameters, usage: execute_tests.sh <build_path> <datasets_path>"
fi

BUILD_PATH=$(realpath $1)
WORKLOADS_PATH=$(realpath $2)


# ARGS="--numa --membind 0 --physcpubind 16" # uncomment to run with numa 
ARGS="" # uncomment to run without numa

SCRIPT_DIR_PATH=$(dirname -- "$( readlink -f -- "$0"; )")


if ! python3 $SCRIPT_DIR_PATH/test.py $ARGS --test corr $WORKLOADS_PATH/corr_test $BUILD_PATH ; then
  echo "[!!] corr_test test failed"
  exit 1
fi
echo "[!!] corr_test test executed successfully"

if ! python3 $SCRIPT_DIR_PATH/test.py $ARGS --test adapt $WORKLOADS_PATH/adapt_test $BUILD_PATH ; then
  echo "[!!] adapt_test test failed"
  exit 1
fi
echo "[!!] adapt_test test executed successfully"

if ! python3 $SCRIPT_DIR_PATH/test.py $ARGS --test fpr $WORKLOADS_PATH/fpr_test $BUILD_PATH ; then
  echo "[!!] fpr_test test failed"
  exit 1
fi
echo "[!!] fpr_test test executed successfully"

if ! python3 $SCRIPT_DIR_PATH/test.py $ARGS --test fpr_real $WORKLOADS_PATH/fpr_real_test $BUILD_PATH ; then
  echo "[!!] fpr_real_test test failed"
  exit 1
fi
echo "[!!] fpr_real_test test executed successfully"

if ! python3 $SCRIPT_DIR_PATH/test.py $ARGS --test true $WORKLOADS_PATH/true_test $BUILD_PATH ; then
  echo "[!!] query_time_test test failed"
  exit 1
fi
echo "[!!] query_time_test test executed successfully"

if ! python3 $SCRIPT_DIR_PATH/test.py $ARGS --test constr_time $WORKLOADS_PATH/constr_time_test $BUILD_PATH ; then
  echo "[!!] constr_time_test test failed"
  exit 1
fi
echo "[!!] constr_time_test test executed successfully"

echo "[!!] success, all tests executed"
