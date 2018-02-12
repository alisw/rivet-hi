#! /usr/bin/env bash

#set -x

export LD_LIBRARY_PATH=$LIBLOCATION:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$LIBLOCATION:$DYLD_LIBRARY_PATH

echo PYTHONPATH=$PYTHONPATH
echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH
echo PATH=$PATH
echo
echo RIVET_ANALYSIS_PATH=$RIVET_ANALYSIS_PATH
echo RIVET_DATA_PATH=$RIVET_DATA_PATH
echo RIVET_REF_PATH=$RIVET_REF_PATH
echo RIVET_INFO_PATH=$RIVET_INFO_PATH
echo
echo PYTHON=$PYTHON

function _clean() {
    rm -f fifo.hepmc
    rm -f file2.hepmc
}

function _setup() {
    _clean
    cp ${RIVET_TESTS_SRC}/testApi.hepmc file2.hepmc
    mkfifo fifo.hepmc
}

function _check() {
    CODE=$?
    if [[ $CODE -ne 0 ]]; then
        _clean
        _exit $CODE
    fi
}
echo "trying to load rivet python module"
$PYTHON -c 'import rivet'  || exit $?
echo "Success"

_setup

echo
rivet --list-analyses > log || exit $?

# this analysis has greek chars in the name
echo
rivet --show-analysis SLD_1999_S37439 > log || exit $?

echo
rivet -a D0_2008_S7554427 ${RIVET_TESTS_SRC}/testApi.hepmc file2.hepmc > log || exit $?
grep -q "20 events" log
_check

echo
cat ${RIVET_TESTS_SRC}/testApi.hepmc | rivet -a D0_2008_S7554427 > log || exit $?
grep -q "10 events" log
_check

echo
cat ${RIVET_TESTS_SRC}/testApi.hepmc > fifo.hepmc &
rivet -a D0_2008_S7554427 fifo.hepmc > log || exit $?
grep -q "10 events" log
_check
_clean
