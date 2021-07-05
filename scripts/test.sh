set -e

test_file="F_document/qq.losslin.all.K6.png"

export TEST_HATBAG_DIR=$1

# Test test_file exists
if test -f "${TEST_HATBAG_DIR}/${test_file}"; then
    echo "test passed; ${test_file} exists"
else
    echo "test failed; ${test_file} does not exist"
    exit 1
fi

# Test significant k-mers
Rscript ./R/test_HATBAG.R
