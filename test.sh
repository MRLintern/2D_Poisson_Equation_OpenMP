#!/bin/bash

echo "Running unit tests..."

PASS=0
FAIL=0

# Test 1: RMS of a known array
test_rms() {

    gcc -o test_rms test_rms.c PSolver.c -lm
    ./test_rms > tmp.out
    result=$(cat tmp.out)
    expected="1.000000"

    if [[ "$result" == "$expected" ]]; then

        echo "test_rms passed"
        ((PASS++))

    else

        echo "test_rms failed: expected $expected, got $result"
        ((FAIL++))

    fi

    rm -f test_rms tmp.out
}

# Test 2: Dirichlet BC
test_bc() {

    gcc -o test_bc test_bc.c PSolver.c -lm
    ./test_bc > tmp.out

    if grep -q "PASS" tmp.out; then

        echo "test_bc passed"
        ((PASS++))

    else

        echo "test_bc failed"

        ((FAIL++))

    fi

    rm -f test_bc tmp.out
}

# Test 3: Exact solution at a known point
test_exact() {

    gcc -o test_exact test_exact.c PSolver.c -lm
    ./test_exact > tmp.out
    result=$(cat tmp.out)
    expected="0.707107"  # sin(pi*0.5*0.5) ~ 0.7071

    if [[ "$result" == "$expected" ]]; then

        echo "test_exact passed"

        ((PASS++))

    else

        echo "test_exact failed: expected $expected, got $result"
        ((FAIL++))

    fi

    rm -f test_exact tmp.out
}

# Run tests
test_rms
test_bc
test_exact

# Report summary
echo ""
echo "Test Summary: $PASS passed, $FAIL failed"

if [[ $FAIL -eq 0 ]]; then

    exit 0

else

    exit 1
fi
