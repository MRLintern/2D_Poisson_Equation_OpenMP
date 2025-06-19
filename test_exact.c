#include <stdio.h>
#include "PSolver.h"

int main() {

    double val = potential_exact_function(0.5, 0.5);

    printf("%.6f\n", val);
}
