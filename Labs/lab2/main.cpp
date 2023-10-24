#include <iostream>

#include <Matrix.hpp>

using namespace dml;

#define OUT(a) std::cout << a << std::endl

int main() {

    Matrix matrix;
    matrix.set({0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1});

    OUT(matrix.isReflectivity());
    OUT(matrix.isAntiReflexivity());
    OUT(matrix.isSymmetry());
    OUT(matrix.isAntiSymmetry());
    OUT(matrix.isAsymmetry());
    OUT(matrix.isTransitivity());

    OUT(matrix.getTransitiveClosureMatrix());
    OUT(matrix.getTransitiveClosureMatrixWarshall());

    OUT(matrix.getSymmetricClosureMatrix().getPower() + matrix.getTransitiveClosureMatrix().getPower() +
        matrix.getReflexiveClosureMatrix().getPower());


    return 0;
}