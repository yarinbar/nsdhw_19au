
#include "matrix.hpp"

Matrix multiply(const Matrix& a, const Matrix& b){

    // check for dim compatibility
    if(a.ncol() != b.nrow())
        throw "sizes dont match!";

    Matrix c(a.nrow(), b.ncol());

    for(unsigned int i = 0; i < a.ncol(); ++i)
        for(unsigned int j = 0; j < a.ncol(); ++j) {
            c(i, j) = 0;
            for(unsigned int k = 0; k < a.ncol(); k++)
                c(i, j) += a(i, k) * b(k, j);
        }
    return c;
}

PYBIND11_MODULE(multiply, mod){
    mod.doc() = "multiply";
    mod.def("multiply", &multiply, "");
}
