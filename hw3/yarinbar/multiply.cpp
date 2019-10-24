
#include "matrix.hpp"

Matrix multiply(const Matrix a, const Matrix b){

    // check for dim compatibility
    if(a.ncol != b.nrows)
        return nullptr;

    Matrix c(a.nrows(), b.ncols());

    for(unsigned int i = 0; i < a.ncol(); ++i)
        for(unsigned int j = 0; j < a.ncol(); ++j) {
            c(i, j) = 0;
            for(unsigned int k = 0; k < a.ncol(); k++)
                c(i, j) += a(i, k) * b(k, j);
        }

}

PYBIND11_MODULE(multiply, mod){
    mod.def("multiply", &multiply, "");
}
