#include <Kokkos_Core.hpp>
#include "flcl-cxx.hpp"


const Scalar SC_ONE = Kokkos::ArithTraits<Scalar>::one();


extern "C" {


void c_stream_sl_kokkos(flcl_ndarray_t *nd_array_y, flcl_ndarray_t *nd_array_x) {

	using flcl::view_from_ndarray;

	auto y = view_from_ndarray<double*>(*nd_array_y);
	auto x = view_from_ndarray<double*>(*nd_array_x);

	const double alpha = 1.0;
	const double beta  = 0.0;  


	KokkosSparse::spmv("N",alpha,mat,x,beta,y)

}

void c_kokkos_create_csr(void** A, int rows, int cols, int* idx_rows, int* idx_cols, double* values){



}

} // extern "C"