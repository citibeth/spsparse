#ifndef SPSPARSE_NETCDF_HPP
#define SPSPARSE_NETCDF_HPP

#include <functional>
#include <ibmisc/ncutil.hpp>
#include <spsparse/array.hpp>

namespace spsparse {

/** @defgroup netcdf netcdf.hpp
@brief Simple way to read/write SpSparse arrays via NetCDF.

@{
*/

template<class ArrayT>
void nc_write(
	netCDF::NcFile *nc,
	std::string const &vname,
	ArrayT const *A);

template<class ArrayT>
void nc_write(
	netCDF::NcFile *nc,
	std::string const &vname,
	ArrayT const *A)
{
	netCDF::NcVar indices_v = nc->getVar(vname + ".indices");
	netCDF::NcVar vals_v = nc->getVar(vname + ".vals");

	std::vector<size_t> startp0 = {0};
	std::vector<size_t> countp0 = {A->size()};
	blitz::Array<typename ArrayT::val_type, 1> vals(A->vals());
	vals_v.putVar(startp0, countp0, &vals(0));

	std::vector<size_t> startp = {0, 0};		// SIZE, RANK
	std::vector<size_t> countp = {1, A->rank};	// Write RANK elements at a time
	for (auto ov = A->begin(); ov != A->end(); ++ov) {
		typename ArrayT::indices_type index = ov.index();

		indices_v.putVar(startp, countp, &index[0]);	// 2-D
		typename ArrayT::val_type val = ov.val();
		vals_v.putVar(startp, countp, &val);	// 1-D

		++startp[0];
	}
}

template<class ArrayT>
 void nc_define(
	netCDF::NcFile &nc, std::string const &vname,
	ArrayT const &A,
	ibmisc::NcWrites &writes);

template<class ArrayT>
 void nc_define(
	netCDF::NcFile &nc,
	std::string const &vname,
	ArrayT const &A,
	ibmisc::NcWrites &writes)
{
	netCDF::NcDim size_d = nc.addDim(vname + ".size", A.size());
	netCDF::NcDim rank_d = nc.addDim(vname + ".rank", A.rank);
	nc.addVar(vname + ".indices", netCDF::ncInt, {size_d, rank_d});
	nc.addVar(vname + ".vals", netCDF::ncDouble, {size_d});

	auto one_d = ibmisc::getOrAddDim(nc, "one", 1);
	auto infoVar = nc.addVar(vname + ".info", netCDF::ncInt, {one_d});
	infoVar.putAtt("shape", netCDF::ncInt64, A.rank, &A.shape[0]);

	writes += std::bind(&nc_write<ArrayT>, &nc, vname, &A);
}



/** @} */

}	// Namespace
#endif	// Guard
