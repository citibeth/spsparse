#ifndef SPSPARSE_NETCDF_HPP
#define SPSPARSE_NETCDF_HPP

#include <spsparse/array.hpp>

namespace spsparse {

/** @defgroup netcdf netcdf.hpp
@brief Simple way to read/write SpSparse arrays via NetCDF.

@{
*/

// These need to move outside of the CooArray class.
// Temporarily comment out NetCDF stuff.
// Not sure if it belongs in the core class.

template<class ArrayT>
 void netcdf_define(
	netCDF::NcFile &nc, std::string const &vname,
	ArrayT const &A,
	std::vector<std::function<void ()>> &writes) const
{
	NcDim size_d = nc.addDim(vname + ".size", this->size());
	NcDim rank_d = nc.addDim(vname + ".rank", RANK);
	nc.add_var(vname + ".indices", ncInt, {size_d, rank_d});
	nc.add_var(vname + ".vals", ncDouble, {size_d});

	one_d = getOrAddDim(nc, "one", 1);
	auto descrVar = nc.add_var(vname + ".descr", ncInt, {one_d});	// TODO: This should be ".info"
	descrVar.putAtt("shape", ncLong, RANK, &shape);

	writes.push_back(&CooArray<IndexT,ValT,RANK>::netcdf_write,
		this, &nc, vname);
}


void netcdf_write(
	netCDF::NcFile &nc, std::string const &vname, const
	ArrayT const &A)
{
	NcVar indices_v = nc.getVar(vname + ".index");
	NcVar vals_v = nc.getVar(vname + ".val");

	vals_v.putVar(val_vec, size());

	std::vector<size_t> startp(2) = {0, 0};		// SIZE, RANK
	std::vector<size_t> countp(2) = {1, RANK};	// Write RANK elements at a time
	for (auto ov = this->begin(); ov != this->end(); ++ov) {
		std::array<size_t> index = index();

		indices_v.putVar(startp, countp, &index[0]);	// 2-D
		vals_v.putVar(startp, countp, &ov.val());	// 1-D

		++startp[0];
	}
}

/** @} */

}	// Namespace
#endif	// Guard
