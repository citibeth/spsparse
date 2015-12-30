#ifndef SPSPARSE_ARRAY_HPP
#define SPSPARSE_ARRAY_HPP

#include <array>
#include <cstdio>
#include <functional>
#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <spsparse/spsparse.hpp>
#include <spsparse/algorithm.hpp>
#include <spsparse/xiter.hpp>
#include <spsparse/accum.hpp>

namespace spsparse {


// -----------------------------------------------------
template<class IndexT, class ValT, int RANK>
class CooArray
{
public:
	static const int rank = RANK;
	typedef IndexT index_type;
	typedef ValT val_type;
	typedef std::array<index_type, rank> indices_type;

	std::array<size_t, RANK> const shape;		// Extent of each dimension

protected:
	typedef CooArray<IndexT, ValT, RANK> ThisCooArrayT;
	std::array<std::vector<IndexT>, RANK> index_vecs;
	std::vector<ValT> val_vec;

	// OPTIONAL
	// If sort_order[0] != -1, this is the index of the beginning of
	// each different value of dimension sort_order[0] (i.e. beginning
	// of each row, if we're sorted row-major).
	bool dim_beginnings_set;
	std::vector<size_t> _dim_beginnings;

public:
	bool edit_mode;		// Are we in edit mode?
	std::array<int,RANK> sort_order;	// Non-negative elements if this is sorted



	CooArray(std::array<size_t, RANK> const &_shape)
	: shape(_shape), edit_mode(true), dim_beginnings_set(false), sort_order() {
		sort_order[0] = -1;
	}

	IndexT &index(int dim, size_t ix)
		{ return index_vecs[dim][ix]; }
	IndexT index(int dim, size_t ix) const
		{ return index_vecs[dim][ix]; }

	ValT &val(size_t ix)
		{ return val_vec[ix]; }
	ValT val(size_t ix) const
		{ return val_vec[ix]; }

	std::array<IndexT, RANK> index(int ix) const {
		std::array<IndexT, RANK> ret;
		for (int k=0; k<RANK; ++k) ret[k] = index(k, ix);
		return ret;
	}
	std::vector<IndexT> index_vec(int ix) const {
		std::vector<IndexT> ret;
		for (int k=0; k<RANK; ++k) ret.push_back(index(k, ix));
		return ret;
	}
	void set_index(int ix, std::array<IndexT, RANK> const &idx)
		{ for (int k=0; k<RANK; ++k) index(k, ix) = idx[k]; }



	blitz::Array<IndexT, 1> indices(int dim)
		{ return vector_to_blitz(index_vecs[dim]); }
	blitz::Array<ValT, 1> vals()
		{ return vector_to_blitz(val_vec); }

	// Move semantics
	CooArray(CooArray &&other) :
		shape(other.shape),
		index_vecs(std::move(other.index_vecs)),
		val_vec(std::move(other.val_vec)),
		dim_beginnings_set(other.dim_beginnings_set),
		_dim_beginnings(std::move(other._dim_beginnings)),
		edit_mode(other.edit_mode),
		sort_order(other.sort_order) {}

	void operator=(ThisCooArrayT &&other) {
		if (shape != other.shape) {
			(*sparse_error)(-1, "Cannot change shape on assignment.");
		}
		index_vecs = std::move(other.index_vecs);
		val_vec = std::move(other.val_vec);
		dim_beginnings_set = other.dim_beginnings_set;
		_dim_beginnings = std::move(other._dim_beginnings);
		edit_mode = other.edit_mode;
		sort_order = other.sort_order;
	}

	// Copy semantics
	CooArray(CooArray const &other) :
		shape(other.shape),
		index_vecs(other.index_vecs),
		val_vec(other.val_vec),
		dim_beginnings_set(other.dim_beginnings_set),
		_dim_beginnings(other._dim_beginnings),
		edit_mode(other.edit_mode),
		sort_order(other.sort_order) {}

	void operator=(ThisCooArrayT const &other) {
		if (shape != other.shape) {
			(*sparse_error)(-1, "Cannot change shape on assignment.");
		}
		index_vecs = other.index_vecs;
		val_vec = other.val_vec;
		dim_beginnings_set = other.dim_beginnings_set;
		_dim_beginnings = other._dim_beginnings;
		edit_mode = other.edit_mode;
		sort_order = other.sort_order;
	}


	// -------------------------------------------------
	size_t size() const
		{ return val_vec.size(); }

	void clear() {
		for (int k=0; k<RANK; ++k) index_vecs[k].clear();
		val_vec.clear();
		dim_beginnings_set = false;
		_dim_beginnings.clear();
		edit_mode = true;
	}

	void reserve(size_t size) {
		for (int k=0; k<RANK; ++k) index_vecs[k].reserve(size);
		val_vec.reserve(size);
	}

	// -------------------------------------------------
	// --------------------------------------------------
	/** Standard STL-type iterator for iterating through a VectorSparseMatrix. */
	class iterator {
	protected:
		ThisCooArrayT * const parent;
		size_t i;
	public:
		static const int rank = RANK;

		typedef IndexT value_type;	// Standard STL: Type we get upon operator*()
		typedef IndexT index_type;	// Our convention
		typedef ValT val_type;	// Extension: Type we get upon val()

		iterator(ThisCooArrayT *p, size_t _i) : parent(p), i(_i) {}

		size_t offset() const { return i; }
		bool operator==(iterator const &rhs) const { return i == rhs.i; }
		bool operator!=(iterator const &rhs) const { return i != rhs.i; }
		bool operator<(iterator const &rhs) const { return i < rhs.i; }
		void operator++() { ++i; }
		IndexT &index(int k) { return parent->index(k,i); }
		IndexT index(int k) const { return parent->index(k,i); }
		void set_index(std::array<IndexT, RANK> const &idx)
			{ parent->set_index(i, idx); }
		std::array<IndexT, RANK> operator*() const { return index(); }
		ValT &val() { return parent->val(i); }
		ValT val() const { return parent->val(i); }


		std::array<IndexT, RANK> index() const
			{ return parent->index(i); }
		std::vector<IndexT> index_vec() const
			{ return parent->index_vec(i); }
	};

	// -------------------------------------------------
	iterator iter(size_t ix) { return iterator(this, ix); }
	iterator begin() { return iter(0); }
	iterator end() { return iter(size()); }

	typedef DimIndexIter<iterator> dim_iterator;
	dim_iterator dim_iter(int dim, size_t ix)
		{ return dim_iterator(dim, iter(ix)); }
	dim_iterator dim_begin(int dim)
		{ return dim_iter(dim, 0); }
	dim_iterator dim_end(int dim)
		{ return dim_iter(dim, size()); }


	// -------------------------------------------------
	class const_iterator {
	protected:
		ThisCooArrayT const * const parent;
		size_t i;
	public:
		static const int rank = RANK;
		typedef IndexT value_type;	// Standard STL: Type we get upon operator*()
		typedef IndexT index_type;	// Our convention
		typedef ValT val_type;	// Extension: Type we get upon val()

		const_iterator(ThisCooArrayT const *p, size_t _i) : parent(p), i(_i) {}

		size_t offset() const { return i; }
		bool operator==(const_iterator const &rhs) const { return i == rhs.i; }
		bool operator!=(const_iterator const &rhs) const { return i != rhs.i; }
		bool operator<(const_iterator const &rhs) const { return i < rhs.i; }
		void operator++() { ++i; }
		IndexT index(int k) const { return parent->index(k,i); }
		ValT val() const { return parent->val(i); }
		std::array<IndexT, RANK> operator*() const { return index(); }


		std::array<IndexT, RANK> index() const
			{ return parent->index(i); }
		std::vector<IndexT> index_vec() const
			{ return parent->index_vec(i); }

	};

	// -------------------------------------------------

	const_iterator iter(size_t ix) const { return const_iterator(this, ix); }
	const_iterator begin() const { return iter(0); }
	const_iterator end() const { return iter(size()); }

	typedef DimIndexIter<const_iterator> const_dim_iterator;
	const_dim_iterator dim_iter(int dim, size_t ix) const
		{ return const_dim_iterator(dim, iter(ix)); }
	const_dim_iterator dim_begin(int dim) const
		{ return dim_iter(dim, 0); }
	const_dim_iterator dim_end(int dim) const
		{ return dim_iter(dim, size()); }

	// -------------------------------------------------





	// -------------------------------------------------
	/** Goes in to add mode: legal to add more things to the vector. */
	void edit()
	{
		edit_mode = true;
	}

	void add(std::array<IndexT, RANK> const index, ValT const val)
	{
		if (!edit_mode) {
			(*sparse_error)(-1, "Must be in edit mode to use CooArray::add()");
		}

		// Check bounds
		for (int i=0; i<RANK; ++i) {
			if (index[i] < 0 || index[i] >= shape[i]) {
				std::ostringstream buf;
				buf << "Sparse index out of bounds: index=(";
				for (int j=0; j<RANK; ++j) {
					buf << index[j];
					buf << " ";
				}
				buf << ") vs. shape=(";
				for (int j=0; j<RANK; ++j) {
					buf << shape[j];
					buf << " ";
				}
				buf << ")";
				(*sparse_error)(-1, buf.str().c_str());
			}
		}

		for (int i=0; i<RANK; ++i) index_vecs[i].push_back(index[i]);
		val_vec.push_back(val);
	}

	/** Mark that this is now in sorted form. */
	void set_sorted(std::array<int,RANK> _sort_order)
	{
		sort_order = _sort_order;
		edit_mode = false;
	}

	// --------------------------------------------------
	// In-place algos
	void consolidate(
		std::array<int, RANK> const &_sort_order,
		DuplicatePolicy duplicate_policy = DuplicatePolicy::ADD,
		bool handle_nan = false)
	{
		// Do nothing if we're already properly consolidated
		if (this->sort_order == _sort_order && !edit_mode) return;

		ThisCooArrayT ret(shape);
		spsparse::consolidate(ret, *this, _sort_order, duplicate_policy, handle_nan);
		*this = std::move(ret);
	}

	void transpose(std::array<int, RANK> const &sort_order)
	{
		OverwriteAccum<iterator> overwrite(begin());
		spsparse::transpose(overwrite, *this, sort_order);
	}

	blitz::Array<ValT, RANK> to_dense()
	{
		blitz::Array<ValT, RANK> ret(0);
		DenseAccum<ThisCooArrayT> accum(ret);
		copy(accum, *this);
		return ret;
	}

	std::vector<size_t> &dim_beginnings()
	{
		// See if we need to compute it; lazy eval
		if (!dim_beginnings_set) {
			this->_dim_beginnings = spsparse::dim_beginnings(*this);
			dim_beginnings_set = true;
		}
		return _dim_beginnings;
	}

	DimBeginningsXiter<ThisCooArrayT> dim_beginnings_xiter()
	{
		auto &db(dim_beginnings());
		int const index_dim = sort_order[0];
		int const val_dim = sort_order[1];
		return DimBeginningsXiter<ThisCooArrayT>(this, index_dim, val_dim, db.begin(), db.end());
	}

};




#if 0
// These need to move outside of the CooArray class.
// Temporarily comment out NetCDF stuff.
// Not sure if it belongs in the core class.

	 void netcdf_define(
		netCDF::NcFile &nc, std::string const &vname,
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


	void netcdf_write(NcFile &nc, std::string const &vname) const
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
#endif


template<class IndexT, class ValT>
using CooMatrix = CooArray<IndexT, ValT, 2>;

template<class IndexT, class ValT>
using CooVector = CooArray<IndexT, ValT, 1>;



}	// Namespace



template<class IndexT, class ValT, int RANK>
std::ostream &operator<<(std::ostream &os, spsparse::CooArray<IndexT, ValT, RANK> const &A);

template<class IndexT, class ValT, int RANK>
std::ostream &operator<<(std::ostream &os, spsparse::CooArray<IndexT, ValT, RANK> const &A)
{
	os << "CooArray<";
	stream(os, &A.shape[0], A.shape.size());
	os << ">(";
	for (auto ii(A.begin()); ii != A.end(); ++ii) {
		os << "(";
		auto idx(ii.index());
		for (int k=0; k<RANK; ++k) os << idx[k] << " ";
		os << ": " << ii.val() << ")";
	}
	os << ")";
	return os;
}


#endif	// Guard
