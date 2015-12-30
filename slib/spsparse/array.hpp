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
#include <spsparse/xiter.hpp>
#include <spsparse/blitz.hpp>

namespace spsparse {

enum class DuplicatePolicy {LEAVE_ALONE, ADD, REPLACE};

// -------------------------------------------------------------
// Values for sort_order formal parameter below
//const std::array<int,2> ROW_MAJOR = {0,1};
//const std::array<int,2> COLUMN_MAJOR = {1,0};



// ====================================================
#if 0
/** Used as an alternate underlying storage for CooArray */
template<class T>
class BlitzCooVec
{
	blitz::Array<T,1> vec;
	size_t size;

	void clear()
		{ size = 0; }
	void push_back(T const &val)
	{
		if (CHECK_BOUNDS && size >= val.extent(0)) {
		  (*sparse_error)(-1, "")
		vec[size++] = val;
	}

};
#endif
// ====================================================
/** Wrapps CooArray::iterator, currying one dimension
so that operator*() produces the index in that dimension. */
template<class IterT>
class DimIndexIter : public WrapForwardValIter<IterT>
{
public:
	typedef typename WrapForwardValIter<IterT>::value_type value_type;
	const int dim;

	typedef WrapForwardValIter<IterT> WrapT;
	// Export as is standard with STL

	DimIndexIter(int _dim, IterT const &_ii) : WrapT(_ii), dim(_dim) {}
	DimIndexIter(int _dim, IterT const &&_ii) : WrapT(std::move(_ii)), dim(_dim) {}

	value_type operator*()
		{ return this->ii.index(dim); }
};

// -----------------------------------------------------
/** A must be sorted properly.
See: CooArray::consolidate() */
template<class CooArrayT>
std::vector<size_t> dim_beginnings(CooArrayT const &A)
{
	const int RANK = CooArrayT::rank;

	std::vector<size_t> abegin;

	// Check that we're sorted by SOME dimension.
	if (A.sort_order[0] < 0) {
		(*sparse_error)(-1, "dim_beginnings() required the CooArray is sorted first.");
	}

	// Get beginning of each row in a (including sentinel at end)
	auto ai(A.begin());
	auto const end(A.end());
	if (ai != end) {		// At least 1 element
		abegin.push_back(ai.offset());			// First item in array is always 0
		int const dim = A.sort_order[0];		// Dimension we're sectioning by
		int last_row = ai.index(dim);

		for (++ai; ; ++ai) {
			if (ai == end) {
				// Add a sential row
				abegin.push_back(ai.offset());
				break;
			}
			if (ai.index(dim) != last_row) {
				// We see a new row!
				abegin.push_back(ai.offset());
				last_row = ai.index(dim);
			}
		}
	}

	return abegin;
}
// ----------------------------------------------------------
/** Iterates over one row/col at a time (in general, one index in a dimension).
The underlying array must be sorted by that dimension. */
template<class CooArrayT>
class DimBeginningsXiter : public STLXiter<std::vector<size_t>::iterator>
{
public:

	typedef typename CooArrayT::index_type IndexT;
	typedef typename CooArrayT::val_type ValueT;
	int const RANK = CooArrayT::rank;
	typedef std::vector<size_t>::iterator DimIterT;

protected:
	CooArrayT *arr;
	int index_dim;	// Dimension corresponding to our "rows"
	int val_dim;	// Dimension corresponding to our "cols"

public:
	DimBeginningsXiter(
		CooArrayT *_arr,
		int _index_dim, int _val_dim,
		DimIterT const &dim_beginnings_begin,
		DimIterT const &dim_beginnings_end)
	: STLXiter<std::vector<size_t>::iterator>(dim_beginnings_begin, dim_beginnings_end),
		arr(_arr), index_dim(_index_dim), val_dim(_val_dim)
	{}

	IndexT operator*()
		{ return arr->index(index_dim, *ii); }

	/** Returns an Xiter along the current row/col in the CooArray. */
	typedef ValSTLXiter<typename CooArrayT::dim_iterator> SubXiterT;
	SubXiterT sub_xiter(int _val_dim = -1)
	{
		if (_val_dim < 0) _val_dim = val_dim;
		size_t a0 = *ii;
		size_t a1 = *(ii + 1);
		return SubXiterT(arr->dim_iter(_val_dim, a0), arr->dim_iter(_val_dim, a1));
	}

	// No val()
};
// -----------------------------------------------------
template<class CooArrayT, class AccumulatorT>
void consolidate(AccumulatorT &ret,
	CooArrayT const &A,
	std::array<int, CooArrayT::rank> const &sort_order,
	DuplicatePolicy duplicate_policy = DuplicatePolicy::ADD,
	bool handle_nan = false)	// TODO: This is currently ignored!
{
	const int RANK = CooArrayT::rank;
	typedef typename CooArrayT::index_type IndexT;
	typedef typename CooArrayT::val_type ValT;

	// Nothing to do for zero-size matrices (or ones with just one element)
	if (A.size() > 1) {

		// Get a sorted permutation
		std::vector<size_t> perm(sorted_permutation(A, sort_order));

		// Scan through to identify duplicates
		auto ii(perm.begin());

		// Skip over initial 0 and NaN
		for (;; ++ii) {
			if (ii == perm.end()) goto finished;		// Nothing more to do
			if (!std::isnan(A.val(*ii)) && A.val(*ii) != 0) break;
		}

		// New temporary entry
		std::array<IndexT,RANK> accum_idx(A.index(*ii));
		ValT accum_val = A.val(*ii);
		++ii;

		for (; ; ++ii) {
			// Skip over initial 0 and NaN
			for (;; ++ii) {
				if (ii == perm.end()) {
					// Write out the last thing we had in our accumulator
					ret.add(accum_idx, accum_val);

					goto finished;		// Nothing more to do
				}
				if (!std::isnan(A.val(*ii)) && A.val(*ii) != 0) break;
			}

			// Test if A.index(*ii) == accum_idx
			auto new_idx = A.index(*ii);
			for (int k=0; k<RANK; ++k) {
				if (new_idx[k] != accum_idx[k]) {
					// They don't match.  Write out our accumulator, and reset to this one
					ret.add(accum_idx, accum_val);
					accum_idx = new_idx;
					accum_val = A.val(*ii);
					goto continue_outer;
				}
			}

			// They do match!  Add it into the accumulator...
			if (duplicate_policy == DuplicatePolicy::ADD)
				accum_val += A.val(*ii);
			else if (duplicate_policy == DuplicatePolicy::REPLACE)
				accum_val = A.val(*ii);

		continue_outer: ;
		}
	}
finished:


	ret.set_sorted(sort_order);
}
// -----------------------------------------------------
template<class IndexT, class ValT, int RANK>
class CooArray
{
public:
	static const int rank = RANK;
	typedef IndexT index_type;
	typedef ValT val_type;

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
		typedef IndexT value_type;	// Standard STL: Type we get upon operator*()
		typedef ValT val_type;	// Extension: Type we get upon val()

		iterator(ThisCooArrayT *p, size_t _i) : parent(p), i(_i) {}

		size_t offset() const { return i; }
		bool operator==(iterator const &rhs) const { return i == rhs.i; }
		bool operator!=(iterator const &rhs) const { return i != rhs.i; }
		bool operator<(iterator const &rhs) const { return i < rhs.i; }
		void operator++() { ++i; }
		IndexT &index(int k) { return parent->index(k,i); }
		IndexT index(int k) const { return parent->index(k,i); }
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
		typedef IndexT value_type;	// Standard STL: Type we get upon operator*()
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
		{ return dim_iterator(dim, size()); }

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
		std::array<int, RANK> const &sort_order,
		DuplicatePolicy duplicate_policy = DuplicatePolicy::ADD,
		bool handle_nan = false)
	{
		ThisCooArrayT ret(shape);
		spsparse::consolidate(ret, *this, sort_order, duplicate_policy, handle_nan);
		*this = std::move(ret);
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


	friend std::ostream &::operator<<(std::ostream &os, CooArray<IndexT, ValT, RANK> const &A);

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

// -----------------------------------------------------------
template<class CooArrayT, class AccumulatorT>
void transpose(AccumulatorT &ret, CooArrayT const &A, std::array<int,CooArrayT::rank> perm)
{
	std::array<int,CooArrayT::rank> idx;
	for (auto ii=A.begin(); ii != A.end(); ++ii) {
		for (int new_k=0; new_k < CooArrayT::rank; ++new_k) {
			int old_k = perm[new_k];
			idx[new_k] = ii.index(old_k);
		}
		ret.add(idx, A.val());
	}
}
// -----------------------------------------------------------
template<class CooArrayT>
struct CmpIndex {
	const int RANK = CooArrayT::rank;

	CooArrayT const *arr;
	std::array<int, CooArrayT::rank> sort_order;

	CmpIndex(CooArrayT const *_arr,
		std::array<int, CooArrayT::rank> const &_sort_order)
	: arr(_arr), sort_order(_sort_order) {}

	bool operator()(int i, int j)
	{
		for (int k=0; k<RANK-1; ++k) {
			int dim = sort_order[k];
			if (arr->index(dim,i) < arr->index(dim,j)) return true;
			if (arr->index(dim,i) > arr->index(dim,j)) return false;
		}
		int dim = sort_order[RANK-1];
		return arr->index(dim,i) < arr->index(dim,j);
	}
};

template<class CooArrayT>
std::vector<size_t> sorted_permutation(CooArrayT const &A,
	std::array<int, CooArrayT::rank> const &sort_order)
{
	// Decide on how we'll sort
	CmpIndex<CooArrayT> cmp(&A, sort_order);

	// Generate a permuatation
	int n = A.size();
	std::vector<size_t> perm; perm.reserve(n);
	for (int i=0; i<n; ++i) perm.push_back(i);

	// Sort the permutation
	std::stable_sort(perm.begin(), perm.end(), cmp);

	return perm;
}


// ---------------------------------------------

#if 0
template<class CooArrayT>
DimBeginningsXiter<CooArrayT> dim_beginnings_xiter(
	CooArrayT *_arr,
	int _index_dim, int _val_dim,
	std::vector<size_t>::iterator const &dim_beginnings_begin,
	std::vector<size_t>::iterator const &dim_beginnings_end)
{
	return DimBeginningsXiter<CooArrayT>(_arr, _index_dim, _val_dim, dim_beginnings_begin, dim_beginnings_end);
}
#endif
// ----------------------------------------------------------
// ----------------------------------------------------------



}	// Namespace

template<class IndexT, class ValT, int RANK>
std::ostream &operator<<(std::ostream &os, spsparse::CooArray<IndexT, ValT, RANK> const &A)
{
	os << "CooArray<" << RANK << ">(";
	for (auto ii(A.begin()); ii != A.end(); ++ii) {
		os << "(";
		auto idx(ii.index());
		for (int k=0; k<RANK; ++k) os << idx[k] << " ";
		os << ": " << ii.val() << ")";
	}
	os << ")";
}


#endif	// Guard
