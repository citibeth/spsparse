#ifndef SPSPARSE_ARRAY_HPP
#define SPSPARSE_ARRAY_HPP

#include <array>
#include <cstdio>
#include <functional>
#include <vector>
#include <cmath>
#include <string>
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
	int dim;

public:
	typedef WrapForwardValIter<IterT> WrapT;

	DimIndexIter(IterT const &_ii, int _dim) : WrapT(_ii), dim(_dim) {}
	DimIndexIter(IterT const &&_ii, int _dim) : WrapT(std::move(_ii)), dim(_dim) {}

	auto operator*() -> decltype(WrapForwardValIter<IterT>::ii.index(dim))
		{ return WrapForwardValIter<IterT>::ii.index(dim); }
};

// -----------------------------------------------------

template<class IndexT, class ValT, int _RANK>
class CooArray
{
public:
	static const int RANK = _RANK;
	static const int rank = RANK;
	typedef IndexT index_type;
	typedef ValT val_type;

	std::array<size_t, RANK> const shape;		// Extent of each dimension


protected:
	typedef CooArray<IndexT, ValT, RANK> ThisCooArrayT;

	std::array<std::vector<IndexT>, RANK> index_vecs;
	std::vector<ValT> val_vec;

	bool in_edit;		// Are we in edit mode?
	std::array<int,RANK> sort_order;	// Non-negative elements if this is sorted

public:
	CooArray(std::array<size_t, RANK> const &_shape)
	: shape(_shape), in_edit(true), sort_order() {
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
printf("ret=%d %d\n", ret[0], ret[1]);
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
		in_edit(other.in_edit),
		sort_order(other.sort_order) {}

	void operator=(ThisCooArrayT &&other) {
		shape = other.shape;
		index_vecs = std::move(other.index_vecs);
		val_vec = std::move(other.val_vec);
		in_edit = other.in_edit;
		sort_order = other.sort_order;
	}

	// Copy semantics
	CooArray(CooArray const &other) :
		shape(other.shape),
		index_vecs(other.index_vecs),
		val_vec(other.val_vec),
		in_edit(other.in_edit),
		sort_order(other.sort_order) {}

	void operator=(ThisCooArrayT const &other) {
		shape = other.shape;
		index_vecs = other.index_vecs;
		val_vec = other.val_vec;
		in_edit = other.in_edit;
		sort_order = other.sort_order;
	}


	// -------------------------------------------------
	size_t size() const
		{ return val_vec.size(); }

	void clear() {
		for (int k=0; k<RANK; ++k) index_vecs[k].clear();
		val_vec.clear();
	}

	void reserve(size_t size) {
		for (int k=0; k<RANK; ++k) index_vecs[k].reserve(size);
		val_vec.reserve(size);
	}

#if 0
Not needed, duplicated below.
	/** Adds an element to the sparse array. */
	void add(std::array<IndexT,RANK> const index, ValT const &val)
	{
		for (int k=0; k<RANK; ++k) index_vecs[k].push_back(index[k]);
		val_vec.push_back(val);
	}
	void add_v(std::vector<IndexT> const &index, ValT const &val)
	{
		for (int k=0; k<RANK; ++k) index_vecs[k].push_back(index[k]);
		val_vec.push_back(val);
	}
#endif


	// -------------------------------------------------
	// --------------------------------------------------
	/** Standard STL-type iterator for iterating through a VectorSparseMatrix. */
	class iterator {
	protected:
		ThisCooArrayT *parent;
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


		std::array<IndexT, RANK> index() const
			{ return parent->index(i); }
		std::vector<IndexT> index_vec() const
			{ return parent->index_vec(i); }

		// Convenience methods
		IndexT &row() const { return index(0); }
		IndexT &col() const { return index(1); }
		ValT &value() { return val(); }
	};

	iterator iter(size_t ix) { return iterator(this, ix); }
	iterator begin() { return iterator(this, 0); }
	iterator end() { return iterator(this, val_vec.size()); }

	typedef DimIndexIter<iterator> dim_iterator;
	iterator dim_iter(size_t ix, int dim)
		{ return dim_iterator(iter(ix), dim); }
	iterator dim_begin(int dim)
		{ return dim_iterator(begin(), dim); }
	iterator dim_end(int dim)
		{ return dim_iterator(end(), dim); }

	// --------------------------------------------------
#if 0
// No const_iterator for now, until we're more settled. 

	class const_iterator {
	protected:
		ThisCooArrayT const *parent;
		size_t i;
	public:
		const_iterator(ThisCooArrayT *p, size_t _i) : parent(p), i(_i) {}
		size_t offset() const { return i; }
		bool operator==(iterator const &rhs) const { return i == rhs.i; }
		bool operator!=(iterator const &rhs) const { return i != rhs.i; }
		bool operator<(iterator const &rhs) const { return i < rhs.i; }
		void operator++() { ++i; }
		std::array<IndexT, RANK> operator*() const { return index(); }
		ValT const &val() { return parent->val(i); }

		IndexT const &index(int k) const { return parent->index(k,i); }
		std::array<IndexT, RANK> index() const {
			std::array<size_t, RANK> ret;
			for (int k=0; k<RANK; ++k) ret[k] = index(k);
			return ret;
		}
		std::vector<IndexT> index_vec() const {
			std::vector<size_t> ret;
			for (int k=0; k<RANK; ++k) ret.push_back(index(k));
			return ret;
		}


		// Convenience methods
		IndexT const &row() const { return index(0); }
		IndexT const &col() const { return index(1); }
		ValT const &value() { return val(); }
	};
	iterator operator[](size_t ix) const { return const_iterator(this, ix); }
	iterator begin() const { return const_iterator(this, 0); }
	iterator end() const { return const_iterator(this, val_vec.size()); }
#endif
	// --------------------------------------------------

	// -------------------------------------------------
	/** Goes in to add mode: legal to add more things to the vector. */
	void edit()
	{
		in_edit = true;
	}

	void add(std::array<IndexT, RANK> const index, ValT const val)
	{

		if (!in_edit) {
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
		in_edit = false;
	}


#if 0
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
};


template<class IndexT, class ValT>
using CooMatrix = CooArray<IndexT, ValT, 2>;

template<class IndexT, class ValT>
using CooVector = CooArray<IndexT, ValT, 1>;






template<class CooArrayT>
struct CmpIndex {
	const int RANK = CooArrayT::RANK;

	CooArrayT const *arr;
	int const *sort_order;

	CmpIndex(CooArrayT const *_arr, int const *_sort_order) :
	arr(_arr), sort_order(_sort_order) {}

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

template<class CooArrayT, class AccumulatorT>
void consolidate(AccumulatorT &ret,
	CooArrayT const A,
	std::array<int, CooArrayT::rank> const &_sort_order,
	DuplicatePolicy duplicate_policy = DuplicatePolicy::ADD,
	bool handle_nan = false)
{
	const int RANK = CooArrayT::rank;
	typedef typename CooArrayT::index_type IndexT;
	typedef typename CooArrayT::val_type ValT;

	// Nothing to do for zero-size matrices (or ones with just one element)
	if (A.size() > 1) {

		// Decide on how we'll sort
		CmpIndex<CooArrayT> cmp(&A, &_sort_order[0]);

		// Generate a permuatation
		int n = A.size();
		std::vector<int> perm; perm.reserve(n);
		for (int i=0; i<n; ++i) perm.push_back(i);
		std::stable_sort(perm.begin(), perm.end(), cmp);

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
			for (;;) {
				if (ii == perm.end()) goto finished;		// Nothing more to do
				if (!std::isnan(A.val(*ii)) && A.val(*ii) != 0) break;
			}

			// Test if A.index(*ii) == accum_idx
			auto new_idx = A.index(*ii);
			for (int k=0; k<RANK; ++k) {
				if (new_idx[k] != accum_idx[k]) {
					// They don't match.  Add our accumulator, and reset to this one
					ret.add(accum_idx, accum_val);
					accum_idx = new_idx;
					accum_val = A.val(*ii);
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
	ret.set_sorted(_sort_order);
}
// ---------------------------------------------



/** A must be sorted properly.
See: CooArray::consolidate() */
template<class CooArrayT, class AccumulatorT>
std::vector<size_t> dim_beginnings(CooArrayT const &A)
{
	const int RANK = CooArrayT::RANK;

	std::vector<int> abegin;

	// Check that we're sorted by SOME dimension.
	if (A.sort_order[0] < 0) {
		(*sparse_error)(-1, "CooArray::get_dim_beginnings() required the CooArray is sorted first.");
	}

	// Get beginning of each row in a (including sentinel at end)
	auto ai(A.begin());
	auto const end(A.end());
	if (ai != end) {		// At least 1 element
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














}	// Namespace

#endif	// Guard
