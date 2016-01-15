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

#include <ibmisc/iter.hpp>

#include <spsparse/spsparse.hpp>
#include <spsparse/algorithm.hpp>
#include <spsparse/xiter.hpp>
#include <spsparse/accum.hpp>
#include <spsparse/blitz.hpp>

namespace spsparse {

/** @defgroup array array.hpp
@brief Basic sparse arrays

@{
*/


// -----------------------------------------------------
/** @brief Select out just one dimension of the index on iteration.

Wraps CooArray::iterator, changing operator*() to produce produces the
index in just one dimension.

Code Example
@code
CooMatrix<int, double> const A;
typedef DimIndexIter<decltype(A)::const_iterator> DIType;
for (DIType ii(1, A.begin()); ii != DIType(1, A.end()); ++ii)
	printf("Element with column %d and value %d\n", *ii, ii.val());
@endcode

@see spsparse::CooMatrix::dim_iter(), spsparse::CooMatrix::dim_begin(), spsparse::CooMatrix::dim_end()
*/
template<class ValueT, class ValT, class IterT>
class DimIndexIter : public ibmisc::forward_iterator<ValueT, DimIndexIter<ValueT, ValT, IterT>>
{
public:
	IterT wrapped;
	const int dim;

	DimIndexIter(int _dim, IterT const &&ii) : wrapped(ii), dim(_dim) {}

	ValueT operator*()
		{ return wrapped.index(dim); }
	ValT &val()
		{ return wrapped.val(); }
	ValT const val() const
		{ return wrapped.val(); }

	DimIndexIter &operator++()
		{ ++wrapped; return *this; }
	bool operator==(const DimIndexIter& rhs) const
		{return wrapped == rhs.wrapped;}
};
// -----------------------------------------------------
template<class IndicesT, class IterIndexT, int RANK, class IterValT, class CollectionT>
class CooIterator
{
protected:
	CollectionT * const parent;
	int i;
public:
	static const int rank = RANK;
	typedef IndicesT indices_type;
	typedef indices_type value_type;	// Standard STL: Type we get upon operator*()
	typedef IterIndexT index_type;	// Our convention
	typedef IterValT val_type;	// Extension: Type we get upon val()

	CooIterator(CollectionT *p, int _i) : parent(p), i(_i) {}


	indices_type operator[](int n)
		{ return parent->index(i+n); }
	indices_type index()
		{ return parent->index(i); }
	indices_type operator*()
		{ return this->operator[](0); }

	CooIterator &operator+=(int n)
		{i += n; return *this; }
	CooIterator& operator--()
		{ return this->operator+=(-1); }
	CooIterator& operator++()
		{ return this->operator+=(1); }
	CooIterator &operator-=(int n)
		{ return this->operator+=(-n); }

	CooIterator operator+(int n) const
		{ return CooIterator(parent, i+n); }
	bool operator==(CooIterator const &rhs) const
		{ return i == rhs.i; }
	bool operator!=(const CooIterator& rhs) const
		{return !this->operator==(rhs); }


	int offset() const { return i; }
	IterIndexT &index(int k)
		{ return parent->index(k,i); }
	void set_index(indices_type const &idx)
		{ parent->set_index(i, idx); }
	IterValT &val() { return parent->val(i); }
};
// -----------------------------------------------------
template<class IndexT, class ValT, int RANK>
class CooArray
{
public:
	static const int rank = RANK;
	typedef IndexT index_type;
	typedef ValT val_type;
	typedef std::array<index_type, rank> indices_type;

	std::array<size_t, RANK> shape;		// Extent of each dimension
	void set_shape(std::array<size_t, RANK> const &_shape) { shape = _shape; }

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

	CooArray();

	CooArray(std::array<size_t, RANK> const &_shape);

	std::unique_ptr<ThisCooArrayT> new_blank() const
		{ return std::unique_ptr<ThisCooArrayT>(new ThisCooArrayT(shape)); }
	ThisCooArrayT make_blank() const
		{ return ThisCooArrayT(shape); }

	IndexT &index(int dim, size_t ix)
		{ return index_vecs[dim][ix]; }
	IndexT const &index(int dim, size_t ix) const
		{ return index_vecs[dim][ix]; }

	ValT &val(size_t ix)
		{ return val_vec[ix]; }
	ValT const &val(size_t ix) const
		{ return val_vec[ix]; }

	std::array<IndexT, RANK> index(int ix) const {
		std::array<IndexT, RANK> index_ret;
		for (int k=0; k<RANK; ++k) index_ret[k] = index(k, ix);
		return index_ret;
	}
	std::vector<IndexT> index_vec(int ix) const {
		std::vector<IndexT> ret;
		for (int k=0; k<RANK; ++k) ret.push_back(index(k, ix));
		return ret;
	}
	void set_index(int ix, std::array<IndexT, RANK> const &idx)
		{ for (int k=0; k<RANK; ++k) index(k, ix) = idx[k]; }



	blitz::Array<IndexT, 1> indices(int dim) const
		{ return vector_to_blitz(index_vecs[dim]); }
	blitz::Array<ValT, 1> vals() const
		{ return vector_to_blitz(val_vec); }

	// Move semantics
	CooArray(CooArray &&other);
	void operator=(ThisCooArrayT &&other);

	// Copy semantics
	CooArray(CooArray const &other);
	void operator=(ThisCooArrayT const &other);


	// -------------------------------------------------
	size_t size() const
		{ return val_vec.size(); }
	void clear();
	void reserve(size_t size);

	// -------------------------------------------------
	// --------------------------------------------------
	/** Standard STL-type iterator for iterating through a VectorSparseMatrix. */

	typedef CooIterator<const std::array<IndexT, RANK>, const IndexT, RANK, const ValT, const ThisCooArrayT> const_iterator;
	typedef CooIterator<std::array<IndexT, RANK>, IndexT, RANK, ValT, ThisCooArrayT> iterator;

	iterator begin(int ix = 0)
		{ return iterator(this, ix); }
	iterator end(int ix = 0)
		{ return iterator(this, size() + ix); }
	const_iterator cbegin(int ix = 0) const
		{ return const_iterator(this, ix); }
	const_iterator cend(int ix = 0) const
		{ return const_iterator(this, size() + ix); }
	const_iterator begin(int ix = 0) const
		{ return const_iterator(this, ix); }
	const_iterator end(int ix = 0) const
		{ return const_iterator(this, size() - ix); }

	// typedef DimIndexIter<IndexT, ValT, iterator> dim_iterator;
	typedef DimIndexIter<const IndexT, const ValT, const_iterator> const_dim_iterator;

	const_dim_iterator dim_iter(int dim, int ix) const
		{ return const_dim_iterator(dim, const_iterator(this, ix)); }
	const_dim_iterator dim_begin(int dim) const
		{ return dim_iter(dim, 0); }
	const_dim_iterator dim_end(int dim) const
		{ return dim_iter(dim, size()); }
	// -------------------------------------------------
	/** Goes in to add mode: legal to add more things to the vector. */
	void edit()
	{
		edit_mode = true;
		sort_order[0] = -1;
	}

	void add(std::array<IndexT, RANK> const index, ValT const val);

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
		bool handle_nan = false);

	void transpose(std::array<int, RANK> const &sort_order)
	{
		OverwriteAccum<iterator> overwrite(begin());
		spsparse::transpose(overwrite, *this, sort_order);
	}

	blitz::Array<ValT, RANK> to_dense();

	// Sets and returns this->_dim_beginnings
	std::vector<size_t> const &dim_beginnings() const;

	DimBeginningsXiter<ThisCooArrayT> dim_beginnings_xiter() const;

	std::ostream &operator<<(std::ostream &out) const;
};



// --------------------------- Method Definitions
template<class IndexT, class ValT, int RANK>
CooArray<IndexT, ValT, RANK>::
	CooArray() : edit_mode(true), dim_beginnings_set(false), sort_order() {
		sort_order[0] = -1;
		for (int k=0; k<RANK; ++k) shape[k] = 0;	// User must set this later
	}

template<class IndexT, class ValT, int RANK>
CooArray<IndexT, ValT, RANK>::
	CooArray(std::array<size_t, RANK> const &_shape)
	: shape(_shape), edit_mode(true), dim_beginnings_set(false), sort_order() {
		sort_order[0] = -1;
	}

template<class IndexT, class ValT, int RANK>
CooArray<IndexT, ValT, RANK>::
	CooArray(CooArray &&other) :
		shape(other.shape),
		index_vecs(std::move(other.index_vecs)),
		val_vec(std::move(other.val_vec)),
		dim_beginnings_set(other.dim_beginnings_set),
		_dim_beginnings(std::move(other._dim_beginnings)),
		edit_mode(other.edit_mode),
		sort_order(other.sort_order) {}

template<class IndexT, class ValT, int RANK>
	void CooArray<IndexT, ValT, RANK>::operator=(ThisCooArrayT &&other) {
		shape = other.shape;
		index_vecs = std::move(other.index_vecs);
		val_vec = std::move(other.val_vec);
		dim_beginnings_set = other.dim_beginnings_set;
		_dim_beginnings = std::move(other._dim_beginnings);
		edit_mode = other.edit_mode;
		sort_order = other.sort_order;
	}

template<class IndexT, class ValT, int RANK>
CooArray<IndexT, ValT, RANK>::
	CooArray(CooArray const &other) :
		shape(other.shape),
		index_vecs(other.index_vecs),
		val_vec(other.val_vec),
		dim_beginnings_set(other.dim_beginnings_set),
		_dim_beginnings(other._dim_beginnings),
		edit_mode(other.edit_mode),
		sort_order(other.sort_order) {}

template<class IndexT, class ValT, int RANK>
	void CooArray<IndexT, ValT, RANK>::operator=(ThisCooArrayT const &other) {
		shape = other.shape;
		index_vecs = other.index_vecs;
		val_vec = other.val_vec;
		dim_beginnings_set = other.dim_beginnings_set;
		_dim_beginnings = other._dim_beginnings;
		edit_mode = other.edit_mode;
		sort_order = other.sort_order;
	}

template<class IndexT, class ValT, int RANK>
void CooArray<IndexT, ValT, RANK>::clear() {
		for (int k=0; k<RANK; ++k) index_vecs[k].clear();
		val_vec.clear();
		dim_beginnings_set = false;
		_dim_beginnings.clear();
		edit_mode = true;
		sort_order[0] = -1;
	}

template<class IndexT, class ValT, int RANK>
void CooArray<IndexT, ValT, RANK>::reserve(size_t size) {
		for (int k=0; k<RANK; ++k) index_vecs[k].reserve(size);
		val_vec.reserve(size);
	}

template<class IndexT, class ValT, int RANK>
	void CooArray<IndexT, ValT, RANK>::add(std::array<IndexT, RANK> const index, ValT const val)
	{
		if (!edit_mode) {
			(*spsparse_error)(-1, "Must be in edit mode to use CooArray::add()");
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
				(*spsparse_error)(-1, buf.str().c_str());
			}
		}

		for (int i=0; i<RANK; ++i) index_vecs[i].push_back(index[i]);
		val_vec.push_back(val);
	}

template<class IndexT, class ValT, int RANK>
void CooArray<IndexT, ValT, RANK>::consolidate(
		std::array<int, RANK> const &_sort_order,
		DuplicatePolicy duplicate_policy,
		bool handle_nan)
	{
		// Do nothing if we're already properly consolidated
		if (this->sort_order == _sort_order && !edit_mode) return;

		ThisCooArrayT ret(shape);
		spsparse::consolidate(ret, *this, _sort_order, duplicate_policy, handle_nan);
		*this = std::move(ret);
	}

template<class IndexT, class ValT, int RANK>
	blitz::Array<ValT, RANK> CooArray<IndexT, ValT, RANK>::to_dense()
	{
		blitz::Array<ValT, RANK> ret(array_to_tiny<int,size_t,rank>(shape));
		ret = 0;
		DenseAccum<ThisCooArrayT> accum(ret);
		copy(accum, *this);
		return ret;
	}

template<class IndexT, class ValT, int RANK>
	// Sets and returns this->_dim_beginnings
	std::vector<size_t> const &CooArray<IndexT, ValT, RANK>::dim_beginnings() const
	{
		// See if we need to compute it; lazy eval
		if (!dim_beginnings_set) {
			// Const cast OK here for lazy eval implementation
			ThisCooArrayT *vthis = const_cast<ThisCooArrayT *>(this);
			vthis->_dim_beginnings = spsparse::dim_beginnings(*this);
			vthis->dim_beginnings_set = true;
		}
		return _dim_beginnings;
	}

template<class IndexT, class ValT, int RANK>
	DimBeginningsXiter <CooArray<IndexT, ValT, RANK>> CooArray<IndexT, ValT, RANK>::dim_beginnings_xiter() const
	{
		auto &db(dim_beginnings());
		int const index_dim = sort_order[0];
		int const val_dim = sort_order[1];
		return DimBeginningsXiter<ThisCooArrayT>(this, index_dim, val_dim, db.begin(), db.end());
	}


// ---------------------------------------------------------------------------
template<class IndexT, class ValT>
using CooMatrix = CooArray<IndexT, ValT, 2>;

template<class IndexT, class ValT>
using CooVector = CooArray<IndexT, ValT, 1>;



/** @} */

}	// Namespace


// ---------------------------------------------------------------------------
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
// ---------------------------------------------------------------------------


#endif	// Guard
