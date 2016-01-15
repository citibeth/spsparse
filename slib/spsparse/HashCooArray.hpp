#ifndef SPSPARSE_VECTOR_COOARRAY_HPP
#define SPSPARSE_VECTOR_COOARRAY_HPP

#include <unordered_map>
#include <spsparse/array.hpp>

namespace spsparse {

template<class IndexT, class ValT, int RANK>
class HashCooArray
{
public:
	static const int rank = RANK;
	typedef IndexT index_type;
	typedef ValT val_type;
	typedef std::array<index_type, rank> indices_type;

	std::array<size_t, RANK> shape;		// Extent of each dimension
	void set_shape(std::array<size_t, RANK> const &_shape) { shape = _shape; }

protected:
	typedef HashCooArray<IndexT, ValT, RANK> ThisHashCooArrayT;
	std::unordered_map<indices_type, val_type> _vals;

public:
	HashCooArray();

	HashCooArray(std::array<size_t, RANK> const &_shape);

	std::unique_ptr<ThisHashCooArrayT> new_blank() const
		{ return std::unique_ptr<ThisHashCooArrayT>(new ThisHashCooArrayT(shape)); }
	ThisHashCooArrayT make_blank() const
		{ return ThisHashCooArrayT(shape); }

	ValT& operator[] (const IndexT& k)
		{ return _vals[k]; }
	ValT& operator[] (IndexT&& k)
		{ return _vals[std::move(k)]; }

#if 0
	// Move semantics
	HashCooArray(HashCooArray &&other);
	void operator=(ThisHashCooArrayT &&other);

	// Copy semantics
	HashCooArray(HashCooArray const &other);
	void operator=(ThisHashCooArrayT const &other);
#endif


	// -------------------------------------------------
	size_t size() const
		{ return _vals.size(); }
	void clear()
		{ _vals.clear(); }
	void reserve(size_t size)
		{}

	// -------------------------------------------------
	// --------------------------------------------------
	/** Standard STL-type iterator for iterating through a VectorSparseMatrix. */
	template <class IterValueT, class IterIndexT, int RANK, class IterValT, class CollectionT>
	class HashIterator : ibmisc::forward_iterator<IterValueT, CollectionT>
	{
		typedef unorderd_map<IterValueT, IterValT>::iterator WrappedT;
	public:
		WrappedIterT wrapped;

		HashIterator(WrappedIterT &&ii) : wrapped(ii) {}

		ValueT &operator*() const
			{ return wrapped->second; }
		SecondIter &operator++()
			{ ++wrapped; return *this; }

		bool operator==(const SecondIter& rhs) const
			{return wrapped == rhs.wrapped;}

	indices_type index()
		{ return wrapped->first; }
	indices_type operator*()
		{ return index(); }

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
	


	typedef CooIterator<const std::array<IndexT, RANK>, const IndexT, RANK, const ValT, const ThisHashCooArrayT> const_iterator;
	typedef CooIterator<std::array<IndexT, RANK>, IndexT, RANK, ValT, ThisHashCooArrayT> iterator;

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
	{}

	void add(std::array<IndexT, RANK> const index, ValT const val);

	blitz::Array<ValT, RANK> to_dense();

	std::ostream &operator<<(std::ostream &out) const;
};



// --------------------------- Method Definitions
template<class IndexT, class ValT, int RANK>
HashCooArray<IndexT, ValT, RANK>::
	HashCooArray() {
		for (int k=0; k<RANK; ++k) shape[k] = 0;	// User must set this later
	}

template<class IndexT, class ValT, int RANK>
HashCooArray<IndexT, ValT, RANK>::
	HashCooArray(std::array<size_t, RANK> const &_shape)
	: shape(_shape) {}

#if 0
template<class IndexT, class ValT, int RANK>
HashCooArray<IndexT, ValT, RANK>::
	HashCooArray(HashCooArray &&other) :
		shape(other.shape),
		index_vecs(std::move(other.index_vecs)),
		val_vec(std::move(other.val_vec)),
		dim_beginnings_set(other.dim_beginnings_set),
		_dim_beginnings(std::move(other._dim_beginnings)),
		edit_mode(other.edit_mode),
		sort_order(other.sort_order) {}

template<class IndexT, class ValT, int RANK>
	void HashCooArray<IndexT, ValT, RANK>::operator=(ThisHashCooArrayT &&other) {
		shape = other.shape;
		index_vecs = std::move(other.index_vecs);
		val_vec = std::move(other.val_vec);
		dim_beginnings_set = other.dim_beginnings_set;
		_dim_beginnings = std::move(other._dim_beginnings);
		edit_mode = other.edit_mode;
		sort_order = other.sort_order;
	}

template<class IndexT, class ValT, int RANK>
HashCooArray<IndexT, ValT, RANK>::
	HashCooArray(HashCooArray const &other) :
		shape(other.shape),
		index_vecs(other.index_vecs),
		val_vec(other.val_vec),
		dim_beginnings_set(other.dim_beginnings_set),
		_dim_beginnings(other._dim_beginnings),
		edit_mode(other.edit_mode),
		sort_order(other.sort_order) {}


template<class IndexT, class ValT, int RANK>
	void HashCooArray<IndexT, ValT, RANK>::operator=(ThisHashCooArrayT const &other) {
		shape = other.shape;
		index_vecs = other.index_vecs;
		val_vec = other.val_vec;
		dim_beginnings_set = other.dim_beginnings_set;
		_dim_beginnings = other._dim_beginnings;
		edit_mode = other.edit_mode;
		sort_order = other.sort_order;
	}

#endif


template<class IndexT, class ValT, int RANK>
	void HashCooArray<IndexT, ValT, RANK>::add(std::array<IndexT, RANK> const index, ValT const val)
	{
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

		auto ii(_vals.find(index));
		if (ii == _vals.end()) {
			_vals.insert(std::make_pair(index, val));
		} else {
			ii->second += val;
		}
	}



// ---------------------------------------------------------------------------
template<class IndexT, class ValT, int RANK>
std::ostream &operator<<(std::ostream &os, spsparse::HashCooArray<IndexT, ValT, RANK> const &A)
	{ return spsparse::_ostream_out_array(os, A); }

template<class IndexT, class ValT>
using HashCooMatrix = HashCooArray<IndexT, ValT, 2>;

template<class IndexT, class ValT>
using HashCooVector = HashCooArray<IndexT, ValT, 1>;


/** @} */

}	// Namespace
#endif	// Guard
