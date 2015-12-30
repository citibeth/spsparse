#ifndef SPSPARSE_ALGORITHM_HPP
#define SPSPARSE_ALGORITHM_HPP

#include <spsparse/spsparse.hpp>
#include <spsparse/xiter.hpp>

namespace spsparse {

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


// -----------------------------------------------------------
template<class CooArrayT, class AccumulatorT>
void copy(AccumulatorT &ret, CooArrayT const &A)
{
	std::array<int,CooArrayT::rank> idx;
	for (auto ii=A.begin(); ii != A.end(); ++ii) {
		ret.add(ii.index(), ii.val());
	}
}
// -----------------------------------------------------------
template<class CooArrayT, class AccumulatorT>
void transpose(AccumulatorT &ret, CooArrayT const &A, std::array<int,CooArrayT::rank> const &perm)
{
	std::array<int,CooArrayT::rank> idx;
	for (auto ii=A.begin(); ii != A.end(); ++ii) {
		for (int new_k=0; new_k < CooArrayT::rank; ++new_k) {
			int old_k = perm[new_k];
			idx[new_k] = ii.index(old_k);
		}
		ret.add(idx, ii.val());
	}
}
// -----------------------------------------------------------
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
	SPSPARSE_LOCAL_TYPES(CooArrayT);
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

	index_type operator*()
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
// --------------------------------------------------------
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


// ----------------------------------------------------------

}

#endif // Guard
