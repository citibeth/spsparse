#ifndef SPSPARSE_ACCUM_HPP
#define SPSPARSE_ACCUM_HPP

#include <cstddef>
#include <array>
#include <vector>
#include <spsparse/spsparse.hpp>

namespace spsparse {


// -----------------------------------------------------------
/** Accumulates by overwriting a CooArray::iterator.  Good for
in-place operations that we KNOW won't change the size of the
CooArray (eg: use with transpose()). */
template<class IterT>
struct OverwriteAccum
{
	SPSPARSE_LOCAL_TYPES(IterT);

	IterT ii;
	OverwriteAccum(IterT &&_ii) : ii(std::move(_ii)) {}

	void add(indices_type const &index, typename IterT::val_type const &val) {
		ii.set_index(index);
		ii.val() = val;
		++ii;
	}
};
// -----------------------------------------------------------

template<int IN_RANK, class AccumulatorT>
struct PermuteAccum
{
	static const int rank = IN_RANK;
	static const int out_rank = AccumulatorT::rank;

	AccumulatorT sub;
	std::vector<int> perm;
	std::array<int, out_rank> out_idx;

	PermuteAccum(AccumulatorT &&_sub, std::array<int, out_rank> const &_perm)
		: sub(std::move(_sub)), perm(_perm) {}

	void add(std::array<int, IN_RANK> const &index, typename AccumulatorT::val_type const &val) {
		for (int i=0; i<IN_RANK; ++i) out_idx[i] = index[perm[i]];
		sub.add(out_idx, val);
	}
};
// -----------------------------------------------------------
template<class CooArrayT>
struct DenseAccum
{
	SPSPARSE_LOCAL_TYPES(CooArrayT);

	DuplicatePolicy duplicate_policy;
	blitz_type dense;
	blitz::TinyVector<index_type, rank> bidx;

	DenseAccum(blitz_type &_dense, DuplicatePolicy _duplicate_policy=DuplicatePolicy::ADD) : dense(_dense), duplicate_policy(_duplicate_policy) {}

	void add(indices_type const &index, val_type const &val)
	{
		array_to_tiny<index_type, rank>(bidx, index);

		val_type &oval(dense(bidx));

		switch(duplicate_policy) {
			case DuplicatePolicy::LEAVE_ALONE :
				if (!std::isnan(oval)) oval = val;
			break;
			case DuplicatePolicy::ADD :
				oval += val;
			break;
			case DuplicatePolicy::REPLACE :
				oval = val;
			break;
		}
	}
};
// -----------------------------------------------------------


// -------------------------------------------------------
template<class IndexT, class ValT, int RANK>
struct ScalarAccumulator {
	ValT val;
	ScalarAccumulator() : val(0) {}

	void add(const std::array<IndexT, RANK> &index, ValT const &_val)
		{ this->val += _val; }
};






}	// Namespace

#endif // Guard
