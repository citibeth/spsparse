#ifndef SPSPARSE_MULTIPLY_SPARSE_HPP
#define SPSPARSE_MULTIPLY_SPARSE_HPP

#include <spsparse/xiter.hpp>
#include <spsparse/array.hpp>

namespace spsparse {

// ----------------------------------------------------------
/** Computes [row] diag(scale) [[col]].
Or: row_j scale_j col_j   (out of overall matrix computation row_ij scale_j col_jk

NOTE: row, col and scale must NOT have duplicate elements.  Consolidate() must have been run before --- which gets ride of nans too, via handle_nan.
*/
template<class RowXiterT, class ScaleXiterT, class ColXiterT, class AccumulatorT>
static void multiply_row_scale_col(
	AccumulatorT &ret,

	RowXiterT &&rowi,			// DimBeginningsXiter
	ScaleXiterT &&scalei,		// Regular Xiter: val() returns value
	ColXiterT &&coli)			// DimBeginningsXiter
{
	for (Join3Xiter<RowXiterT, ScaleXiterT, ColXiterT> join_xiter(
		std::move(rowi), std::move(scalei), std::move(coli));
	!join_xiter.eof(); ++join_xiter) {

		auto term = join_xiter.i1.val() * join_xiter.i2.val() * join_xiter.i3.val();
		ret.add({*join_xiter.i1}, term);	// 1-D accumulator
	}
}

#if 0

// -------------------------------------------------------------
/** NOTES:
A must be consolidated ROW_MAJOR, B must be consolidated COLUMN_MAJOR
*/
template<class IndexT, class ValueT, class AccumulatorT>
static void multiply(
	AccumulatorT &ret,
	double C,		// Multiply everything by this
	CooVector<IndexT, ValueT> const &scalei,
	CooMatrix<IndexT, ValueT> const &A,
	CooVector<IndexT, ValueT> const &scalej,
	CooMatrix<IndexT, ValueT> const &B,
	CooVector<IndexT, ValueT> const &scalek)
{
Compute this in consolidate...?
Check A and B are consolidated properly...
Check dimensions...?
	if (C == 0) return;		// Everything will be zero, so give up now!

	// Get beginnings of each row in A and each col in B
	std::vector<size_t> abegin(get_rowcol_beginnings(A, 0));
	std::vector<size_t> bbegin(get_rowcol_beginnings(B, 1));

	typedef ValSTLXiter<typename CooVector<IndexT, ValueT>::dim_iterator> ScaleXiterT;
	typedef DimBeginningsXiter<CooMatrix<IndexT, ValueT>> MatrixXiterT;


	// Multiply each row by each column
	// ------ Outer loop: Rows in A
	Join2Xiter<ScaleXiterT, MatrixXiterT> join_a(
		ScaleXiterT(scalei.dim_begin(0), scalei.dim_end(0)),
		MatrixXiterT(&A, 0, 1, abegin.begin(), abegin.end()));

	ScaleXiter &scalei_ii(join_a.i1);
	MatrixXiter &a_ii(join_a.i2);

	for (; !join_a.eof(); ++join_a) {
		double sival = scalei_ii.val();
		if (sival == 0) continue;	// shortcut if zero in the vector


		// -------- Inner loop: Cols in B
		Join2Xiter<MatrixXiterT, ScaleXiterT> join_b(
			MatrixXiterT(&B, 1, 0, bbegin.begin(), bbegin.end()),
			ScaleXiterT(scalek.dim_begin(0), scalek.dim_end(0)));
		ScaleXiter &b_ii(join_b.i1);
		MatrixXiter &scalek_ii(join_b.i2);

		for (; !join_b.eof(); ++join_b) {
			double skval = scalek.vals[scalek_ii.offset()];
			if (skval == 0) continue;	// shortcut if zero in the vector

			// Multiply A row by B column
			multiply_row_scale_col(rcval,
				a_ii.sub_xiter(),
				scalej.dim_iter(0, 0),
				b_ii.sub_xiter());

			double val = C * sival * rcval.value() * skval;
			if (val == 0) continue;

			// *a_ii == index of current row in A
			ret.add({*a_ii, *b_ii}, val);

		}
	}
}

#endif

#if 0


/** NOTES:
A must be consolidated ROW_MAJOR, B must be consolidated COLUMN_MAJOR
*/
template<class IndexT, class ValueT, class AccumulatorT>
static double multiply(
	AccumulatorT &ret,
	CooVector<IndexT, ValueT> const *scalei,
	CooMatrix<IndexT, ValueT> const &A,
	CooVector<IndexT, ValueT> const &b)
{
	std::vector<size_t> abegin(get_rowcol_beginnings(A, 0));
	std::vector<size_t> bbegin(get_rowcol_beginnings(B, 1));

	// Multiply each row by each column

	// ---------- Outer loop: rows in A
	size_t ai = 0;
	size_t sii = 0;		// Index into sparse scalei representation
	IndexT next_match_a = A.indices[0][ai];
	if (scalei) next_match_a = std::max(next_match_a, scalei->indices[0][sii]);
	for (;; ++ai, ++sii) {
		// ---- Increment forward to the next row in A w/ matching scalei
		if (!scalei) {
			if (ai >= abegin.size()-1) goto break_outer;
		} else {
			for (;;++ai) {
				if (ai >= abegin.size()-1) goto break_outer;
				if (A.indices[0][ai] >= netxt_match_a) break;
			}
			next_match_a = A.indices[0][ai];

			for (;;++sii) {
				if (si >= scalei->size()) goto break_outer;
				if (scalei->indices[0][sii] >= next_match_a) break;
			}
			next_match_a = scalei->indices[0][sii];
		}

		ValueT sival;
		if (!scalei) {
			sival = 1.0;
		} else {
			sival = scalei->indices[0][sii];
			if (sival == 0) continue;
		}

		// ------------ Just one column of b

		// Multiply a row by a column
		IndexT const a0 = abegin[ai];
		ScalarAccumulator<IndexT, ValueT, RANK> rcval;
		multiply_row_col(rcval,
			&A.indices[0][a0], &A.vals[a0], abegin[ai+1]-a0,
			0, 0, 0, 0,
			&b.indices[0][0], &b.vals[0], b.size())
			handle_nan);

		double val = sival * rcval.value();
		if (val == 0) continue;

		IndexT aRow = a.indices[0][a0];
		ret.add({{aRow}}, val);

	}
	break_outer: ;

}






template<class IndexT, class ValueT, class AccumulatorT>
static double multiply(
	AccumulatorT &ret,
	CooVector<IndexT, ValueT> const &A,
	CooVector<IndexT, ValueT> const &B,
	bool handle_nan = false)
{
	multiply_row_col(ret,
		&a.indices[0][0], &a.vals[0], a.size(),
		0, 0, 0, 0,
		&b.indices[0][0], &b.vals[0], b.size())
		handle_nan);
}



/** Copy a to b while transposing.  Does not clear b. */
template<class CooMatrixT>
CooMatrixT transpose(
	AccumulatorT &ret,
	CooMatrixT &A)
{
	CooMatrixT ret;
	for (auto ii = A.begin(); ii != A.end(); ++ii) {
		ret.add({ii.col(), ii.row()}, ii.val(), dups);
	}
	return ret;
}

template<class AccumulatorT, class IndexT, class ValueT, int RANK>
void invert(
	AccumulatorT &ret,
	CooArray<IndexT, ValueT, RANK> &A)
{
	size_t size = A.size();
	ret.reserve(size);
	for (auto ii = A.begin(); ii != A.end(); ++ii) {
		std::array<ValueT, RANK> index;
		for (int k=0; i<RANK; ++k) index[k] = A.indices[k][i];
		ret.add(index, A.vals[i]);
	}
}


template<class AccumulatorT, class IndexT, class ValueT, int RANK>
void invert(
	AccumulatorT &ret,
	CooArray<IndexT, ValueT, RANK> &A)


// --------------------------------------------------------------


/** Accumulates into a dense array */
template<class IndexT, class ValueT, int RANK, int DUPLICATE_POLICY>
class DenseAccumulator : public BlankAccumulator {
	blitz::Array<ValueT, RANK> &dense;
	DuplicatePolicy duplicate_policy;
public:
	int rank() { return RANK; }
	DenseAccumulator(blitz::Array<ValueT, RANK> &_dense, DuplicatePolicy _duplicate_policy)
		: dense(_dense), duplicate_policy(_duplicate_policy) {}


	void add(TinyVector<int,RANK> &index, ValueT const &val)
	{
		ValueT &oval(dense(bindex));

		switch(DUPLICATE_POLICY) {
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

	void add(std::array<IndexT, RANK> const index, ValueT const &val)
	{
		TinyVector<int,RANK> bindex;
		for (int k=0; k<RANK; ++k) bindex(k) = index[k];
		add(bindex, val);
	}

	void add(std::vector<IndexT> const index, ValueT const &val)
	{
		TinyVector<int,RANK> bindex;
		for (int k=0; k<RANK; ++k) bindex(k) = index[k];
		add(bindex, val);
	}


	blitz::Array<ValueT, RANK> &value()
		{ return dense; }
};

template<class IndexT, class ValueT, class LesserAccumulatorT>
class CollapseAccumulator : public BlankAccumulator
{
	LesserAccumulatorT sub;
	std::vector<int> keep_dims;		// Keep these dimensions; sum over the rest

	std::vector<IndexT> sub_index;	// Temporary
public:
	int rank() { return sub.rank(); }
	CollapseAccumulator(LesserAccumulatorT &&_sub, std::vector<int> &&_keep_dims)
		: sub(std::move(_sub)), keep_dims(std::move(_keep_dims))
	{
		if (keep_dims.size() != sub.rank()) {
			(*sparse_error)(-1, "Rank mismatch: %d vs %d\n", keep_dims.size(), sub.rank());
		}
		index.reserve(sub.rank());
	}

	void reserve(size_t size) { sub.reserve(size); }
	void add(std::array<IndexT, RANK> const index, ValueT const &val)
	{
		for (int k=0; k<sub.rank(); ++k) sub_index[k] = index[keep_dims[i]];
		sub.add(sub_index, val);
	}

	LesserAccumulatorT &value()
		{ return sub; }
}

#endif	// #if 0


} // namespace


#endif	// Guard
