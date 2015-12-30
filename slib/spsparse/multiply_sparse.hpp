#ifndef SPSPARSE_MULTIPLY_SPARSE_HPP
#define SPSPARSE_MULTIPLY_SPARSE_HPP

#include <spsparse/xiter.hpp>
#include <spsparse/array.hpp>

namespace spsparse {

// ----------------------------------------------------------
/** Computes [row] diag(scale) [[col]].

Or: row_j scale_j col_j (out of overall matrix computation row_ij
scale_j col_jk

NOTE: row, col and scale must NOT have duplicate elements.
Consolidate() must have been run before --- which gets ride of nans
too, via handle_nan.
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


// -------------------------------------------------------------
#if 1
/** NOTES:
A must be consolidated ROW_MAJOR, B must be consolidated COLUMN_MAJOR
*/
template<class ScaleIT, class MatAT, class ScaleJT, class MatBT, class ScaleKT, class AccumulatorT>
static void multiply(
	AccumulatorT &ret,
	double C,		// Multiply everything by this
	ScaleIT const &scalei,		// Vector
	MatAT   &A,			// Matrix
	ScaleJT const &scalej,		// Vector
	MatBT   &B,			// Matrix
	ScaleKT const &scalek)		// Vector
{

//Check dimensions...?

	if (C == 0) return;		// Everything will be zero, so give up now!

	// This will do nothing if they're already consolidated
	A.consolidate({0,1});
	B.consolidate({1,0});

	// Multiply each row by each column
	// ------ Outer loop: Rows in A
	for (auto join_a(join2_xiter(
		make_val_xiter(scalei.dim_begin(0), scalei.dim_end(0)),
		A.dim_beginnings_xiter()));
		!join_a.eof(); ++join_a)
	{

		// -------- Inner loop: Cols in B
		for (auto join_b(join2_xiter(
			B.dim_beginnings_xiter(),
			make_val_xiter(scalek.dim_begin(0), scalek.dim_end(0))));
			!join_b.eof(); ++join_b)
		{

			// ---------- Inside of loop
			// Multiply A row by B column
			ScalarAccumulator<typename MatAT::index_type, typename MatAT::val_type, 2> rcval;
			multiply_row_scale_col(rcval,
				join_a.i2.sub_xiter(),
				make_val_xiter(scalej.dim_begin(0), scalej.dim_end(0)),
				join_b.i1.sub_xiter());

			double val = C * join_a.i1.val() * rcval.val * join_b.i2.val();
			if (val == 0) continue;

			// *a_ii == index of current row in A
			ret.add({*join_a.i2, *join_b.i1}, val);
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


#endif


} // namespace


#endif	// Guard
