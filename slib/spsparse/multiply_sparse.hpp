#ifndef SPSPARSE_MULTIPLY_SPARSE_HPP
#define SPSPARSE_MULTIPLY_SPARSE_HPP

#include <spsparse/xiter.hpp>
#include <spsparse/array.hpp>

namespace spsparse {

// -------------------------------------------------------------
/** Superclass to chose whether we do/do not multiply our matrix by a scaling diag. matrix */
template<class MatT>
struct MultXiter {
	virtual typename MatT::index_type index() = 0;
	virtual bool eof() = 0;
	virtual void operator++() = 0;
	virtual typename MatT::val_type scale_val() = 0;
	virtual typename DimBeginningsXiter<MatT>::sub_xiter_type sub_xiter() = 0;
};


template<class MatT>
class SimpleMultXiter : public MultXiter<MatT>
{
	// SPSPARSE_LOCAL_TYPES(MatT);
	DimBeginningsXiter<MatT> ii;
public:
	SimpleMultXiter(MatT const &A)
		: ii(A.dim_beginnings_xiter())
	{}

	typename MatT::index_type index() { return *ii; }
	bool eof() { return ii.eof(); }
	void operator++() { ++ii; }
	typename MatT::val_type scale_val() { return 1; }
	typename DimBeginningsXiter<MatT>::sub_xiter_type sub_xiter() { return ii.sub_xiter(); }
};



/** Used to iterate along row/columns for multiply, along with a
diagonal scale matrix. */
template<class MatT, class ScaleT>
class ScaledMultXiter : public MultXiter<MatT>
{
	// SPSPARSE_LOCAL_TYPES(MatT);

	Join2Xiter<
		DimBeginningsXiter<MatT>,
		ValSTLXiter<typename ScaleT::const_dim_iterator>> ii;
public:
	ScaledMultXiter(MatT const &A, ScaleT const &scale)
		: ii(A.dim_beginnings_xiter(),
			make_val_xiter(scale.dim_begin(0), scale.dim_end(0)))
	{}

	typename MatT::index_type index() { return *ii.i1; }
	bool eof() { return ii.eof(); }
	void operator++() { ++ii; }
	typename MatT::val_type scale_val() { return ii.i2.val(); }
	typename DimBeginningsXiter<MatT>::sub_xiter_type sub_xiter() { return ii.i1.sub_xiter(); }
};

template<class MatT, class ScaleT>
std::unique_ptr<MultXiter<MatT>> new_mult_xiter(MatT const &A, ScaleT const *scale)
{
	if (scale) {
		return std::unique_ptr<MultXiter<MatT>>(
			new ScaledMultXiter<MatT, ScaleT>(A, *scale));
	} else {
		return std::unique_ptr<MultXiter<MatT>>(
			new SimpleMultXiter<MatT>(A));
	}

};
// -------------------------------------------------------------
/** NOTES:
A must be consolidated ROW_MAJOR, B must be consolidated COLUMN_MAJOR
*/
template<class ScaleIT, class MatAT, class ScaleJT, class MatBT, class ScaleKT, class AccumulatorT>
static void multiply(
	AccumulatorT &ret,
	double C,		// Multiply everything by this
	ScaleIT const *scalei,
	MatAT const &A,
	char transpose_A,			// 'T' for transpose, '.' otherwise
	ScaleJT const *scalej,		// Vector
	MatBT const &B,
	char transpose_B,			// 'T' for transpose, '.' otherwise
	ScaleKT const *scalek,
	DuplicatePolicy duplicate_policy = DuplicatePolicy::ADD,
	bool zero_nan = false)
{
	// Set dimensions of output, even if we store nothing in it.
	std::array<int,2> const &a_sort_order(transpose_A == 'T' ? COL_MAJOR : ROW_MAJOR);
	std::array<int,2> const &b_sort_order(transpose_B == 'T' ? ROW_MAJOR : COL_MAJOR);
	ret.set_shape({A.shape[a_sort_order[0]], B.shape[b_sort_order[0]]});

	// Short-circuit return on empty output
	// (DimBeginningXiter doesn't like size()==0)
	if ((isnone(C)) 
		|| (scalei && scalei->size() == 0)
		|| (A.size() == 0)
		|| (scalej && scalej->size() == 0)
		|| (B.size() == 0)
		|| (scalek && scalek->size() == 0))
	{ return; }

	// --------- Consolidate the matrices if needed
	Consolidate<MatAT> Acon(&A, a_sort_order, duplicate_policy, zero_nan);
	Consolidate<MatBT> Bcon(&B, b_sort_order, duplicate_policy, zero_nan);

//	ret.set_shape({Acon().shape[Acon().sort_order[0]], Bcon().shape[Bcon().sort_order[0]]});

	// Multiply each row by each column
	// ------ Loop 1: Rows in A	
	for (auto join_a(new_mult_xiter(Acon(), scalei));
		!join_a->eof(); ++(*join_a))
	{
		if (isnone(join_a->scale_val())) continue;
		auto aix(join_a->index());
		auto a_scale(join_a->scale_val());

#if 0
printf("Starting row %d:", aix);
for (auto ii=join_a->sub_xiter(); !ii.eof(); ++ii) {
	printf(" (%d : %g)", *ii, ii.val());
	}
printf("\n");
#endif

		// ------ Loop 2: Cols in B
		for (auto join_b(new_mult_xiter(Bcon(), scalek));
			!join_b->eof(); ++(*join_b))
		{
			if (isnone(join_b->scale_val())) continue;
			auto bix(join_b->index());
			auto b_scale(join_b->scale_val());
#if 0
printf("  Starting col %d\n", bix);
#endif

			// ---------- Loop 3: Multiply A row by B column
			typename AccumulatorT::val_type sum = 0;

			if (scalej) {
				// ----- With scalej
				for (auto ab(join3_xiter(
					join_a->sub_xiter(),
					make_val_xiter(scalej->dim_begin(0), scalej->dim_end(0)),
					join_b->sub_xiter()));
					!ab.eof(); ++ab)
				{ sum += ab.i1.val() * ab.i2.val() * ab.i3.val(); }
			} else {
				// ----- Without scalej
				for (auto ab(join2_xiter(
					join_a->sub_xiter(),
					join_b->sub_xiter()));
					!ab.eof(); ++ab)
				{ sum += ab.i1.val() * ab.i2.val(); }
			}

			if (!isnone(sum)) {
#if 0
printf("set: (%d %d : %g)\n", aix, bix, sum * C * a_scale * b_scale);
#endif
				ret.add({aix, bix}, sum * C * a_scale * b_scale);
			}

		}
	}

}

#if 0
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
printf("answer add %d %d: %g\n", *join_a.i2, *join_b.i1, val);
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
