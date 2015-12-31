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
/** Matrix-matrix multiply */
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

	// Check inner dimensions
	if (A.shape[a_sort_order[1]] != B.shape[b_sort_order[1]]) {
		(*sparse_error)(-1, "Inner dimensions for A (%ld) and B (%ld) must match!", A.shape[a_sort_order[1]], B.shape[b_sort_order[1]]);
	}

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
// ------------------------------------------------------------
/** Matrix-vector multiply */
template<class ScaleIT, class MatAT, class ScaleJT, class VecT, class AccumulatorT>
static void multiply(
	AccumulatorT &ret,
	double C,		// Multiply everything by this
	ScaleIT const *scalei,
	MatAT const &A,
	char transpose_A,			// 'T' for transpose, '.' otherwise
	ScaleJT const *scalej,		// Vector
	VecT const &V,
	DuplicatePolicy duplicate_policy = DuplicatePolicy::ADD,
	bool zero_nan = false)
{
	// Set dimensions of output, even if we store nothing in it.
	std::array<int,2> const &a_sort_order(transpose_A == 'T' ? COL_MAJOR : ROW_MAJOR);
	ret.set_shape({A.shape[a_sort_order[0]]});

	// Check inner dimensions
	if (A.shape[a_sort_order[1]] != V.shape[0]) {
		(*sparse_error)(-1, "Inner dimensions for A (%ld) and V (%ld) must match!", A.shape[a_sort_order[1]], V.shape[0]);
	}

	// Short-circuit return on empty output
	// (DimBeginningXiter doesn't like size()==0)
	if ((isnone(C)) 
		|| (scalei && scalei->size() == 0)
		|| (A.size() == 0)
		|| (scalej && scalej->size() == 0)
		|| (V.size() == 0))
	{ return; }

	// --------- Consolidate the matrices if needed
	Consolidate<MatAT> Acon(&A, a_sort_order, duplicate_policy, zero_nan);
	Consolidate<VecT> Vcon(&V, {0}, duplicate_policy, zero_nan);

//std::cout << "Vcon: " << Vcon() << std::endl;

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
		// ---------- Loop 3: Multiply A row by the single B column
		typename AccumulatorT::val_type sum = 0;

		if (scalej) {
			// ----- With scalej
			for (auto ab(join3_xiter(
				join_a->sub_xiter(),
				make_val_xiter(scalej->dim_begin(0), scalej->dim_end(0)),
				make_val_xiter(Vcon().dim_begin(0), Vcon().dim_end(0))));
				!ab.eof(); ++ab)
			{ sum += ab.i1.val() * ab.i2.val() * ab.i3.val(); }
		} else {
			// ----- Without scalej
			for (auto ab(join2_xiter(
				join_a->sub_xiter(),
				make_val_xiter(Vcon().dim_begin(0), Vcon().dim_end(0))));
				!ab.eof(); ++ab)
			{
//printf("    triplet %d=?%d: %g*%g = %g\n", *ab.i1, *ab.i2, ab.i1.val(), ab.i2.val(), ab.i1.val()*ab.i2.val());
				sum += ab.i1.val() * ab.i2.val();
			}
		}

		if (!isnone(sum)) {
#if 0
printf("    set: (%d : %g)\n", aix, sum * C * a_scale);
#endif
			ret.add({aix}, sum * C * a_scale);
		}

	}

}





} // namespace


#endif	// Guard
