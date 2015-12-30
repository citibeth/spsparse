// https://github.com/google/googletest/blob/master/googletest/docs/Primer.md

#include <gtest/gtest.h>
#include <spsparse/array.hpp>
#include <spsparse/multiply_sparse.hpp>
#include <iostream>
#ifdef USE_EVERYTRACE
#include <everytrace.h>
#endif

using namespace spsparse;

// The fixture for testing class Foo.
class SpSparseTest : public ::testing::Test {
protected:

	// You can do set-up work for each test here.
	SpSparseTest() {}

	// You can do clean-up work that doesn't throw exceptions here.
	virtual ~SpSparseTest() {}

	// If the constructor and destructor are not enough for setting up
	// and cleaning up each test, you can define the following methods:

	// Code here will be called immediately after the constructor (right
	// before each test).
	virtual void SetUp() {}

	// Code here will be called immediately after each test (right
	// before the destructor).
	virtual void TearDown() {}

//	  // The mock bar library shaed by all tests
//	  MockBar m_bar;
};


TEST_F(SpSparseTest, CooArray)
{
	typedef CooArray<int, double, 2> CooArrayT;

	CooMatrix<int, double> row({2,10});
	row.add({0,8}, 5.);
	row.add({0,4}, 4.);
	row.add({0,0}, 2.);
	row.add({0,3}, 3.);
	row.add({1,7}, 3.);
	row.consolidate({0,1});
	auto rowdbxi(row.dim_beginnings_xiter());	// Iterate over rows

	CooVector<int, double> scale({10});
	scale.add({0}, 2.);
//	scale.add({3}, 3.);
	scale.add({4}, 4.);
	scale.add({8}, 5.);

	CooMatrix<int, double> col({10,1});
	col.add({0,0}, 2.);
	col.add({3,0}, 3.);
//	col.add({4,0}, 4.);
	col.add({8,0}, 5.);
	col.consolidate({1,0});
	auto coldbxi(col.dim_beginnings_xiter());	// Iterate over cols

	CooArray<int, double, 1> ret1({10});
	multiply_row_scale_col
#if 0
		<decltype(rowdbxi.sub_xiter()),
		decltype(make_val_xiter(scale.dim_begin(0), scale.dim_end(0))),
		decltype(coldbxi.sub_xiter())>
#endif
	(ret1,
		rowdbxi.sub_xiter(),
		make_val_xiter(scale.dim_begin(0), scale.dim_end(0)),
		coldbxi.sub_xiter());

	// std::cout << ret1 << std::endl;

	EXPECT_EQ(std::vector<int>({0,8}), blitz_to_vector(ret1.indices(0)));
	EXPECT_EQ(std::vector<double>({8., 125.}), blitz_to_vector(ret1.vals()));

	double sum1 = 0;
	for (auto ii = ret1.begin(); ii != ret1.end(); ++ii) sum1 += ii.val();

	// ------------------
	CooVector<int, double> eye({10});
	for (int i=0; i<10; ++i) eye.add({i}, 1.);

	CooArray<int, double, 2> ret2({1,1});
	multiply(ret2, 1.0,
		eye, row, eye, col, eye);

	EXPECT_EQ(1, ret2.size());
	EXPECT_EQ(0, ret2.index(0,0));	// dim, ix
	EXPECT_EQ(0, ret2.index(1,0));
	EXPECT_EQ(sum1, ret2.val(0));
}


int main(int argc, char **argv) {
#ifdef USE_EVERYTRACE
	everytrace_init();
#endif
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
