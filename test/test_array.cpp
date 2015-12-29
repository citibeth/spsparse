// https://github.com/google/googletest/blob/master/googletest/docs/Primer.md

#include <gtest/gtest.h>
#include <spsparse/array.hpp>
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


TEST_F(SpSparseTest, CooArray) {
	// Make a simple CooArray
	CooArray<int, double, 1> arr1({4});
	arr1.add({1}, 2.);
	arr1.add({3}, 6.);
	EXPECT_EQ(2, arr1.size());
	EXPECT_EQ(1, arr1.index(0, 0));
	EXPECT_EQ(3, arr1.index(0,1));
	EXPECT_EQ(2., arr1.val(0));

	// Test bounds checking
	try {
		arr1.add({17}, 4.);
		FAIL() << "Excpected spsparse::Exception";
	} catch(spsparse::Exception const &err) {
	} catch(...) {
		FAIL() << "Excpected spsparse::Exception";
	}
}

TEST_F(SpSparseTest, consolidate) {
	// 2-D CooArray; test consolidate
	CooArray<int, double, 2> arr2({2,4});
	arr2.add({1,3}, 5.);
	arr2.add({1,2}, 3.);
	arr2.add({0,3}, 17.);
	arr2.add({1,2}, 15.);

	CooArray<int, double, 2> arr3(arr2.shape);
printf("-------------------- consolidate\n");
	consolidate(arr3, arr2, {1,0});

	EXPECT_EQ(3, arr3.size());

printf("AA1\n");
	EXPECT_EQ(std::vector<int>({0,1,1}), blitz_to_vector(arr2.indices(0)));
	EXPECT_EQ(std::vector<int>({3,2,3}), blitz_to_vector(arr2.indices(1)));	// j
	EXPECT_EQ(std::vector<double>({17., 18., 5.}), blitz_to_vector(arr2.vals()));

	// Try out iterators a bit
	auto ii(arr2.begin());
	EXPECT_NE(ii, arr2.end());
	EXPECT_EQ(0, ii.index(0));
	EXPECT_EQ(3, ii.index(1));
	EXPECT_EQ(17., ii.val());
	std::array<int,2> arr({0,3});
	EXPECT_EQ(arr, *ii);
	++ii;
	EXPECT_NE(ii, arr2.end());
	EXPECT_EQ(1, ii.index(0));
	EXPECT_EQ(2, ii.index(1));
	EXPECT_EQ(18., ii.val());

}

int main(int argc, char **argv) {
#ifdef USE_EVERYTRACE
	everytrace_init();
#endif
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
