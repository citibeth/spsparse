// https://github.com/google/googletest/blob/master/googletest/docs/Primer.md

#include <gtest/gtest.h>
#include <sparsesparse.hpp>
#include <iostream>

using namespace sparsesparse;

// The fixture for testing class Foo.
class SparseSparseTest : public ::testing::Test {
protected:

	// You can do set-up work for each test here.
	SparseSparseTest() {}

	// You can do clean-up work that doesn't throw exceptions here.
	virtual ~SparseSparseTest() {}

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


TEST_F(SparseSparseTest, CooArray) {
	// Make a simple CooArray
	CooArray<int, double, 1> arr1({{4}}, {{"dim0"}});
	arr1.add({1}, 2.);
	arr1.add({3}, 6.);
	EXPECT_EQ(2, arr1.size());
	EXPECT_EQ(1, arr1.indices[0][0]);
	EXPECT_EQ(3, arr1.indices[0][1]);
	EXPECT_EQ(2., arr1.vals[0]);

	// Test bounds checking
	try {
		arr1.add({17}, 4.);
		FAIL() << "Excpected sparsesparse::Exception";
	} catch(sparsesparse::Exception const &err) {
	} catch(...) {
		FAIL() << "Excpected sparsesparse::Exception";
	}
}


TEST_F(SparseSparseTest, consolidate) {
	// 2-D CooArray; test consolidate
	CooArray<int, double, 2> arr2({2,4});
	arr2.add({1,3}, 5.);
	arr2.add({1,2}, 3.);
	arr2.add({0,3}, 17.);
	arr2.add({1,2}, 15.);

	arr2.consolidate({1,0});
	EXPECT_EQ(3, arr2.size());
	EXPECT_EQ(std::vector<int>({0,1,1}), arr2.indices[0]);	// i
	EXPECT_EQ(std::vector<int>({3,2,3}), arr2.indices[1]);	// j
	EXPECT_EQ(std::vector<double>({17., 18., 5.}), arr2.vals);

	// Try out iterators a bit
	auto ii(arr2.begin());
	EXPECT_NE(ii, arr2.end());
	EXPECT_EQ(0, ii.index(0));
	EXPECT_EQ(3, ii.index(1));
	EXPECT_EQ(17., ii.val());
	EXPECT_EQ(17., *ii);
	++ii;
	EXPECT_NE(ii, arr2.end());
	EXPECT_EQ(1, ii.index(0));
	EXPECT_EQ(2, ii.index(1));
	EXPECT_EQ(18., *ii);


}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
