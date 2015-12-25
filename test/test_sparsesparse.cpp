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


// Tests that STLSparseSparse replicates standard
// STL iterator with new interface
TEST_F(SparseSparseTest, CooArray) {
	CooArray<int, double, 1> arr1;
	arr1.add({{1}}, 2.);
	arr1.add({{3}}, 6.);
	EXPECT_EQ(2, arr1.size());
}


int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
