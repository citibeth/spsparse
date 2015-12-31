#ifndef SPSPARSE_SPSPARSE_H
#define SPSPARSE_SPSPARSE_H

#include <exception>
#include <iostream>

#include <spsparse/blitz.hpp>

namespace spsparse {

enum class DuplicatePolicy {LEAVE_ALONE, ADD, REPLACE};


class Exception : public std::exception
{
public:
	virtual ~Exception()
		{}

	virtual const char* what() const noexcept
		{ return "spsparse::Exception()"; }
};


// This is not possible in C++11
//extern std::function<void(int, const char *, ...)> error;

// Use this instead.
// http://www.thecodingforums.com/threads/function-pointers-to-printf.317925/
typedef void (*error_ptr) (int retcode, char const *str, ...);
extern error_ptr sparse_error;

#define SPSPARSE_LOCAL_TYPES(ArrayOrIterT) \
	static const int rank = ArrayOrIterT::rank; \
	typedef typename ArrayOrIterT::index_type index_type; \
	typedef typename ArrayOrIterT::val_type val_type; \
	typedef std::array<index_type, rank> indices_type; \
	typedef blitz::Array<val_type, rank> blitz_type

// -------------------------------------------------------------
// Values for sort_order formal parameter below
extern const std::array<int,2> ROW_MAJOR;
extern const std::array<int,2> COL_MAJOR;

// The expression "std::isnan(n) || (n == 0)" for different data types.
// Use template specilization here...
template<class NumT>
inline bool isnone(NumT const n, bool const zero_nan=false)
{
	if (zero_nan) {
		return std::isnan(n) || (n == 0);
	} else {
		return (n == 0);
	}
}

// ----------------------------------------------------------
} // Namespace

// -------------------------------------------------------------
/** Hack to write std::array to ostream. */
template<class T>
std::ostream &stream(std::ostream &os, T const * const a, int RANK);

template<class T>
std::ostream &stream(std::ostream &os, T const * const a, int RANK)
{
	if (RANK == 0) {
		os << "{}";
	} else {
		os << "{";
		for (int k=0; ; ) {
			os << a[k];
			++k;
			if (k == RANK) break;
			os << ", ";
		}
		os << "}";
	}
	return os;
}


#endif // Guard
