#ifndef XITER_HPP
#define XITER_HPP

namespace spsparse {

/** Convert standard STL iterator into our XIter. */
template<class STLIter>
class STLXiter {
	STLIter const begin;
	STLIter ii;
	STLIter const end;
public:
	typedef typename STLIter::value_type value_type;

	STLXiter(STLIter const &_begin, STLIter const &_end) :
		begin(_begin), ii(_begin), end(_end) {}

	bool eof()
		{ return (ii == end); }

	size_t offset()
		{ return ii - begin; }

	void operator++()
		{ ++ii; }

	auto operator*() -> decltype(*ii)
		{ return *ii; }
};
// --------------------------------------------------------

#if 0
// http://stackoverflow.com/questions/984394/why-not-infer-template-parameter-from-constructor
template<class STLIter>
STLXiter<STLIter> make_stl_xiter(typename STLIter::value_type const &_begin, typename STLIter::value_type const &_end)
	{ return STLXiter<STLIter>(_begin, _end); }

#endif

// -------------------------------------------------------





// -------------------------------------------------------
/** Join three (sorted, non-repeating) iterators, producing output when the two match. */
template<class Xiter1T, class Xiter2T, class Xiter3T>
class Join3Xiter
{
	typename Xiter1T::value_type next_match;
	bool _eof;

public:
	// Allow user to access underlying Xiters, to get useful stuff out of them.
	Xiter1T i1;
	Xiter2T i2;
	Xiter3T i3;
	int total_in_use;

	bool eof()
	{
		return _eof;
	}

	Join3Xiter(Xiter1T &&_i1, Xiter2T &&_i2, Xiter3T &&_i3) :
		i1(std::move(_i1)),
		i2(std::move(_i2)),
		i3(std::move(_i3))
	{
		_eof = i1.eof();
		if (_eof) return;
		next_match = *i1;
		next_noincr();
	}

private:
	void next_noincr()
	{
#define JOIN_RANK 3
#include "next_noincr_body.hpp"
#undef JOIN_RANK
	}

public:
	void operator++()
	{
		++i1;
		++i2;
		++i3;
		// We now don't know, i1 or i2 might now be EOF
		next_noincr();
	}

};

// ========================================================
template<class Xiter1T, class Xiter2T>
class Join2Xiter
{
	typename Xiter1T::value_type next_match;
	bool _eof;

public:
	// Allow user to access underlying Xiters, to get useful stuff out of them.
	Xiter1T i1;
	Xiter2T i2;
	int total_in_use;

	bool eof()
		{ return _eof; }

	Join2Xiter(Xiter1T &&_i1, Xiter2T &&_i2) :
		i1(std::move(_i1)),
		i2(std::move(_i2))
	{
		_eof = i1.eof();
		if (_eof) return;
		next_match = *i1;
		next_noincr();
	}

private:
	void next_noincr()
	{
#define JOIN_RANK 2
#include "next_noincr_body.hpp"
#undef JOIN_RANK
	}

public:
	void operator++()
	{
		++i1;
		++i2;
		// We now don't know, i1 or i2 might now be EOF
		next_noincr();
	}

};

// ========================================================

}

#endif
