#ifndef XITER_HPP
#define XITER_HPP

namespace spsparse {


// -----------------------------------------------
template<class IterT>
class WrapForwardValIter {
public:
	typedef typename IterT::value_type value_type;

	// http://www.cplusplus.com/reference/iterator/ForwardIterator/
	IterT ii;

	WrapForwardValIter(IterT const &_ii) : ii(_ii) {}
	WrapForwardValIter(IterT const &&_ii) : ii(std::move(_ii)) {}

	bool operator==(WrapForwardValIter<IterT> const &other)
		{ return ii == other.ii; }
	bool operator!=(WrapForwardValIter<IterT> const &other)
		{ return ii != other.ii; }

	auto operator*() -> decltype(*ii)
		{ return *ii; }
//	auto operator->*() -> decltype(*ii)
//		{ return ii.operator->*(); }

	auto val() -> decltype(ii.val())
		{ return ii.val(); }

	void operator++()
		{ ++ii; }

//	auto operator++(int) -> delctype(ii++)	// postfix
//		{ return ii++; }

};

// -----------------------------------------------
/** Convert standard STL iterator into our XIter. */
template<class STLIter>
class STLXiter {
public:
	STLIter const begin;
	STLIter ii;
	STLIter const end;

	typedef typename STLIter::value_type value_type;

	STLXiter(STLIter const &_begin, STLIter const &_end) :
		begin(_begin), ii(_begin), end(_end) {}

	bool eof() { return (ii == end); }
	size_t offset() { return ii - begin; }
	void operator++() { ++ii; }
	auto operator*() -> decltype(*ii) { return *ii; }
};

/** Works for iterators with .val() method */
template<class ValSTLIter>
class ValSTLXiter : public STLXiter<ValSTLIter>
{
public:

	ValSTLXiter(ValSTLIter const &_begin, ValSTLIter const &_end) :
		STLXiter<ValSTLIter>(_begin, _end) {}

	auto val() -> decltype(STLXiter<ValSTLIter>::ii.val())
		{ return this->ii.val(); }
};
// --------------------------------------------------------

// --------------------------------------------------------

// http://stackoverflow.com/questions/984394/why-not-infer-template-parameter-from-constructor

template<class STLIter>
STLXiter<STLIter> make_xiter(
	STLIter **_begin, STLIter &&_end)
{ return STLXiter<STLIter>(std::move(_begin), std::move(_end)); }

template<class ValSTLIter>
ValSTLXiter<ValSTLIter> make_val_xiter(
	ValSTLIter &&_begin, ValSTLIter &&_end)
{ return ValSTLXiter<ValSTLIter>(std::move(_begin), std::move(_end)); }


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

	bool eof() { return _eof; }

	Join3Xiter(Xiter1T &&_i1, Xiter2T &&_i2, Xiter3T &&_i3) :
		next_match(0),
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

template<class Xiter1T, class Xiter2T, class Xiter3T>
class Join3Xiter<Xiter1T, Xiter2T, Xiter3T> join3_xiter(Xiter1T &&i1, Xiter2T &&i2, Xiter3T &&i3)
	{ return Join3Xiter<Xiter1T, Xiter2T, Xiter3T>(std::move(i1), std::move(i2), std::move(i3)); }

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

template<class Xiter1T, class Xiter2T>
class Join2Xiter<Xiter1T, Xiter2T> join2_xiter(Xiter1T &&i1, Xiter2T &&i2)
	{ return Join2Xiter<Xiter1T, Xiter2T>(std::move(i1), std::move(i2)); }

// ========================================================

}

#endif
