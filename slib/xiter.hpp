#ifndef XITER_HPP
#define XITER_HPP

namespace sparsesparse {

/** Convert standard STL iterator into our XIter. */
template<class STLIter>
class STLXiter {
	STLIter ii;
	STLIter end;
public:
	typedef typename STLIter::value_type value_type;

	STLXiter(STLIter const &_begin, STLIter const &_end) {
		ii = _begin;
		end = _end;
	}

	bool eof()
		{ return (ii == end); }

	void operator++() {
		++ii;
	}

	typename STLIter::value_type operator*()
		{ return *ii; }
};

#if 0
// http://stackoverflow.com/questions/984394/why-not-infer-template-parameter-from-constructor
template<class STLIter>
STLXiter<STLIter> make_stl_xiter(typename STLIter::value_type const &_begin, typename STLIter::value_type const &_end)
	{ return STLXiter<STLIter>(_begin, _end); }

#endif

// -------------------------------------------------------
/** Join two (sorted, non-repeating) iterators, producing output when the two match. */
template<class Xiter1T, class Xiter2T>
class Join2Xiter
{
	typename Xiter1T::value_type next_match;
	bool _eof;

public:
	// Allow user to access underlying Xiters, to get useful stuff out of them.
	Xiter1T i1;
	Xiter2T i2;

	bool eof()
	{
		return _eof;
	}

	Join2Xiter(Xiter1T &&_i1, Xiter2T &&_i2) :
		i1(std::move(_i1)),
		i2(std::move(_i2))
	{
		if (!i1.eof() && !i2.eof()) {
			_eof = false;
			next_match = std::max(*i1, *i2);
			next_noincr();
		} else {
			_eof = true;
		}
	}

private:
	void next_noincr()
	{
		do {
			// Scan forward first iterator
			for (;;++i1) {
				if (i1.eof()) {
					_eof = true;
					return;
				}
				if (*i1 >= next_match) break;
			}

			next_match = *i1;

			// Scan forward second iterator
			for (;;++i2) {
				if (i2.eof()) {
					_eof = true;
					return;
				}
				if (*i2 >= next_match) break;
			}
			next_match = *i2;
		} while (*i1 != *i2);
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

	bool eof()
	{
		return _eof;
	}

	Join3Xiter(Xiter1T &&_i1, Xiter2T &&_i2, Xiter3T &&_i3) :
		i1(std::move(_i1)),
		i2(std::move(_i2)),
		i3(std::move(_i3))
	{
		if (!i1.eof() && !i2.eof() && !i3.eof()) {
			_eof = false;
			next_match = std::max(std::max(*i1, *i2), *i3);
			next_noincr();
		} else {
			_eof = true;
		}
	}

private:
	void next_noincr()
	{
		do {
			// Scan forward first iterator
			for (;;++i1) {
				if (i1.eof()) {
					_eof = true;
					return;
				}
				if (*i1 >= next_match) break;
			}
			next_match = *i1;

			// Scan forward second iterator
			for (;;++i2) {
				if (i2.eof()) {
					_eof = true;
					return;
				}
				if (*i2 >= next_match) break;
			}
			next_match = *i2;

			// Scan forward second iterator
			for (;;++i3) {
				if (i3.eof()) {
					_eof = true;
					return;
				}
				if (*i3 >= next_match) break;
			}
			next_match = *i3;


		} while ((*i1 != *i2) || (*i1 != *i3));
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

}

#endif
