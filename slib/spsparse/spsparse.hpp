#include <exception>

namespace spsparse {

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


}
