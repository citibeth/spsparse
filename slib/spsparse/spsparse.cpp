#include <functional>
#include <cstdlib>
#include <spsparse/spsparse.hpp>
#include <exception>

namespace spsparse {

void default_error(int retcode, const char *format, ...)
{
	va_list arglist;

	va_start(arglist, format);
	fprintf(stderr, format, arglist);
	va_end(arglist);
	fprintf(stderr, "\n");

	throw spsparse::Exception();
//	exit(-1);
}

error_ptr sparse_error = &default_error;

const std::array<int,2> ROW_MAJOR = {0,1};
const std::array<int,2> COLUMN_MAJOR = {1,0};

} 	// Namespace
