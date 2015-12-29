#include <functional>
#include <cstdlib>
#include "spsparse.hpp"
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

} 	// Namespace
