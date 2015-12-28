#include <functional>
#include <cstdlib>
#include "sparsesparse.hpp"
#include <exception>

namespace sparsesparse {

void default_error(int retcode, const char *format, ...)
{
	va_list arglist;

	va_start(arglist, format);
	fprintf(stderr, format, arglist);
	va_end(arglist);
	fprintf(stderr, "\n");

	throw sparsesparse::Exception();
//	exit(-1);
}

error_ptr sparse_error = &default_error;

} 	// Namespace
