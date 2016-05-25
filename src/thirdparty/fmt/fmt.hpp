#ifndef FMT_TPL_H_
#define FMT_TPL_H_

// customize the format library for the ASC toolkit project

// Don't use exceptions
#undef MUST_UNDEF_FMT_EXCEPTIONS
#ifndef FMT_EXCEPTIONS
  #define FMT_EXCEPTIONS 0
  #define MUST_UNDEF_FMT_EXCEPTIONS
#endif

// Use the library as header-only
#undef MUST_UNDEF_FMT_HEADER_ONLY
#ifndef FMT_HEADER_ONLY
  #define FMT_HEADER_ONLY 1
  #define MUST_UNDEF_FMT_HEADER_ONLY
#endif


#include "detail/format.h"


// Undef the defines above
#ifdef MUST_UNDEF_FMT_EXCEPTIONS
  #undef MUST_UNDEF_FMT_EXCEPTIONS
  #undef FMT_EXCEPTIONS
#endif

#ifdef FMT_HEADER_ONLY
  #undef FMT_HEADER_ONLY
  #undef MUST_UNDEF_FMT_HEADER_ONLY
#endif

#endif  // FMT_TPL_H_
