/**************************************************************************
Copyright (C) 2021 Alexander Shaduri
License: 0BSD (Zero-Clause BSD)
Copyright (C) 2025 Anna DESILLES
License: 0BSD (Zero-Clause BSD)
***************************************************************************/
#ifndef VIABLAB_EXPORT_H
#define VIABLAB_EXPORT_H


/*
Preconditions:
Define VIABLAB_LIBRARY_STATIC when building and using as a static library.
Define VIABLAB_LIBRARY_BUILD when building the library.

Mark public symbols with VIABLAB_LIBRARY_EXPORT.
*/
#ifndef VIABLAB_LIBRARY_STATIC
	/* It's a dynamic library.
	The public symbols must be marked as "exported" when building the library,
	and "imported" when using the library.
	*/
	#ifdef VIABLAB_LIBRARY_BUILD
		/* Building the library */
		#ifdef _WIN32
			#define VIABLAB_LIBRARY_EXPORT __declspec(dllexport)
		#elif __GNUC__ >= 4
			#define VIABLAB_LIBRARY_EXPORT __attribute__((visibility("default")))
		#else
			#define VIABLAB_LIBRARY_EXPORT 
		#endif
	#else
		/* Using the library */
		#ifdef _WIN32
			#define VIABLAB_LIBRARY_EXPORT __declspec(dllimport)
		#else
			#define VIABLAB_LIBRARY_EXPORT
		#endif
	#endif
#endif

#ifndef VIABLAB_LIBRARY_EXPORT
	/* It's a static library, no need to import/export anything */
	#define VIABLAB_LIBRARY_EXPORT
#endif

#endif