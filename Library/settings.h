
/* ------------------------------------------------------------
 * This is the file "settings.h" of the KIPS package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

/** @file settings.h
 *  @author Steffen B&ouml;rm
 */

#ifndef SETTINGS_H
#define SETTINGS_H

/** @defgroup settings settings
 *  @brief Fundamental types and macros.
 *  @{ */

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#ifdef USE_COMPLEX
#include <complex.h>
#endif

/* ------------------------------------------------------------
 * Compilation settings
 * ------------------------------------------------------------ */

/** @brief Prefix for inline functions. */
#ifdef __cplusplus
#define INLINE_PREFIX static extern "C"
#else
#define INLINE_PREFIX static
#endif

/** @brief Prefix for function declarations. */
#ifdef __cplusplus
#define HEADER_PREFIX extern "C"
#else
#define HEADER_PREFIX
#endif

/** @brief Prefix for external declarations. */
#ifdef __cplusplus
#define IMPORT_PREFIX extern "C"
#else
#define IMPORT_PREFIX
#endif

/* ------------------------------------------------------------
 * Approximation of the spectral norm
 * ------------------------------------------------------------ */

/** @brief Number of steps of the power iteration used to
 *  approximate the spectral norm. */
#define NORM_STEPS 20

/* ------------------------------------------------------------
 * Types
 * ------------------------------------------------------------ */

/** @brief Unsigned integer type.
 *
 *  This type is mostly used to access components of arrays,
 *  vectors and matrices. */
typedef uint32_t uint;

/** @brief Unsigned long type.
 *
 *  This type is used to access components of particularly large
 *  arrays, e.g., matrices in column-major array representation. */
typedef uint_least64_t longindex;

/** @brief @ref real floating point type.
 *
 *  This type is used, e.g., for geometric coordinates, norms
 *  and diagonal elements of self-adjoint matrices. */
#ifdef USE_FLOAT
typedef float real;
#else
typedef double real;
#endif

/** Relative tolerance for run-time checks. */
#ifdef USE_FLOAT
#define KIPS_CHECK_TOLERANCE 1.0e-6
#else
#define KIPS_CHECK_TOLERANCE 1.0e-12
#endif

/** Bound for determining when a number is essentially zero. */
#ifdef USE_FLOAT
#define KIPS_ALMOST_ZERO 1e-30
#else
#define KIPS_ALMOST_ZERO 1e-300
#endif

/** @brief Pointer to @ref real array. */
typedef real *preal;

/** @brief Pointer to constant @ref real array. */
typedef const real *pcreal;

/** @brief Field type.
 *
 *  This type is used in the linear algebra modules to represent
 *  the coefficients of matrices and vectors. */
#ifdef USE_FLOAT
#  ifdef USE_COMPLEX
typedef float _Complex field;
#  else
typedef float field;
#  endif
#else
#  ifdef USE_COMPLEX
typedef double _Complex field;
#  else
typedef double field;
#  endif
#endif

/** @brief Pointer to @ref field array. */
typedef field *pfield;

/** @brief Pointer to constant @ref field array. */
typedef const field *pcfield;

/** @} */

#endif
