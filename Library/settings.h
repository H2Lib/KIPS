
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
 * Types
 * ------------------------------------------------------------ */

/** @brief Unsigned integer type.
 *
 *  This type is mostly used to access components of arrays,
 *  vectors and matrices. */
typedef uint32_t uint;

/** @brief Signed integer constant zero. */
extern const int i_zero;

/** @brief Signed integer constant one. */
extern const int i_one;

/** @brief Unsigned integer constant zero. */
extern const uint u_zero;

/** @brief Unsigned integer constant one. */
extern const uint u_one;

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
#define H2_CHECK_TOLERANCE 1.0e-6
#else
#define H2_CHECK_TOLERANCE 1.0e-12
#endif

/** Bound for determining when a number is essentially zero. */
#ifdef USE_FLOAT
#define H2_ALMOST_ZERO 1e-30
#else
#define H2_ALMOST_ZERO 1e-300
#endif

/** @brief Pointer to @ref real array. */
typedef real *preal;

/** @brief Pointer to constant @ref real array. */
typedef const real *pcreal;

/** @brief @ref real constant zero */
extern const real r_zero;

/** @brief @ref real constant one. */
extern const real r_one;

/** @brief @ref real constant minus one. */
extern const real r_minusone;

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

/** @brief @ref field constant zero */
extern const field f_zero;

/** @brief @ref field constant one. */
extern const field f_one;

/** @brief @ref field constant minus one. */
extern const field f_minusone;

#ifdef USE_COMPLEX
/** @brief @ref field constant for the imaginary number. */
extern const field f_i;
#endif

#ifndef USE_COMPLEX
#  ifdef I
#  undef I
#  endif
#define I 0.0
#endif

/** @} */

#endif
