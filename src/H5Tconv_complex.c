/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the COPYING file, which can be found at the root of the source code       *
 * distribution tree, or in https://www.hdfgroup.org/licenses.               *
 * If you do not have access to either file, you may request a copy from     *
 * help@hdfgroup.org.                                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * Purpose: Datatype conversion functions for complex number datatypes
 */

/****************/
/* Module Setup */
/****************/
#include "H5Tmodule.h" /* This source code file is part of the H5T module */

/***********/
/* Headers */
/***********/
#include "H5Eprivate.h"
#include "H5Tconv.h"
#include "H5Tconv_macros.h"
#include "H5Tconv_complex.h"

/******************/
/* Local Typedefs */
/******************/

/********************/
/* Local Prototypes */
/********************/

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_complex
 *
 * Purpose:     Convert one complex number type to another. This is the
 *              catch-all function for complex number conversions and is
 *              probably not particularly fast.
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_complex(const H5T_t *src, const H5T_t *dst, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                  size_t nelmts, size_t buf_stride, size_t bkg_stride, void *_buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_complex_struct
 *
 * Purpose:     Convert a complex number type to the equivalent compound
 *              type representation. A compound type must match one of the
 *              the following representations exactly to be considered
 *              equivalent.
 *
 *              H5T_COMPOUND {
 *                  <float_type> "r"; OFFSET 0
 *                  <float_type> "i"; OFFSET SIZEOF("r")
 *              }
 *
 *              H5T_COMPOUND {
 *                  <float_type> "re"; OFFSET 0
 *                  <float_type> "im"; OFFSET SIZEOF("re")
 *              }
 *
 *              H5T_COMPOUND {
 *                  <float_type> "real"; OFFSET 0
 *                  <float_type> "imag"; OFFSET SIZEOF("real")
 *              }
 *
 *              H5T_COMPOUND {
 *                  <float_type> "real";      OFFSET 0
 *                  <float_type> "imaginary"; OFFSET SIZEOF("real")
 *              }
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_complex_struct(const H5T_t *src, const H5T_t *dst, H5T_cdata_t *cdata,
                         const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride, size_t bkg_stride,
                         void *_buf, void *bkg)
{
    return FAIL; /* TODO */
}

#ifdef H5_HAVE_COMPLEX_NUMBERS
/*-------------------------------------------------------------------------
 * Function:    H5T__conv_fcomplex_schar
 *
 * Purpose:     Converts `float _Complex' to `signed char'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_fcomplex_schar(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                         size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */ /* H5T_CONV_su(FLOAT_COMPLEX, SCHAR, float _Complex, signed char, -, -); */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_fcomplex_uchar
 *
 * Purpose:     Converts `float _Complex' to `unsigned char'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_fcomplex_uchar(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                         size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_fcomplex_short
 *
 * Purpose:     Converts `float _Complex' to `short'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_fcomplex_short(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                         size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_fcomplex_ushort
 *
 * Purpose:     Converts `float _Complex' to `unsigned short'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_fcomplex_ushort(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata,
                          const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride, size_t bkg_stride,
                          void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_fcomplex_int
 *
 * Purpose:     Converts `float _Complex' to `int'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_fcomplex_int(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                       size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_fcomplex_uint
 *
 * Purpose:     Converts `float _Complex' to `unsigned int'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_fcomplex_uint(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                        size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_fcomplex_long
 *
 * Purpose:     Converts `float _Complex' to `long'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_fcomplex_long(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                        size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_fcomplex_ulong
 *
 * Purpose:     Converts `float _Complex' to `unsigned long'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_fcomplex_ulong(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                         size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_fcomplex_llong
 *
 * Purpose:     Converts `float _Complex' to `long long'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_fcomplex_llong(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                         size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_fcomplex_ullong
 *
 * Purpose:     Converts `float _Complex' to `unsigned long long'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_fcomplex_ullong(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata,
                          const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride, size_t bkg_stride,
                          void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

#ifdef H5_HAVE__FLOAT16
/*-------------------------------------------------------------------------
 * Function:    H5T__conv_fcomplex__Float16
 *
 * Purpose:     Converts `float _Complex' to `_Float16'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_fcomplex__Float16(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata,
                            const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                            size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}
#endif

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_fcomplex_float
 *
 * Purpose:     Converts `float _Complex' to `float'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_fcomplex_float(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                         size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_fcomplex_double
 *
 * Purpose:     Converts `float _Complex' to `double'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_fcomplex_double(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata,
                          const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride, size_t bkg_stride,
                          void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_fcomplex_ldouble
 *
 * Purpose:     Converts `float _Complex' to `long double'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_fcomplex_ldouble(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata,
                           const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                           size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_dcomplex_schar
 *
 * Purpose:     Converts `double _Complex' to `signed char'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_dcomplex_schar(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                         size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_dcomplex_uchar
 *
 * Purpose:     Converts `double _Complex' to `unsigned char'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_dcomplex_uchar(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                         size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_dcomplex_short
 *
 * Purpose:     Converts `double _Complex' to `short'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_dcomplex_short(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                         size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_dcomplex_ushort
 *
 * Purpose:     Converts `double _Complex' to `unsigned short'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_dcomplex_ushort(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata,
                          const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride, size_t bkg_stride,
                          void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_dcomplex_int
 *
 * Purpose:     Converts `double _Complex' to `int'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_dcomplex_int(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                       size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_dcomplex_uint
 *
 * Purpose:     Converts `double _Complex' to `unsigned int'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_dcomplex_uint(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                        size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_dcomplex_long
 *
 * Purpose:     Converts `double _Complex' to `long'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_dcomplex_long(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                        size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_dcomplex_ulong
 *
 * Purpose:     Converts `double _Complex' to `unsigned long'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_dcomplex_ulong(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                         size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_dcomplex_llong
 *
 * Purpose:     Converts `double _Complex' to `long long'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_dcomplex_llong(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                         size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_dcomplex_ullong
 *
 * Purpose:     Converts `double _Complex' to `unsigned long long'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_dcomplex_ullong(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata,
                          const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride, size_t bkg_stride,
                          void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

#ifdef H5_HAVE__FLOAT16
/*-------------------------------------------------------------------------
 * Function:    H5T__conv_dcomplex__Float16
 *
 * Purpose:     Converts `double _Complex' to `_Float16'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_dcomplex__Float16(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata,
                            const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                            size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}
#endif

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_dcomplex_float
 *
 * Purpose:     Converts `double _Complex' to `float'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_dcomplex_float(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                         size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_dcomplex_double
 *
 * Purpose:     Converts `double _Complex' to `double'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_dcomplex_double(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata,
                          const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride, size_t bkg_stride,
                          void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_dcomplex_ldouble
 *
 * Purpose:     Converts `double _Complex' to `long double'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_dcomplex_ldouble(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata,
                           const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                           size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_lcomplex_schar
 *
 * Purpose:     Converts `long double _Complex' to `signed char'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_lcomplex_schar(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                         size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_lcomplex_uchar
 *
 * Purpose:     Converts `long double _Complex' to `unsigned char'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_lcomplex_uchar(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                         size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_lcomplex_short
 *
 * Purpose:     Converts `long double _Complex' to `short'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_lcomplex_short(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                         size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_lcomplex_ushort
 *
 * Purpose:     Converts `long double _Complex' to `unsigned short'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_lcomplex_ushort(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata,
                          const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride, size_t bkg_stride,
                          void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_lcomplex_int
 *
 * Purpose:     Converts `long double _Complex' to `int'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_lcomplex_int(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                       size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_lcomplex_uint
 *
 * Purpose:     Converts `long double _Complex' to `unsigned int'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_lcomplex_uint(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                        size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_lcomplex_long
 *
 * Purpose:     Converts `long double _Complex' to `long'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_lcomplex_long(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                        size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_lcomplex_ulong
 *
 * Purpose:     Converts `long double _Complex' to `unsigned long'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_lcomplex_ulong(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                         size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_lcomplex_llong
 *
 * Purpose:     Converts `long double _Complex' to `long long'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_lcomplex_llong(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                         size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_lcomplex_ullong
 *
 * Purpose:     Converts `long double _Complex' to `unsigned long long'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_lcomplex_ullong(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata,
                          const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride, size_t bkg_stride,
                          void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

#ifdef H5_HAVE__FLOAT16
/*-------------------------------------------------------------------------
 * Function:    H5T__conv_lcomplex__Float16
 *
 * Purpose:     Converts `long double _Complex' to `_Float16'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_lcomplex__Float16(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata,
                            const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                            size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}
#endif

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_lcomplex_float
 *
 * Purpose:     Converts `long double _Complex' to `float'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_lcomplex_float(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                         size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_lcomplex_double
 *
 * Purpose:     Converts `long double _Complex' to `double'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_lcomplex_double(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata,
                          const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride, size_t bkg_stride,
                          void *buf, void *bkg)
{
    return FAIL; /* TODO */
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_lcomplex_ldouble
 *
 * Purpose:     Converts `long double _Complex' to `long double'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_lcomplex_ldouble(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata,
                           const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                           size_t bkg_stride, void *buf, void *bkg)
{
    return FAIL; /* TODO */
}
#endif /* H5_HAVE_COMPLEX_NUMBERS */
