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
#include "H5Tconv_integer.h"
#include "H5Tconv_float.h"

/******************/
/* Local Typedefs */
/******************/

/********************/
/* Local Prototypes */
/********************/

static herr_t H5T__conv_complex_loop(const H5T_t *src_p, const H5T_t *dst_p, const H5T_conv_ctx_t *conv_ctx,
                                     size_t nelmts, size_t buf_stride, void *buf);

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
H5T__conv_complex(const H5T_t *src_p, const H5T_t *dst_p, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                  size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                  void H5_ATTR_UNUSED *bkg)
{
    herr_t ret_value = SUCCEED;

    FUNC_ENTER_PACKAGE

    switch (cdata->command) {
        case H5T_CONV_INIT: {
            H5T_atomic_t src_atomic; /* source datatype atomic info      */
            H5T_atomic_t dst_atomic; /* destination datatype atomic info */

            if (!src_p || !dst_p)
                HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a datatype");
            if (!H5T_IS_ATOMIC(src_p->shared->parent->shared))
                HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "invalid source complex number datatype");
            if (!H5T_IS_ATOMIC(dst_p->shared->parent->shared))
                HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "invalid destination complex number datatype");
            if (!src_p->shared->u.cplx.homogeneous || !dst_p->shared->u.cplx.homogeneous)
                HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
                            "source and destination complex number datatypes must be homogeneous");
            src_atomic = src_p->shared->parent->shared->u.atomic;
            dst_atomic = dst_p->shared->parent->shared->u.atomic;
            if (H5T_ORDER_LE != src_atomic.order && H5T_ORDER_BE != src_atomic.order &&
                H5T_ORDER_VAX != src_atomic.order)
                HGOTO_ERROR(H5E_DATATYPE, H5E_UNSUPPORTED, FAIL,
                            "unsupported byte order for source datatype");
            if (H5T_ORDER_LE != dst_atomic.order && H5T_ORDER_BE != dst_atomic.order &&
                H5T_ORDER_VAX != src_atomic.order)
                HGOTO_ERROR(H5E_DATATYPE, H5E_UNSUPPORTED, FAIL,
                            "unsupported byte order for destination datatype");
            if (dst_p->shared->size > 2 * TEMP_FLOAT_CONV_BUFFER_SIZE)
                HGOTO_ERROR(H5E_DATATYPE, H5E_UNSUPPORTED, FAIL, "destination datatype size is too large");
            if (8 * sizeof(int64_t) - 1 < src_atomic.u.f.esize ||
                8 * sizeof(int64_t) - 1 < dst_atomic.u.f.esize)
                HGOTO_ERROR(H5E_DATATYPE, H5E_UNSUPPORTED, FAIL, "exponent field is too large");
            cdata->need_bkg = H5T_BKG_NO;

            break;
        }

        case H5T_CONV_FREE:
            break;

        case H5T_CONV_CONV:
            if (!src_p || !dst_p)
                HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a datatype");
            if (NULL == conv_ctx)
                HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "invalid datatype conversion context pointer");

            if (H5T__conv_complex_loop(src_p, dst_p, conv_ctx, nelmts, buf_stride, buf) < 0)
                HGOTO_ERROR(H5E_DATATYPE, H5E_CANTCONVERT, FAIL, "unable to convert data values");

            break;

        default:
            HGOTO_ERROR(H5E_DATATYPE, H5E_UNSUPPORTED, FAIL, "unknown conversion command");
    }

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5T__conv_complex() */

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_complex_loop
 *
 * Purpose:     Implements the body of the conversion loop when converting
 *              complex number values to another complex number type.
 *
 * NOTE:        The conversion logic in this function is essentially
 *              identical to the logic in the H5T__conv_f_f_loop function.
 *              However, conversion has to be performed on both the real
 *              and imaginary parts of each complex number element. Since
 *              complex numbers have the same representation as an array
 *              of two elements of the base floating-point type, this could
 *              be simulated in some cases with the H5T__conv_f_f_loop
 *              function by doubling the number of elements to be converted
 *              and halving the sizes involved. However, overlapping
 *              elements or a non-zero `buf_stride` value would complicate
 *              the buffer pointer advancements since each part of the
 *              complex number value has to be processed before advancing
 *              the buffer pointer. Application conversion exception
 *              callbacks also pose a problem since they would expect to
 *              receive an entire complex number rather than part of one.
 *              Therefore, the H5T__conv_f_f_loop logic is mostly
 *              duplicated here and fixes to one function should be made to
 *              the other, if appropriate.
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
/* TODO: cleanup and organize */
static herr_t
H5T__conv_complex_loop(const H5T_t *src_p, const H5T_t *dst_p, const H5T_conv_ctx_t *conv_ctx, size_t nelmts,
                       size_t buf_stride, void *buf)
{
    H5T_conv_float_specval_t realval_type; /* floating-point value type (regular, +/-Inf, +/-0, NaN) */
    H5T_conv_float_specval_t imagval_type; /* floating-point value type (regular, +/-Inf, +/-0, NaN) */
    H5T_atomic_t             src_atomic;   /* source datatype atomic info                            */
    H5T_atomic_t             dst_atomic;   /* destination datatype atomic info                       */
    hssize_t                 expo_max;     /* maximum possible dst exponent                          */
    ssize_t  src_delta, dst_delta;         /* source & destination stride                            */
    uint8_t *s, *sp, *d, *dp;              /* source and dest traversal ptrs                         */
    uint8_t *src_rev = NULL;               /* order-reversed source buffer                           */
    uint8_t  dbuf[2 * TEMP_FLOAT_CONV_BUFFER_SIZE]; /* temp destination buffer */
    size_t   src_part_size; /* size of each complex number part                       */
    size_t   dst_part_size; /* size of each complex number part                       */
    size_t   olap;          /* num overlapping elements                               */
    int      direction;     /* forward or backward traversal                          */
    herr_t   ret_value = SUCCEED;

    assert(src_p);
    assert(src_p->shared->type == H5T_COMPLEX);
    assert(src_p->shared->u.cplx.homogeneous);
    assert(dst_p);
    assert(dst_p->shared->type == H5T_COMPLEX);
    assert(dst_p->shared->u.cplx.homogeneous);
    assert(conv_ctx);
    assert(buf);

    FUNC_ENTER_PACKAGE

    src_atomic    = src_p->shared->parent->shared->u.atomic;
    dst_atomic    = dst_p->shared->parent->shared->u.atomic;
    src_part_size = src_p->shared->size / 2;
    dst_part_size = dst_p->shared->size / 2;
    expo_max      = ((hssize_t)1 << dst_atomic.u.f.esize) - 1;

    /*
     * Do we process the values from beginning to end or vice versa? Also,
     * how many of the elements have the source and destination areas
     * overlapping?
     */
    if (src_p->shared->size == dst_p->shared->size || buf_stride) {
        sp = dp   = (uint8_t *)buf;
        direction = 1;
        olap      = nelmts;
    }
    else if (src_p->shared->size >= dst_p->shared->size) {
        double olap_d =
            ceil((double)(dst_p->shared->size) / (double)(src_p->shared->size - dst_p->shared->size));
        olap = (size_t)olap_d;
        sp = dp   = (uint8_t *)buf;
        direction = 1;
    }
    else {
        double olap_d =
            ceil((double)(src_p->shared->size) / (double)(dst_p->shared->size - src_p->shared->size));
        olap      = (size_t)olap_d;
        sp        = (uint8_t *)buf + (nelmts - 1) * src_p->shared->size;
        dp        = (uint8_t *)buf + (nelmts - 1) * dst_p->shared->size;
        direction = -1;
    }

    /* Direction & size of buffer traversal */
    H5_CHECK_OVERFLOW(buf_stride, size_t, ssize_t);
    H5_CHECK_OVERFLOW(src_p->shared->size, size_t, ssize_t);
    H5_CHECK_OVERFLOW(dst_p->shared->size, size_t, ssize_t);
    src_delta = (ssize_t)direction * (ssize_t)(buf_stride ? buf_stride : src_p->shared->size);
    dst_delta = (ssize_t)direction * (ssize_t)(buf_stride ? buf_stride : dst_p->shared->size);

    /* Allocate space for order-reversed source buffer */
    if (conv_ctx->u.conv.cb_struct.func)
        if (NULL == (src_rev = H5MM_calloc(src_p->shared->size)))
            HGOTO_ERROR(H5E_DATATYPE, H5E_CANTALLOC, FAIL, "couldn't allocate temporary buffer");

    /* The conversion loop */
    for (size_t elmtno = 0; elmtno < nelmts; elmtno++) {
        H5T_conv_ret_t except_ret = H5T_CONV_UNHANDLED; /* return of conversion exception callback function */
        ssize_t        real_mant_msb = 0;               /* most significant bit set in mantissa             */
        ssize_t        imag_mant_msb = 0;               /* most significant bit set in mantissa             */
        int64_t        real_expo, imag_expo;            /* exponent                                         */
        size_t         real_implied;                    /* destination implied bits                         */
        size_t         imag_implied;                    /* destination implied bits                         */
        size_t         mpos;                            /* offset to useful mant in src                     */
        size_t         real_msize = 0, imag_msize = 0;  /* useful size of mantissa in src                   */
        size_t         real_mrsh, imag_mrsh;            /* amount to right shift mantissa                   */
        bool           reverse           = true;        /* if reversed the order of destination             */
        bool           real_denormalized = false;       /* is either source or destination denormalized?    */
        bool           imag_denormalized = false;       /* is either source or destination denormalized?    */
        bool           real_carry        = false;       /* carry after rounding mantissa                    */
        bool           imag_carry        = false;       /* carry after rounding mantissa                    */
        bool           real_zero         = false;       /* if real part is +/-0                             */
        bool           imag_zero         = false;       /* if imaginary part is +/-0                        */
        bool           real_except       = false;       /* if an exception happened for the real part       */
        bool           imag_except       = false;       /* if an exception happened for the imaginary part  */

        /*
         * If the source and destination buffers overlap then use a
         * temporary buffer for the destination.
         */
        s = sp;
        if (direction > 0)
            d = elmtno < olap ? dbuf : dp;
        else
            d = elmtno + olap >= nelmts ? dbuf : dp;
        if (d == dbuf)
            memset(dbuf, 0, sizeof(dbuf));

#ifndef NDEBUG
        if (d == dbuf) {
            assert((dp >= sp && dp < sp + src_p->shared->size) ||
                   (sp >= dp && sp < dp + dst_p->shared->size));
        }
        else {
            assert((dp < sp && dp + dst_p->shared->size <= sp) ||
                   (sp < dp && sp + src_p->shared->size <= dp));
        }
#endif

        /*
         * Put the data in little endian order so our loops aren't so
         * complicated. We'll do all the conversion stuff assuming
         * little endian and then we'll fix the order at the end.
         */
        if (H5T_ORDER_BE == src_atomic.order) {
            uint8_t *cur_part = s;
            /* Swap real part of complex number element */
            for (size_t j = 0; j < src_part_size / 2; j++)
                H5_SWAP_BYTES(cur_part, j, src_part_size - (j + 1));
            /* Swap imaginary part of complex number element */
            cur_part += src_part_size;
            for (size_t j = 0; j < src_part_size / 2; j++)
                H5_SWAP_BYTES(cur_part, j, src_part_size - (j + 1));
        }
        else if (H5T_ORDER_VAX == src_atomic.order)
            HGOTO_ERROR(H5E_DATATYPE, H5E_UNSUPPORTED, FAIL,
                        "VAX byte ordering is unsupported for complex number type conversions");

        /* Check for special cases: +0, -0, +Inf, -Inf, NaN */
        realval_type = H5T__conv_float_find_special(s, &src_atomic, NULL);
        imagval_type = H5T__conv_float_find_special(s + (src_p->shared->size / 2), &src_atomic, NULL);

        real_zero = (realval_type == H5T_CONV_FLOAT_SPECVAL_POSZERO ||
                     realval_type == H5T_CONV_FLOAT_SPECVAL_NEGZERO);
        imag_zero = (imagval_type == H5T_CONV_FLOAT_SPECVAL_POSZERO ||
                     imagval_type == H5T_CONV_FLOAT_SPECVAL_NEGZERO);
        real_except =
            (realval_type == H5T_CONV_FLOAT_SPECVAL_POSINF || realval_type == H5T_CONV_FLOAT_SPECVAL_NEGINF ||
             realval_type == H5T_CONV_FLOAT_SPECVAL_NAN);
        imag_except =
            (imagval_type == H5T_CONV_FLOAT_SPECVAL_POSINF || imagval_type == H5T_CONV_FLOAT_SPECVAL_NEGINF ||
             imagval_type == H5T_CONV_FLOAT_SPECVAL_NAN);

        /* A complex number is zero if both parts are +/-0 */
        if (real_zero && imag_zero) {
            H5T__bit_copy(d, dst_atomic.u.f.sign, s, src_atomic.u.f.sign, (size_t)1);
            H5T__bit_copy(d + dst_part_size, dst_atomic.u.f.sign, s + src_part_size, src_atomic.u.f.sign,
                          (size_t)1);
            H5T__bit_set(d, dst_atomic.u.f.epos, dst_atomic.u.f.esize, false);
            H5T__bit_set(d + dst_part_size, dst_atomic.u.f.epos, dst_atomic.u.f.esize, false);
            H5T__bit_set(d, dst_atomic.u.f.mpos, dst_atomic.u.f.msize, false);
            H5T__bit_set(d + dst_part_size, dst_atomic.u.f.mpos, dst_atomic.u.f.msize, false);
            goto padding;
        }
        else if (real_except || imag_except) {
            /* If user's exception handler is present, use it */
            if (conv_ctx->u.conv.cb_struct.func) {
                H5T_conv_except_t except_type; /* type of conversion exception that occurred */

                /* reverse source buffer order first */
                H5T__reverse_order(src_rev, s, src_p);

                /*
                 * A complex number is infinity if either part is infinity,
                 * even if the other part is NaN. If a part is infinity,
                 * since we can only throw one type of conversion exception,
                 * arbitrarily choose the exception type to throw based
                 * on the infinity type for the real part (if it's infinity),
                 * followed by the infinity type for the imaginary part. For
                 * now, it will be assumed that the conversion exception
                 * callback will inspect and handle both parts of the complex
                 * number value.
                 */
                if (realval_type == H5T_CONV_FLOAT_SPECVAL_POSINF)
                    except_type = H5T_CONV_EXCEPT_PINF;
                else if (realval_type == H5T_CONV_FLOAT_SPECVAL_NEGINF)
                    except_type = H5T_CONV_EXCEPT_NINF;
                else if (imagval_type == H5T_CONV_FLOAT_SPECVAL_POSINF)
                    except_type = H5T_CONV_EXCEPT_PINF;
                else if (imagval_type == H5T_CONV_FLOAT_SPECVAL_NEGINF)
                    except_type = H5T_CONV_EXCEPT_NINF;
                else {
                    assert(realval_type == H5T_CONV_FLOAT_SPECVAL_NAN ||
                           imagval_type == H5T_CONV_FLOAT_SPECVAL_NAN);
                    except_type = H5T_CONV_EXCEPT_NAN;
                }

                except_ret = (conv_ctx->u.conv.cb_struct.func)(except_type, conv_ctx->u.conv.src_type_id,
                                                               conv_ctx->u.conv.dst_type_id, src_rev, d,
                                                               conv_ctx->u.conv.cb_struct.user_data);
            }

            if (except_ret == H5T_CONV_UNHANDLED) {
                if (realval_type == H5T_CONV_FLOAT_SPECVAL_POSINF ||
                    realval_type == H5T_CONV_FLOAT_SPECVAL_NEGINF ||
                    imagval_type == H5T_CONV_FLOAT_SPECVAL_POSINF ||
                    imagval_type == H5T_CONV_FLOAT_SPECVAL_NEGINF) {
                    H5T__bit_copy(d, dst_atomic.u.f.sign, s, src_atomic.u.f.sign, (size_t)1);
                    H5T__bit_copy(d + dst_part_size, dst_atomic.u.f.sign, s + src_part_size,
                                  src_atomic.u.f.sign, (size_t)1);

                    if (realval_type == H5T_CONV_FLOAT_SPECVAL_POSINF ||
                        realval_type == H5T_CONV_FLOAT_SPECVAL_NEGINF) {
                        H5T__bit_set(d, dst_atomic.u.f.epos, dst_atomic.u.f.esize, true);
                        H5T__bit_set(d, dst_atomic.u.f.mpos, dst_atomic.u.f.msize, false);
                        /* If the destination has no implied mantissa bit, we'll need to set
                         * the 1st bit of mantissa to 1. The Intel-Linux "long double" is
                         * this case. */
                        if (H5T_NORM_NONE == dst_atomic.u.f.norm)
                            H5T__bit_set(d, dst_atomic.u.f.mpos + dst_atomic.u.f.msize - 1, (size_t)1, true);
                    }
                    if (imagval_type == H5T_CONV_FLOAT_SPECVAL_POSINF ||
                        imagval_type == H5T_CONV_FLOAT_SPECVAL_NEGINF) {
                        H5T__bit_set(d + dst_part_size, dst_atomic.u.f.epos, dst_atomic.u.f.esize, true);
                        H5T__bit_set(d + dst_part_size, dst_atomic.u.f.mpos, dst_atomic.u.f.msize, false);
                        /* If the destination has no implied mantissa bit, we'll need to set
                         * the 1st bit of mantissa to 1. The Intel-Linux "long double" is
                         * this case. */
                        if (H5T_NORM_NONE == dst_atomic.u.f.norm)
                            H5T__bit_set(d + dst_part_size, dst_atomic.u.f.mpos + dst_atomic.u.f.msize - 1,
                                         (size_t)1, true);
                    }
                }
                else {
                    /* There are many NaN values, so we just set all bits of the significand. */
                    if (realval_type == H5T_CONV_FLOAT_SPECVAL_NAN) {
                        H5T__bit_copy(d, dst_atomic.u.f.sign, s, src_atomic.u.f.sign, (size_t)1);
                        H5T__bit_set(d, dst_atomic.u.f.epos, dst_atomic.u.f.esize, true);
                        H5T__bit_set(d, dst_atomic.u.f.mpos, dst_atomic.u.f.msize, true);
                    }
                    if (imagval_type == H5T_CONV_FLOAT_SPECVAL_NAN) {
                        H5T__bit_copy(d + dst_part_size, dst_atomic.u.f.sign, s + src_part_size,
                                      src_atomic.u.f.sign, (size_t)1);
                        H5T__bit_set(d + dst_part_size, dst_atomic.u.f.epos, dst_atomic.u.f.esize, true);
                        H5T__bit_set(d + dst_part_size, dst_atomic.u.f.mpos, dst_atomic.u.f.msize, true);
                    }
                }
            }
            else if (except_ret == H5T_CONV_HANDLED) {
                /* No need to reverse the order of destination because user handles it */
                reverse = false;
                goto next;
            }
            else if (except_ret == H5T_CONV_ABORT)
                HGOTO_ERROR(H5E_DATATYPE, H5E_CANTCONVERT, FAIL, "can't handle conversion exception");

            goto padding;
        }

        if (real_zero) {
            H5T__bit_copy(d, dst_atomic.u.f.sign, s, src_atomic.u.f.sign, (size_t)1);
            H5T__bit_set(d, dst_atomic.u.f.epos, dst_atomic.u.f.esize, false);
            H5T__bit_set(d, dst_atomic.u.f.mpos, dst_atomic.u.f.msize, false);
        }
        else {
            /*
             * Get the exponent as an unsigned quantity from the section of
             * the source bit field where it's located. Don't worry about
             * the exponent bias yet.
             */
            real_expo = (int64_t)H5T__bit_get_d(s, src_atomic.u.f.epos, src_atomic.u.f.esize);

            if (real_expo == 0)
                real_denormalized = true;

            /* Determine size of mantissa for real and imaginary parts */
            if (0 == real_expo || H5T_NORM_NONE == src_atomic.u.f.norm) {
                if ((real_mant_msb =
                         H5T__bit_find(s, src_atomic.u.f.mpos, src_atomic.u.f.msize, H5T_BIT_MSB, true)) > 0)
                    real_msize = (size_t)real_mant_msb;
                else if (0 == real_mant_msb) {
                    real_msize = 1;
                    H5T__bit_set(s, src_atomic.u.f.mpos, (size_t)1, false);
                }
            }
            else if (H5T_NORM_IMPLIED == src_atomic.u.f.norm)
                real_msize = src_atomic.u.f.msize;
            else
                HGOTO_ERROR(H5E_DATATYPE, H5E_CANTCONVERT, FAIL, "normalization method not implemented yet");

            /*
             * The sign for the destination is the same as the sign for the
             * source in all cases.
             */
            H5T__bit_copy(d, dst_atomic.u.f.sign, s, src_atomic.u.f.sign, (size_t)1);

            /*
             * Calculate the true source exponent by adjusting according to
             * the source exponent bias.
             */
            if (0 == real_expo || H5T_NORM_NONE == src_atomic.u.f.norm) {
                assert(real_mant_msb >= 0);
                real_expo -=
                    (int64_t)((src_atomic.u.f.ebias - 1) + (src_atomic.u.f.msize - (size_t)real_mant_msb));
            }
            else if (H5T_NORM_IMPLIED == src_atomic.u.f.norm)
                real_expo -= (int64_t)src_atomic.u.f.ebias;
            else
                HGOTO_ERROR(H5E_DATATYPE, H5E_CANTCONVERT, FAIL, "normalization method not implemented yet");

            /*
             * If the destination is not normalized then right shift the
             * mantissa by one.
             */
            real_mrsh = 0;
            if (H5T_NORM_NONE == dst_atomic.u.f.norm)
                real_mrsh++;

            /*
             * Calculate the destination exponent by adding the destination
             * bias and clipping by the minimum and maximum possible
             * destination exponent values.
             */
            real_expo += (int64_t)dst_atomic.u.f.ebias;

            if (real_expo < -(hssize_t)(dst_atomic.u.f.msize)) {
                /* The exponent is way too small. Result is zero. */
                real_expo = 0;
                H5T__bit_set(d, dst_atomic.u.f.mpos, dst_atomic.u.f.msize, false);
                real_msize = 0;
            }
            else if (real_expo <= 0) {
                /*
                 * The exponent is too small to fit in the exponent field,
                 * but by shifting the mantissa to the right we can
                 * accommodate that value. The mantissa of course is no
                 * longer normalized.
                 */
                real_mrsh += (size_t)(1 - real_expo);
                real_expo         = 0;
                real_denormalized = true;
            }
            else if (real_expo >= expo_max) {
                /*
                 * The exponent is too large to fit in the available region
                 * or it results in the maximum possible value. Use positive
                 * or negative infinity instead unless the application
                 * specifies something else. Before calling the overflow
                 * handler make sure the source buffer we hand it is in the
                 * original byte order.
                 */
                if (conv_ctx->u.conv.cb_struct.func) { /* If user's exception handler is present, use it */
                    /* reverse source buffer order first */
                    H5T__reverse_order(src_rev, s, src_p);

                    except_ret = (conv_ctx->u.conv.cb_struct.func)(
                        H5T_CONV_EXCEPT_RANGE_HI, conv_ctx->u.conv.src_type_id, conv_ctx->u.conv.dst_type_id,
                        src_rev, d, conv_ctx->u.conv.cb_struct.user_data);
                }

                if (except_ret == H5T_CONV_UNHANDLED) {
                    real_expo = expo_max;
                    H5T__bit_set(d, dst_atomic.u.f.mpos, dst_atomic.u.f.msize, false);
                    real_msize = 0;
                }
                else if (except_ret == H5T_CONV_HANDLED) {
                    reverse = false;
                    goto next;
                }
                else if (except_ret == H5T_CONV_ABORT)
                    HGOTO_ERROR(H5E_DATATYPE, H5E_CANTCONVERT, FAIL, "can't handle conversion exception");
            }

            /*
             * If the destination mantissa is smaller than the source
             * mantissa then round the source mantissa. Rounding may cause a
             * carry in which case the exponent has to be re-evaluated for
             * overflow. That is, if `carry' is clear then the implied
             * mantissa bit is `1', else it is `10' binary.
             */
            real_implied = 1;
            mpos         = src_atomic.u.f.mpos;
            if (real_msize > 0 && real_mrsh <= dst_atomic.u.f.msize &&
                real_mrsh + real_msize > dst_atomic.u.f.msize) {
                real_mant_msb = (ssize_t)(real_mrsh + real_msize - dst_atomic.u.f.msize);
                assert(real_mant_msb >= 0 && (size_t)real_mant_msb <= real_msize);
                /* If the 1st bit being cut off is set and source isn't denormalized. */
                if (H5T__bit_get_d(s, (mpos + (size_t)real_mant_msb) - 1, (size_t)1) && !real_denormalized) {
                    /* Don't do rounding if exponent is 111...110 and mantissa is 111...11.
                     * To do rounding and increment exponent in this case will create an infinity value. */
                    if ((H5T__bit_find(s, mpos + (size_t)real_mant_msb, real_msize - (size_t)real_mant_msb,
                                       H5T_BIT_LSB, false) >= 0 ||
                         real_expo < expo_max - 1)) {
                        real_carry = H5T__bit_inc(s, mpos + (size_t)real_mant_msb - 1,
                                                  1 + real_msize - (size_t)real_mant_msb);
                        if (real_carry)
                            real_implied = 2;
                    }
                }
                else if (H5T__bit_get_d(s, (mpos + (size_t)real_mant_msb) - 1, (size_t)1) &&
                         real_denormalized)
                    /* For either source or destination, denormalized value doesn't increment carry. */
                    H5T__bit_inc(s, mpos + (size_t)real_mant_msb - 1, 1 + real_msize - (size_t)real_mant_msb);
            }
            else
                real_carry = false;

            /* Write the mantissa to the destination */
            if (real_mrsh > dst_atomic.u.f.msize + 1) {
                H5T__bit_set(d, dst_atomic.u.f.mpos, dst_atomic.u.f.msize, false);
            }
            else if (real_mrsh == dst_atomic.u.f.msize + 1) {
                H5T__bit_set(d, dst_atomic.u.f.mpos + 1, dst_atomic.u.f.msize - 1, false);
                H5T__bit_set(d, dst_atomic.u.f.mpos, (size_t)1, true);
            }
            else if (real_mrsh == dst_atomic.u.f.msize) {
                H5T__bit_set(d, dst_atomic.u.f.mpos, dst_atomic.u.f.msize, false);
                H5T__bit_set_d(d, dst_atomic.u.f.mpos, MIN(2, dst_atomic.u.f.msize), (hsize_t)real_implied);
            }
            else {
                if (real_mrsh > 0) {
                    H5T__bit_set(d, dst_atomic.u.f.mpos + dst_atomic.u.f.msize - real_mrsh, real_mrsh, false);
                    H5T__bit_set_d(d, dst_atomic.u.f.mpos + dst_atomic.u.f.msize - real_mrsh, (size_t)2,
                                   (hsize_t)real_implied);
                }
                if (real_mrsh + real_msize >= dst_atomic.u.f.msize) {
                    H5T__bit_copy(d, dst_atomic.u.f.mpos, s,
                                  (mpos + real_msize + real_mrsh - dst_atomic.u.f.msize),
                                  dst_atomic.u.f.msize - real_mrsh);
                }
                else {
                    H5T__bit_copy(d, dst_atomic.u.f.mpos + dst_atomic.u.f.msize - (real_mrsh + real_msize), s,
                                  mpos, real_msize);
                    H5T__bit_set(d, dst_atomic.u.f.mpos, dst_atomic.u.f.msize - (real_mrsh + real_msize),
                                 false);
                }
            }

            /* Write the exponent */
            if (real_carry) {
                real_expo++;
                if (real_expo >= expo_max) {
                    /*
                     * The exponent is too large to fit in the available
                     * region or it results in the maximum possible value.
                     * Use positive or negative infinity instead unless the
                     * application specifies something else. Before calling
                     * the overflow handler make sure the source buffer we
                     * hand it is in the original byte order.
                     */
                    if (conv_ctx->u.conv.cb_struct
                            .func) { /* If user's exception handler is present, use it */
                        /* reverse source buffer order first */
                        H5T__reverse_order(src_rev, s, src_p);

                        except_ret = (conv_ctx->u.conv.cb_struct.func)(
                            H5T_CONV_EXCEPT_RANGE_HI, conv_ctx->u.conv.src_type_id,
                            conv_ctx->u.conv.dst_type_id, src_rev, d, conv_ctx->u.conv.cb_struct.user_data);
                    }

                    if (except_ret == H5T_CONV_UNHANDLED) {
                        real_expo = expo_max;
                        H5T__bit_set(d, dst_atomic.u.f.mpos, dst_atomic.u.f.msize, false);
                    }
                    else if (except_ret == H5T_CONV_HANDLED) {
                        reverse = false;
                        goto next;
                    }
                    else if (except_ret == H5T_CONV_ABORT)
                        HGOTO_ERROR(H5E_DATATYPE, H5E_CANTCONVERT, FAIL, "can't handle conversion exception");
                }
            }

            real_carry = false;

            H5_CHECK_OVERFLOW(real_expo, hssize_t, hsize_t);
            H5T__bit_set_d(d, dst_atomic.u.f.epos, dst_atomic.u.f.esize, (hsize_t)real_expo);
        }

        if (imag_zero) {
            H5T__bit_copy(d + dst_part_size, dst_atomic.u.f.sign, s + src_part_size, src_atomic.u.f.sign,
                          (size_t)1);
            H5T__bit_set(d + dst_part_size, dst_atomic.u.f.epos, dst_atomic.u.f.esize, false);
            H5T__bit_set(d + dst_part_size, dst_atomic.u.f.mpos, dst_atomic.u.f.msize, false);
        }
        else {
            /*
             * Get the exponent as an unsigned quantity from the section of
             * the source bit field where it's located. Don't worry about
             * the exponent bias yet.
             */
            imag_expo = (int64_t)H5T__bit_get_d(s + src_part_size, src_atomic.u.f.epos, src_atomic.u.f.esize);

            if (imag_expo == 0)
                imag_denormalized = true;

            if (0 == imag_expo || H5T_NORM_NONE == src_atomic.u.f.norm) {
                if ((imag_mant_msb = H5T__bit_find(s + src_part_size, src_atomic.u.f.mpos,
                                                   src_atomic.u.f.msize, H5T_BIT_MSB, true)) > 0)
                    imag_msize = (size_t)imag_mant_msb;
                else if (0 == imag_mant_msb) {
                    imag_msize = 1;
                    H5T__bit_set(s + src_part_size, src_atomic.u.f.mpos, (size_t)1, false);
                }
            }
            else if (H5T_NORM_IMPLIED == src_atomic.u.f.norm)
                imag_msize = src_atomic.u.f.msize;
            else
                HGOTO_ERROR(H5E_DATATYPE, H5E_CANTCONVERT, FAIL, "normalization method not implemented yet");

            /*
             * The sign for the destination is the same as the sign for the
             * source in all cases.
             */
            H5T__bit_copy(d + dst_part_size, dst_atomic.u.f.sign, s + src_part_size, src_atomic.u.f.sign,
                          (size_t)1);

            /*
             * Calculate the true source exponent by adjusting according to
             * the source exponent bias.
             */
            if (0 == imag_expo || H5T_NORM_NONE == src_atomic.u.f.norm) {
                assert(imag_mant_msb >= 0);
                imag_expo -=
                    (int64_t)((src_atomic.u.f.ebias - 1) + (src_atomic.u.f.msize - (size_t)imag_mant_msb));
            }
            else if (H5T_NORM_IMPLIED == src_atomic.u.f.norm)
                imag_expo -= (int64_t)src_atomic.u.f.ebias;
            else
                HGOTO_ERROR(H5E_DATATYPE, H5E_CANTCONVERT, FAIL, "normalization method not implemented yet");

            /*
             * If the destination is not normalized then right shift the
             * mantissa by one.
             */
            imag_mrsh = 0;
            if (H5T_NORM_NONE == dst_atomic.u.f.norm)
                imag_mrsh++;

            /*
             * Calculate the destination exponent by adding the destination
             * bias and clipping by the minimum and maximum possible
             * destination exponent values.
             */
            imag_expo += (int64_t)dst_atomic.u.f.ebias;

            if (imag_expo < -(hssize_t)(dst_atomic.u.f.msize)) {
                /* The exponent is way too small. Result is zero. */
                imag_expo = 0;
                H5T__bit_set(d + dst_part_size, dst_atomic.u.f.mpos, dst_atomic.u.f.msize, false);
                imag_msize = 0;
            }
            else if (imag_expo <= 0) {
                /*
                 * The exponent is too small to fit in the exponent field,
                 * but by shifting the mantissa to the right we can
                 * accommodate that value. The mantissa of course is no
                 * longer normalized.
                 */
                imag_mrsh += (size_t)(1 - imag_expo);
                imag_expo         = 0;
                imag_denormalized = true;
            }
            else if (imag_expo >= expo_max) {
                /*
                 * The exponent is too large to fit in the available region
                 * or it results in the maximum possible value. Use positive
                 * or negative infinity instead unless the application
                 * specifies something else. Before calling the overflow
                 * handler make sure the source buffer we hand it is in the
                 * original byte order.
                 */
                if (conv_ctx->u.conv.cb_struct.func) { /* If user's exception handler is present, use it */
                    /* reverse source buffer order first */
                    H5T__reverse_order(src_rev, s, src_p);

                    except_ret = (conv_ctx->u.conv.cb_struct.func)(
                        H5T_CONV_EXCEPT_RANGE_HI, conv_ctx->u.conv.src_type_id, conv_ctx->u.conv.dst_type_id,
                        src_rev, d, conv_ctx->u.conv.cb_struct.user_data);
                }

                if (except_ret == H5T_CONV_UNHANDLED) {
                    imag_expo = expo_max;
                    H5T__bit_set(d + dst_part_size, dst_atomic.u.f.mpos, dst_atomic.u.f.msize, false);
                    imag_msize = 0;
                }
                else if (except_ret == H5T_CONV_HANDLED) {
                    reverse = false;
                    goto next;
                }
                else if (except_ret == H5T_CONV_ABORT)
                    HGOTO_ERROR(H5E_DATATYPE, H5E_CANTCONVERT, FAIL, "can't handle conversion exception");
            }

            /*
             * If the destination mantissa is smaller than the source
             * mantissa then round the source mantissa. Rounding may cause a
             * carry in which case the exponent has to be re-evaluated for
             * overflow. That is, if `carry' is clear then the implied
             * mantissa bit is `1', else it is `10' binary.
             */
            imag_implied = 1;
            mpos         = src_atomic.u.f.mpos;
            if (imag_msize > 0 && imag_mrsh <= dst_atomic.u.f.msize &&
                imag_mrsh + imag_msize > dst_atomic.u.f.msize) {
                imag_mant_msb = (ssize_t)(imag_mrsh + imag_msize - dst_atomic.u.f.msize);
                assert(imag_mant_msb >= 0 && (size_t)imag_mant_msb <= imag_msize);
                /* If the 1st bit being cut off is set and source isn't denormalized. */
                if (H5T__bit_get_d(s + src_part_size, (mpos + (size_t)imag_mant_msb) - 1, (size_t)1) &&
                    !imag_denormalized) {
                    /* Don't do rounding if exponent is 111...110 and mantissa is 111...11.
                     * To do rounding and increment exponent in this case will create an infinity value. */
                    if ((H5T__bit_find(s + src_part_size, mpos + (size_t)imag_mant_msb,
                                       imag_msize - (size_t)imag_mant_msb, H5T_BIT_LSB, false) >= 0 ||
                         imag_expo < expo_max - 1)) {
                        imag_carry = H5T__bit_inc(s + src_part_size, mpos + (size_t)imag_mant_msb - 1,
                                                  1 + imag_msize - (size_t)imag_mant_msb);
                        if (imag_carry)
                            imag_implied = 2;
                    }
                }
                else if (H5T__bit_get_d(s + src_part_size, (mpos + (size_t)imag_mant_msb) - 1, (size_t)1) &&
                         imag_denormalized)
                    /* For either source or destination, denormalized value doesn't increment carry. */
                    H5T__bit_inc(s + src_part_size, mpos + (size_t)imag_mant_msb - 1,
                                 1 + imag_msize - (size_t)imag_mant_msb);
            }
            else
                imag_carry = false;

            /* Write the mantissa to the destination */
            if (imag_mrsh > dst_atomic.u.f.msize + 1) {
                H5T__bit_set(d + dst_part_size, dst_atomic.u.f.mpos, dst_atomic.u.f.msize, false);
            }
            else if (imag_mrsh == dst_atomic.u.f.msize + 1) {
                H5T__bit_set(d + dst_part_size, dst_atomic.u.f.mpos + 1, dst_atomic.u.f.msize - 1, false);
                H5T__bit_set(d + dst_part_size, dst_atomic.u.f.mpos, (size_t)1, true);
            }
            else if (imag_mrsh == dst_atomic.u.f.msize) {
                H5T__bit_set(d + dst_part_size, dst_atomic.u.f.mpos, dst_atomic.u.f.msize, false);
                H5T__bit_set_d(d + dst_part_size, dst_atomic.u.f.mpos, MIN(2, dst_atomic.u.f.msize),
                               (hsize_t)imag_implied);
            }
            else {
                if (imag_mrsh > 0) {
                    H5T__bit_set(d + dst_part_size, dst_atomic.u.f.mpos + dst_atomic.u.f.msize - imag_mrsh,
                                 imag_mrsh, false);
                    H5T__bit_set_d(d + dst_part_size, dst_atomic.u.f.mpos + dst_atomic.u.f.msize - imag_mrsh,
                                   (size_t)2, (hsize_t)imag_implied);
                }
                if (imag_mrsh + imag_msize >= dst_atomic.u.f.msize) {
                    H5T__bit_copy(d + dst_part_size, dst_atomic.u.f.mpos, s + src_part_size,
                                  (mpos + imag_msize + imag_mrsh - dst_atomic.u.f.msize),
                                  dst_atomic.u.f.msize - imag_mrsh);
                }
                else {
                    H5T__bit_copy(d + dst_part_size,
                                  dst_atomic.u.f.mpos + dst_atomic.u.f.msize - (imag_mrsh + imag_msize),
                                  s + src_part_size, mpos, imag_msize);
                    H5T__bit_set(d + dst_part_size, dst_atomic.u.f.mpos,
                                 dst_atomic.u.f.msize - (imag_mrsh + imag_msize), false);
                }
            }

            /* Write the exponent */
            if (imag_carry) {
                imag_expo++;
                if (imag_expo >= expo_max) {
                    /*
                     * The exponent is too large to fit in the available
                     * region or it results in the maximum possible value.
                     * Use positive or negative infinity instead unless the
                     * application specifies something else. Before calling
                     * the overflow handler make sure the source buffer we
                     * hand it is in the original byte order.
                     */
                    if (conv_ctx->u.conv.cb_struct
                            .func) { /* If user's exception handler is present, use it */
                        /* reverse source buffer order first */
                        H5T__reverse_order(src_rev, s, src_p);

                        except_ret = (conv_ctx->u.conv.cb_struct.func)(
                            H5T_CONV_EXCEPT_RANGE_HI, conv_ctx->u.conv.src_type_id,
                            conv_ctx->u.conv.dst_type_id, src_rev, d, conv_ctx->u.conv.cb_struct.user_data);
                    }

                    if (except_ret == H5T_CONV_UNHANDLED) {
                        imag_expo = expo_max;
                        H5T__bit_set(d + dst_part_size, dst_atomic.u.f.mpos, dst_atomic.u.f.msize, false);
                    }
                    else if (except_ret == H5T_CONV_HANDLED) {
                        reverse = false;
                        goto next;
                    }
                    else if (except_ret == H5T_CONV_ABORT)
                        HGOTO_ERROR(H5E_DATATYPE, H5E_CANTCONVERT, FAIL, "can't handle conversion exception");
                }
            }

            imag_carry = false;

            H5_CHECK_OVERFLOW(imag_expo, hssize_t, hsize_t);
            H5T__bit_set_d(d + dst_part_size, dst_atomic.u.f.epos, dst_atomic.u.f.esize, (hsize_t)imag_expo);
        }

padding:
        /*
         * Set external padding areas
         */
        if (dst_atomic.offset > 0) {
            assert(H5T_PAD_ZERO == dst_atomic.lsb_pad || H5T_PAD_ONE == dst_atomic.lsb_pad);
            H5T__bit_set(d, (size_t)0, dst_atomic.offset, (bool)(H5T_PAD_ONE == dst_atomic.lsb_pad));
            H5T__bit_set(d + dst_part_size, (size_t)0, dst_atomic.offset,
                         (bool)(H5T_PAD_ONE == dst_atomic.lsb_pad));
        }
        {
            size_t type_size = dst_p->shared->parent->shared->size;

            if (dst_atomic.offset + dst_atomic.prec != 8 * type_size) {
                assert(H5T_PAD_ZERO == dst_atomic.msb_pad || H5T_PAD_ONE == dst_atomic.msb_pad);
                H5T__bit_set(d, dst_atomic.offset + dst_atomic.prec,
                             8 * type_size - (dst_atomic.offset + dst_atomic.prec),
                             (bool)(H5T_PAD_ONE == dst_atomic.msb_pad));
                H5T__bit_set(d + dst_part_size, dst_atomic.offset + dst_atomic.prec,
                             8 * type_size - (dst_atomic.offset + dst_atomic.prec),
                             (bool)(H5T_PAD_ONE == dst_atomic.msb_pad));
            }
        }

        /* Put the destination in the correct byte order. See note at beginning of loop. */
        if (H5T_ORDER_BE == dst_atomic.order && reverse) {
            uint8_t *cur_part = d;
            /* Swap real part of complex number element */
            for (size_t j = 0; j < dst_part_size / 2; j++)
                H5_SWAP_BYTES(cur_part, j, dst_part_size - (j + 1));
            /* Swap imaginary part of complex number element */
            cur_part += dst_part_size;
            for (size_t j = 0; j < dst_part_size / 2; j++)
                H5_SWAP_BYTES(cur_part, j, dst_part_size - (j + 1));
        }
        else if (H5T_ORDER_VAX == dst_atomic.order)
            HGOTO_ERROR(H5E_DATATYPE, H5E_UNSUPPORTED, FAIL,
                        "VAX byte ordering is unsupported for complex number type conversions");

next:
        /*
         * If we had used a temporary buffer for the destination then we
         * should copy the value to the true destination buffer.
         */
        if (d == dbuf)
            H5MM_memcpy(dp, d, dst_p->shared->size);

        /* Advance source & destination pointers by delta amounts */
        sp += src_delta;
        dp += dst_delta;
    } /* end conversion loop */

done:
    H5MM_free(src_rev);

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5T__conv_complex_loop() */

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_complex_i
 *
 * Purpose:     Convert complex number values to integer values. This is
 *              the catch-all function for complex number -> integer
 *              conversions and is probably not particularly fast.
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_complex_i(const H5T_t *src_p, const H5T_t *dst_p, H5T_cdata_t *cdata,
                    const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                    size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    herr_t ret_value = SUCCEED;

    FUNC_ENTER_PACKAGE

    switch (cdata->command) {
        case H5T_CONV_INIT: {
            H5T_atomic_t src_atomic; /* source datatype atomic info      */
            H5T_atomic_t dst_atomic; /* destination datatype atomic info */

            if (!src_p || !dst_p)
                HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a datatype");
            if (!H5T_IS_ATOMIC(src_p->shared->parent->shared))
                HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "invalid complex number datatype");
            src_atomic = src_p->shared->parent->shared->u.atomic;
            dst_atomic = dst_p->shared->u.atomic;
            if (H5T_ORDER_LE != src_atomic.order && H5T_ORDER_BE != src_atomic.order &&
                H5T_ORDER_VAX != src_atomic.order)
                HGOTO_ERROR(H5E_DATATYPE, H5E_UNSUPPORTED, FAIL,
                            "unsupported byte order for source datatype");
            if (H5T_ORDER_LE != dst_atomic.order && H5T_ORDER_BE != dst_atomic.order &&
                H5T_ORDER_VAX != dst_atomic.order)
                HGOTO_ERROR(H5E_DATATYPE, H5E_UNSUPPORTED, FAIL,
                            "unsupported byte order for destination datatype");
            if (dst_p->shared->size > TEMP_INT_CONV_BUFFER_SIZE)
                HGOTO_ERROR(H5E_DATATYPE, H5E_UNSUPPORTED, FAIL, "destination datatype size is too large");
            if (8 * sizeof(hssize_t) - 1 < src_atomic.u.f.esize)
                HGOTO_ERROR(H5E_DATATYPE, H5E_UNSUPPORTED, FAIL, "exponent field is too large");
            cdata->need_bkg = H5T_BKG_NO;

            break;
        }

        case H5T_CONV_FREE:
            break;

        case H5T_CONV_CONV:
            if (!src_p || !dst_p)
                HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a datatype");
            if (NULL == conv_ctx)
                HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "invalid datatype conversion context pointer");

            if (H5T__conv_f_i_loop(src_p, dst_p, conv_ctx, nelmts, buf_stride, buf) < 0)
                HGOTO_ERROR(H5E_DATATYPE, H5E_CANTCONVERT, FAIL, "unable to convert data values");

            break;

        default:
            HGOTO_ERROR(H5E_DATATYPE, H5E_UNSUPPORTED, FAIL, "unknown conversion command");
    }

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5T__conv_complex_i() */

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_complex_f
 *
 * Purpose:     Convert complex number values to floating-point values.
 *              This is the catch-all function for complex number -> float
 *              conversions and is probably not particularly fast.
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_complex_f(const H5T_t *src_p, const H5T_t *dst_p, H5T_cdata_t *cdata,
                    const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                    size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    bool   equal_cplx_conv = false; /* if converting between complex and matching float */
    herr_t ret_value       = SUCCEED;

    FUNC_ENTER_PACKAGE

    switch (cdata->command) {
        case H5T_CONV_INIT: {
            H5T_atomic_t src_atomic; /* source datatype atomic info                      */
            H5T_atomic_t dst_atomic; /* destination datatype atomic info                 */

            if (!src_p || !dst_p)
                HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a datatype");
            if (!H5T_IS_ATOMIC(src_p->shared->parent->shared))
                HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "invalid complex number datatype");
            src_atomic = src_p->shared->parent->shared->u.atomic;
            dst_atomic = dst_p->shared->u.atomic;
            if (H5T_ORDER_LE != src_atomic.order && H5T_ORDER_BE != src_atomic.order &&
                H5T_ORDER_VAX != src_atomic.order)
                HGOTO_ERROR(H5E_DATATYPE, H5E_UNSUPPORTED, FAIL,
                            "unsupported byte order for source datatype");
            if (H5T_ORDER_LE != dst_atomic.order && H5T_ORDER_BE != dst_atomic.order &&
                H5T_ORDER_VAX != dst_atomic.order)
                HGOTO_ERROR(H5E_DATATYPE, H5E_UNSUPPORTED, FAIL,
                            "unsupported byte order for destination datatype");
            if (dst_p->shared->size > TEMP_FLOAT_CONV_BUFFER_SIZE)
                HGOTO_ERROR(H5E_DATATYPE, H5E_UNSUPPORTED, FAIL, "destination datatype size is too large");
            if (8 * sizeof(int64_t) - 1 < src_atomic.u.f.esize ||
                8 * sizeof(int64_t) - 1 < dst_atomic.u.f.esize)
                HGOTO_ERROR(H5E_DATATYPE, H5E_UNSUPPORTED, FAIL, "exponent field is too large");
            cdata->need_bkg = H5T_BKG_NO;

            break;
        }

        case H5T_CONV_FREE:
            break;

        case H5T_CONV_CONV:
            if (!src_p || !dst_p)
                HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a datatype");
            if (NULL == conv_ctx)
                HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "invalid datatype conversion context pointer");

            /* Are we converting between a floating-point type and a complex number
             * type consisting of the same floating-point type?
             */
            /* TODO: cache result */
            equal_cplx_conv = (0 == H5T_cmp(src_p->shared->parent, dst_p, false));
            if (!equal_cplx_conv) {
                /* If floating-point types differ, use generic f_f loop */
                if (H5T__conv_f_f_loop(src_p, dst_p, conv_ctx, nelmts, buf_stride, buf) < 0)
                    HGOTO_ERROR(H5E_DATATYPE, H5E_CANTCONVERT, FAIL, "unable to convert data values");
            }
            else {
                /* If floating-point types are the same, use specialized loop */
                if (H5T__conv_complex_f_matched(src_p, dst_p, conv_ctx, nelmts, buf_stride, buf) < 0)
                    HGOTO_ERROR(H5E_DATATYPE, H5E_CANTCONVERT, FAIL, "unable to convert data values");
            }

            break;

        default:
            HGOTO_ERROR(H5E_DATATYPE, H5E_UNSUPPORTED, FAIL, "unknown conversion command");
    }

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5T__conv_complex_f() */

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_complex_f_matched
 *
 * Purpose:     Implements the body of the conversion loop when converting
 *              between a floating-point type and a complex number type
 *              consisting of the same floating-point type. Encapsulates
 *              common code that is shared between the H5T__conv_complex_f
 *              and H5T__conv_f_complex functions. Values can be directly
 *              converted between the types after checking for conversion
 *              exceptions.
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_complex_f_matched(const H5T_t *src_p, const H5T_t *dst_p, const H5T_conv_ctx_t *conv_ctx,
                            size_t nelmts, size_t buf_stride, void *buf)
{
    H5T_conv_float_specval_t specval_type;      /* floating-point value type (regular, +/-Inf, +/-0, NaN) */
    H5T_atomic_t             src_atomic;        /* source datatype atomic info                            */
    H5T_atomic_t             dst_atomic;        /* destination datatype atomic info                       */
    ssize_t  src_delta, dst_delta;              /* source & destination stride                            */
    uint8_t *s, *sp, *d, *dp;                   /* source and dest traversal ptrs                         */
    uint8_t *src_rev = NULL;                    /* order-reversed source buffer                           */
    uint8_t  dbuf[TEMP_FLOAT_CONV_BUFFER_SIZE]; /* temp destination buffer                                */
    size_t   olap;                              /* num overlapping elements                               */
    int      direction;                         /* forward or backward traversal                          */
    herr_t   ret_value = SUCCEED;

    assert(src_p);
    assert(src_p->shared->type == H5T_FLOAT || src_p->shared->type == H5T_COMPLEX);
    assert(dst_p);
    assert(dst_p->shared->type == H5T_FLOAT || dst_p->shared->type == H5T_COMPLEX);
    assert(conv_ctx);
    assert(buf);

    FUNC_ENTER_PACKAGE

    if (src_p->shared->type == H5T_COMPLEX)
        src_atomic = src_p->shared->parent->shared->u.atomic;
    else
        src_atomic = src_p->shared->u.atomic;
    if (dst_p->shared->type == H5T_COMPLEX)
        dst_atomic = dst_p->shared->parent->shared->u.atomic;
    else
        dst_atomic = dst_p->shared->u.atomic;

#ifndef NDEBUG
    {
        /* Make sure the floating-point types match */
        const H5T_t *src_base = (src_p->shared->type == H5T_FLOAT) ? src_p : src_p->shared->parent;
        const H5T_t *dst_base = (dst_p->shared->type == H5T_FLOAT) ? dst_p : dst_p->shared->parent;
        assert(0 == (H5T_cmp(src_base, dst_base, false)));
    }
#endif

    /*
     * Do we process the values from beginning to end or vice versa? Also,
     * how many of the elements have the source and destination areas
     * overlapping?
     */
    if (src_p->shared->size == dst_p->shared->size || buf_stride) {
        sp = dp   = (uint8_t *)buf;
        direction = 1;
        olap      = nelmts;
    }
    else if (src_p->shared->size >= dst_p->shared->size) {
        double olap_d =
            ceil((double)(dst_p->shared->size) / (double)(src_p->shared->size - dst_p->shared->size));
        olap = (size_t)olap_d;
        sp = dp   = (uint8_t *)buf;
        direction = 1;
    }
    else {
        double olap_d =
            ceil((double)(src_p->shared->size) / (double)(dst_p->shared->size - src_p->shared->size));
        olap      = (size_t)olap_d;
        sp        = (uint8_t *)buf + (nelmts - 1) * src_p->shared->size;
        dp        = (uint8_t *)buf + (nelmts - 1) * dst_p->shared->size;
        direction = -1;
    }

    /* Direction & size of buffer traversal */
    H5_CHECK_OVERFLOW(buf_stride, size_t, ssize_t);
    H5_CHECK_OVERFLOW(src_p->shared->size, size_t, ssize_t);
    H5_CHECK_OVERFLOW(dst_p->shared->size, size_t, ssize_t);
    src_delta = (ssize_t)direction * (ssize_t)(buf_stride ? buf_stride : src_p->shared->size);
    dst_delta = (ssize_t)direction * (ssize_t)(buf_stride ? buf_stride : dst_p->shared->size);

    /* Allocate space for order-reversed source buffer */
    if (conv_ctx->u.conv.cb_struct.func)
        if (NULL == (src_rev = H5MM_calloc(src_p->shared->size)))
            HGOTO_ERROR(H5E_DATATYPE, H5E_CANTALLOC, FAIL, "couldn't allocate temporary buffer");

    /* The conversion loop */
    for (size_t elmtno = 0; elmtno < nelmts; elmtno++) {
        H5T_conv_ret_t except_ret = H5T_CONV_UNHANDLED; /* return of conversion exception callback function */
        bool           reverse    = true;               /* if reversed the order of destination             */

        /*
         * If the source and destination buffers overlap then use a
         * temporary buffer for the destination.
         */
        s = sp;
        if (direction > 0)
            d = elmtno < olap ? dbuf : dp;
        else
            d = elmtno + olap >= nelmts ? dbuf : dp;
        if (d == dbuf)
            memset(dbuf, 0, sizeof(dbuf));

#ifndef NDEBUG
        if (d == dbuf) {
            assert((dp >= sp && dp < sp + src_p->shared->size) ||
                   (sp >= dp && sp < dp + dst_p->shared->size));
        }
        else {
            assert((dp < sp && dp + dst_p->shared->size <= sp) ||
                   (sp < dp && sp + src_p->shared->size <= dp));
        }
#endif

        /*
         * Put the data in little endian order so our loops aren't so
         * complicated. We'll do all the conversion stuff assuming
         * little endian and then we'll fix the order at the end.
         */
        if (H5T_ORDER_BE == src_atomic.order) {
            size_t half_size = src_p->shared->size / 2;

            if (H5T_FLOAT == src_p->shared->type) {
                for (size_t j = 0; j < half_size; j++)
                    H5_SWAP_BYTES(s, j, src_p->shared->size - (j + 1));
            }
            else {
                uint8_t *cur_part = s;
                /* Swap real part of complex number element */
                for (size_t j = 0; j < half_size / 2; j++)
                    H5_SWAP_BYTES(cur_part, j, half_size - (j + 1));
                /* Swap imaginary part of complex number element */
                cur_part += half_size;
                for (size_t j = 0; j < half_size / 2; j++)
                    H5_SWAP_BYTES(cur_part, j, half_size - (j + 1));
            }
        }
        else if (H5T_ORDER_VAX == src_atomic.order) {
            if (H5T_FLOAT == src_p->shared->type) {
                uint8_t tmp1, tmp2;
                size_t  tsize = src_p->shared->size;
                assert(0 == tsize % 2);

                for (size_t i = 0; i < tsize; i += 4) {
                    tmp1 = s[i];
                    tmp2 = s[i + 1];

                    s[i]     = s[(tsize - 2) - i];
                    s[i + 1] = s[(tsize - 1) - i];

                    s[(tsize - 2) - i] = tmp1;
                    s[(tsize - 1) - i] = tmp2;
                }
            }
            else
                HGOTO_ERROR(H5E_DATATYPE, H5E_UNSUPPORTED, FAIL,
                            "VAX byte ordering is unsupported for complex number type conversions");
        }

        /* Check for special cases: +0, -0, +Inf, -Inf, NaN */
        specval_type = H5T__conv_float_find_special(s, &src_atomic, NULL);
        if (specval_type == H5T_CONV_FLOAT_SPECVAL_POSZERO ||
            specval_type == H5T_CONV_FLOAT_SPECVAL_NEGZERO) {
            H5T__bit_copy(d, dst_atomic.u.f.sign, s, src_atomic.u.f.sign, (size_t)1);
            H5T__bit_set(d, dst_atomic.u.f.epos, dst_atomic.u.f.esize, false);
            H5T__bit_set(d, dst_atomic.u.f.mpos, dst_atomic.u.f.msize, false);
            goto padding;
        }
        else if (specval_type != H5T_CONV_FLOAT_SPECVAL_REGULAR) {
            /* If user's exception handler is present, use it */
            if (conv_ctx->u.conv.cb_struct.func) {
                H5T_conv_except_t except_type; /* type of conversion exception that occurred */

                /* reverse source buffer order first */
                H5T__reverse_order(src_rev, s, src_p);

                if (specval_type == H5T_CONV_FLOAT_SPECVAL_POSINF)
                    except_type = H5T_CONV_EXCEPT_PINF;
                else if (specval_type == H5T_CONV_FLOAT_SPECVAL_NEGINF)
                    except_type = H5T_CONV_EXCEPT_NINF;
                else
                    except_type = H5T_CONV_EXCEPT_NAN;

                except_ret = (conv_ctx->u.conv.cb_struct.func)(except_type, conv_ctx->u.conv.src_type_id,
                                                               conv_ctx->u.conv.dst_type_id, src_rev, d,
                                                               conv_ctx->u.conv.cb_struct.user_data);
            }

            if (except_ret == H5T_CONV_UNHANDLED) {
                H5T__bit_copy(d, dst_atomic.u.f.sign, s, src_atomic.u.f.sign, (size_t)1);
                H5T__bit_set(d, dst_atomic.u.f.epos, dst_atomic.u.f.esize, true);
                if (specval_type == H5T_CONV_FLOAT_SPECVAL_NAN)
                    /* There are many NaN values, so we just set all bits of the significand. */
                    H5T__bit_set(d, dst_atomic.u.f.mpos, dst_atomic.u.f.msize, true);
                else {
                    /* +/-Inf */
                    H5T__bit_set(d, dst_atomic.u.f.mpos, dst_atomic.u.f.msize, false);
                    /* If the destination has no implied mantissa bit, we'll need to set
                     * the 1st bit of mantissa to 1. The Intel-Linux "long double" is
                     * this case. */
                    if (H5T_NORM_NONE == dst_atomic.u.f.norm)
                        H5T__bit_set(d, dst_atomic.u.f.mpos + dst_atomic.u.f.msize - 1, (size_t)1, true);
                }
            }
            else if (except_ret == H5T_CONV_HANDLED) {
                /* No need to reverse the order of destination because user handles it */
                reverse = false;
                goto next;
            }
            else if (except_ret == H5T_CONV_ABORT)
                HGOTO_ERROR(H5E_DATATYPE, H5E_CANTCONVERT, FAIL, "can't handle conversion exception");

            goto padding;
        }

        /* Direct copy between complex number and floating-point type */
        if (H5T_FLOAT == src_p->shared->type)
            memcpy(d, s, src_p->shared->size);
        else
            memcpy(d, s, src_p->shared->size / 2);

padding:
        /*
         * Set external padding areas
         */
        if (dst_atomic.offset > 0) {
            assert(H5T_PAD_ZERO == dst_atomic.lsb_pad || H5T_PAD_ONE == dst_atomic.lsb_pad);
            H5T__bit_set(d, (size_t)0, dst_atomic.offset, (bool)(H5T_PAD_ONE == dst_atomic.lsb_pad));
        }
        {
            size_t type_size;

            if (dst_p->shared->type == H5T_FLOAT)
                type_size = dst_p->shared->size;
            else
                type_size = dst_p->shared->parent->shared->size;

            if (dst_atomic.offset + dst_atomic.prec != 8 * type_size) {
                assert(H5T_PAD_ZERO == dst_atomic.msb_pad || H5T_PAD_ONE == dst_atomic.msb_pad);
                H5T__bit_set(d, dst_atomic.offset + dst_atomic.prec,
                             8 * type_size - (dst_atomic.offset + dst_atomic.prec),
                             (bool)(H5T_PAD_ONE == dst_atomic.msb_pad));
            }
        }

        /*
         * Put the destination in the correct byte order. See note at
         * beginning of loop. Only the "real" part of a complex number
         * element is swapped. By the C standard, the "imaginary" part
         * should just be zeroed when converting a real value to a
         * complex value.
         */
        if (H5T_ORDER_BE == dst_atomic.order && reverse) {
            size_t half_size = dst_p->shared->size / 2;

            if (H5T_FLOAT == dst_p->shared->type) {
                for (size_t j = 0; j < half_size; j++)
                    H5_SWAP_BYTES(d, j, dst_p->shared->size - (j + 1));
            }
            else {
                for (size_t j = 0; j < half_size / 2; j++)
                    H5_SWAP_BYTES(d, j, half_size - (j + 1));
            }
        }
        else if (H5T_ORDER_VAX == dst_atomic.order && reverse) {
            if (H5T_FLOAT == dst_p->shared->type) {
                uint8_t tmp1, tmp2;
                size_t  tsize = dst_p->shared->size / 2;
                assert(0 == tsize % 2);

                for (size_t i = 0; i < tsize; i += 4) {
                    tmp1 = d[i];
                    tmp2 = d[i + 1];

                    d[i]     = d[(tsize - 2) - i];
                    d[i + 1] = d[(tsize - 1) - i];

                    d[(tsize - 2) - i] = tmp1;
                    d[(tsize - 1) - i] = tmp2;
                }
            }
            else
                HGOTO_ERROR(H5E_DATATYPE, H5E_UNSUPPORTED, FAIL,
                            "VAX byte ordering is unsupported for complex number type conversions");
        }

next:
        /*
         * If we had used a temporary buffer for the destination then we
         * should copy the value to the true destination buffer.
         */
        if (d == dbuf) {
            if (H5T_FLOAT == dst_p->shared->type)
                H5MM_memcpy(dp, d, dst_p->shared->size);
            else
                H5MM_memcpy(dp, d, dst_p->shared->size / 2);
        }

        /* Ensure imaginary part of complex number is zeroed */
        if (H5T_COMPLEX == dst_p->shared->type)
            memset(dp + (dst_p->shared->size / 2), 0, dst_p->shared->size / 2);

        /* Advance source & destination pointers by delta amounts */
        sp += src_delta;
        dp += dst_delta;
    } /* end conversion loop */

done:
    H5MM_free(src_rev);

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5T__conv_complex_f_matched() */

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
    (void)src;
    (void)dst;
    (void)cdata;
    (void)conv_ctx;
    (void)nelmts;
    (void)buf_stride;
    (void)bkg_stride;
    (void)_buf;
    (void)bkg;
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
                         size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                         void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(FLOAT_COMPLEX, SCHAR, float _Complex, signed char, SCHAR_MIN, SCHAR_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                         size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                         void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(FLOAT_COMPLEX, UCHAR, float _Complex, unsigned char, 0, UCHAR_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                         size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                         void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(FLOAT_COMPLEX, SHORT, float _Complex, short, SHRT_MIN, SHRT_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                          const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                          size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(FLOAT_COMPLEX, USHORT, float _Complex, unsigned short, 0, USHRT_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                       size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                       void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(FLOAT_COMPLEX, INT, float _Complex, int, INT_MIN, INT_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                        size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                        void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(FLOAT_COMPLEX, UINT, float _Complex, unsigned int, 0, UINT_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                        size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                        void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(FLOAT_COMPLEX, LONG, float _Complex, long, LONG_MIN, LONG_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                         size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                         void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(FLOAT_COMPLEX, ULONG, float _Complex, unsigned long, 0, ULONG_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                         size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                         void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(FLOAT_COMPLEX, LLONG, float _Complex, long long, LLONG_MIN, LLONG_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                          const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                          size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(FLOAT_COMPLEX, ULLONG, float _Complex, unsigned long long, 0, ULLONG_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                            size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    /* Suppress warning about non-standard floating-point literal suffix */
    H5_GCC_CLANG_DIAG_OFF("pedantic")
    H5T_CONV_Zf(FLOAT_COMPLEX, FLOAT16, float _Complex, H5__Float16, -FLT16_MAX, FLT16_MAX);
    H5_GCC_CLANG_DIAG_ON("pedantic")
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
                         size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                         void H5_ATTR_UNUSED *bkg)
{
    H5T_CONV_zf(FLOAT_COMPLEX, FLOAT, float _Complex, float, -, -);
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
                          const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                          size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    H5T_CONV_zF(FLOAT_COMPLEX, DOUBLE, float _Complex, double, -, -);
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
                           size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    H5T_CONV_zF(FLOAT_COMPLEX, LDOUBLE, float _Complex, long double, -, -);
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_fcomplex_dcomplex
 *
 * Purpose:     Converts `float _Complex' to `double _Complex'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_fcomplex_dcomplex(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata,
                            const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                            size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    H5T_CONV_zZ(FLOAT_COMPLEX, DOUBLE_COMPLEX, float _Complex, double _Complex, -, -);
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_fcomplex_lcomplex
 *
 * Purpose:     Converts `float _Complex' to `long double _Complex'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_fcomplex_lcomplex(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata,
                            const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                            size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    H5T_CONV_zZ(FLOAT_COMPLEX, LDOUBLE_COMPLEX, float _Complex, long double _Complex, -, -);
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
                         size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                         void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(DOUBLE_COMPLEX, SCHAR, double _Complex, signed char, SCHAR_MIN, SCHAR_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                         size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                         void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(DOUBLE_COMPLEX, UCHAR, double _Complex, unsigned char, 0, UCHAR_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                         size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                         void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(DOUBLE_COMPLEX, SHORT, double _Complex, short, SHRT_MIN, SHRT_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                          const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                          size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(DOUBLE_COMPLEX, USHORT, double _Complex, unsigned short, 0, USHRT_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                       size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                       void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(DOUBLE_COMPLEX, INT, double _Complex, int, INT_MIN, INT_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                        size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                        void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(DOUBLE_COMPLEX, UINT, double _Complex, unsigned int, 0, UINT_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                        size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                        void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(DOUBLE_COMPLEX, LONG, double _Complex, long, LONG_MIN, LONG_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                         size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                         void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(DOUBLE_COMPLEX, ULONG, double _Complex, unsigned long, 0, ULONG_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                         size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                         void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(DOUBLE_COMPLEX, LLONG, double _Complex, long long, LLONG_MIN, LLONG_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                          const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                          size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(DOUBLE_COMPLEX, ULLONG, double _Complex, unsigned long long, 0, ULLONG_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                            size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    /* Suppress warning about non-standard floating-point literal suffix */
    H5_GCC_CLANG_DIAG_OFF("pedantic")
    H5T_CONV_Zf(DOUBLE_COMPLEX, FLOAT16, double _Complex, H5__Float16, -FLT16_MAX, FLT16_MAX);
    H5_GCC_CLANG_DIAG_ON("pedantic")
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
                         size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                         void H5_ATTR_UNUSED *bkg)
{
    H5T_CONV_Zf(DOUBLE_COMPLEX, FLOAT, double _Complex, float, -FLT_MAX, FLT_MAX);
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
                          const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                          size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    H5T_CONV_zf(DOUBLE_COMPLEX, DOUBLE, double _Complex, double, -, -);
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
                           size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    H5T_CONV_zF(DOUBLE_COMPLEX, LDOUBLE, double _Complex, long double, -, -);
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_dcomplex_fcomplex
 *
 * Purpose:     Converts `double _Complex' to `float _Complex'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_dcomplex_fcomplex(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata,
                            const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                            size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    H5T_CONV_Zz(DOUBLE_COMPLEX, FLOAT_COMPLEX, double _Complex, float _Complex, -FLT_MAX, FLT_MAX);
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_dcomplex_lcomplex
 *
 * Purpose:     Converts `double _Complex' to `long double _Complex'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_dcomplex_lcomplex(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata,
                            const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                            size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    H5T_CONV_zZ(DOUBLE_COMPLEX, LDOUBLE_COMPLEX, double _Complex, long double _Complex, -, -);
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
                         size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                         void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(LDOUBLE_COMPLEX, SCHAR, long double _Complex, signed char, SCHAR_MIN, SCHAR_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                         size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                         void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(LDOUBLE_COMPLEX, UCHAR, long double _Complex, unsigned char, 0, UCHAR_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                         size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                         void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(LDOUBLE_COMPLEX, SHORT, long double _Complex, short, SHRT_MIN, SHRT_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                          const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                          size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(LDOUBLE_COMPLEX, USHORT, long double _Complex, unsigned short, 0, USHRT_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                       size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                       void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(LDOUBLE_COMPLEX, INT, long double _Complex, int, INT_MIN, INT_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                        size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                        void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(LDOUBLE_COMPLEX, UINT, long double _Complex, unsigned int, 0, UINT_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                        size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                        void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(LDOUBLE_COMPLEX, LONG, long double _Complex, long, LONG_MIN, LONG_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
                         size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                         void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(LDOUBLE_COMPLEX, ULONG, long double _Complex, unsigned long, 0, ULONG_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
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
#ifdef H5T_CONV_INTERNAL_LDOUBLE_LLONG
herr_t
H5T__conv_lcomplex_llong(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata, const H5T_conv_ctx_t *conv_ctx,
                         size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                         void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(LDOUBLE_COMPLEX, LLONG, long double _Complex, long long, LLONG_MIN, LLONG_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
}
#endif /* H5T_CONV_INTERNAL_LDOUBLE_LLONG */

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_lcomplex_ullong
 *
 * Purpose:     Converts `long double _Complex' to `unsigned long long'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
#ifdef H5T_CONV_INTERNAL_LDOUBLE_ULLONG
herr_t
H5T__conv_lcomplex_ullong(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata,
                          const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                          size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    H5_GCC_CLANG_DIAG_OFF("float-equal")
    H5T_CONV_Zx(LDOUBLE_COMPLEX, ULLONG, long double _Complex, unsigned long long, 0, ULLONG_MAX);
    H5_GCC_CLANG_DIAG_ON("float-equal")
}
#endif /* H5T_CONV_INTERNAL_LDOUBLE_ULLONG */

#ifdef H5_HAVE__FLOAT16
#ifdef H5T_CONV_INTERNAL_LDOUBLE_FLOAT16
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
                            size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    /* Suppress warning about non-standard floating-point literal suffix */
    H5_GCC_CLANG_DIAG_OFF("pedantic")
    H5T_CONV_Zf(LDOUBLE_COMPLEX, FLOAT16, long double _Complex, H5__Float16, -FLT16_MAX, FLT16_MAX);
    H5_GCC_CLANG_DIAG_ON("pedantic")
}
#endif
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
                         size_t nelmts, size_t buf_stride, size_t H5_ATTR_UNUSED bkg_stride, void *buf,
                         void H5_ATTR_UNUSED *bkg)
{
    H5T_CONV_Zf(LDOUBLE_COMPLEX, FLOAT, long double _Complex, float, -FLT_MAX, FLT_MAX);
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
                          const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                          size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    H5T_CONV_Zf(LDOUBLE_COMPLEX, DOUBLE, long double _Complex, double, -DBL_MAX, DBL_MAX);
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
                           size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    H5T_CONV_zf(LDOUBLE_COMPLEX, LDOUBLE, long double _Complex, long double, -, -);
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_lcomplex_fcomplex
 *
 * Purpose:     Converts `long double _Complex' to `float _Complex'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_lcomplex_fcomplex(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata,
                            const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                            size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    H5T_CONV_Zz(LDOUBLE_COMPLEX, FLOAT_COMPLEX, long double _Complex, float _Complex, -FLT_MAX, FLT_MAX);
}

/*-------------------------------------------------------------------------
 * Function:    H5T__conv_lcomplex_dcomplex
 *
 * Purpose:     Converts `long double _Complex' to `double _Complex'
 *
 * Return:      Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5T__conv_lcomplex_dcomplex(const H5T_t *st, const H5T_t *dt, H5T_cdata_t *cdata,
                            const H5T_conv_ctx_t *conv_ctx, size_t nelmts, size_t buf_stride,
                            size_t H5_ATTR_UNUSED bkg_stride, void *buf, void H5_ATTR_UNUSED *bkg)
{
    H5T_CONV_Zz(LDOUBLE_COMPLEX, DOUBLE_COMPLEX, long double _Complex, double _Complex, -DBL_MAX, DBL_MAX);
}
#endif /* H5_HAVE_COMPLEX_NUMBERS */
