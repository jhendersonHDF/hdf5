/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the LICENSE file, which can be found at the root of the source code       *
 * distribution tree, or in https://www.hdfgroup.org/licenses.               *
 * If you do not have access to either file, you may request a copy from     *
 * help@hdfgroup.org.                                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef H5bcursor_H
#define H5bcursor_H

#include "H5private.h"

/* TODO: docs */
typedef struct H5_bcursor_t {
    const uint8_t *buf_ptr;
    const uint8_t *buf_end;
    size_t         buf_bytes_left;
} H5_bcursor_t;

/* TODO: docs */
static inline herr_t
H5_bcursor_init(const uint8_t *buf, size_t buf_size, H5_bcursor_t *byte_cursor) {
    assert(buf);
    assert(byte_cursor);

    byte_cursor->buf_ptr        = buf;
    byte_cursor->buf_bytes_left = buf_size;

    /* TODO: account for SIZE_MAX */
    byte_cursor->buf_end = buf;
    if (buf_size > 0)
        byte_cursor->buf_end = buf + buf_size - 1;

    return SUCCEED;
}

/* TODO: docs */
static inline bool
H5_bcursor_has_num_bytes(H5_bcursor_t *byte_cursor, size_t num_bytes)
{
    return (num_bytes <= byte_cursor->buf_bytes_left);
}

/* TODO: docs */
/* XXX: Beware endianness */
static inline herr_t
H5_bcursor_read(H5_bcursor_t *byte_cursor, void *buf_out, size_t size)
{
    bool enough_space;

    assert(byte_cursor && byte_cursor->buf_ptr);
    assert(buf_out);

    enough_space = (size <= byte_cursor->buf_bytes_left);
    if (!enough_space)
        return FAIL;

    memcpy(buf_out, byte_cursor->buf_ptr, size);
    byte_cursor->buf_ptr        += size; /* TODO: handle pointer wraparound with large size */
    byte_cursor->buf_bytes_left -= size;

    return SUCCEED;
}

/* TODO: docs */
static inline herr_t
H5_bcursor_decode_file_addr(H5_bcursor_t *byte_cursor, const H5F_t *f, haddr_t *file_addr_out)
{
    size_t file_addr_size = (size_t)H5F_sizeof_addr(f);
    bool   enough_space;

    assert(byte_cursor && byte_cursor->buf_ptr);
    assert(f);
    assert(file_addr_out);

    enough_space = (file_addr_size <= byte_cursor->buf_bytes_left);
    if (!enough_space)
        return FAIL;

    H5F_addr_decode(f, &byte_cursor->buf_ptr, file_addr_out);
    byte_cursor->buf_bytes_left -= file_addr_size;

    return SUCCEED;
}

/* TODO: docs */
static inline herr_t
H5_bcursor_decode_file_hsize(H5_bcursor_t *byte_cursor, const H5F_t *f, hsize_t *file_hsize_out)
{
    size_t file_hsize_size = (size_t)H5F_sizeof_size(f);
    bool   enough_space;

    assert(byte_cursor && byte_cursor->buf_ptr);
    assert(f);
    assert(file_hsize_out);

    enough_space = (file_hsize_size <= byte_cursor->buf_bytes_left);
    if (!enough_space)
        return FAIL;

    H5F_DECODE_LENGTH(f, byte_cursor->buf_ptr, *file_hsize_out);
    byte_cursor->buf_bytes_left -= file_hsize_size;

    return SUCCEED;
}

/* TODO: docs */
static inline herr_t
H5_bcursor_decode_u16(H5_bcursor_t *byte_cursor, uint16_t *u16_out)
{
    bool enough_space;

    assert(byte_cursor && byte_cursor->buf_ptr);
    assert(u16_out);

    enough_space = (sizeof(uint16_t) <= byte_cursor->buf_bytes_left);
    if (!enough_space)
        return FAIL;

    UINT16DECODE(byte_cursor->buf_ptr, *u16_out);
    byte_cursor->buf_bytes_left -= sizeof(uint16_t);

    return SUCCEED;
}

/* TODO: docs */
static inline herr_t
H5_bcursor_decode_u32(H5_bcursor_t *byte_cursor, uint32_t *u32_out)
{
    bool enough_space;

    assert(byte_cursor && byte_cursor->buf_ptr);
    assert(u32_out);

    enough_space = (sizeof(uint32_t) <= byte_cursor->buf_bytes_left);
    if (!enough_space)
        return FAIL;

    UINT32DECODE(byte_cursor->buf_ptr, *u32_out);
    byte_cursor->buf_bytes_left -= sizeof(uint32_t);

    return SUCCEED;
}

/* TODO: docs */
static inline herr_t
H5_bcursor_decode_var_u64(H5_bcursor_t *byte_cursor, size_t size, uint64_t *u64_out)
{
    bool enough_space;

    assert(byte_cursor && byte_cursor->buf_ptr);
    assert(u64_out);

    enough_space = (size <= byte_cursor->buf_bytes_left);
    if (!enough_space)
        return FAIL;

    UINT64DECODE_VAR(byte_cursor->buf_ptr, *u64_out, size);
    byte_cursor->buf_bytes_left -= size;

    return SUCCEED;
}

/* TODO: docs */
static inline herr_t
H5_bcursor_skip(H5_bcursor_t *byte_cursor, size_t size)
{
    bool enough_space;

    assert(byte_cursor && byte_cursor->buf_ptr);

    enough_space = (size <= byte_cursor->buf_bytes_left);
    if (!enough_space)
        return FAIL;

    byte_cursor->buf_ptr        += size; /* TODO: handle pointer wraparound with large size */
    byte_cursor->buf_bytes_left -= size;

    return SUCCEED;
}

#endif /* H5bcursor_H */
