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
    size_t         buf_size;
    size_t         buf_bytes_left;
} H5_bcursor_t;

/* TODO: docs */
H5_DLL herr_t H5_bcursor_init(const uint8_t *buf, size_t buf_size, H5_bcursor_t *byte_cursor);

/* TODO: docs */
H5_DLL const uint8_t *H5_bcursor_get_buf_ptr(H5_bcursor_t *byte_cursor);

/* TODO: docs */
H5_DLL size_t H5_bcursor_get_num_bytes(H5_bcursor_t *byte_cursor);

/* TODO: docs */
H5_DLL size_t H5_bcursor_get_position(H5_bcursor_t *byte_cursor);

/* TODO: docs */
H5_DLL bool H5_bcursor_has_num_bytes(H5_bcursor_t *byte_cursor, size_t num_bytes);

/* TODO: docs */
/* XXX: Beware endianness */
H5_DLL herr_t H5_bcursor_read(H5_bcursor_t *byte_cursor, size_t num_bytes, void *buf_out);

/* TODO: docs */
H5_DLL herr_t H5_bcursor_memcmp(H5_bcursor_t *byte_cursor, const void *cmp_buf, size_t cmp_size, int *cmp_ret);

/* TODO: docs */
H5_DLL size_t H5_bcursor_strnlen(H5_bcursor_t *byte_cursor, size_t maxlen);

/* TODO: docs */
H5_DLL herr_t H5_bcursor_strncpy(H5_bcursor_t *byte_cursor, char *dest, size_t size);

/* TODO: H5_bcursor_memcmp_and_skip? */

/* TODO: docs */
H5_DLL herr_t H5_bcursor_skip(H5_bcursor_t *byte_cursor, size_t num_bytes);

/* TODO: docs */
H5_DLL herr_t H5_bcursor_decode_u16(H5_bcursor_t *byte_cursor, uint16_t *u16_out);

/* TODO: docs */
H5_DLL herr_t H5_bcursor_decode_u32(H5_bcursor_t *byte_cursor, uint32_t *u32_out);

/* TODO: docs */
H5_DLL herr_t H5_bcursor_decode_u64(H5_bcursor_t *byte_cursor, uint64_t *u64_out);

/* TODO: docs */
H5_DLL herr_t H5_bcursor_decode_var_u64(H5_bcursor_t *byte_cursor, size_t size, uint64_t *u64_out);

/* TODO: docs */
H5_DLL herr_t H5_bcursor_decode_i64(H5_bcursor_t *byte_cursor, int64_t *i64_out);

/* TODO: docs */
H5_DLL herr_t H5_bcursor_decode_file_addr(H5_bcursor_t *byte_cursor, const H5F_t *f, haddr_t *file_addr_out);

/* TODO: docs */
H5_DLL herr_t H5_bcursor_decode_file_addr_len(H5_bcursor_t *byte_cursor, size_t sizeof_addr, haddr_t *file_addr_out);

/* TODO: docs */
H5_DLL herr_t H5_bcursor_decode_file_hsize(H5_bcursor_t *byte_cursor, const H5F_t *f, hsize_t *file_hsize_out);

/* TODO: docs */
H5_DLL herr_t H5_bcursor_decode_len(H5_bcursor_t *byte_cursor, size_t sizeof_len, size_t *len_out);


#endif /* H5bcursor_H */
