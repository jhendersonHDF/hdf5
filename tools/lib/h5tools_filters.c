/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the COPYING file, which can be found at the root of the source code       *
 * distribution tree, or in https://www.hdfgroup.org/licenses.               *
 * If you do not have access to either file, you may request a copy from     *
 * help@hdfgroup.org.                                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "H5private.h"
#include "h5tools.h"

#define H5TOOLS_UD_FILTER_ID_PARM       0
#define H5TOOLS_UD_FILTER_FLAG_PARM     1
#define H5TOOLS_UD_FILTER_CD_COUNT_PARM 2

/*-------------------------------------------------------------------------
 * print a warning message
 *-------------------------------------------------------------------------
 */
static void
print_filter_warning(const char *dname, const char *fname)
{
    HDfprintf(stderr, "Warning: dataset <%s> cannot be read, %s filter is not available\n", dname, fname);
}

/*-------------------------------------------------------------------------
 * Function: h5tools_canreadf
 *
 * Purpose:  check if the dataset creation property list has filters that
 *           are not registered in the current configuration
 *               1) the external filters GZIP and SZIP might not be available
 *               2) the internal filters might be turned off
 *
 * Return:
 *           1 can read,
 *           0 cannot,
 *           -1 error
 *-------------------------------------------------------------------------
 */
int
h5tools_canreadf(const char *name, /* object name, serves also as boolean print */
                 hid_t       dcpl_id)    /* dataset creation property list */
{
    int          nfilters;       /* number of filters */
    H5Z_filter_t filtn;          /* filter identification number */
    int          i;              /* index */
    int          udfilter_avail; /* index */
    int          ret_value = 1;

    /* get information about filters */
    if ((nfilters = H5Pget_nfilters(dcpl_id)) < 0)
        H5TOOLS_GOTO_ERROR(FAIL, "H5Pget_nfilters failed");

    /* if we do not have filters, we can read the dataset safely */
    if (!nfilters)
        H5TOOLS_GOTO_DONE(1);

    /* check availability of filters */
    for (i = 0; i < nfilters; i++) {
        if ((filtn = H5Pget_filter2(dcpl_id, (unsigned)i, 0, 0, 0, (size_t)0, 0, NULL)) < 0)
            H5TOOLS_GOTO_ERROR(FAIL, "H5Pget_filter2 failed");

        switch (filtn) {
            /*-------------------------------------------------------------------------
             * user defined filter
             *-------------------------------------------------------------------------
             */
            default:
                if ((udfilter_avail = H5Zfilter_avail(filtn)) < 0) {
                    H5TOOLS_GOTO_ERROR(FAIL, "H5Zfilter_avail failed");
                }
                else if (!udfilter_avail) {
                    if (name)
                        print_filter_warning(name, "user defined");
                    ret_value = 0;
                }
                break;

                /*-------------------------------------------------------------------------
                 * H5Z_FILTER_DEFLATE      1 , deflation like gzip
                 *-------------------------------------------------------------------------
                 */
            case H5Z_FILTER_DEFLATE:
#ifndef H5_HAVE_FILTER_DEFLATE
                if (name)
                    print_filter_warning(name, "deflate");
                ret_value = 0;
#endif
                break;
                /*-------------------------------------------------------------------------
                 * H5Z_FILTER_SZIP       4 , szip compression
                 *-------------------------------------------------------------------------
                 */
            case H5Z_FILTER_SZIP:
#ifndef H5_HAVE_FILTER_SZIP
                if (name)
                    print_filter_warning(name, "SZIP");
                ret_value = 0;
#endif
                break;
                /*-------------------------------------------------------------------------
                 * H5Z_FILTER_SHUFFLE    2 , shuffle the data
                 *-------------------------------------------------------------------------
                 */
            case H5Z_FILTER_SHUFFLE:
                break;
                /*-------------------------------------------------------------------------
                 * H5Z_FILTER_FLETCHER32 3 , fletcher32 checksum of EDC
                 *-------------------------------------------------------------------------
                 */
            case H5Z_FILTER_FLETCHER32:
                break;
                /*-------------------------------------------------------------------------
                 * H5Z_FILTER_NBIT
                 *-------------------------------------------------------------------------
                 */
            case H5Z_FILTER_NBIT:
                break;
                /*-------------------------------------------------------------------------
                 * H5Z_FILTER_SCALEOFFSET
                 *-------------------------------------------------------------------------
                 */
            case H5Z_FILTER_SCALEOFFSET:
                break;
        } /*switch*/
    }     /*for*/

done:
    return ret_value;
}

/*-------------------------------------------------------------------------
 * Function: h5tools_canwritef
 *
 * Purpose:  check if the filter is available and can write data.
 *
 * Return:   1 can write,
 *           0 cannot,
 *           -1 error
 *-------------------------------------------------------------------------
 */
H5_ATTR_CONST int
h5tools_can_encode(H5Z_filter_t filtn)
{
    int ret_value = 1;

    switch (filtn) {
        /* user defined filter     */
        default:
            H5TOOLS_GOTO_DONE(0);
            break;
        case H5Z_FILTER_DEFLATE:
#ifndef H5_HAVE_FILTER_DEFLATE
            H5TOOLS_GOTO_DONE(0);
#endif
            break;

        case H5Z_FILTER_SZIP:
#ifndef H5_HAVE_FILTER_SZIP
            H5TOOLS_GOTO_DONE(0);
#else
        {
            unsigned int filter_config_flags;

            if (H5Zget_filter_info(filtn, &filter_config_flags) < 0)
                H5TOOLS_GOTO_ERROR(FAIL, "H5Zget_filter_info failed");
            if ((filter_config_flags &
                 (H5Z_FILTER_CONFIG_ENCODE_ENABLED | H5Z_FILTER_CONFIG_DECODE_ENABLED)) == 0) {
                /* filter present but neither encode nor decode is supported (???) */
                H5TOOLS_GOTO_ERROR(FAIL, "neither encode nor decode is supported");
            }
            else if ((filter_config_flags &
                      (H5Z_FILTER_CONFIG_ENCODE_ENABLED | H5Z_FILTER_CONFIG_DECODE_ENABLED)) ==
                     H5Z_FILTER_CONFIG_DECODE_ENABLED) {
                /* decoder only: read but not write */
                H5TOOLS_GOTO_DONE(0);
            }
            else if ((filter_config_flags &
                      (H5Z_FILTER_CONFIG_ENCODE_ENABLED | H5Z_FILTER_CONFIG_DECODE_ENABLED)) ==
                     H5Z_FILTER_CONFIG_ENCODE_ENABLED) {
                /* encoder only: write but not read (???) */
                H5TOOLS_GOTO_ERROR(FAIL, "encoder only: write but not read");
            }
            else if ((filter_config_flags &
                      (H5Z_FILTER_CONFIG_ENCODE_ENABLED | H5Z_FILTER_CONFIG_DECODE_ENABLED)) ==
                     (H5Z_FILTER_CONFIG_ENCODE_ENABLED | H5Z_FILTER_CONFIG_DECODE_ENABLED)) {
                H5TOOLS_GOTO_DONE(1);
            }
        }
#endif
            break;

        case H5Z_FILTER_SHUFFLE:
            break;

        case H5Z_FILTER_FLETCHER32:
            break;

        case H5Z_FILTER_NBIT:
            break;

        case H5Z_FILTER_SCALEOFFSET:
            break;
    } /*switch*/

done:
    return ret_value;
}

/*-------------------------------------------------------------------------
 * Function: h5tools_filter_name_to_id
 *
 * Purpose:  Given the name of a filter, returns the H5Z_filter_t ID for
 *           that filter.
 *
 * Return:   Filter ID if name matches a known filter
 *           H5Z_FILTER_RESERVED for a user-defined filter
 *           H5Z_FILTER_ERROR if name doesn't match a known filter
 *
 *-------------------------------------------------------------------------
 */
H5Z_filter_t
h5tools_filter_name_to_id(const char *filter_name)
{
    H5Z_filter_t ret_value = H5Z_FILTER_ERROR;
    size_t       name_len;
    char *       name_copy = NULL;

    HDassert(filter_name);

    name_len = HDstrlen(filter_name);

    /*
     * Make a copy of the filter name and convert it to
     * all lowercase letters to simplify comparisons
     */
    if (NULL == (name_copy = HDstrdup(filter_name)))
        H5TOOLS_GOTO_ERROR(FAIL, "couldn't copy filter name");

    for (size_t i = 0; i < name_len; i++)
        if (HDisalpha((int)name_copy[i]))
            name_copy[i] = (char)tolower((int)name_copy[i]);

    if (!HDstrcmp(name_copy, "gzip")) {
        ret_value = H5Z_FILTER_DEFLATE;
    }
    else if (!HDstrcmp(name_copy, "szip")) {
        ret_value = H5Z_FILTER_SZIP;
    }
    else if (!HDstrcmp(name_copy, "shuffle") || !HDstrcmp(name_copy, "shuf")) {
        ret_value = H5Z_FILTER_SHUFFLE;
    }
    else if (!HDstrcmp(name_copy, "fletcher32") || !HDstrcmp(name_copy, "flet")) {
        ret_value = H5Z_FILTER_FLETCHER32;
    }
    else if (!HDstrcmp(name_copy, "nbit")) {
        ret_value = H5Z_FILTER_NBIT;
    }
    else if (!HDstrcmp(name_copy, "scaleoffset") || !HDstrcmp(name_copy, "soff")) {
        ret_value = H5Z_FILTER_SCALEOFFSET;
    }
    else if (!HDstrcmp(name_copy, "user-defined") || !HDstrcmp(name_copy, "ud")) {
        ret_value = H5Z_FILTER_RESERVED;
    }
    else if (!HDstrcmp(name_copy, "none")) {
        ret_value = H5Z_FILTER_NONE;
    }

done:
    HDfree(name_copy);

    return ret_value;
}

/*-------------------------------------------------------------------------
 * Function: h5tools_filter_id_to_name
 *
 * Purpose:  Given a filter ID, returns the name of the filter
 *
 * Return:   Name of the filter if filter is known
 *           "User-defined (Filter ID)" for user-defined filters
 *           "H5Z_FILTER_ERROR" for H5Z_FILTER_ERROR ID value
 *           "H5Z_FILTER_RESERVED" for H5Z_FILTER_RESERVED ID value
 *           "Unknown" otherwise
 *
 *-------------------------------------------------------------------------
 */
const char *
h5tools_filter_id_to_name(H5Z_filter_t filter_id)
{
    static char ud_filter_buf[64];

    switch (filter_id) {
        case H5Z_FILTER_NONE:
            return "None";
        case H5Z_FILTER_DEFLATE:
            return "GZIP/Deflate";
        case H5Z_FILTER_SHUFFLE:
            return "Shuffle";
        case H5Z_FILTER_FLETCHER32:
            return "Fletcher32";
        case H5Z_FILTER_SZIP:
            return "SZIP";
        case H5Z_FILTER_NBIT:
            return "Nbit";
        case H5Z_FILTER_SCALEOFFSET:
            return "ScaleOffset";
        case H5Z_FILTER_ERROR:
            return "H5Z_FILTER_ERROR";
        case H5Z_FILTER_RESERVED:
            return "H5Z_FILTER_RESERVED";
        default:
            if ((filter_id >= H5Z_FILTER_RESERVED) && (filter_id <= H5Z_FILTER_MAX)) {
                HDsnprintf(ud_filter_buf, sizeof(ud_filter_buf), "User-defined (ID: %lld)", (long long)filter_id);
                return ud_filter_buf;
            }
            else
                return "Unknown";
    }
}

/*-------------------------------------------------------------------------
 * Function: h5tools_parse_filter
 *
 * Purpose:  Parses a filter string containing the name of the filter and
 *           any parameters for the filter. The format for the filter
 *           string is the following:
 *
 *           Filter Name:{Filter Parameters}
 *
 *           where 'Filter Name' is the name of the filter, e.g. GZIP, and
 *           'Filter Parameters' is an optional comma-separated list of
 *           auxiliary data for the filter. For example, "GZIP:{7}" would
 *           signify that the GZIP/Deflate filter should be used with a
 *           compression level of 7. "SZIP:{NN, 32}" would signify that
 *           the SZIP filter should be used with the H5_SZIP_NN_OPTION_MASK
 *           flag for the filter's coding method and a Pixels-Per-block
 *           setting of 32. If the filter takes no extra parameters, only
 *           the filter name should be specified in the filter string and
 *           everything else should be omitted, e.g. "Shuffle" is a valid
 *           filter string, but "Shuffle:" and "Shuffle:{}" should be
 *           considered malformed. User-defined filters should specify
 *           "User-defined" for the filter name.
 *
 *           Note that auxiliary data for filters is passed as unsigned
 *           integer values. This routine may interpret special strings
 *           like 'NN' for filters it knows about, but parameters should
 *           otherwise be parseable as unsigned integer values. Also note
 *           that the order these parameters are specified in is important.
 *           The following is a list of known filters which accept filter
 *           parameters. The parameters are listed in the order they are
 *           accepted.
 *
 *            GZIP/Deflate - GZIP compression
 *              - Compression Level - An integer value between 1 (less
 *                                    compression) - 9 (more compression)
 *
 *            SZIP - SZIP compression
 *              - Coding Method    - 'EC' or 'NN'
 *              - Pixels-per-block - An even number between 2 - 32
 *
 *            SOFF/ScaleOffset - Scale/Offset filter
 *              - Scale Type   - 'IN' or 'DS'
 *              - Scale Factor - An integer value
 *
 *            UD/User-defined - User-defined filter
 *              - Filter ID      - Required - ID of the user-defined filter
 *              - Filter Flags   - Required - 1 to mark the filter as
 *                                            optional
 *                                            0 to mark the filter as
 *                                            mandatory
 *              - CD Value Count - Required - Number of client data values
 *                                            the filter expects and
 *                                            accepts
 *              - CD Values      - Optional - List of client data values
 *                                            for the filter
 *
 *           The parsed filter's info is returned through the `filter_info`
 *           parameter.
 *
 * Return:   Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
h5tools_parse_filter(const char *filter_str, h5tools_filter_info_t *filter_info)
{
    const char * delims = ":{},";
    H5Z_filter_t filter_id;
    unsigned     filter_flag = H5Z_FLAG_MANDATORY;
    unsigned     cd_values[MAX_CD_VALUES] = { 0 };
    hbool_t      user_defined             = FALSE;
    size_t       num_parsed_params        = 0;
    size_t       expected_cd_nvalues      = 0;
    char *       filter_str_copy          = NULL;
    char *       filter_name              = NULL;
    char *       token                    = NULL;
    herr_t       ret_value                = SUCCEED;

    HDassert(filter_str);
    HDassert(filter_info);

    /* Initialize filter info */
    HDmemset(filter_info, 0, sizeof(h5tools_filter_info_t));

    /* Make a copy of the filter string for tokenizing */
    if (NULL == (filter_str_copy = HDstrdup(filter_str)))
        H5TOOLS_GOTO_ERROR(FAIL, "couldn't copy filter string for tokenizing");

    /* Parse filter name and set filter ID */
    if (NULL == (filter_name = HDstrtok(filter_str_copy, delims)))
        H5TOOLS_GOTO_ERROR(FAIL, "couldn't parse filter name in <%s>", filter_str);

    if (H5Z_FILTER_ERROR == (filter_id = h5tools_filter_name_to_id(filter_name)))
        H5TOOLS_GOTO_ERROR(FAIL, "unknown filter name in <%s>", filter_name);

    /* User-defined filters should have the filter ID as part of the filter string */
    user_defined = (H5Z_FILTER_RESERVED == filter_id);

    /*
     * Preemptively set the expected number of client data values based on the filter.
     * For user-defined filters, this will be parsed out of the filter string instead.
     */
    switch (filter_id) {
        case H5Z_FILTER_DEFLATE:
            expected_cd_nvalues = 1;
            break;
        case H5Z_FILTER_SZIP:
            expected_cd_nvalues = H5Z_SZIP_USER_NPARMS;
            break;
        case H5Z_FILTER_SHUFFLE:
            expected_cd_nvalues = H5Z_SHUFFLE_USER_NPARMS;
            break;
        case H5Z_FILTER_FLETCHER32:
            expected_cd_nvalues = 0;
            break;
        case H5Z_FILTER_NBIT:
            expected_cd_nvalues = H5Z_NBIT_USER_NPARMS;
            break;
        case H5Z_FILTER_SCALEOFFSET:
            expected_cd_nvalues = H5Z_SCALEOFFSET_USER_NPARMS;
            break;
        case H5Z_FILTER_NONE:
        default:
            expected_cd_nvalues = 0;
            break;
    }

    /* Parse any filter parameters given */
    while (NULL != (token = HDstrtok(NULL, delims))) {
        unsigned long token_value;

        if (num_parsed_params >= MAX_CD_VALUES)
            H5TOOLS_GOTO_ERROR(FAIL, "number of filter client data values exceeded maximum (%d)", MAX_CD_VALUES);

        /*
         * Do any special processing or validation based on filter
         * type and parameter index
         */
        if ((filter_id == H5Z_FILTER_SZIP) && (num_parsed_params == H5Z_SZIP_PARM_MASK)) {
            /* Parse SZIP coding method string */
            if (!HDstrcmp(token, "NN"))
                token_value = H5_SZIP_NN_OPTION_MASK;
            else if (!HDstrcmp(token, "EC"))
                token_value = H5_SZIP_EC_OPTION_MASK;
            else
                H5TOOLS_GOTO_ERROR(FAIL, "SZIP coding method must be 'NN' or 'EC' in <%s>", filter_str);
        }
        else if ((filter_id == H5Z_FILTER_SCALEOFFSET) && (num_parsed_params == 0)) {
            /* Parse ScaleOffset scale type string */
            if (!HDstrcmp(token, "IN"))
                token_value = H5Z_SO_INT;
            else if (!HDstrcmp(token, "DS"))
                token_value = H5Z_SO_FLOAT_DSCALE;
            else
                H5TOOLS_GOTO_ERROR(FAIL, "ScaleOffset scale type must be 'IN' or 'DS' in <%s>", filter_str);
        }
        else {
            errno = 0;

            /* Parse token as an unsigned integer value */
            token_value = HDstrtoul(token, NULL, 0);

            if (errno)
                H5TOOLS_GOTO_ERROR(FAIL, "couldn't parse client data value from token <%s>", token);
        }

        /* Make sure the value will fit in an unsigned int */
        if (token_value > UINT_MAX)
            H5TOOLS_GOTO_ERROR(FAIL, "client data value <%lu> was larger than an unsigned int for token <%s>",
                    token_value, token);

        /*
         * For known filters, set the value in the client data values array.
         * For user-defined filters, the parsed value might be the filter ID,
         * filter flag, client data value count or a regular client data value,
         * so we have to do some extra processing.
         */
        if (!user_defined)
            cd_values[num_parsed_params] = (unsigned)token_value;
        else {
            if (num_parsed_params == H5TOOLS_UD_FILTER_ID_PARM) {
                if (token_value <= H5Z_FILTER_RESERVED)
                    H5TOOLS_GOTO_ERROR(FAIL, "invalid user-defined filter ID in <%s>", token);
                if (token_value > H5Z_FILTER_MAX)
                    H5TOOLS_GOTO_ERROR(FAIL, "user-defined filter ID too large in <%s>", token);

                filter_id = (H5Z_filter_t)token_value;
            }
            else if (num_parsed_params == H5TOOLS_UD_FILTER_FLAG_PARM) {
                if (token_value != H5Z_FLAG_MANDATORY && token_value != H5Z_FLAG_OPTIONAL)
                    H5TOOLS_GOTO_ERROR(FAIL, "invalid user-defined filter flag in <%s>; must be MANDATORY (%u) or OPTIONAL (%u)",
                            token, H5Z_FLAG_MANDATORY, H5Z_FLAG_OPTIONAL);

                filter_flag = (unsigned)token_value;
            }
            else if (num_parsed_params == H5TOOLS_UD_FILTER_CD_COUNT_PARM) {
                if (token_value > MAX_CD_VALUES)
                    H5TOOLS_GOTO_ERROR(FAIL, "number of user-defined filter client data values (%lu) exceeds maximum (%d)",
                            token_value, MAX_CD_VALUES);

                expected_cd_nvalues = (size_t)token_value;
            }
            else {
                /* Make sure we processed the filter ID, filter flag and filter CD value count */
                HDassert(num_parsed_params >= 3);

                cd_values[num_parsed_params - 3] = (unsigned)token_value;
            }
        }

        num_parsed_params++;
    }

    /*
     * Since user-defined filter strings carry the filter ID, filter flag and
     * client data value count along with any client data values, subtract
     * these from the parsed parameters count before comparing against the
     * expected number of real client data values for the user-defined filter.
     */
    if (user_defined)
        num_parsed_params -= 3;

    /*-------------------------------------------------------------------------
     * validate filter parameters
     *-------------------------------------------------------------------------
     */

    if (num_parsed_params != expected_cd_nvalues)
        H5TOOLS_GOTO_ERROR(FAIL,
                "expected number of filter client data values (%zu) didn't match number of parsed values (%zu) in <%s>",
                expected_cd_nvalues, num_parsed_params, filter_str);

    switch (filter_id) {
        case H5Z_FILTER_DEFLATE:
            if (cd_values[0] > 9)
                H5TOOLS_GOTO_ERROR(FAIL, "invalid compression parameter in <%s>", filter_str);

            break;

        case H5Z_FILTER_SZIP: {
            unsigned pixels_per_block = cd_values[0];

            if ((pixels_per_block % 2) == 1)
                H5TOOLS_GOTO_ERROR(FAIL, "pixels_per_block is not even in <%s>", filter_str);
            if (pixels_per_block > H5_SZIP_MAX_PIXELS_PER_BLOCK)
                H5TOOLS_GOTO_ERROR(FAIL, "pixels_per_block is too large in <%s>", filter_str);

            break;
        }

        default:
            break;
    }

    /*-------------------------------------------------------------------------
     * set filter info
     *-------------------------------------------------------------------------
     */

    filter_info->filter_id   = filter_id;
    filter_info->filter_flag = filter_flag;
    filter_info->cd_nelmts   = num_parsed_params;
    HDmemcpy(filter_info->cd_values, cd_values, sizeof(filter_info->cd_values));

done:
    HDfree(filter_str_copy);

    return ret_value;
}

/*-------------------------------------------------------------------------
 * Function: h5tools_apply_filters
 *
 * Purpose:  Given a Dataset Creation Property List, applies all the
 *           filters in the given filter info list in the order that they
 *           appear.
 *
 * Return:   Non-negative on success/Negative on failure
 *
 *-------------------------------------------------------------------------
 */
herr_t
h5tools_apply_filters(hid_t dcpl_id, h5tools_filter_info_t *filter_info_list, size_t num_filters)
{
    herr_t ret_value = SUCCEED;

    HDassert(dcpl_id >= 0);
    HDassert(filter_info_list);

    for (size_t i = 0; i < num_filters; i++) {
        h5tools_filter_info_t filter_info = filter_info_list[i];

        if (filter_info.filter_id < 0 || filter_info.filter_id > H5Z_FILTER_MAX) {
            if (filter_info.filter_flag == H5Z_FLAG_MANDATORY)
                H5TOOLS_GOTO_ERROR(FAIL, "invalid filter in filter list");
            else
                continue;
        }

        switch (filter_info.filter_id) {
            case H5Z_FILTER_DEFLATE:
                HDassert(filter_info.cd_nelmts == 1);
                if (H5Pset_deflate(dcpl_id, filter_info.cd_values[0]) < 0)
                    H5TOOLS_GOTO_ERROR(FAIL, "H5Pset_deflate failed");
                break;

            case H5Z_FILTER_SHUFFLE:
                HDassert(filter_info.cd_nelmts == H5Z_SHUFFLE_USER_NPARMS);
                if (H5Pset_shuffle(dcpl_id) < 0)
                    H5TOOLS_GOTO_ERROR(FAIL, "H5Pset_shuffle failed");
                break;

            case H5Z_FILTER_FLETCHER32:
                HDassert(filter_info.cd_nelmts == 0);
                if (H5Pset_fletcher32(dcpl_id) < 0)
                    H5TOOLS_GOTO_ERROR(FAIL, "H5Pset_fletcher32 failed");
                break;

            case H5Z_FILTER_SZIP:
                HDassert(filter_info.cd_nelmts == H5Z_SZIP_USER_NPARMS);
                if (H5Pset_szip(dcpl_id, filter_info.cd_values[H5Z_SZIP_PARM_MASK],
                        filter_info.cd_values[H5Z_SZIP_PARM_PPB]) < 0)
                    H5TOOLS_GOTO_ERROR(FAIL, "H5Pset_szip failed");
                break;

            case H5Z_FILTER_NBIT:
                HDassert(filter_info.cd_nelmts == H5Z_NBIT_USER_NPARMS);
                if (H5Pset_nbit(dcpl_id) < 0)
                    H5TOOLS_GOTO_ERROR(FAIL, "H5Pset_nbit failed");
                break;

            case H5Z_FILTER_SCALEOFFSET:
                HDassert(filter_info.cd_nelmts == H5Z_SCALEOFFSET_USER_NPARMS);
                if (filter_info.cd_values[1] > INT_MAX)
                    H5TOOLS_GOTO_ERROR(FAIL, "Scale Offset scale_factor too large");

                if (H5Pset_scaleoffset(dcpl_id, filter_info.cd_values[0],
                        (int)filter_info.cd_values[1]) < 0)
                    H5TOOLS_GOTO_ERROR(FAIL, "H5Pset_scaleoffset failed");
                break;

            default:
                if (H5Pset_filter(dcpl_id, filter_info.filter_id, filter_info.filter_flag,
                        filter_info.cd_nelmts, filter_info.cd_values) < 0)
                    H5TOOLS_GOTO_ERROR(FAIL, "H5Pset_filter failed");
                break;
        }
    }

done:
    return ret_value;
}

/*-------------------------------------------------------------------------
 * Function: h5tools_dump_filter_list
 *
 * Purpose:  Given a list of h5tools_filter_info_t structures, prints out
 *           a comma-separated list of the filters, along with the settings
 *           for each filter.
 *
 *-------------------------------------------------------------------------
 */
void
h5tools_dump_filter_list(FILE *stream, h5tools_filter_info_t *filter_info_list, size_t num_filters)
{
    const char *filter_entry_begin = "{";
    const char *filter_param_begin = "{";
    const char *filter_param_end   = "}";
    const char *filter_entry_end   = "}";

    HDassert(stream);
    HDassert(filter_info_list);

    for (size_t i = 0; i < num_filters; i++) {
        h5tools_filter_info_t filter_info = filter_info_list[i];
        const char *filter_name;

        if (i > 0)
            HDfprintf(stream, ", ");

        HDfprintf(stream, "%s ", filter_entry_begin);

        if (filter_info.filter_id < 0 || filter_info.filter_id > H5Z_FILTER_MAX) {
            HDfprintf(stream, "Invalid Filter }");
            continue;
        }

        filter_name = h5tools_filter_id_to_name(filter_info.filter_id);
        HDassert(filter_name);

        HDfprintf(stream, "%s (%s)", filter_name,
                filter_info.filter_flag == H5Z_FLAG_MANDATORY ? "MANDATORY" : "OPTIONAL");

        if (filter_info.cd_nelmts > 0) {
            HDfprintf(stream, " - %s ", filter_param_begin);

            if (filter_info.filter_id < H5Z_FILTER_RESERVED) {
                /* Library filter */
                switch (filter_info.filter_id) {
                    case H5Z_FILTER_DEFLATE:
                        HDassert(filter_info.cd_nelmts == 1);
                        HDfprintf(stream, "Compression Level: %u", filter_info.cd_values[0]);
                        break;
                    case H5Z_FILTER_SZIP:
                        HDassert(filter_info.cd_nelmts == H5Z_SZIP_USER_NPARMS);
                        HDfprintf(stream, "Coding Method: %s",
                                filter_info.cd_values[0] == H5_SZIP_NN_OPTION_MASK ? "NN" : "EC");
                        HDfprintf(stream, ", Pixels-per-block: %u", filter_info.cd_values[1]);
                        break;
                    case H5Z_FILTER_SCALEOFFSET:
                        HDassert(filter_info.cd_nelmts == H5Z_SCALEOFFSET_USER_NPARMS);
                        HDfprintf(stream, "Scale Type: %s",
                                filter_info.cd_values[0] == H5Z_SO_INT ? "IN" : "DS");
                        HDfprintf(stream, ", Scale Factor: %u", filter_info.cd_values[1]);
                        break;
                    default:
                        break;
                }
            }
            else {
                /* User-defined filter */
                for (size_t j = 0; j < filter_info.cd_nelmts; j++) {
                    if (j > 0)
                        HDfprintf(stream, ", ");

                    HDfprintf(stream, "%u", filter_info.cd_values[j]);
                }
            }

            HDfprintf(stream, " %s", filter_param_end);
        }

        HDfprintf(stream, " %s", filter_entry_end);
    }
}
