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

#include "H5_api_attribute_test_parallel.h"

static void print_attribute_test_header(void);

static void
print_attribute_test_header(void)
{
    if (MAINPROCESS) {
        printf("\n");
        printf("**********************************************\n");
        printf("*                                            *\n");
        printf("*        API Parallel Attribute Tests        *\n");
        printf("*                                            *\n");
        printf("**********************************************\n\n");
    }
}

void
H5_api_attribute_test_parallel_add(void)
{
    /* Add a fake test to print out a header to distinguish different test interfaces */
    AddTest("print_attribute_test_header", print_attribute_test_header, NULL, "Prints header for attribute tests", NULL);

    /* No tests yet */
}
