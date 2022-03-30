/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* All rights reserved.                                              */
/*                                                                   */
/* Redistribution and use in source and binary forms, with or        */
/* without modification, are permitted provided that the following   */
/* conditions are met:                                               */
/*                                                                   */
/*   1. Redistributions of source code must retain the above         */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer.                                                  */
/*                                                                   */
/*   2. Redistributions in binary form must reproduce the above      */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer in the documentation and/or other materials       */
/*      provided with the distribution.                              */
/*                                                                   */
/*    THIS  SOFTWARE IS PROVIDED  BY THE  UNIVERSITY OF  TEXAS AT    */
/*    AUSTIN  ``AS IS''  AND ANY  EXPRESS OR  IMPLIED WARRANTIES,    */
/*    INCLUDING, BUT  NOT LIMITED  TO, THE IMPLIED  WARRANTIES OF    */
/*    MERCHANTABILITY  AND FITNESS FOR  A PARTICULAR  PURPOSE ARE    */
/*    DISCLAIMED.  IN  NO EVENT SHALL THE UNIVERSITY  OF TEXAS AT    */
/*    AUSTIN OR CONTRIBUTORS BE  LIABLE FOR ANY DIRECT, INDIRECT,    */
/*    INCIDENTAL,  SPECIAL, EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES    */
/*    (INCLUDING, BUT  NOT LIMITED TO,  PROCUREMENT OF SUBSTITUTE    */
/*    GOODS  OR  SERVICES; LOSS  OF  USE,  DATA,  OR PROFITS;  OR    */
/*    BUSINESS INTERRUPTION) HOWEVER CAUSED  AND ON ANY THEORY OF    */
/*    LIABILITY, WHETHER  IN CONTRACT, STRICT  LIABILITY, OR TORT    */
/*    (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY WAY OUT    */
/*    OF  THE  USE OF  THIS  SOFTWARE,  EVEN  IF ADVISED  OF  THE    */
/*    POSSIBILITY OF SUCH DAMAGE.                                    */
/*                                                                   */
/* The views and conclusions contained in the software and           */
/* documentation are those of the authors and should not be          */
/* interpreted as representing official policies, either expressed   */
/* or implied, of The University of Texas at Austin.                 */
/*********************************************************************/

#include <stdio.h>
#include "common.h"
#include <arm_sve.h>

// TODO: write in assembly with proper unrolling of inner loop
int CNAME(BLASLONG m, BLASLONG n, IFLOAT *a, BLASLONG lda, IFLOAT *b){

    BLASLONG j;
    IFLOAT *a_offset, *a_offset1, *a_offset2, *a_offset3, *a_offset4;
    IFLOAT *b_offset, *b_offset1, *b_offset2, *b_offset3, *b_offset4;

    svfloat64_t a_vec, b_vec, c_vec, d_vec;

    svint64_t lda_vec = svindex_s64(0LL, lda);
    uint64_t sve_size = svcntd();

    a_offset = a;
    b_offset = b;

    j = 0;
    svbool_t pg = svwhilelt_b64(j, n);
    uint64_t active = svcntp_b64(svptrue_b64(), pg);
    do {

        a_offset1 = a_offset;
        a_offset2 = a_offset + 1;
        a_offset3 = a_offset + 2;
        a_offset4 = a_offset + 3;

        b_offset1 = b_offset;
        b_offset2 = b_offset1 + active;
        b_offset3 = b_offset2 + active;
        b_offset4 = b_offset3 + active;

        uint64_t i_cnt;

        for(i_cnt=m; i_cnt>3; i_cnt-=4){ //for 4 loop
            a_vec = svld1_gather_index(pg, (double *) a_offset1, lda_vec);
            b_vec = svld1_gather_index(pg, (double *) a_offset2, lda_vec);
            c_vec = svld1_gather_index(pg, (double *) a_offset3, lda_vec);
            d_vec = svld1_gather_index(pg, (double *) a_offset4, lda_vec);
            a_offset1 += 4;
            a_offset2 += 4;
            a_offset3 += 4;
            a_offset4 += 4;
            svst1_f64(pg, (double *) b_offset1, a_vec);
            svst1_f64(pg, (double *) b_offset2, b_vec);
            svst1_f64(pg, (double *) b_offset3, c_vec);
            svst1_f64(pg, (double *) b_offset4, d_vec);
            b_offset1 += active * 4;
            b_offset2 += active * 4;
            b_offset3 += active * 4;
            b_offset4 += active * 4;
            b_offset += active * 4;
        }

        if (i_cnt == 3) { // for left 3
            a_vec = svld1_gather_index(pg, (double *) a_offset1, lda_vec);
            b_vec = svld1_gather_index(pg, (double *) a_offset2, lda_vec);
            c_vec = svld1_gather_index(pg, (double *) a_offset3, lda_vec);
            svst1_f64(pg, (double *) b_offset1, a_vec);
            svst1_f64(pg, (double *) b_offset2, b_vec);
            svst1_f64(pg, (double *) b_offset3, c_vec);
            b_offset += active * 3;
        }

        if (i_cnt == 2) { // for left 2
            a_vec = svld1_gather_index(pg, (double *) a_offset1, lda_vec);
            b_vec = svld1_gather_index(pg, (double *) a_offset2, lda_vec);
            svst1_f64(pg, (double *) b_offset1, a_vec);
            svst1_f64(pg, (double *) b_offset2, b_vec);
            b_offset += active * 2;
        }

        if (i_cnt == 1) { // for left 1
            a_vec = svld1_gather_index(pg, (double *) a_offset1, lda_vec);
            svst1_f64(pg, (double *) b_offset1, a_vec);
            b_offset += active;
        }

        a_offset += sve_size * lda;

        j += svcntd();
        pg = svwhilelt_b64(j, n);
        active = svcntp_b64(svptrue_b64(), pg);


    } while (svptest_any(svptrue_b64(), pg));

    return 0;
}

