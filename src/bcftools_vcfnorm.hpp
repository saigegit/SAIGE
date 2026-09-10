/*  vcfnorm.c -- Left-align and normalize indels.

    Copyright (C) 2013-2026 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.  */

#include <htslib/vcf.h>

namespace bcftools {

typedef struct
{
    bcf1_t *rec;
    int pass;
}
line_t;

typedef struct
{
    bcf1_t  **tmp_lines;
    line_t *lines;
    int ntmp_lines, mtmp_lines;
    char **als;
    uint8_t *tmp_arr1;
    int ntmp_arr1;
    bcf_hdr_t *hdr, *out_hdr;
    int force, force_warned, keep_sum_ad;
    int ma_use_ref_allele;
}
args_t;

void split_multiallelic_to_biallelics(args_t *args, bcf1_t *line);
void init_data(args_t *args);
void destroy_data(args_t *args);

} // namespace