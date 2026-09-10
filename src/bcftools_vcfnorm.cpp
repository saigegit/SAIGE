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

/*
This is a stripped down version of the original file to split multiallelic variants into biallelic ones.
- Everything is moved into a bcftools namespace
- Functions (and fields of arg_t) not needed by split_multiallelic_to_biallelics have been removed
- split_multiallelic_to_biallelics has been altered to strip a common suffix (GAGA CAGA,G is now changed into G C and GAGA G instead of GACA CAGA and GACA G)
- The constructor and destructor no longer allocate/destroy out_hdr
- Removed static keyword from functions so they get exported
*/

#include "bcftools.hpp"
#include "bcftools_vcfnorm.hpp"
#include <string>


namespace bcftools {

static void split_info_numeric(args_t *args, bcf1_t *src, bcf_info_t *info, int ialt, bcf1_t *dst)
{
    #define BRANCH_NUMERIC(type,type_t,is_vector_end,is_missing) \
    { \
        const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,info->key); \
        int ntmp = args->ntmp_arr1 / sizeof(type_t); \
        int ret = bcf_get_info_##type(args->hdr,src,tag,&args->tmp_arr1,&ntmp); \
        args->ntmp_arr1 = ntmp * sizeof(type_t); \
        assert( ret>0 ); \
        type_t *vals = (type_t*) args->tmp_arr1; \
        int len = bcf_hdr_id2length(args->hdr,BCF_HL_INFO,info->key); \
        if ( len==BCF_VL_A ) \
        { \
            if ( ret!=src->n_allele-1 ) \
            { \
                if ( args->force && !args->force_warned ) \
                { \
                    fprintf(stderr, \
                        "Warning: wrong number of fields in INFO/%s at %s:%"PRId64", expected %d, found %d\n" \
                        "         (This warning is printed only once.)\n", \
                        tag,bcf_seqname(args->hdr,src),(int64_t) src->pos+1,src->n_allele-1,ret); \
                    args->force_warned = 1; \
                } \
                if ( args->force ) \
                { \
                    bcf_update_info_##type(args->out_hdr,dst,tag,NULL,0); \
                    return; \
                } \
                error("Error: wrong number of fields in INFO/%s at %s:%"PRId64", expected %d, found %d. Use --force to proceed anyway.\n", \
                        tag,bcf_seqname(args->hdr,src),(int64_t) src->pos+1,src->n_allele-1,ret); \
            } \
            bcf_update_info_##type(args->out_hdr,dst,tag,vals+ialt,1); \
        } \
        else if ( len==BCF_VL_R ) \
        { \
            if ( ret!=src->n_allele ) \
            { \
                if ( args->force && !args->force_warned ) \
                { \
                    fprintf(stderr, \
                        "Warning: wrong number of fields in INFO/%s at %s:%"PRId64", expected %d, found %d\n" \
                        "         (This warning is printed only once.)\n", \
                        tag,bcf_seqname(args->hdr,src),(int64_t) src->pos+1,src->n_allele,ret); \
                    args->force_warned = 1; \
                } \
                if ( args->force ) \
                { \
                    bcf_update_info_##type(args->out_hdr,dst,tag,NULL,0); \
                    return; \
                } \
                error("Error: wrong number of fields in INFO/%s at %s:%"PRId64", expected %d, found %d. Use --force to proceed anyway.\n", \
                        tag,bcf_seqname(args->hdr,src),(int64_t) src->pos+1,src->n_allele,ret); \
            } \
            if ( args->keep_sum_ad >= 0 && args->keep_sum_ad==info->key ) \
            { \
                int j; \
                for (j=1; j<info->len; j++) \
                    if ( j!=ialt+1 && !(is_missing) && !(is_vector_end) ) vals[0] += vals[j]; \
                vals[1] = vals[ialt+1]; \
            } \
            else \
            { \
                if ( ialt!=0 ) vals[1] = vals[ialt+1]; \
            } \
            bcf_update_info_##type(args->out_hdr,dst,tag,vals,2); \
        } \
        else if ( len==BCF_VL_G ) \
        { \
            if ( ret!=src->n_allele*(src->n_allele+1)/2 ) \
            { \
                if ( args->force && !args->force_warned ) \
                { \
                    fprintf(stderr, \
                        "Warning: wrong number of fields in INFO/%s at %s:%"PRId64", expected %d, found %d\n" \
                        "         (This warning is printed only once.)\n", \
                        tag,bcf_seqname(args->hdr,src),(int64_t) src->pos+1,src->n_allele*(src->n_allele+1)/2,ret); \
                    args->force_warned = 1; \
                } \
                if ( args->force ) \
                { \
                    bcf_update_info_##type(args->out_hdr,dst,tag,NULL,0); \
                    return; \
                } \
                error("Error: wrong number of fields in INFO/%s at %s:%"PRId64", expected %d, found %d. Use --force to proceed anyway.\n", \
                        tag,bcf_seqname(args->hdr,src),(int64_t) src->pos+1,src->n_allele*(src->n_allele+1)/2,ret); \
            } \
            if ( ialt!=0 ) \
            { \
                vals[1] = vals[bcf_alleles2gt(0,ialt+1)]; \
                vals[2] = vals[bcf_alleles2gt(ialt+1,ialt+1)]; \
            } \
            bcf_update_info_##type(args->out_hdr,dst,tag,vals,3); \
        } \
        else \
            bcf_update_info_##type(args->out_hdr,dst,tag,vals,ret); \
    }
    switch (bcf_hdr_id2type(args->hdr,BCF_HL_INFO,info->key))
    {
        case BCF_HT_INT:  BRANCH_NUMERIC(int32, int32_t, vals[j]==bcf_int32_vector_end, vals[j]==bcf_int32_missing); break;
        case BCF_HT_REAL: BRANCH_NUMERIC(float, float, bcf_float_is_vector_end(vals[j]), bcf_float_is_missing(vals[j])); break;
    }
    #undef BRANCH_NUMERIC
}
// Find nth field in a comma-separated list in src and move it to dst.
// The dst and src memory areas may overlap, end points just after the last valid src
// character.
// On success returns pointer to the end of the parsed field and increments ndst by the
// number of memmoved characters. Returns NULL if the field was not found.
static inline char *string_move_nth(char *dst, char *src, char *end, int nth, size_t *ndst)
{
    if ( src>=end ) return NULL;
    char *ss = src, *se = src;
    int j = 0;
    while ( *se && se<(end) )
    {
        if ( *se==',' )
        {
            if ( j==nth ) break;
            j++;
            ss = se+1;
        }
        se++;
    }
    if ( j!=nth ) return NULL;
    if ( ss>=end ) return NULL;
    if ( !*ss ) return NULL;

    int n = se - ss;
    memmove((dst),ss,n);
    *ndst += n;
    return se;
}

static void split_info_string(args_t *args, bcf1_t *src, bcf_info_t *info, int ialt, bcf1_t *dst)
{
    const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,info->key);
    int ret = bcf_get_info_string(args->hdr,src,tag,&args->tmp_arr1,&args->ntmp_arr1);
    assert( ret>0 );

    kstring_t str;
    str.m = args->ntmp_arr1;
    str.l = ret;
    str.s = (char*) args->tmp_arr1;

    int len = bcf_hdr_id2length(args->hdr,BCF_HL_INFO,info->key);
    if ( len==BCF_VL_A )
    {
        char *end = str.s + str.l;
        char *tmp = str.s;
        str.l = 0;
        tmp = string_move_nth(str.s,tmp,end,ialt,&str.l);
        if ( !tmp ) str.l = 1, str.s[0] = '.';
        kputc_(0,&str);
        bcf_update_info_string(args->out_hdr,dst,tag,str.s);
    }
    else if ( len==BCF_VL_R )
    {
        char *end = str.s + str.l;
        char *tmp = str.s;
        str.l = 0;
        tmp = string_move_nth(str.s,tmp,end,0,&str.l);
        if ( tmp )
        {
            kputc_(',',&str);
            tmp = string_move_nth(str.s+str.l,tmp+1,end,ialt,&str.l); // ialt is 0-based index to ALT
        }
        if ( !tmp ) str.l = 1, str.s[0] = '.';
        kputc_(0,&str);
        bcf_update_info_string(args->out_hdr,dst,tag,str.s);
    }
    else if ( len==BCF_VL_G )
    {
        int i0a = bcf_alleles2gt(0,ialt+1), iaa = bcf_alleles2gt(ialt+1,ialt+1);
        char *end = str.s + str.l;
        char *tmp = str.s;
        str.l = 0;
        tmp = string_move_nth(str.s,tmp,end,0,&str.l);
        if ( tmp )
        {
            kputc_(',',&str);
            tmp = string_move_nth(str.s+str.l,tmp+1,end,i0a-1,&str.l);
        }
        if ( tmp )
        {
            kputc_(',',&str);
            tmp = string_move_nth(str.s+str.l,tmp+1,end,iaa-i0a-1,&str.l);
        }
        if ( !tmp ) str.l = 1, str.s[0] = '.';
        kputc_(0,&str);
        bcf_update_info_string(args->out_hdr,dst,tag,str.s);
    }
    else
        bcf_update_info_string(args->out_hdr,dst,tag,str.s);
    if ( args->ntmp_arr1 < str.m )
    {
        args->ntmp_arr1 = str.m;
        args->tmp_arr1 = (uint8_t*)str.s;
    }
}
static void split_info_flag(args_t *args, bcf1_t *src, bcf_info_t *info, int ialt, bcf1_t *dst)
{
    const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,info->key);
    int ret = bcf_get_info_flag(args->hdr,src,tag,&args->tmp_arr1,&args->ntmp_arr1);
    bcf_update_info_flag(args->out_hdr,dst,tag,NULL,ret);
}

static void split_format_genotype(args_t *args, bcf1_t *src, bcf_fmt_t *fmt, int ialt, bcf1_t *dst)
{
    int ntmp = args->ntmp_arr1 / 4;
    int ngts = bcf_get_genotypes(args->hdr,src,&args->tmp_arr1,&ntmp);
    args->ntmp_arr1 = ntmp * 4;
    assert( ngts >0 );

    int32_t *gt = (int32_t*) args->tmp_arr1;
    int i, j, nsmpl = bcf_hdr_nsamples(args->hdr);
    ngts /= nsmpl;
    for (i=0; i<nsmpl; i++)
    {
        for (j=0; j<ngts; j++)
        {
            if ( gt[j]==bcf_int32_vector_end ) break;
            if ( bcf_gt_is_missing(gt[j]) ) continue; // missing allele: leave as is
            if ( bcf_gt_allele(gt[j])==0 ) continue; // ref && `--multi-overlaps 0`: leave as is
            if ( bcf_gt_allele(gt[j])==ialt+1 )
                gt[j] = bcf_gt_unphased(1) | bcf_gt_is_phased(gt[j]); // set to first ALT
            else if ( args->ma_use_ref_allele )
                gt[j] = bcf_gt_unphased(0) | bcf_gt_is_phased(gt[j]); // set to REF
            else
                gt[j] = bcf_gt_missing | bcf_gt_is_phased(gt[j]);     // set to missing
        }
        gt += ngts;
    }
    bcf_update_genotypes(args->out_hdr,dst,args->tmp_arr1,ngts*nsmpl);
}
static void split_format_numeric(args_t *args, bcf1_t *src, bcf_fmt_t *fmt, int ialt, bcf1_t *dst)
{
    #define BRANCH_NUMERIC(type,type_t,is_vector_end,is_missing,set_vector_end,set_missing) \
    { \
        const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,fmt->id); \
        int ntmp = args->ntmp_arr1 / sizeof(type_t); \
        int nvals = bcf_get_format_##type(args->hdr,src,tag,&args->tmp_arr1,&ntmp); \
        args->ntmp_arr1 = ntmp * sizeof(type_t); \
        assert( nvals>0 ); \
        type_t *vals = (type_t *) args->tmp_arr1; \
        int len = bcf_hdr_id2length(args->hdr,BCF_HL_FMT,fmt->id); \
        int i,j, nsmpl = bcf_hdr_nsamples(args->hdr); \
        if ( nvals==nsmpl ) /* all values are missing */ \
        { \
            bcf_update_format_##type(args->out_hdr,dst,tag,vals,nsmpl); \
            return; \
        } \
        if ( len==BCF_VL_A ) \
        { \
            if ( nvals!=(src->n_allele-1)*nsmpl ) \
            { \
                if ( args->force && !args->force_warned ) \
                { \
                    fprintf(stderr, \
                        "Warning: wrong number of fields in FMT/%s at %s:%"PRId64", expected %d, found %d. Removing the field.\n" \
                        "         (This warning is printed only once.)\n", \
                        tag,bcf_seqname(args->hdr,src),(int64_t) src->pos+1,(src->n_allele-1)*nsmpl,nvals); \
                    args->force_warned = 1; \
                } \
                if ( args->force ) \
                { \
                    bcf_update_format_##type(args->out_hdr,dst,tag,NULL,0); \
                    return; \
                } \
                error("Error: wrong number of fields in FMT/%s at %s:%"PRId64", expected %d, found %d. Use --force to proceed anyway.\n", \
                    tag,bcf_seqname(args->hdr,src),(int64_t) src->pos+1,(src->n_allele-1)*nsmpl,nvals); \
            } \
            nvals /= nsmpl; \
            type_t *src_vals = vals, *dst_vals = vals; \
            for (i=0; i<nsmpl; i++) \
            { \
                int idst = 0; \
                int isrc = ialt; \
                if ( is_missing || is_vector_end ) set_missing; \
                else dst_vals[idst] = src_vals[isrc]; \
                dst_vals += 1; \
                src_vals += nvals; \
            } \
            bcf_update_format_##type(args->out_hdr,dst,tag,vals,nsmpl); \
        } \
        else if ( len==BCF_VL_R ) \
        { \
            if ( nvals!=src->n_allele*nsmpl ) \
            { \
                if ( args->force && !args->force_warned ) \
                { \
                    fprintf(stderr, \
                        "Warning: wrong number of fields in FMT/%s at %s:%"PRId64", expected %d, found %d. Removing the field.\n" \
                        "         (This warning is printed only once.)\n", \
                        tag,bcf_seqname(args->hdr,src),(int64_t) src->pos+1,(src->n_allele-1)*nsmpl,nvals); \
                    args->force_warned = 1; \
                } \
                if ( args->force ) \
                { \
                    bcf_update_format_##type(args->out_hdr,dst,tag,NULL,0); \
                    return; \
                } \
                error("Error: wrong number of fields in FMT/%s at %s:%"PRId64", expected %d, found %d. Use --force to proceed anyway.\n", \
                    tag,bcf_seqname(args->hdr,src),(int64_t) src->pos+1,src->n_allele*nsmpl,nvals); \
            } \
            nvals /= nsmpl; \
            type_t *src_vals = vals, *dst_vals = vals; \
            if ( args->keep_sum_ad >= 0 && args->keep_sum_ad==fmt->id ) \
            { \
                for (i=0; i<nsmpl; i++) \
                { \
                    dst_vals[0] = src_vals[0]; \
                    for (j=1; j<nvals; j++) \
                    { \
                        int isrc = j; \
                        if ( j!=ialt+1 && !(is_missing) && !(is_vector_end) ) dst_vals[0] += src_vals[j]; \
                    } \
                    int isrc = ialt + 1; \
                    int idst = 1; \
                    if ( is_vector_end ) set_missing; \
                    else dst_vals[idst] = src_vals[isrc]; \
                    dst_vals += 2; \
                    src_vals += nvals; \
                } \
            } \
            else \
            { \
                for (i=0; i<nsmpl; i++) \
                { \
                    dst_vals[0] = src_vals[0]; \
                    int isrc = ialt + 1; \
                    int idst = 1; \
                    if ( is_vector_end ) set_missing; \
                    else dst_vals[idst] = src_vals[isrc]; \
                    dst_vals += 2; \
                    src_vals += nvals; \
                } \
            } \
            bcf_update_format_##type(args->out_hdr,dst,tag,vals,nsmpl*2); \
        } \
        else if ( len==BCF_VL_G ) \
        { \
            if ( nvals!=src->n_allele*(src->n_allele+1)/2*nsmpl && nvals!=src->n_allele*nsmpl ) \
            { \
                if ( args->force && !args->force_warned ) \
                { \
                    fprintf(stderr, \
                        "Warning: wrong number of fields in FMT/%s at %s:%"PRId64", expected %d, found %d. Removing the field.\n" \
                        "         (This warning is printed only once.)\n", \
                        tag,bcf_seqname(args->hdr,src),(int64_t) src->pos+1,(src->n_allele-1)*nsmpl,nvals); \
                    args->force_warned = 1; \
                } \
                if ( args->force ) \
                { \
                    bcf_update_format_##type(args->out_hdr,dst,tag,NULL,0); \
                    return; \
                } \
                error("Error at %s:%"PRId64", the tag %s has wrong number of fields. Use --force to proceed anyway.\n", \
                    bcf_seqname(args->hdr,src),(int64_t) src->pos+1,bcf_hdr_int2id(args->hdr,BCF_DT_ID,fmt->id)); \
            } \
            nvals /= nsmpl; \
            int all_haploid = nvals==src->n_allele ? 1 : 0; \
            type_t *src_vals = vals, *dst_vals = vals; \
            for (i=0; i<nsmpl; i++) \
            { \
                int haploid = all_haploid; \
                if ( !haploid ) \
                { \
                    int j; \
                    for (j=0; j<nvals; j++) \
                    { \
                        int isrc = j; \
                        if ( is_vector_end ) break; \
                    } \
                    if ( j!=nvals ) haploid = 1; \
                } \
                dst_vals[0] = src_vals[0]; \
                if ( haploid ) \
                { \
                    dst_vals[1] = src_vals[ialt+1]; \
                    if ( !all_haploid ) { int idst = 2; set_vector_end; } \
                } \
                else \
                { \
                    dst_vals[1] = src_vals[bcf_alleles2gt(0,ialt+1)]; \
                    dst_vals[2] = src_vals[bcf_alleles2gt(ialt+1,ialt+1)]; \
                } \
                dst_vals += all_haploid ? 2 : 3; \
                src_vals += nvals; \
            } \
            bcf_update_format_##type(args->out_hdr,dst,tag,vals,all_haploid ? nsmpl*2 : nsmpl*3); \
        } \
        else \
            bcf_update_format_##type(args->out_hdr,dst,tag,vals,nvals); \
    }
    switch (bcf_hdr_id2type(args->hdr,BCF_HL_FMT,fmt->id))
    {
        case BCF_HT_INT:  BRANCH_NUMERIC(int32, int32_t, src_vals[isrc]==bcf_int32_vector_end, src_vals[isrc]==bcf_int32_missing, dst_vals[idst]=bcf_int32_vector_end, dst_vals[idst]=bcf_int32_missing); break;
        case BCF_HT_REAL: BRANCH_NUMERIC(float, float, bcf_float_is_vector_end(src_vals[isrc]), bcf_float_is_missing(src_vals[isrc]), bcf_float_set_vector_end(dst_vals[idst]), bcf_float_set_missing(dst_vals[idst])); break;
    }
    #undef BRANCH_NUMERIC
}
static void squeeze_format_char(char *str, int src_blen, int dst_blen, int n)
{
    int i, isrc = 0, idst = 0;
    for (i=0; i<n; i++)
    {
        memmove(str+idst,str+isrc,dst_blen);
        idst += dst_blen;
        isrc += src_blen;
    }
}
// ialt is 0-based index to ALT
static void split_format_string(args_t *args, bcf1_t *src, bcf_fmt_t *fmt, int ialt, bcf1_t *dst)
{
    const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,fmt->id);
    int ret = bcf_get_format_char(args->hdr,src,tag,&args->tmp_arr1,&args->ntmp_arr1);
    if ( !ret ) return; // all values can be empty, leave out the tag, no need to panic
    assert( ret>0 );

    int nsmpl = bcf_hdr_nsamples(args->hdr);
    int blen = ret/nsmpl;   // per-sample field length
    assert( blen>0 );

    kstring_t str;
    str.m = args->ntmp_arr1;
    str.s = (char*) args->tmp_arr1;
    str.l = ret;

    int tag_len = bcf_hdr_id2length(args->hdr,BCF_HL_FMT,fmt->id);
    if ( tag_len==BCF_VL_A )
    {
        int i, maxlen = 0;
        char *ptr = str.s;
        for (i=0; i<nsmpl; i++)
        {
            char *tmp = ptr;
            char *end = ptr + blen;
            size_t len = 0;
            tmp = string_move_nth(ptr,tmp,end,ialt,&len);
            if ( !tmp ) ptr[0] = '.', len = 1;
            if ( maxlen < len ) maxlen = len;
            while (len<blen) ptr[len++] = 0;
            ptr += blen;
        }
        if ( maxlen<blen ) squeeze_format_char(str.s,blen,maxlen,nsmpl);
        bcf_update_format_char(args->out_hdr,dst,tag,str.s,nsmpl*maxlen);
    }
    else if ( tag_len==BCF_VL_R )
    {
        int i, maxlen = 0;
        char *ptr = str.s;
        for (i=0; i<nsmpl; i++)
        {
            char *tmp = ptr;
            char *end = ptr + blen;
            size_t len = 0;
            tmp = string_move_nth(ptr,tmp,end,0,&len);
            if ( tmp )
            {
                ptr[len++] = ',';
                tmp = string_move_nth(ptr+len,tmp+1,end,ialt,&len);
            }
            if ( !tmp ) ptr[0] = '.', len = 1;
            if ( maxlen < len ) maxlen = len;
            while (len<blen) ptr[len++] = 0;
            ptr += blen;
        }
        if ( maxlen<blen ) squeeze_format_char(str.s,blen,maxlen,nsmpl);
        bcf_update_format_char(args->out_hdr,dst,tag,str.s,nsmpl*maxlen);
    }
    else if ( tag_len==BCF_VL_G )
    {
        int i, maxlen = 0, i0a = bcf_alleles2gt(0,ialt+1), iaa = bcf_alleles2gt(ialt+1,ialt+1);
        char *ptr = str.s;
        for (i=0; i<nsmpl; i++)
        {
            char *se = ptr, *sx = ptr+blen;
            int nfields = 1;
            while ( *se && se<sx )
            {
                if ( *se==',' ) nfields++;
                se++;
            }
            if ( nfields==1 && se-ptr==1 && *ptr=='.' ) continue;   // missing value
            if ( nfields!=src->n_allele*(src->n_allele+1)/2 && nfields!=src->n_allele )
            {
                if ( args->force && !args->force_warned )
                {
                    fprintf(stderr,
                            "Warning: wrong number of fields in FMT/%s at %s:%"PRId64", expected %d or %d, found %d. Removing the field.\n"
                            "         (This warning is printed only once.)\n",
                            tag,bcf_seqname(args->hdr,src),(int64_t)src->pos+1,src->n_allele*(src->n_allele+1)/2,src->n_allele,nfields);
                    args->force_warned = 1;
                }
                if ( args->force )
                {
                    bcf_update_format_char(args->out_hdr,dst,tag,NULL,0);
                    return;
                }
                error("Error: wrong number of fields in FMT/%s at %s:%"PRId64", expected %d or %d, found %d. Use --force to proceed anyway.\n",
                        tag,bcf_seqname(args->hdr,src),(int64_t) src->pos+1,src->n_allele*(src->n_allele+1)/2,src->n_allele,nfields);
            }

            char *tmp = ptr;
            char *end = ptr + blen;
            size_t len = 0;
            tmp = string_move_nth(ptr,tmp,end,0,&len);
            if ( nfields==src->n_allele )   // haploid
            {
                if ( tmp )
                {
                    ptr[len++] = ',';
                    tmp = string_move_nth(ptr+len,tmp+1,end,ialt,&len);
                }
            }
            else    // diploid
            {
                if ( tmp )
                {
                    ptr[len++] = ',';
                    tmp = string_move_nth(ptr+len,tmp+1,end,i0a-1,&len);
                }
                if ( tmp )
                {
                    ptr[len++] = ',';
                    tmp = string_move_nth(ptr+len,tmp+1,end,iaa-i0a-1,&len);
                }
            }
            if ( !tmp ) ptr[0] = '.', len = 1;
            if ( maxlen < len ) maxlen = len;
            while (len<blen) ptr[len++] = 0;
            ptr += blen;
        }
        if ( maxlen<blen ) squeeze_format_char(str.s,blen,maxlen,nsmpl);
        bcf_update_format_char(args->out_hdr,dst,tag,str.s,nsmpl*maxlen);
    }
    else
        bcf_update_format_char(args->out_hdr,dst,tag,str.s,str.l);
}

void split_multiallelic_to_biallelics(args_t *args, bcf1_t *line)
{
    int i;

    bcf_unpack(line, BCF_UN_ALL);

    // Init the target biallelic lines
    args->ntmp_lines = line->n_allele-1;
    if ( args->mtmp_lines < args->ntmp_lines )
    {
        args->tmp_lines = (bcf1_t **)realloc(args->tmp_lines,sizeof(bcf1_t*)*args->ntmp_lines);
        for (i=args->mtmp_lines; i<args->ntmp_lines; i++)
            args->tmp_lines[i] = NULL;
        args->mtmp_lines = args->ntmp_lines;
    }
    kstring_t tmp = {0,0,0};
    kputs(line->d.allele[0], &tmp);
    kputc(',', &tmp);
    int rlen  = tmp.l;
    int gt_id = bcf_hdr_id2int(args->hdr,BCF_DT_ID,"GT");
    for (i=0; i<args->ntmp_lines; i++)  // for each ALT allele
    {
        if ( !args->tmp_lines[i] ) args->tmp_lines[i] = bcf_init1();
        bcf1_t *dst = args->tmp_lines[i];
        bcf_clear(dst);

        dst->rid  = line->rid;
        dst->pos  = line->pos;
        dst->qual = line->qual;

        // Not quite sure how to handle IDs, they can be assigned to a specific
        // ALT.  For now we leave the ID unchanged for all.
        bcf_update_id(args->out_hdr, dst, line->d.id ? line->d.id : ".");

        tmp.l = rlen;
        kputs(line->d.allele[i+1],&tmp);
        int n1 = strlen(line->d.allele[0]);
        int n2 = strlen(line->d.allele[i+1]);
        if(n1>1 && n2>1) {
            // remove common suffix, e.g., shorten "GAGA-CAGA to "G-C"
            std::string ref(line->d.allele[0]);
            std::string alt(line->d.allele[i+1]);
            int common_suffix_length = 0;
            while(common_suffix_length < n1 && common_suffix_length < n2 && ref[n1 - common_suffix_length - 1] == alt[n2 - common_suffix_length - 1]) common_suffix_length++;
            std::string alleles_str = ref.substr(0, n1 - common_suffix_length) + "," + alt.substr(0, n2 - common_suffix_length);
            bcf_update_alleles_str(args->out_hdr, dst, alleles_str.c_str());
        } else {
            bcf_update_alleles_str(args->out_hdr,dst,tmp.s);
        }

        if ( line->d.n_flt ) bcf_update_filter(args->hdr, dst, line->d.flt, line->d.n_flt);

        int j;
        for (j=0; j<line->n_info; j++)
        {
            bcf_info_t *info = &line->d.info[j];
            int type = bcf_hdr_id2type(args->hdr,BCF_HL_INFO,info->key);
            if ( type==BCF_HT_INT || type==BCF_HT_REAL ) split_info_numeric(args, line, info, i, dst);
            else if ( type==BCF_HT_FLAG ) split_info_flag(args, line, info, i, dst);
            else split_info_string(args, line, info, i, dst);
        }

        dst->n_sample = line->n_sample;
        for (j=0; j<line->n_fmt; j++)
        {
            bcf_fmt_t *fmt = &line->d.fmt[j];
            int type = bcf_hdr_id2type(args->hdr,BCF_HL_FMT,fmt->id);
            if ( fmt->id==gt_id ) split_format_genotype(args, line, fmt, i, dst);
            else if ( type==BCF_HT_INT || type==BCF_HT_REAL ) split_format_numeric(args, line, fmt, i, dst);
            else split_format_string(args, line, fmt, i, dst);
        }
    }
    free(tmp.s);
}

void init_data(args_t *args)
{
    //args->hdr = args->files->readers[0].header;
    if ( args->keep_sum_ad )
    {
        args->keep_sum_ad = bcf_hdr_id2int(args->hdr,BCF_DT_ID,"AD");
        if ( args->keep_sum_ad < 0 ) error("Error: --keep-sum-ad requested but the tag AD is not present\n");
    }
    else
        args->keep_sum_ad = -1;
}

void destroy_data(args_t *args)
{
    int i;
    for (i=0; i<args->mtmp_lines; i++)
        if ( args->tmp_lines[i] ) bcf_destroy1(args->tmp_lines[i]);
    free(args->tmp_lines);
    free(args->als);
    free(args->tmp_arr1);
}


} // namespace