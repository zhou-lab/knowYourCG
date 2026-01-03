#ifndef _CDATA_H
#define _CDATA_H

#include <R.h>
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <limits.h>
#include <inttypes.h>
#include "khash.h"
#include "wzmisc.h"
#include "wzmisc.h"
#include "wzbed.h"
#include "bgzf.h"

#define CDSIG 266563789635

typedef struct f2_aux_t {
  uint64_t nk;                  // num keys
  char **keys;                  // pointer to keys, doesn't own memory
  uint8_t *data;                // pointer to data, doesn't own memory
} f2_aux_t;

/** The header design, 17 bytes
    uint64_t: signature, used for validation
    uint8_t: format (0=vec; 1=rle)
    uint64_t: length (n_cs or n_bytes for rle)
**/
typedef struct cdata_t {
  uint8_t *s;
  uint64_t n; /* number of bytes, except for fmt 0, which is sub-byte you need the actual length */
  int compressed;
  char fmt;
  uint8_t unit; // how many bytes is needed for each decompressed data unit, use 0 for format 0,1,6
  void *aux;
} cdata_t;

static inline uint64_t cdata_nbytes(cdata_t *c) {
  uint64_t n = 0;
  switch(c->fmt) {
  case '0': n = ((c->n+7)>>3); break;
  case '6': n = ((c->n+3)>>2); break;
  default: n = c->n;
  }
  return n;
}

static inline void free_cdata(cdata_t *c) {
  if (c->s) free(c->s);
  if (c->fmt == '2' && c->aux) {
    free(((f2_aux_t*) c->aux)->keys);
    free(c->aux);
  }
  if (c->fmt == '7' && c->aux) free(c->aux);
  c->s = NULL;
}

static inline size_t bit_count(cdata_t c) {

  /* create a look-up table */
  int byte2cnt[256]; int p;
  for (p=0; p<256; ++p) {
    unsigned char q = p; int ii, cnt = 0;
    for (ii=0; ii<8; ++ii) { if (q&1) cnt++; q>>=1; }
    byte2cnt[p] = cnt;
  }
  
  size_t i,k,m = 0;
  for (i=0; i<(c.n>>3); ++i) m += byte2cnt[c.s[i]];   // full bytes
  for (k=0; k<(c.n&0x7); ++k) m += (c.s[i]>>k) & 0x1; // last byte
  return m;
}

void fmta_tryBinary2byteRLE_ifsmaller(cdata_t *c);

void decompress(cdata_t *c, cdata_t *expanded);
void decompress2(cdata_t *c);
void cdata_compress(cdata_t *c);

static inline uint64_t cdata_n(cdata_t *c) {
  if (!c->compressed) return c->n;
  cdata_t c2 = {0};
  decompress(c, &c2);
  uint64_t n = c2.n;
  free(c2.s);
  return n;
}

void fmt0_decompress(cdata_t *c, cdata_t *inflated);
void convertToFmt0(cdata_t *c);
#define FMT0_IN_SET(c, i) ((c).s[(i)>>3] & (1<<((i)&0x7)))
#define FMT0_SET(c, i) (c.s[(i)>>3] |= (1<<((i)&0x7)))

void fmt1_compress(cdata_t *c);
void fmt1_decompress(cdata_t *c, cdata_t *inflated);

// ----- format 2 (state data) ----
// key section + data section
// The key section and data section are separated by an extra '\0'.
// The key section is made of multiple c-strings concatenated by '\0'.
// The data section is either an RLE (compressed) or a integer vector (inflated).
// When compressed, the RLE is made of a value part and a length part.
// The value part size is defined by a uint8_t that leads the data section.
// The length part is always 2 bytes in size.
void fmt2_compress(cdata_t *c);
void fmt2_decompress(cdata_t *c, cdata_t *inflated);
void fmt2_set_aux(cdata_t *c);
uint8_t* fmt2_get_data(cdata_t *c);
uint64_t fmt2_get_keys_n(cdata_t *c);
uint64_t f2_get_uint64(cdata_t *c, uint64_t i);
char* f2_get_string(cdata_t *c, uint64_t i);

void fmt3_compress(cdata_t *c);
void fmt3_decompress(cdata_t *c, cdata_t *inflated);
void f3_set_mu(cdata_t *c, uint64_t i, uint64_t M, uint64_t U);
uint64_t f3_get_mu(cdata_t *c, uint64_t i);
#define MU2beta(mu) (double) ((mu)>>32) / (((mu)>>32) + ((mu)&0xffffffff))
#define MU2cov(mu) (((mu)>>32) + ((mu)&0xffffffff))

void fmt4_compress(cdata_t *c);
void fmt4_decompress(cdata_t *c, cdata_t *inflated);

void fmt5_compress(cdata_t *c);
void fmt5_decompress(cdata_t *c, cdata_t *inflated);

void fmt6_compress(cdata_t *c);
void fmt6_decompress(cdata_t *c, cdata_t *inflated);
#define FMT6_IN_SET(c, i) ((c).s[i>>2] & (1<<((i&0x3)*2)))
#define FMT6_IN_UNI(c, i) ((c).s[i>>2] & (1<<((i&0x3)*2+1)))
#define FMT6_SET0(c, i) ((c).s[i>>2] = ((c).s[i>>2] & ~(3<<((i&0x3)*2))) | (2<<((i&0x3)*2))) // 10
#define FMT6_SET1(c, i) ((c).s[i>>2] |= (3<<((i&0x3)*2))) // 11
#define FMT6_SET_NA(c, i) ((c).s[i>>2] &= (~(3<<((i&0x3)*2)))) // 00

int fmt7_next_bed(cdata_t *c);
uint64_t fmt7_data_length(cdata_t *c);
cdata_t fmt7_sliceToBlock(cdata_t *cr, uint64_t beg, uint64_t end);
cdata_t fmt7_sliceToIndices(cdata_t *cr, int64_t *row_indices, int64_t n_indices);
cdata_t fmt7_sliceToMask(cdata_t *cr, cdata_t *c_mask);

static inline void slice(cdata_t *c, uint64_t beg, uint64_t end, cdata_t *c_sliced) {

  if (c->compressed) {
    REprintf("[%s:%d] Cannot slice compressed data.\n", __func__, __LINE__);
    error("Abort.");
  }
  if (end > c->n-1) end = c->n-1;
  if (end < beg) wzfatal("Slicing negative span.");

  c_sliced->s = realloc(c_sliced->s, (end-beg+1)*c->unit);
  memcpy(c_sliced->s, c->s+beg*c->unit, (end-beg+1)*c->unit);
  c_sliced->n = end-beg+1;
  c_sliced->compressed = 0;
  c_sliced->fmt = c->fmt;
  c_sliced->unit = c->unit;
}

typedef struct row_reader_t {
  uint64_t index;
  char *chrm;                   // on cdata_t.s
  uint64_t loc;                 // on cdata_t.s
  uint64_t value;
} row_reader_t;

KHASH_MAP_INIT_STR(str2int, uint64_t) // Initialize a hashmap with keys as strings and values as uint64_t

typedef struct chromosome_t {
  uint64_t *locs;
  uint64_t *vals;
  uint64_t *inds;
  uint64_t n;
} chromosome_t;

typedef struct row_finder_t {
  chromosome_t *chrms;
  int n;
  khash_t(str2int) *h; // chromosome string > chromosome_t
} row_finder_t;

static inline void free_row_finder(row_finder_t *fdr) {
  for (int i=0; i<fdr->n; ++i) {
    free(fdr->chrms[i].locs);
    free(fdr->chrms[i].vals);
    free(fdr->chrms[i].inds);
  }
  free(fdr->chrms);
  kh_destroy(str2int, fdr->h);
}

int row_reader_next_loc(row_reader_t *rdr, cdata_t *c);
row_finder_t init_finder(cdata_t *cr);
uint64_t row_finder_search(char *chrm, uint64_t beg1, row_finder_t *fdr, cdata_t *cr);

#endif /* _CDATA_H */
