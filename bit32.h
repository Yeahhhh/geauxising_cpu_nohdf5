#ifndef BIT32_H
#define BIT32_H

#include <stdint.h>

/*
  LENG -           length of an integer
  NBETA -          number of parallel betas
  NBETA_MAX -      align NBETA into "LENG" boundaries
  NBETA_PER_WORD - number of betas combined into a word
  NWORD -          number of words for all betas

  must garantee NBETA <= LENG,
  this restriction can be overcomed by distributing temperature replicas on multiple words

  6 neighbors
  left     (xa,y ,z )    J0
  right    (xb,y ,z )    J1
  up       (x ,ya,z )    J2
  down     (x ,yb,z )    J3
  front    (x ,y ,za)    J4
  back     (x ,y ,zb)    J5
*/

// ACMSC-1 is buggy
#define ACMSC_FORMAT 0

typedef int32_t MSC_DATATYPE;
#define LENG 32
#define NBETA_MAX LENG
#define MASK_A  0xffffffff

#define NBETA 24
#define NBETA_PER_WORD 24
#define NWORD 1
#define NBIT_PER_SEG 4


#if ACMSC_FORMAT == 0
#define NSEG_PER_WORD 6
#define NBETA_PER_SEG 4
#define MASK_J  0xfc000000
#define MASK_J0 0x04000000
#define MASK_J1 0x08000000
#define MASK_J2 0x10000000
#define MASK_J3 0x20000000
#define MASK_J4 0x40000000
#define MASK_J5 0x80000000
#define SHIFT_J0 26
#define SHIFT_J1 27
#define SHIFT_J2 28
#define SHIFT_J3 29
#define SHIFT_J4 30
#define SHIFT_J5 31
#define MASK_S  0x00111111
#define MASK_S0 0x00ffffff
#define MASK_E  0xf
#define SHIFT_MAX 23
#endif
/*
  MASK_J  1111 11-- 0000 0000 0000 0000 0000 0000
  MASK_J0 0000 01-- 0000 0000 0000 0000 0000 0000
  MASK_J1 0000 10-- 0000 0000 0000 0000 0000 0000
  MASK_J2 0001 00-- 0000 0000 0000 0000 0000 0000
  MASK_J3 0010 00-- 0000 0000 0000 0000 0000 0000
  MASK_J4 0100 00-- 0000 0000 0000 0000 0000 0000
  MASK_J5 1000 00-- 0000 0000 0000 0000 0000 0000
  MASK_S  ---- ---- 0001 0001 0001 0001 0001 0001
  MASK_S0 ---- ---- 1111 1111 1111 1111 1111 1111
  iter0                *    *    *    *    *    *
  iter1               *    *    *    *    *    *
  iter2              *    *    *    *    *    *
  iter3             *    *    *    *    *    *
*/


#if ACMSC_FORMAT == 1
#define NSEG_PER_WORD 8
#define NBETA_PER_SEG 3
#define MASK_J  0x00888888
#define MASK_J0 0x00000008
#define MASK_J1 0x00000080
#define MASK_J2 0x00000800
#define MASK_J3 0x00008000
#define MASK_J4 0x00080000
#define MASK_J5 0x00800000
#define SHIFT_J0 3
#define SHIFT_J1 7
#define SHIFT_J2 11
#define SHIFT_J3 15
#define SHIFT_J4 19
#define SHIFT_J5 23
#define MASK_S  0x11111111
#define MASK_S0 0x77777777
#define MASK_E  0xf
#define SHIFT_MAX 31
#endif
/*
  MASK_J  0000 0000 1000 1000 1000 1000 1000 1000
  MASK_J0 0000 0000 0000 0000 0000 0000 0000 1000
  MASK_J1 0000 0000 0000 0000 0000 0000 1000 0000
  MASK_J2 0000 0000 0000 0000 0000 1000 0000 0000
  MASK_J3 0000 0000 0000 0000 1000 0000 0000 0000
  MASK_J4 0000 0000 0000 1000 0000 0000 0000 0000
  MASK_J5 0000 0000 1000 0000 0000 0000 0000 0000
  MASK_S  0001 0001 0001 0001 0001 0001 0001 0001
  MASK_S0 0111 0111 0111 0111 0111 0111 0111 0111
  iter0      *    *    *    *    *    *    *    *
  iter1     *    *    *    *    *    *    *    *
  iter2    *    *    *    *    *    *    *    *
*/



#endif /* BIT32_H */
