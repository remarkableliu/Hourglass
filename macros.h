
#ifndef MACROS_H
#define MACROS_H

#ifdef __cplusplus
extern "C" {
#endif

#define MASK(i) ((1ULL << i) - 1)
#define MASK128(i) (((__uint128_t)1 << i) - 1)
#define POW(i) (1ULL << i)
#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))

#ifdef __cplusplus
}
#endif

#endif
