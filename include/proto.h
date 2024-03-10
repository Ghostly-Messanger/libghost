#ifndef PROTO_H
#define PROTO_H

#include <gmp.h>

/**
 * ECDH - Elliptic curve Diffie-Hallman key exchange
 */

#define P "7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFED" // 2^255-19

/**
 * struct to present 2D point
 */
typedef struct {
    mpz_t x;
    mpz_t y;
} point;


/**
 * Check if point lies on curve
 * @param p point to check
 * @return 1 if point on curve, 0 otherwise
 */
int on_curve(point* p);

/**
 * Check if point is infinite (when x = 0 and y = 0)
 * @param p point to check
 * @return 1 if point is infinite, 0 otherwise
 */
int is_infinite(point* p);

/**
 * initialize point
 */
point* point_init();

/**
 * initialize point and copy from another variable
 */
point* point_init_set(point*);

/**
 * deinitialize point
 */
void point_close(point*);

/**
 * set X coordinate of point
 */
void point_set_x(point*,mpz_t);

/**
 * set Y coordinate of point
 */
void point_set_y(point*,mpz_t);

/**
 * set X coordinate of point from string (char*) with base
 */
void point_set_x_str(point*,char*,int);

/**
 * set Y coordinate of point from string (char*) with base
 */
void point_set_y_str(point*,char*,int);

/**
 * adds points
 */
point* point_add(point*, point*);

/**
 * multiply point by scalar
 */
point* point_mul(point*, mpz_t);

/**
 * print point in base n
 */
void point_print(point*,int n);

#endif

