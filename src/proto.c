#include <stdio.h>
#include <malloc.h>
#include "proto.h"

int on_curve(point* p){
    if(!is_infinite(p)){
        // F -> y^2 = x^3 + 486662*x^2 + x
        mpz_t ys,xt,xs,xr;
        // y^2
        mpz_init_set(ys,p->y);
        mpz_pow_ui(ys,ys,2);
        // x^3
        mpz_init_set(xt,p->x);
        mpz_pow_ui(xt,xt,3);
        // 48666*x^2
        mpz_init_set(xs,p->x);
        mpz_pow_ui(xs,xs,2);
        mpz_mul_ui(xs,xs,486662);

        mpz_init(xr); // x^3 + 486662*x^2 + x
        mpz_add(xr,xr,xt);
        mpz_clear(xt);
        mpz_add(xr,xr,xs);
        mpz_clear(xs);
        mpz_add(xr,xr,p->x);

        mpz_t r;
        mpz_init(r);
        mpz_sub(r,ys,xr); // y^2 - (x^3 + 486662*x^2 + x)
        mpz_clear(xr);
        mpz_clear(ys);
        mpz_t prime;
        mpz_init_set_str(prime,P,16);
        mpz_mod(r,r,prime);
        int res = mpz_cmp_ui(r,0) != 0;
        mpz_clear(r);
        mpz_clear(prime);
        return res;
    }else{
        // p is infinity
        return 1;
    }
}

int is_infinite(point* p){
    return mpz_cmp_ui(p->x,0) == 0 && mpz_cmp_ui(p->y,0) == 0;
}

point* point_init(){
    point *p = malloc(sizeof (point));
    mpz_init(p->x);
    mpz_init(p->y);
    return p;
}

point* point_init_set(point* p){
    point *pt = point_init();
    mpz_set(pt->x,p->x);
    mpz_set(pt->y,p->y);
    return pt;
}

void point_close(point* p){
    mpz_clear(p->x);
    mpz_clear(p->y);
    free(p);
}

void point_set_x(point* p,mpz_t t){
    mpz_set(p->x,t);
}

void point_set_y(point* p,mpz_t t){
    mpz_set(p->x,t);
}

void point_set_x_str(point* p ,char *s,int b){
    mpz_set_str(p->x,s,b);
}

void point_set_y_str(point* p,char *s,int b){
    mpz_set_str(p->y,s,b);
}

point* point_add(point* a, point* b){
    point *p = point_init();
    if(!on_curve(a)){
        fprintf(stderr,"Error : point a not relies on curve\n");
        return p;
    }
    if(!on_curve(b)){
        fprintf(stderr,"Error : point b not relies on curve\n");
        return p ;
    }
    if(is_infinite(a)){
        return b;
    }else if(is_infinite(b)){
        return a;
    }else{
        mpz_t prime;
        mpz_init_set_str(prime,P,16);
        if(mpz_cmp(a->x,b->x) == 0 && mpz_cmp(a->y,b->y) != 0){
            return p;
        }
        mpz_t m,t;
        mpz_init(m);
        if(mpz_cmp(a->x,b->x) == 0){
            // inverse_mod(2 * y1, prime)
            mpz_init_set(t,a->y);    // t = y1
            mpz_mul_ui(t,t,2);      // t = 2 * y1
            mpz_invert(t,t,prime);  // t = inverse_mod(2 * y1, prime)
            mpz_set(m,t);
            // (3 * x1 * x1 + curve.a)
            mpz_set(t,a->x);         // t = x1
            mpz_pow_ui(t,t,2);      // t = x1^2
            mpz_mul_ui(t,t,3);      // t = 3 * x1^2
            mpz_add_ui(t,t,486662); // t = 3 * x1^2 + 486662
            mpz_mul(m,m,t); //  m = (3 * x1 * x1 + 486662) * inverse_mod(2 * y1, prime)
            mpz_clear(t);
        }else{
            mpz_init(t);
            mpz_sub(t,a->x,b->x);     //x1 - x2
            mpz_invert(t,t,prime);  // inverse_mod(x1 - x2, prime)
            mpz_set(m,t);
            mpz_sub(t,a->y,b->y);     // y1 - y2
            mpz_mul(m,m,t);         // m = (y1 - y2) * inverse_mod(x1 - x2, curve.p)
            mpz_clear(t);
        }
        mpz_t x,y;
        mpz_inits(t,x,y,NULL);
        //x = m * m - x1 - x2
        //y = y1 + m * (x - x1)

        mpz_pow_ui(x,m,2);  // m * m
        mpz_sub(t,a->x,b->x); // x1 - x2
        mpz_sub(x,x,t);     // m * m - x1 - x2

        mpz_mod(x,x,prime);

        mpz_sub(t,x,a->x);   // x - x1
        mpz_mul(t,t,m);     // m * (x - x1)
        mpz_add(y,t,a->y);   // y1 + m * (x - x1)

        mpz_neg(y,y);
        mpz_mod(y,y,prime);

        mpz_clear(t);
        point_set_x(p,x);
        mpz_clear(x);
        point_set_x(p,y);
        mpz_clear(y);
        if(on_curve(p)) {
            return p;
        }else{
            fprintf(stderr,"Error : point not relies on curve\n");
            return point_init();
        }
    }
}

point* point_mul(point* a,mpz_t b){
    point *p = point_init();
    if(mpz_cmp_ui(b,0) == 0){
        return p;
    }else {
        if (on_curve(a)) {
            point *t = point_init_set(a);
            int s = mpz_sizeinbase(b,2)+2;
            char binary[s];
            mpz_get_str(binary,2,b);
            for(int i = 0; i < s;i++){
                if(binary[i] == '1'){
                    p = point_add(p,t);
                }
                t = point_add(t,t);
            }
        } else {
            fprintf(stderr, "Error: point not relies on curve\n");
        }
    }
    return p;
}

void point_print(point* p,int n){
    int s = mpz_sizeinbase(p->x,n) + 2;
    char* temp = malloc(sizeof (char) * s);
    mpz_get_str(temp,n,p->x);
    printf("X : %s\n",temp);
    s = mpz_sizeinbase(p->y,n) + 2;
    temp = reallocarray(temp,sizeof (char),s);
    mpz_get_str(temp,n,p->y);
    printf("Y : %s\n",temp);
    free(temp);
}