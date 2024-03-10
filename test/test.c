#include <proto.h>
#include <stdio.h>

int main(){
    point *p = point_init();
    point_set_x_str(p,"469",10);
    point_set_y_str(p,"4562",10);
    point_print(p,10);
    int a = on_curve(p);
    printf("a = %i\n",a);
    mpz_t t;
    mpz_init_set_ui(t,2);
    point *q = point_mul(p,t);
    point_print(q,10);
    int b = on_curve(q);
    printf("b = %i\n",b);
    point_close(p);
    return 0;
}