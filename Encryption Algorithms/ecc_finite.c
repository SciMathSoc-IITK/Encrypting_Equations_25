// #include <stdio.h>
// #include <gmp.h>
// #include <stdlib.h>

// #define READ_INPUT(varname, prompt)                   \
//         printf("Enter %s: ", prompt);                     \
//         scanf("%1023s", input_str);                       \
//         if (mpz_set_str(varname, input_str, 10) != 0) {   \
//             printf("Invalid input for %s.\n", prompt);    \
//             return 1;                                     \
//         }

// struct point{
//     mpz_t x, y;
// };

// void cube(mpz_t x, mpz_t x3, mpz_t m){
//     mpz_pow_ui(x3, x, 3);
//     mpz_mod(x3, x3, m);
// }
// void y_2(mpz_t y_2, mpz_t y, mpz_t m) {
//     // Handle trivial cases
//     if (mpz_cmp_ui(y_2, 0) == 0) {
//         mpz_set_ui(y, 0);
//         return;
//     }
//     // Check Legendre symbol (y_2|m); must be 1 for a root to exist
//     if (mpz_legendre(y_2, m) != 1) {
//         mpz_set_ui(y, 0);
//         return;
//     }

//     // Special case m ≡ 3 (mod 4): root = y_2^((m+1)/4) mod m
//     if (mpz_tstbit(m, 0) && mpz_tstbit(m, 1)) {  // m % 4 == 3
//         mpz_t exp;
//         mpz_init(exp);
//         // exp = (m + 1) / 4
//         mpz_add_ui(exp, m, 1);
//         mpz_fdiv_q_ui(exp, exp, 4);
//         mpz_powm(y, y_2, exp, m);
//         mpz_clear(exp);
//         return;
//     }

//     // Tonelli–Shanks for general odd prime
//     mpz_t q, s, z, c, t, r, b, tmp;
//     unsigned long int i, e;

//     mpz_inits(q, s, z, c, t, r, b, tmp, NULL);

//     // Factor out powers of 2 from m-1: m-1 = q * 2^s with q odd
//     mpz_sub_ui(q, m, 1);
//     mpz_set_ui(s, 0);
//     while (mpz_even_p(q)) {
//         mpz_fdiv_q_2exp(q, q, 1);
//         mpz_add_ui(s, s, 1);
//     }

//     // Find a quadratic non-residue z (Legendre(z,m) == -1)
//     mpz_set_ui(z, 2);
//     while (mpz_legendre(z, m) != -1) {
//         mpz_add_ui(z, z, 1);
//     }

//     // c = z^q mod m
//     mpz_powm(c, z, q, m);
//     // r = y_2^((q+1)/2) mod m
//     mpz_add_ui(tmp, q, 1);
//     mpz_fdiv_q_ui(tmp, tmp, 2);
//     mpz_powm(r, y_2, tmp, m);
//     // t = y_2^q mod m
//     mpz_powm(t, y_2, q, m);

//     // m = s (we'll reuse variable s for the exponent count)
//     e = mpz_get_ui(s);
//     while (mpz_cmp_ui(t, 1) != 0) {
//         // find smallest 0 < i < e such that t^(2^i) ≡ 1
//         mpz_set(tmp, t);
//         for (i = 1; i < e; i++) {
//             mpz_powm_ui(tmp, tmp, 2, m);
//             if (mpz_cmp_ui(tmp, 1) == 0) break;
//         }
//         // b = c^(2^(e-i-1)) mod m
//         mpz_powm_ui(b, c, 1UL << (e - i - 1), m);
//         // update
//         mpz_mul(r, r, b); mpz_mod(r, r, m);
//         mpz_mul(c, b, b); mpz_mod(c, c, m);
//         mpz_mul(t, t, c); mpz_mod(t, t, m);
//         e = i;
//     }

//     mpz_mod(r, r, m);
//     mpz_set(y, r);

//     mpz_clears(q, s, z, c, t, r, b, tmp, NULL);
// }

// void find_y(struct point* p, mpz_t a, mpz_t b, mpz_t m){
//     mpz_t x_3, a_x , b_m;
//     mpz_inits(x_3, a_x, b_m, NULL);

//     cube(p->x, x_3, m);

//     mpz_mul(a_x, p->x, a);
//     mpz_mod(a_x, a_x, m);
    
//     mpz_mod(b_m, b, m);

//     mpz_add(p->y, x_3, a_x);
//     mpz_add(p->y, p->y, b_m);

//     mpz_mod(p->y, p->y, m);

//     y_2(p->y, p->y, m);

//     return;
// }



// void find_slope(struct point* p, struct point* q, mpz_t slope,mpz_t a, mpz_t m)
// {
//     mpz_t num, den, inv;
//     mpz_inits(num, den, inv, NULL);

//     if (mpz_cmp(p->x, q->x) != 0) {
//         // chord: (y_q - y_p) / (x_q - x_p)
//         mpz_sub(num, q->y, p->y);
//         mpz_sub(den, q->x, p->x);
//     } else {
//         // tangent: (3*x_p^2 + a) / (2*y_p)
//         mpz_mul(num, p->x, p->x);        // x_p^2
//         mpz_mul_ui(num, num, 3);         // 3*x_p^2
//         mpz_add(num, num, a);            // 3*x_p^2 + a

//         mpz_mul_ui(den, p->y, 2);        // 2*y_p
//     }

//     // den^{-1} mod m
//     if (mpz_invert(inv, den, m) == 0) {
//         // no inverse ⇒ slope undefined (den ≡ 0 mod m)
//         mpz_set_ui(slope, 0);
//     } else {
//         mpz_mul(slope, num, inv);
//         mpz_mod(slope, slope, m);
//     }

//     mpz_clears(num, den, inv, NULL);
// }
// void find_const(struct point* p, struct point* q, mpz_t v, mpz_t a, mpz_t b, mpz_t m)
// {
//     mpz_t num, den, inv, tmp;
//     mpz_inits(num, den, inv, tmp, NULL);

//     if (mpz_cmp(p->x, q->x) != 0) {
//         // chord: v = (y_p*x_q − x_p*y_q) / (x_q − x_p)
//         mpz_mul(num,  p->y,  q->x);  // y_p * x_q
//         mpz_mul(tmp,  p->x,  q->y);  // x_p * y_q
//         mpz_sub(num,  num,         tmp);
//         mpz_sub(den,  q->x,        p->x);
//     } else {
//         // tangent: v = (a*x_p + b − x_p^3) / (2*y_p)
//         // tmp = x_p^3
//         mpz_mul(tmp,  p->x,  p->x); 
//         mpz_mul(tmp,  tmp,   p->x);

//         // num = a*x_p + b − tmp
//         mpz_mul(num,  a,     p->x);
//         mpz_add(num,  num,   b);
//         mpz_sub(num,  num,   tmp);

//         // den = 2*y_p
//         mpz_mul_ui(den, p->y, 2);
//     }

//     // inv = den^{-1} mod m
//     // since m is prime, inverse always exists for den ≠ 0
//     if (mpz_invert(inv, den, m) == 0) {
//         // (shouldn't happen if den != 0 mod m)
//         mpz_set_ui(v, 0);
//     } else {
//         // v = num * inv mod m
//         mpz_mul(v, num, inv);
//         mpz_mod(v, v, m);
//     }

//     mpz_clears(num, den, inv, tmp, NULL);
// }
// void add_points(struct point *p, struct point *q, struct point *p2, mpz_t a, mpz_t b, mpz_t m){
// mpz_t slope, constant, slope_2, slope_3, k;
// mpz_inits(slope, constant, slope_2, slope_3, k, NULL);

// find_slope(p, q, slope, a, m);
// find_const(p, q, constant, a, b, m);

// mpz_pow_ui(slope_2, slope, 2);
// mpz_pow_ui(slope_3, slope, 3);

// mpz_sub(p2->x, slope_2, p->x);
// mpz_sub(p2->x, p2->x, q->x);
// mpz_mod(p2->x, p2->x, m);

// mpz_add(p2->y, p->x, q->x);
// mpz_mul(p2->y, p2->y, slope);
// mpz_sub(p2->y, p2->y, slope_3);
// mpz_sub(p2->y, p2->y, constant);
// mpz_mod(p2->y, p2->y, m);

// return;
// }
// int main(){
//     char input_str[1024];

//     struct point p1, p2;
//     mpz_t a, b, m;

//     // Initialize all mpz_t variables
//     mpz_inits(p1.x, p2.x, p1.y, p2.y, a, b, m, NULL);

//     // Output the values
//     READ_INPUT(p1.x, "p1.x");
//     READ_INPUT(p2.x, "p2.x");
//     READ_INPUT(a, "a");
//     READ_INPUT(b, "b");
//     READ_INPUT(m, "m");

//     struct point p3;
//     mpz_t slope, constant;
//     mpz_inits(p3.x, p3.y, slope, constant, NULL);

//     find_y(&p1, a, b, m);
//     find_y(&p2, a, b, m);

//     gmp_printf("\nthe value of coordinates (x1, y1) and (x2, y2) is given by: (%Zd, %Zd), (%Zd, %Zd)\n", p1.x, p1.y, p2.x, p2.y);

//     add_points(&p1, &p2, &p3, a, b, m);

//     gmp_printf("the value of point 3 is: (%Zd, %Zd)", p3.x, p3.y);
//     // Clear memory
//     mpz_clears(p1.x, p2.x, a, b, m, NULL);


//     return 0; 
// }

#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>

#define READ_INPUT(varname, prompt)                     \
    printf("Enter %s: ", prompt);                       \
    scanf("%1023s", input_str);                         \
    if (mpz_set_str(varname, input_str, 10) != 0) {     \
        printf("Invalid input for %s.\n", prompt);      \
        return 1;                                       \
    }

// forward‐declare
int y_2(mpz_t y2, mpz_t y, mpz_t m);

struct point {
    mpz_t x, y;
};

void cube(mpz_t x, mpz_t x3, mpz_t m){
    mpz_pow_ui(x3, x, 3);
    mpz_mod(x3, x3, m);
}

// now returns 1 on success, 0 if no square‐root exists
int y_2(mpz_t y2, mpz_t y, mpz_t m) {
    if (mpz_cmp_ui(y2, 0) == 0) {
        mpz_set_ui(y, 0);
        return 1;
    }
    if (mpz_legendre(y2, m) != 1) {
        return 0;      // no root exists
    }

    // fast path m ≡ 3 mod 4
    if (mpz_tstbit(m, 0) && mpz_tstbit(m, 1)) {
        mpz_t exp;
        mpz_init(exp);
        mpz_add_ui(exp, m, 1);
        mpz_fdiv_q_ui(exp, exp, 4);
        mpz_powm(y, y2, exp, m);
        mpz_clear(exp);
        return 1;
    }

    // Tonelli–Shanks...
    mpz_t q, s, z, c, t, r, b, tmp;
    unsigned long int i, e;
    mpz_inits(q, s, z, c, t, r, b, tmp, NULL);

    mpz_sub_ui(q, m, 1);
    mpz_set_ui(s, 0);
    while (mpz_even_p(q)) {
        mpz_fdiv_q_2exp(q, q, 1);
        mpz_add_ui(s, s, 1);
    }

    mpz_set_ui(z, 2);
    while (mpz_legendre(z, m) != -1) {
        mpz_add_ui(z, z, 1);
    }

    mpz_powm(c, z, q, m);
    mpz_add_ui(tmp, q, 1);
    mpz_fdiv_q_ui(tmp, tmp, 2);
    mpz_powm(r, y2, tmp, m);
    mpz_powm(t, y2, q, m);

    e = mpz_get_ui(s);
    while (mpz_cmp_ui(t, 1) != 0) {
        mpz_set(tmp, t);
        for (i = 1; i < e; i++) {
            mpz_powm_ui(tmp, tmp, 2, m);
            if (mpz_cmp_ui(tmp, 1) == 0) break;
        }
        mpz_powm_ui(b, c, 1UL << (e - i - 1), m);
        mpz_mul(r, r, b); mpz_mod(r, r, m);
        mpz_mul(c, b, b); mpz_mod(c, c, m);
        mpz_mul(t, t, c); mpz_mod(t, t, m);
        e = i;
    }

    mpz_mod(r, r, m);
    mpz_set(y, r);
    mpz_clears(q, s, z, c, t, r, b, tmp, NULL);
    return 1;
}

void find_y(struct point* p, mpz_t a, mpz_t b, mpz_t m){
    mpz_t x3, ax, bm;
    mpz_inits(x3, ax, bm, NULL);

    cube(p->x, x3, m);
    mpz_mul(ax, p->x, a);  mpz_mod(ax, ax, m);
    mpz_mod(bm, b, m);

    mpz_add(p->y, x3, ax);
    mpz_add(p->y, p->y, bm);
    mpz_mod(p->y, p->y, m);

    // if no sqrt exists, error out
    if (! y_2(p->y, p->y, m)) {
        fprintf(stderr, "Invalid x input: no square root for y^2 mod m\n");
        exit(1);
    }

    mpz_clears(x3, ax, bm, NULL);
}

void find_slope(struct point* p, struct point* q, mpz_t slope,mpz_t a, mpz_t m)
{
    mpz_t num, den, inv;
    mpz_inits(num, den, inv, NULL);

    if (mpz_cmp(p->x, q->x) != 0) {
        // chord: (y_q - y_p) / (x_q - x_p)
        mpz_sub(num, q->y, p->y);
        mpz_sub(den, q->x, p->x);
    } else {
        // tangent: (3*x_p^2 + a) / (2*y_p)
        mpz_mul(num, p->x, p->x);        // x_p^2
        mpz_mul_ui(num, num, 3);         // 3*x_p^2
        mpz_add(num, num, a);            // 3*x_p^2 + a

        mpz_mul_ui(den, p->y, 2);        // 2*y_p
    }

    // den^{-1} mod m
    if (mpz_invert(inv, den, m) == 0) {
        // no inverse ⇒ slope undefined (den ≡ 0 mod m)
        mpz_set_ui(slope, 0);
    } else {
        mpz_mul(slope, num, inv);
        mpz_mod(slope, slope, m);
    }

    mpz_clears(num, den, inv, NULL);
}
void find_const(struct point* p, struct point* q, mpz_t v, mpz_t a, mpz_t b, mpz_t m)
{
    mpz_t num, den, inv, tmp;
    mpz_inits(num, den, inv, tmp, NULL);

    if (mpz_cmp(p->x, q->x) != 0) {
        // chord: v = (y_p*x_q − x_p*y_q) / (x_q − x_p)
        mpz_mul(num,  p->y,  q->x);  // y_p * x_q
        mpz_mul(tmp,  p->x,  q->y);  // x_p * y_q
        mpz_sub(num,  num,         tmp);
        mpz_sub(den,  q->x,        p->x);
    } else {
        // tangent: v = (a*x_p + b − x_p^3) / (2*y_p)
        // tmp = x_p^3
        mpz_mul(tmp,  p->x,  p->x); 
        mpz_mul(tmp,  tmp,   p->x);

        // num = a*x_p + b − tmp
        mpz_mul(num,  a,     p->x);
        mpz_add(num,  num,   b);
        mpz_sub(num,  num,   tmp);

        // den = 2*y_p
        mpz_mul_ui(den, p->y, 2);
    }

    // inv = den^{-1} mod m
    // since m is prime, inverse always exists for den ≠ 0
    if (mpz_invert(inv, den, m) == 0) {
        // (shouldn't happen if den != 0 mod m)
        mpz_set_ui(v, 0);
    } else {
        // v = num * inv mod m
        mpz_mul(v, num, inv);
        mpz_mod(v, v, m);
    }

    mpz_clears(num, den, inv, tmp, NULL);
}
void add_points(struct point *p, struct point *q, struct point *p2, mpz_t a, mpz_t b, mpz_t m){
mpz_t slope, constant, slope_2, slope_3, k;
mpz_inits(slope, constant, slope_2, slope_3, k, NULL);

find_slope(p, q, slope, a, m);
find_const(p, q, constant, a, b, m);

mpz_pow_ui(slope_2, slope, 2);
mpz_pow_ui(slope_3, slope, 3);

mpz_sub(p2->x, slope_2, p->x);
mpz_sub(p2->x, p2->x, q->x);
mpz_mod(p2->x, p2->x, m);

mpz_add(p2->y, p->x, q->x);
mpz_mul(p2->y, p2->y, slope);
mpz_sub(p2->y, p2->y, slope_3);
mpz_sub(p2->y, p2->y, constant);
mpz_mod(p2->y, p2->y, m);

return;
}

int main(){
    char input_str[1024];
    struct point p1, p2;
    mpz_t a, b, m;

    mpz_inits(p1.x, p1.y, p2.x, p2.y, a, b, m, NULL);

    READ_INPUT(p1.x, "p1.x");
    READ_INPUT(p2.x, "p2.x");
    READ_INPUT(a, "a");
    READ_INPUT(b, "b");
    READ_INPUT(m, "m");

    find_y(&p1, a, b, m);
    find_y(&p2, a, b, m);

    struct point p3;
    mpz_inits(p3.x, p3.y, NULL);
    gmp_printf("\nthe value of coordinates (x1, y1) and (x2, y2) is given by: (%Zd, %Zd), (%Zd, %Zd)\n", p1.x, p1.y, p2.x, p2.y);

    add_points(&p1, &p2, &p3, a, b, m);

    gmp_printf("the value of point 3 is: (%Zd, %Zd)", p3.x, p3.y);
   
    mpz_clears(p1.x, p1.y, p2.x, p2.y, a, b, m, NULL);
    return 0;
}
