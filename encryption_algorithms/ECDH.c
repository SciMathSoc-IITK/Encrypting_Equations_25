#include <stdio.h>
#include <openssl/rand.h>
#include <gmp.h>
#include <sys/types.h>

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

// R = P + Q on E: y^2 = x^3 + a x + b  over Fₘ
// We reserve (0,0) as the “point at infinity.”
// R = P + Q on E: y^2 = x^3 + a x + b  over Fₘ
// We reserve (0,0) as the “point at infinity.”

void add_points(struct point *P, struct point *Q, struct point *R, mpz_t a, mpz_t b, mpz_t m)
{
    // 1) Identity cases: (0,0) is ∞
    if (mpz_cmp_ui(P->x, 0)==0 && mpz_cmp_ui(P->y,0)==0) {
        mpz_set(R->x, Q->x);
        mpz_set(R->y, Q->y);
        return;
    }
    if (mpz_cmp_ui(Q->x, 0)==0 && mpz_cmp_ui(Q->y,0)==0) {
        mpz_set(R->x, P->x);
        mpz_set(R->y, P->y);
        return;
    }

    // 2) P == –Q ?  then P+Q = ∞
    {
        mpz_t sumy;
        mpz_init(sumy);
        mpz_add(sumy, P->y, Q->y);
        mpz_mod(sumy, sumy, m);
        if (mpz_cmp(P->x, Q->x)==0 && mpz_cmp_ui(sumy,0)==0) {
            mpz_set_ui(R->x, 0);
            mpz_set_ui(R->y, 0);
            mpz_clear(sumy);
            return;
        }
        mpz_clear(sumy);
    }

    // 3) Compute slope λ
    mpz_t lam, inv, tmp;
    mpz_inits(lam, inv, tmp, NULL);

    if (mpz_cmp(P->x, Q->x)==0 && mpz_cmp(P->y, Q->y)==0) {
        // tangent: λ = (3*x₁² + a) / (2*y₁)
        mpz_mul(lam, P->x, P->x);       // x₁²
        mpz_mul_ui(lam, lam, 3);        // 3*x₁²
        mpz_add(lam, lam, a);           // 3*x₁² + a
        mpz_mul_ui(tmp, P->y, 2);       // 2*y₁
        mpz_invert(inv, tmp, m);        // inv = (2*y₁)⁻¹
    } else {
        // chord:   λ = (y₂ - y₁) / (x₂ - x₁)
        mpz_sub(lam, Q->y, P->y);       // y₂ - y₁
        mpz_sub(tmp, Q->x, P->x);       // x₂ - x₁
        mpz_invert(inv, tmp, m);        // inv = (x₂ - x₁)⁻¹
    }
    mpz_mul(lam, lam, inv);
    mpz_mod(lam, lam, m);

    // 4) x₃ = λ² - x₁ - x₂   (mod m)
    mpz_mul(tmp, lam, lam);
    mpz_sub(tmp, tmp, P->x);
    mpz_sub(tmp, tmp, Q->x);
    mpz_mod(R->x, tmp, m);

    // 5) y₃ = λ*(x₁ - x₃) - y₁   (mod m)
    mpz_sub(tmp, P->x, R->x);        // x₁ - x₃
    mpz_mul(tmp, lam, tmp);          // λ*(x₁ - x₃)
    mpz_sub(tmp, tmp, P->y);         // λ*(x₁ - x₃) - y₁
    mpz_mod(R->y, tmp, m);

    // 6) Cleanup
    mpz_clears(lam, inv, tmp, NULL);
}


void point_multiplier(struct point *G, struct point *out, mpz_t a, mpz_t b, mpz_t m,mpz_t n)
{
    // Prepare result = identity, base = G
    struct point result, base, tmp;
    mpz_inits(result.x, result.y,
              base.x, base.y,
              tmp.x, tmp.y, NULL);

    // identity = (0,0)
    mpz_set_ui(result.x, 0);
    mpz_set_ui(result.y, 0);

    mpz_set(base.x, G->x);
    mpz_set(base.y, G->y);

    // Scan bits of n from LSB to MSB
    size_t n_bits = mpz_sizeinbase(n, 2);
    for (size_t i = 0; i < n_bits; ++i) {
        if (mpz_tstbit(n, i)) {
            // result = result + base
            add_points(&result, &base, &tmp, a, b, m);
            mpz_set(result.x, tmp.x);
            mpz_set(result.y, tmp.y);
        }
        // base = base + base (point‐doubling)
        add_points(&base, &base, &tmp, a, b, m);
        mpz_set(base.x, tmp.x);
        mpz_set(base.y, tmp.y);
    }

    // Copy final result into out
    mpz_set(out->x, result.x);
    mpz_set(out->y, result.y);

    // Clean up
    mpz_clears(result.x, result.y,
               base.x, base.y,
               tmp.x, tmp.y, NULL);
}


void generate_secure_mpz(mpz_t result, int num_bits) {
    int num_bytes = (num_bits + 7) / 8;  // convert bits to bytes
    unsigned char buffer[64];  // up to 512 bits
    if (num_bytes > sizeof(buffer)) {
        fprintf(stderr, "Too many bits requested.\n");
        return;
    }

    // Generate random bytes
    if (RAND_bytes(buffer, num_bytes) != 1) {
        fprintf(stderr, "OpenSSL RAND_bytes failed\n");
        return;
    }

    // Import random bytes into GMP integer
    mpz_import(result, num_bytes, 1, 1, 0, 0, buffer);

    // Make sure it's within desired bit length
    if (mpz_sizeinbase(result, 2) > num_bits) {
        mpz_fdiv_r_2exp(result, result, num_bits);  // truncate to num_bits
    }
}

int main() {
    mpz_t priv_key_a, priv_key_b;
    mpz_inits(priv_key_a, priv_key_b, NULL);

    generate_secure_mpz(priv_key_a, 256); 
    generate_secure_mpz(priv_key_b, 256); // generate a 256-bit random number

    gmp_printf("Your Private keys are: %Zx, %Zx\n", priv_key_a, priv_key_b);  // print in hex

    printf("\nThe private keys are 256-bit Randomly generated numbers. \nNote that we are printing the private keys here only for presentation.\nIn a real application, you would never print or expose your private keys.\n\n");

    //curve parameters:
    mpz_t a, b, m;
    mpz_inits(a, b, m, NULL);
    struct point g;
    mpz_inits(g.x, g.y, NULL);
    mpz_set_d(g.x, 3);
    mpz_set_d(g.y, 47926);
    mpz_set_d(a, 5);
    mpz_set_d(b, 87);
    mpz_set_d(m, 524287);

    //generate public keys:
    struct point ag, bg, abg, bag;
    mpz_inits(ag.x, ag.y, bg.x, bg.y, abg.x, abg.y, bag.x, bag.y, NULL);

    point_multiplier(&g, &ag, a, b, m, priv_key_a);
    point_multiplier(&g, &bg, a, b, m, priv_key_b);

    gmp_printf("the elliptic curve is given by: y^2 = x^3 + %Zd x + %Zd (mod %Zd)\n", a, b, m);
    gmp_printf("the value of the generator is: g = (%Zx, %Zx)\n", g.x, g.y);
    gmp_printf("the values of our public keys are: ag = (%Zx, %Zx), bg = (%Zx,%Zx)\n", ag.x, ag.y, bg.x, bg.y);

    point_multiplier(&ag, &bag, a, b, m, priv_key_b);
    point_multiplier(&bg, &abg, a, b, m, priv_key_a);
    
    gmp_printf("the value of the points abg and bag is: abg = (%Zx, %Zx), bag = (%Zx, %Zx)\n", abg.x, abg.y, bag.x, bag.y);
    if (mpz_cmp(abg.x, bag.x) == 0 && mpz_cmp(abg.y, bag.y) == 0)
    {
        printf("Success: a*b*g == b*a*g, shared secret matches\n");
    }else {
        printf("Error: a*b*g != b*a*g, something went wrong\n");
    }
    mpz_clears(priv_key_a, priv_key_b, a, b, m, g.x, g.y, ag.x, ag.y, abg.x, abg.y, bag.x, bag.y, NULL);
    
return 0;
}
