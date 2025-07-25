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

void find_y(struct point* p, mpz_t a, mpz_t b, mpz_t m, mpz_t count){
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
        mpz_add_ui(p->x, p->x, 1);
        mpz_add_ui(count, count, 1); // try next x
        find_y(p, a, b, m, count);
        return;
    }

    mpz_clears(x3, ax, bm, NULL);
    return;
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


void generate_random_mpz(mpz_t r, const mpz_t p) {
    mpz_t upper_limit, candidate;
    mpz_inits(upper_limit, candidate, NULL);

    // upper_limit = p - 2
    mpz_sub_ui(upper_limit, p, 2);

    size_t nbytes = (mpz_sizeinbase(upper_limit, 2) + 7) / 8;
    unsigned char *buffer = malloc(nbytes);

    do {
        if (RAND_bytes(buffer, nbytes) != 1) {
            fprintf(stderr, "Error: RAND_bytes failed\n");
            exit(1);
        }

        mpz_import(candidate, nbytes, 1, 1, 0, 0, buffer);
        // candidate = candidate + 1 to ensure r ∈ [1, p-2]
        mpz_add_ui(candidate, candidate, 1);
    } while (mpz_cmp(candidate, upper_limit) > 0);

    mpz_set(r, candidate);

    free(buffer);
    mpz_clears(upper_limit, candidate, NULL);
    return;
}

int main() {

// private key and public key generation

mpz_t priv_key, prime, a, b;
mpz_inits(priv_key, prime, a, b, NULL);
mpz_set_d(prime, 524287);
mpz_set_d(a, 5);
mpz_set_d(b, 87);

generate_random_mpz(priv_key, prime);


gmp_printf("Your private key is given by: %Zx\n", priv_key);

struct point g;
mpz_inits(g.x, g.y, NULL);

mpz_set_d(g.x, 3);
mpz_set_d(g.y, 47926);

struct point h;
mpz_inits(h.x, h.y, NULL);

point_multiplier(&g, &h, a, b, prime, priv_key);

gmp_printf("your public key is given in the format (p, g, h): (%Zd, (%Zd, %Zd), (%Zd, %Zd))", prime, g.x, g.y, h.x, h.y);

// part 1: Encryption
mpz_t x, y, msg, count;
mpz_inits(x, y, msg, count, NULL);
mpz_set_d(x, 3);
mpz_set_d(y, 47926);
mpz_set_d(count, 0);
gmp_printf("\n Please enter your message to encrypt: ");
gmp_scanf("%Zd", &msg);

struct point message, c1, c2;
mpz_inits(message.x, message.y, c1.x, c1.y, c2.x, c2.y, NULL);
mpz_set(message.x, msg);
find_y(&message, a, b, prime, count);

gmp_printf("The point corresponding to your message is (%Zd, %Zd)\n", message.x, message.y);
gmp_printf("your message x-coordinate shifted by %Zd\n", count);

generate_random_mpz(y, prime);
point_multiplier(&g, &c1, a, b, prime, y);
gmp_printf("The first ciphertext is given by: (%Zd, %Zd)\n", c1.x, c1.y);

point_multiplier(&h, &c2, a, b, prime, y);
add_points(&message, &c2, &c2, a, b, prime);
gmp_printf("The second ciphertext is given by: (%Zd, %Zd)\n", c2.x, c2.y);

//part 2: Decryption
struct point decrypted_message, s;
mpz_inits(decrypted_message.x, decrypted_message.y,s.x, s.y, NULL);
point_multiplier(&c1, &s, a, b, prime, priv_key);

gmp_printf("The value of s is given by: (%Zd, %Zd)\n", s.x, s.y);
mpz_neg(s.y, s.y); // negate y-coordinate for decryption
gmp_printf("the value of s_inv is (%Zd, %Zd)\n", s.x, s.y); // negate y-coordinate for decryption
add_points(&c2, &s, &decrypted_message, a, b, prime);
gmp_printf("The decrypted message is given by: (%Zd, %Zd)\n", decrypted_message.x, decrypted_message.y);

mpz_sub(decrypted_message.x, decrypted_message.x, count);
gmp_printf("The decrypted message x-coordinate shifted back to original by (count) %Zd is: %Zd\n", count, decrypted_message.x);

return 0;
}


