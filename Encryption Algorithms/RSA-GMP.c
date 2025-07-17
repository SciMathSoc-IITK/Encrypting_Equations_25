#include <stdio.h>
#include <gmp.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>

int main(){
    mpz_t p, q, candidate;
    gmp_randstate_t state;
    int runs = 25;
    int bits_p = 166;
    int bits_q = 173;
    mpz_init(p);
    mpz_init(q);
    mpz_init(candidate);
    gmp_randinit_mt(state);
    

    gmp_randseed_ui(state, time(NULL));

    mpz_urandomb(candidate, state, bits_p);
    mpz_setbit(candidate, bits_p-1);
    mpz_setbit(candidate, 0);
    mpz_nextprime(p, candidate);

    mpz_urandomb(candidate, state, bits_q);
    mpz_setbit(candidate, bits_q-1);
    mpz_setbit(candidate, 0);
    mpz_nextprime(q, candidate);

    if (mpz_probab_prime_p(p, runs) < 1 ||
        mpz_probab_prime_p(q, runs) < 1) {
        fprintf(stderr, "Primality test failed—try more MR rounds\n");
        return 1;
    }
    else{
        fprintf(stdout, "Primality test passed I guess.\n");

    }

    gmp_printf("p = %Zx\nq = %Zx\n",p , q);



    char *str_p = mpz_get_str(NULL, 10 , p);
    char *str_q = mpz_get_str(NULL, 10 , q);

    printf("\nString_p: %s\nString_q: %s \n", str_p, str_q);

    mpz_t n, p_n, e, d, p_1, q_1;
    mpz_inits(n, p_n, e, d, p_1, q_1, NULL);

    mpz_mul(n, p, q);

    gmp_printf("the value of n is equal to: %Zd", n);

    mpz_sub_ui(p_1, p, 1);
    mpz_sub_ui(q_1, q, 1);

    mpz_mul(p_n, p_1, q_1);
    gmp_printf("\nphi(n) is given by: %Zx\n", p_n);

    mpz_t gcd;
    mpz_init(gcd);

    mpz_ui_pow_ui(e, 10, 20);

    mpz_gcd(gcd, e, p_n);

    while(mpz_cmp_ui(gcd, 1) != 0){
        mpz_add_ui(e, e, 1);
        mpz_gcd(gcd, e, p_n);
    }
    gmp_printf("e is equal to: %Zd\n", e);
    gmp_printf("gcd of phi(n) and e is equal to: %Zx\n", gcd);

    mpz_invert(d, e, p_n);

    gmp_printf("d is equal to: %Zd\n", d);

     int size;
    printf("Enter maximum string length (make sure len(msg)<n): ");
    if (scanf("%d", &size) != 1 || size <= 0) {
        fprintf(stderr, "Invalid size\n");
        return 1;
    }

    // allocate space for 'size' chars + 1 for the terminating '\0'
    char *str = malloc((size + 1) * sizeof(char));
    if (!str) {
        perror("malloc");
        return 1;
    }

    // consume leftover newline from the scanf
    getchar();

    printf("Enter your string (up to %d characters): ", size);
    if (!fgets(str, size + 1, stdin)) {
        perror("fgets");
        free(str);
        return 1;
    }

    // strip trailing newline if present
    char *nl = str;
    while (*nl) {
        if (*nl == '\n') { *nl = '\0'; break; }
        nl++;
    }

    printf("You entered: \"%s\"\n", str);

    mpz_t msg;
    mpz_init(msg);

    mpz_set_str(msg, str, 10);

    mpz_t enc_msg;
    mpz_init(enc_msg);
    
    mpz_powm(enc_msg, msg, e, n);

    gmp_printf("\nyour encrypted message is: %Zx\n", enc_msg);

    mpz_t dec_msg;
    mpz_init(dec_msg);

    mpz_powm(dec_msg, enc_msg, d, n);
    gmp_printf("your decrypted message is: %Zd\n", dec_msg);
   free(str);
    mpz_clears(p, q, candidate,
               n, p_n, e, d, p_1, q_1, gcd,
               msg, enc_msg,dec_msg,NULL);
               
    gmp_randclear(state);
    return 0;



}
