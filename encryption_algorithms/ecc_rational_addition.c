#include <stdio.h>
#include <gmp.h>
#include <stdbool.h>

// A struct representing a point on an elliptic curve over the rationals.
// Uses GMP's mpq_t (rational) type for x and y coordinates.
struct point_rat{
    mpq_t x;  // x-coordinate as a rational number
    mpq_t y;  // y-coordinate as a rational number
};

// Computes the cube of a rational number x and stores the result in x_3.
// x_3 = x * x * x
void cube_rat(mpq_t x, mpq_t x_3){
    mpq_mul(x_3, x, x);      // x_3 = x * x
    mpq_mul(x_3, x, x_3);    // x_3 = x * (x^2)
    return;
}

// Checks whether y_2 is a perfect square in Q, and if so, sets y to its square root.
// y_2: value to test (rational squared value)
// y: output parameter, set to sqrt(y_2) if integer root exists
bool check_update_y_rat(mpq_t y_2, mpq_t y){
    // Get numerator and denominator of y_2
    mpz_ptr num = mpq_numref(y_2);
    mpz_ptr den = mpq_denref(y_2);

    // If both numerator and denominator are perfect squares, take their root
    if(!(mpz_perfect_square_p(num) && mpz_perfect_square_p(den))) {
        return false;  // Not a square in Q
    } else {
        mpz_t sqrt_num, sqrt_den;
        mpz_inits(sqrt_num, sqrt_den, NULL);
        mpz_sqrt(sqrt_num, num);  // sqrt numerator
        mpz_sqrt(sqrt_den, den);  // sqrt denominator

        // Set y = sqrt_num / sqrt_den
        mpq_set_z(y, sqrt_num);
        mpq_set_den(y, sqrt_den);
        mpq_canonicalize(y);
        return true;     
    }
}

// Finds a rational y-coordinate for a given x on the elliptic curve y^2 = x^3 + A*x + B.
// If no rational y exists, increments x by 1 and recurses.
// p: pointer to point_rat containing x (input) and y (output)
// a, b: curve parameters (rational)
void find_y_rat(struct point_rat *p, mpq_t a, mpq_t b){
    mpq_t temp_x3;
    mpq_t A_x;
    mpq_inits(temp_x3, A_x, NULL);

    cube_rat(p->x, temp_x3);      // temp_x3 = x^3
    mpq_mul(A_x, p->x, a);        // A_x = A * x
    
    mpq_t y_2;
    mpq_init(y_2);
    mpq_add(y_2, temp_x3, A_x);   // y_2 = x^3 + A*x
    mpq_add(y_2, b, y_2);         // y_2 = x^3 + A*x + B

    // Try to compute rational y from y_2
    if(check_update_y_rat(y_2, p->y)){
        // Found valid y, clean up and return
        mpq_clear(temp_x3);
        mpq_clear(A_x);
        return;
    } else {
        // No rational y, increment x by 1 and retry
        mpq_clear(temp_x3);
        mpq_clear(A_x);
        
        mpq_t one;
        mpq_init(one);
        mpq_set_ui(one, 1, 1);      // one = 1
        mpq_add(p->x, p->x, one);    // x = x + 1
        mpq_clear(one);
        find_y_rat(p, a ,b);         // recurse with new x
    }
}

// Computes the slope (lambda) of the line joining p1 and p2 on the elliptic curve.
// p1, p2: input points
// slope: output rational
// a: curve parameter A
void lambda_rat(struct point_rat *p1, struct point_rat *p2, mpq_t slope, mpq_t a){
    // Declare temporaries
    mpq_t slope_num, slope_den, k;
    mpq_inits(slope_num, slope_den, k, NULL);

    if (mpq_cmp(p1->x, p2->x) != 0) {
        // Distinct x-coordinates: slope = (y2 - y1) / (x2 - x1)
        mpq_sub(slope_num, p2->y, p1->y);
        mpq_sub(slope_den, p2->x, p1->x);
    } else {
        // Same x (tangent): slope = (3*x1^2 + A) / (2*y1)
        mpq_mul(slope_num, p1->x, p1->x);      // x1^2
        mpq_set_si(k, 3, 1);
        mpq_mul(slope_num, slope_num, k);      // 3*x1^2
        mpq_add(slope_num, slope_num, a);      // 3*x1^2 + A

        mpq_set_si(k, 2, 1);
        mpq_mul(slope_den, p1->y, k);          // 2*y1
    }

    // Compute slope = numerator / denominator
    mpq_div(slope, slope_num, slope_den);

    // Clean up temporaries
    mpq_clears(slope_num, slope_den, k, NULL);
}

// Computes the y-intercept (constant) of the line through p1 and p2.
// v: output rational constant
// a, b: curve parameters (unused here but signature consistent)
void const_rat(struct point_rat *p1, struct point_rat *p2, mpq_t v, mpq_t a, mpq_t b){
    mpq_t const_num, const_den, k;
    mpq_inits(const_num, const_den, k, NULL);

    if(mpq_cmp(p1->x, p2->x) != 0){
        // Standard line: v = (y1*x2 - x1*y2) / (x2 - x1)
        mpq_mul(const_num, p1->y, p2->x);
        mpq_mul(k, p1->x, p2->y);
        mpq_sub(const_num, const_num, k);
        mpq_sub(const_den, p2->x, p1->x);
    }
    else{
        // Tangent line: v = (A*x1 + B - x1^3) / (2*y1)
        cube_rat(p1->x, const_num);           // x1^3
        mpq_mul(k, a, p1->x);                 // A*x1
        mpq_add(k, k, b);                     // A*x1 + B
        mpq_sub(const_num, k, const_num);     // A*x1 + B - x1^3

        mpq_set_si(k, 2, 1);
        mpq_mul(const_den, p1->y, k);        // 2*y1
    }
    mpq_div(v, const_num, const_den);        // v = num/den
    mpq_clears(const_num, const_den, k, NULL);
}

// Adds two rational points p1 and p2 on the elliptic curve and stores result in p3.
// a, b: curve parameters
void add_rat(struct point_rat *p1, struct point_rat *p2, struct point_rat *p3, mpq_t a, mpq_t b){
   mpq_t slope, constant, slope_2, slope_3, k;
   mpq_inits(slope, constant, slope_2, slope_3, k, NULL);

   // Compute slope and constant for the line through p1,p2
   lambda_rat(p1, p2, slope, a);
   const_rat(p1, p2, constant, a, b);

   // Initialize p3 coordinates
   mpq_inits(p3->x, p3->y, NULL);

   // p3.x = slope^2 - x1 - x2
   mpq_mul(slope_2, slope, slope);
   cube_rat(slope, slope_3);  // slope^3 for later use
   mpq_sub(p3->x, slope_2, p1->x);
   mpq_sub(p3->x, p3->x, p2->x);

   // p3.y = slope*(x1 + x2) - slope^3 - constant
   mpq_add(k, p1->x, p2->x);
   mpq_mul(p3->y, slope, k);
   mpq_sub(p3->y, p3->y, slope_3);
   mpq_sub(p3->y, p3->y, constant);

   return;
}

int main(){
    long long int x_1, x_2, A, B;

    // Prompt user for initial x-coordinates and curve parameters (integers)
    printf("please give your input for the values in the order x1, x2, A, B (make sure they are integers):");
    scanf("%lld %lld %lld %lld",&x_1, &x_2, &A, &B);

    // Temporary for converting to GMP types
    mpz_t temp;
    mpz_init(temp);

    // Declare points and curve parameters as rationals
    struct point_rat p1;
    struct point_rat p2;
    mpq_t a, b;
    mpq_inits(p1.x, p1.y, p2.x, p2.y, a,b , NULL);

    // Set p1.x and p2.x from user input
    mpz_set_si(temp, x_1);
    mpq_set_z(p1.x, temp);
    mpz_set_si(temp, x_2);
    mpq_set_z(p2.x, temp);

    // Set curve parameters a and b
    mpz_set_si(temp, A);
    mpq_set_z(a, temp);
    mpz_set_si(temp, B);
    mpq_set_z(b, temp);

    mpz_clear(temp);

    // Compute y-coordinates for p1 and p2
    find_y_rat(&p1, a, b);
    find_y_rat(&p2, a, b);

    // Print the two points
    gmp_printf("the coordinates are: (%Qd,%Qd), (%Qd,%Qd)\n", p1.x,p1.y, p2.x, p2.y);

    // Compute slope and constant of line joining p1,p2
    mpq_t slope, constant;
    struct point_rat p3;
    mpq_inits(slope, constant, p3.x, p3.y, NULL);

    lambda_rat(&p1, &p2, slope, a);
    const_rat(&p1, &p2, constant, a, b);
    gmp_printf("the slope and the constant of the two points is given by: %Qd, %Qd\n", slope, constant);

    // Add the points to get p3
    add_rat(&p1, &p2, &p3, a, b);
    gmp_printf("the value of P+Q is equal to: (%Qd, %Qd)\n",p3.x, p3.y);
    
    return 0;
}




