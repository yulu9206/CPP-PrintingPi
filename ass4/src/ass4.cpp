#include <iostream>
#include <iomanip>
#include <mpir.h>
#include <stdlib.h>
#include <string.h>
#include <chrono>

using namespace std;

const int MAX_ITERATIONS = 100;
const int PLACES         = 10000;        // desired decimal places
const int PRECISION      = PLACES + 1;  // +1 for the digit 3 before the decimal

const int BASE       = 10;  // base 10 numbers
const int BIT_COUNT  = 8;   // bits per machine word

const int BLOCK_SIZE = 10;                // print digits in blocks
const int LINE_SIZE  = 100;               // digits to print per line
const int LINE_COUNT = PLACES/LINE_SIZE;  // lines to print
const int GROUP_SIZE = 5;                 // line grouping size

/**
 * Compute the cube root of a positive integer.
 * @param x where to store the result.
 * @param a the number whose cube root to compute.
 */
void cube_root(mpf_t& x, mpf_t a);

void nonic_algo(mpf_t& pi);
/**
 * Print the decimal places of a multiple-precision number x.
 * @param pi the multiple-precision number to print.
 */

void print(const mpf_t& pi);


/**
 * The main.
 */
int main()
{
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();
//	printf("111");
	cout << endl;
	mpf_set_default_prec(BIT_COUNT*PRECISION);  // precision in bits
	mpf_t pi; mpf_init(pi);
	nonic_algo(pi);
//	print pi
	cout << "calculated first ten thousand digits of pi is: " << endl;
	print(pi);
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	std::cout << "Total elapsed time: " << elapsed_seconds.count() << " seconds.";
	std::cout << std::endl;
//	printf("222");
	return 0;
}

void cube_root(mpf_t& x, mpf_t a)
{
//	std::chrono::time_point<std::chrono::system_clock> start, end;
//	start = std::chrono::system_clock::now();

	//multiple-precision variables
	mpf_t x_prev; mpf_init(x_prev);
	mpf_t temp1; mpf_init(temp1);
	mpf_t temp2; mpf_init(temp2);
	mpf_t two_a; mpf_init(two_a);
	mpf_t x_cubed; mpf_init(x_cubed);

	//constant 3
	mpf_t three; mpf_init(three); mpf_set_str(three, "3", BASE);

	//set an initial estimate for x
	mpf_div(x, a, three);  // x = a/3

	int n = 0; //  iteration counter

	 // loop until two consecutive values are equal or up to MAX_ITERRATION times.
	do
	{
		mpf_set(x_prev, x);

		mpf_mul(x_cubed, x, x);
		mpf_mul(x_cubed, x_cubed, x); //x_cubed = x^3
		mpf_add(two_a, a, a); //two_a = 2a
		mpf_add(temp1, x_cubed, two_a); //temp1 = x^3 + 2a
		mpf_add(temp2, x_cubed, x_cubed); //temp2 = 2x^3
		mpf_add(temp2, temp2, a); //temp2 = 2x^3 + a
		mpf_div(temp1, temp1, temp2); //temp1 = (x^3 + 2a)/(2x^3 + a);
		mpf_mul(x, x, temp1); //x = x((x^3 + 2a)/(2x^3 + a));

		n++;
	} while ((mpf_cmp(x, x_prev) != 0)&&(n < MAX_ITERATIONS));
//	end = std::chrono::system_clock::now();
//	std::chrono::duration<double> elapsed_seconds = end-start;
//	std::cout << "Calculating cube root took " << elapsed_seconds.count() << " seconds.";
//	std::cout << std::endl;
}


//The nonic algorithm.

void nonic_algo(mpf_t& pi)
{
//	multiple-precision constants.
	mpf_t one, two, three, nine, twenty_seven, one_third;
	mpf_init_set_str(one, "1", BASE);
	mpf_init_set_str(two, "2", BASE);
	mpf_init_set_str(three, "3", BASE);
	mpf_init_set_str(nine, "9", BASE);
	mpf_init_set_str(twenty_seven, "27", BASE);
	mpf_init(one_third);
	mpf_div(one_third, one, three);

//	multiple-precision variables
	mpf_t a, r, s, t, u, v, w, power3, prev_a;
	mpf_init(a);
	mpf_init(r);
	mpf_init(s);
	mpf_init(t);
	mpf_init(u);
	mpf_init(v);
	mpf_init(w);
	mpf_init(power3);
	mpf_init(prev_a);

//	temporaties
	mpf_t temp1, temp2;
	mpf_init(temp1);
	mpf_init(temp2);

//	initializitions

	//a = 1/3
	mpf_set(a, one_third);
	// r = (sqrt(3.0) - 1.0)/2.0;
	mpf_sqrt(temp1, three);
	mpf_sub(temp1, temp1, one);
	mpf_div(r, temp1, two);
	// s = cbrt(1 - r^3)
	mpf_mul(temp1, r,r);
	mpf_mul(temp1, temp1, r);
	mpf_sub(temp1, one, temp1);
	cube_root(s, temp1);
	// power3
	mpf_set(power3, one_third);

//	loop until two consecutive values are equal or up to MAX_ITERRATION times.
	int n = 0;

	do
	{
		std::chrono::time_point<std::chrono::system_clock> t1, t2;
		// start time count
		t1 = std::chrono::system_clock::now();
//		high_resolution_clock::time_point t1 = high_resolution_clock::now();

		// prev_a = a
		mpf_set(prev_a, a);
		// t = 1 + 2r
		mpf_add(temp1, r, r);
		mpf_add(t, temp1, one);
		// u = cbrt(9r(1 + r + r^2))
		mpf_mul(temp2, r, r);
		mpf_add(temp1, one, r);
		mpf_add(temp1, temp1, temp2);
		mpf_mul(temp1, temp1, nine);
		mpf_mul(temp1, temp1, r);
		cube_root(u, temp1);
		//  v = t^2 + tu + u^2
		mpf_mul(temp1, t, t);
		mpf_mul(temp2, t, u);
		mpf_add(temp1, temp1, temp2);
		mpf_mul(temp2, u, u);
		mpf_add(v, temp1, temp2);
		// w = (27(1 + s + s^2))/v
		mpf_add(temp1, one, s);
		mpf_mul(temp2, s, s);
		mpf_add(temp1, temp1, temp2);
		mpf_mul(temp1, temp1, twenty_seven);
		mpf_div(w, temp1, v);
		  // a = wa + (3^(2n-1))(1 - w)
		mpf_mul(temp1, w, a);
		mpf_sub(temp2, one, w);
		mpf_mul(temp2, power3, temp2);
		mpf_add(a, temp1, temp2);
		// s = ((1 - r)^3)/((t + 2u)v)
		mpf_sub(temp2, one, r);
		mpf_mul(temp1, temp2, temp2);
		mpf_mul(temp1, temp1, temp2);
		mpf_add(temp2, t, u);
		mpf_add(temp2, temp2, u);
		mpf_mul(temp2, temp2, v);
		mpf_div(s, temp1, temp2);
		// r = (1 - s^3)^(1/3)
		mpf_mul(temp1, s, s);
		mpf_mul(temp1, temp1, s);
		mpf_sub(temp1, one, temp1);
		cube_root(r, temp1);
		// power3 = 3^(2n-1)
		mpf_mul(power3, power3, nine);

		// end time count
		t2 = std::chrono::system_clock::now();

		// calculate elapsed time
		std::chrono::duration<double> time_span = t2 - t1;

		// print elapsed time
		std::cout << "Iteration " << n << " took " << time_span.count() << " seconds.";
		std::cout << std::endl;

		n++;
	} while (((n < 2) || (mpf_eq(a, prev_a, PRECISION) == 0)) && (n < MAX_ITERATIONS));

	// pi = 1/a
	mpf_div(pi, one, a);
}

void print(const mpf_t& pi)
{
	    mp_exp_t exp;  // exponent (not used)

	    // Convert the multiple-precision number x to a C string.
	    char *str = NULL;
	    char *s = mpf_get_str(str, &exp, BASE, PRECISION, pi);
	    char *p = s+1;  // skip the 3 before the decimal point

	    cout << endl;
	    cout << "3.";

	    char block[BLOCK_SIZE + 1];  // 1 extra for the ending \0
	    // Loop for each line.
	    for (int i = 1; i <= LINE_COUNT; i++)
	    {
	        // Loop to print blocks of digits in each line.
	        for (int j = 0; j < LINE_SIZE; j += BLOCK_SIZE)
	        {
	            strncpy(block, p+j, BLOCK_SIZE);
	            block[BLOCK_SIZE] = '\0';
	            cout << block << " ";
	        }

	        cout << endl << "  ";

	        // Print a blank line for grouping.
	        if (i%GROUP_SIZE == 0) cout << endl << "  ";

	        p += LINE_SIZE;
	    }

	    free(s);
	    cout << endl;

}
