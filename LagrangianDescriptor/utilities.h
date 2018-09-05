/*Class for using some utilities functions*/

#pragma once
#include <stdlib.h>				//size_t
#include <math.h>				//sqrt
#include <vector>				//vector	

//Constant definitions
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0 / IM1)
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1 + IMM1 / NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)

class utilities
{
public:
	utilities();
	double powerd(double, int);
	float RandomFloat(float, float);
	float ran2(long*);
	double** create_matrix(size_t, size_t);
	std::vector<double> box_muller(int nstep, long *ran);
	double absol(double val);
	double squareRt(double x);

};

///<sumary>	
/// Utilities functions.
///</sumary>
utilities::utilities(){}

///<sumary>	
/// Random float betwen 2 Floats
///</sumary>
float utilities::RandomFloat(float a, float b)
{
	float random = ((float)rand()) / (float)RAND_MAX;
	float diff = b - a;
	float r = random * diff;
	return a + r;
}

///<sumary>
/// Long period(> 2 × 10 18) random number generator of L’Ecuyer with Bays - Durham shuffle
/// and added safeguards.Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
/// the endpoint values).Call with idum a negative integer to initialize; thereafter, do not alter
///	idum between successive deviates in a sequence.RNMX should approximate the largest floating
///	value that is less than 1.
///</sumary>
float utilities::ran2(long *idum)
{
	int j;
	long k;
	static long idum2 = 123456789;
	static long iy = 0;
	static long iv[NTAB];
	double temp;
	if (*idum <= 0) {								//Initialize.
		if (-(*idum) < 1) *idum = 1;				//Be sure to prevent idum = 0.
		else *idum = -(*idum);
		idum2 = (*idum);
		for (j = NTAB + 7; j >= 0; j--) {			//Load the shuffle table(after 8 warm - ups).
			k = (*idum) / IQ1;
			*idum = IA1 * (*idum - k * IQ1) - k * IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum) / IQ1;								//Start here when not initializing.
	*idum = IA1 * (*idum - k * IQ1) - k * IR1;			//Compute idum = (IA1*idum) % IM1 without
	if (*idum < 0) *idum += IM1;					//overflows by Schrage’s method.
	k = idum2 / IQ2;
	idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;			//Compute idum2 = (IA2*idum) % IM2 likewise.
	if (idum2 < 0) idum2 += IM2;
	j = iy / NDIV;									//Will be in the range 0..NTAB - 1.
	iy = iv[j] - idum2;								//Here idum is shuffled, idum and idum2 are
	iv[j] = *idum;									//	combined to generate output.
	if (iy < 1) iy += IMM1;							//Because users don’t expect endpoint values.
	if ((temp = AM * iy) > RNMX) return RNMX;
	else return temp;
}

///<sumary>	
/// Extended version of power function that can work
///for double base and negative exponent.
///</sumary>
double utilities::powerd(double base, int exponent)
{
	double temp;		//temporal variable
	if (exponent == 0)
		return 1;
	temp = powerd(base, exponent / 2);
	if ((exponent % 2) == 0) {
		return temp * temp;
	}
	else {
		if (exponent > 0)
			return base * temp * temp;
		else
			return (temp * temp) / base;
	}
}

///<sumary>	
/// Creates a Matrix of size = size1 X size2.
///</sumary>
double** utilities::create_matrix(size_t size1, size_t size2)
{
	double** m = new double*[size1];
	for (size_t i = 0; i < size1; ++i)
	{
		m[i] = new double[size2];
	}
	return m;
}

///<sumary>	
/// Box Muller Random Numbers
/// The Box-Muller transform takes two random variables,
/// evenly distributed in the interval (0,1) and transforms them to two independent deviates,
/// which are sampled from a Gaussian distribution. 
///</sumary>
std::vector<double> utilities::box_muller(int nstep, long *ran)
{
	double u1, u2;		//Two random Varialbes
	double term;
	std::vector<double> gasdev(nstep, 0.);
	
	for (size_t i = 0; i < nstep; i += 2)
	{
		u1 = ((double)rand() / (RAND_MAX));
		u2 = ((double)rand() / (RAND_MAX));

		u1 = ran2(ran);
		u2 = ran2(ran + 1);

		term = sqrt(-2 * log(u1));
		gasdev[i] = term * cos(2 * M_PIl * u2);
		gasdev[i + 1] = term * sin(2 * M_PIl * u2);
	}
	return gasdev;
}

///<sumary>	
/// Returns the absolute value of a double number.
///</sumary>
double utilities::absol(double value)
{
	if (value < 0.) { return value * -1; }
	else { return value; }
}

double utilities::squareRt(double x)
{
	if (x <= 0)
		return 0;       // if negative number throw an exception?
	int exp = 0;
	x = frexp(x, &exp); // extract binary exponent from x
	if (exp & 1) {      // we want exponent to be even
		exp--;
		x *= 2;
	}
	double y = (1 + x) / 2; // first approximation
	double z = 0;
	while (y != z) {    // yes, we CAN compare doubles here!
		z = y;
		y = (y + x / y) / 2;
	}
	return ldexp(y, exp / 2); // multiply answer by 2^(exp/2)
}
