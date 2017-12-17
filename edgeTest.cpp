#include "edgeTest.h"
using namespace std;

#ifndef FALSE
#define FALSE 0
#endif /* !FALSE */

#ifndef TRUE
#define TRUE 1
#endif /* !TRUE */

/* PI */
#ifndef M_PI
#define M_PI   3.14159265358979323846
#endif /* !M_PI */

/* fatal error, print a message to standard error and exit*/
void error(char * msg)
{
	fprintf(stderr, "error: %s\n", msg);
	exit(EXIT_FAILURE);
}

/* memory allocation, print an error and exit if fail*/
void * xmalloc(size_t size)
{
	void * p;
	if (size == 0) error("xmalloc: zero size");
	p = malloc(size);
	if (p == NULL) error("xmalloc: out of memory");
	return p;
}

/* compute a > b considering the rounding errors due to the representation of double numbers*/
int greater(double a, double b)
{
	if (a <= b) return FALSE;  /* trivial case, return as soon as possible */
	if ((a - b) < 1000 * DBL_EPSILON) return FALSE;
	return TRUE; /* greater */
}

/* Euclidean distance between x1,y1 and x2,y2*/
double dist(double x1, double y1, double x2, double y2)
{
	return sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
}

/* compute a Gaussian kernel of length n, standard deviation sigma,and centered at value mean.
The size of the kernel is selected to guarantee that the first discarded term is at least 10^prec times
smaller than the central value.For that,the half size of the kernel must be larger than x, with£º
e^(-x^2/2sigma^2) = 1/10^prec
Then,x = sigma * sqrt( 2 * prec * ln(10) ).

for example, if mean=0.5, the Gaussian will be centered in the middle point between values kernel[0] and kernel[1].
kernel must be allocated to a size n.
*/
void gaussian_kernel(double * pfGaussFilter, int n, double sigma, double mean)
{
	double sum = 0.0;
	double val;
	int i;
	if (pfGaussFilter == NULL) error("gaussian_kernel: kernel not allocated");
	if (sigma <= 0.0) error("gaussian_kernel: sigma must be positive");

	// compute Gaussian kernel
	for (i = 0; i<n; i++)
	{
		val = ((double)i - mean) / sigma;
		pfGaussFilter[i] = exp(-0.5 * val * val);
		sum += pfGaussFilter[i];
	}

	// normalization
	if (sum > 0.0) for (i = 0; i < n; i++)
		pfGaussFilter[i] /= sum;
}

/* filter an image with a Gaussian kernel of parameter sigma. return a pointer
to a newly allocated filtered image, of the same size as the input image. */
void gaussian_filter(uchar* image, uchar* out, int X, int Y, double sigma)
{
	int x, y, offset, i, j, nx2, ny2, n;
	double * kernel;
	double * tmp;
	double val, prec;

	if (sigma <= 0.0) error("gaussian_filter: sigma must be positive");
	if (image == NULL || X < 1 || Y < 1) error("gaussian_filter: invalid image");


	tmp = (double *)xmalloc(X * Y * sizeof(double));
	prec = 3.0;
	offset = (int)ceil(sigma * sqrt(2.0 * prec * log(10.0)));
	n = 1 + 2 * offset; //kernel size
	kernel = (double *)xmalloc(n * sizeof(double));
	gaussian_kernel(kernel, n, sigma, (double)offset);

	// auxiliary variables for the double of the image size
	nx2 = 2 * X;
	ny2 = 2 * Y;

	// x axis convolution
	for (x = 0; x<X; x++)
	for (y = 0; y<Y; y++)
	{
		val = 0.0;
		for (i = 0; i<n; i++)
		{
			j = x - offset + i;

			// symmetry boundary condition
			while (j<0) 
				j += nx2;
			while (j >= nx2) 
				j -= nx2;
			if (j >= X) 
				j = nx2 - 1 - j;

			val += (double)image[j + y*X] * kernel[i];
		}
		tmp[x + y*X] = val;
	}

	// y axis convolution
	for (x = 0; x<X; x++)
	for (y = 0; y<Y; y++)
	{
		val = 0.0;
		for (i = 0; i<n; i++)
		{
			j = y - offset + i;

			// symmetry boundary condition 
			while (j<0) 
				j += ny2;
			while (j >= ny2) 
				j -= ny2;
			if (j >= Y) 
				j = ny2 - 1 - j;

			val += tmp[x + j*X] * kernel[i];
		}
		out[x + y*X] = (uchar)val;
	}

	free((void *)kernel);
	free((void *)tmp);
}

/* return a score for chaining pixels 'from' to 'to', favoring closet point:
	= 0.0 invalid chaining;
	> 0.0 valid forward chaining; the larger the value, the better the chaining;
	< 0.0 valid backward chaining; the smaller the value, the better the chaining;

input:
	from, to: the two pixel IDs to evaluate their potential chaining
	Ex[i], Ey[i]: the sub-pixel position of point i, if i is an edge point;
	they take values -1,-1 if i is not an edge point;
	Gx[i], Gy[i]: the image gradient at pixel i;
	X, Y: the size of the image;
*/
double chain(int from, int to, double * Ex, double * Ey,double * Gx, double * Gy, int X, int Y)
{
	double dx, dy;
	if (Ex == NULL || Ey == NULL || Gx == NULL || Gy == NULL)
		error("chain: invalid input");
	if (from < 0 || to < 0 || from >= X*Y || to >= X*Y)
		error("chain: one of the points is out the image");

	//check that the points are different and valid edge points,otherwise return invalid chaining
	if (from == to) 
		return 0.0; // same pixel, not a valid chaining
	if (Ex[from] < 0.0 || Ey[from] < 0.0 || Ex[to] < 0.0 || Ey[to] < 0.0)
		return 0.0; // one of them is not an edge point, not a valid chaining

	/* in a good chaining, the gradient should be roughly orthogonal
	to the line joining the two points to be chained:
	when Gy * dx - Gx * dy > 0, it corresponds to a forward chaining,
	when Gy * dx - Gx * dy < 0, it corresponds to a backward chaining.

	first check that the gradient at both points to be chained agree
	in one direction, otherwise return invalid chaining. */
	dx = Ex[to] - Ex[from];
	dy = Ey[to] - Ey[from];
	if ((Gy[from] * dx - Gx[from] * dy) * (Gy[to] * dx - Gx[to] * dy) <= 0.0)
		return 0.0; /* incompatible gradient angles, not a valid chaining */

	/* return the chaining score: positive for forward chaining,negative for backwards. 
	the score is the inverse of the distance to the chaining point, to give preference to closer points */
	if ((Gy[from] * dx - Gx[from] * dy) >= 0.0)
		return  1.0 / dist(Ex[from], Ey[from], Ex[to], Ey[to]); /* forward chaining  */
	else
		return -1.0 / dist(Ex[from], Ey[from], Ex[to], Ey[to]); /* backward chaining */
}

/* compute the image gradient, giving its x and y components as well as the
modulus. Gx, Gy, and modG must be already allocated. */
void compute_gradient(double * Gx, double * Gy, double * modG,uchar * image, int X, int Y)
{
	int x, y;

	if (Gx == NULL || Gy == NULL || modG == NULL || image == NULL)
		error("compute_gradient: invalid input");

	// approximate image gradient using centered differences
	for (x = 1; x<(X - 1); x++)
	for (y = 1; y<(Y - 1); y++)
	{
		Gx[x + y*X] = (double)image[(x + 1) + y*X] - (double)image[(x - 1) + y*X];
		Gy[x + y*X] = (double)image[x + (y + 1)*X] - (double)image[x + (y - 1)*X];
		modG[x + y*X] = sqrt(Gx[x + y*X] * Gx[x + y*X] + Gy[x + y*X] * Gy[x + y*X]);
	}
}

/* compute sub-pixel edge points using adapted Canny and Devernay methods.
input: 
	Gx, Gy, and modG are the x and y components and modulus of the image gradient. X,Y is the image size.
output: 
	Ex and Ey :will have the x and y sub-pixel coordinates of the edge points found, 
	or -1 and -1 when not an edge point. Ex and Ey must be already allocated.

Devernay correction works very well and is very consistent when used to interpolate only along horizontal
or vertical direction. this modified version requires that a pixel, to be an edge point, must be a local maximum
horizontally or vertically, depending on the gradient orientation: if the x component of the gradient is larger 
than the y component, Gx > Gy, this means that the gradient is roughly horizontal and a horizontal maximum is
required to be an edge point. */
void compute_edge_points(double * Ex, double * Ey, double * modG,
	double * Gx, double * Gy, int X, int Y)
{
	int x, y, i;

	if (Ex == NULL || Ey == NULL || modG == NULL || Gx == NULL || Gy == NULL)
		error("compute_edge_points: invalid input");
	/* initialize Ex and Ey as non-edge points for all pixels */
	for (i = 0; i<X*Y; i++) Ex[i] = Ey[i] = -1.0;
	/* explore pixels inside a 2 pixel margin (so modG[x,y +/- 1,1] is defined) */
	for (x = 2; x<(X - 2); x++)
	for (y = 2; y<(Y - 2); y++)
	{
		int Dx = 0;                     /* interpolation is along Dx,Dy		*/
		int Dy = 0;                     /* which will be selected below		*/
		double mod = modG[x + y*X];     /* modG at pixel					*/
		double L = modG[x - 1 + y*X];   /* modG at pixel on the left			*/
		double R = modG[x + 1 + y*X];   /* modG at pixel on the right		*/
		double U = modG[x + (y + 1)*X]; /* modG at pixel up					*/
		double D = modG[x + (y - 1)*X]; /* modG at pixel below				*/
		double gx = fabs(Gx[x + y*X]);  /* absolute value of Gx				*/
		double gy = fabs(Gy[x + y*X]);  /* absolute value of Gy				*/
		/* when local horizontal maxima of the gradient modulus and the gradient direction 
		is more horizontal (|Gx| >= |Gy|),=> a "horizontal" (H) edge found else, 
		if local vertical maxima of the gradient modulus and the gradient direction is more
		vertical (|Gx| <= |Gy|),=> a "vertical" (V) edge found */

		/* it can happen that two neighbor pixels have equal value and are both	maxima, for example 
		when the edge is exactly between both pixels. in such cases, as an arbitrary convention, 
		the edge is marked on the left one when an horizontal max or below when a vertical max.
		for	this the conditions are L < mod >= R and D < mod >= U,respectively. the comparisons are
		done using the function greater() instead of the operators > or >= so numbers differing only
		due to rounding errors are considered equal */
		if (greater(mod, L) && !greater(R, mod) && gx >= gy) Dx = 1; /* H */
		else if (greater(mod, D) && !greater(U, mod) && gx <= gy) Dy = 1; /* V */

		/* Devernay sub-pixel correction

		the edge point position is selected as the one of the maximum of a quadratic interpolation of the magnitude of
		the gradient along a unidimensional direction. the pixel must be a local maximum. so we	have the values:
		
		the x position of the maximum of the parabola passing through(-1,a), (0,b), and (1,c) is
		offset = (a - c) / 2(a - 2b + c),and because b >= a and b >= c, -0.5 <= offset <= 0.5	*/
		if (Dx > 0 || Dy > 0)
		{
			/* offset value is in [-0.5, 0.5] */
			double a = modG[x - Dx + (y - Dy) * X];
			double b = modG[x + y     * X];
			double c = modG[x + Dx + (y + Dy) * X];
			double offset = 0.5 * (a - c) / (a - b - b + c);

			/* store edge point */
			Ex[x + y*X] = x + offset * Dx;
			Ey[x + y*X] = y + offset * Dy;
		}
	}
}

/* chain edge points
input: 
Ex/Ey:the sub-pixel coordinates when an edge point is present or -1,-1 otherwise. 
Gx/Gy/modG:the x and y components and the modulus of the image gradient. X,Y is the image size.

output: 
next and prev:contain the number of next and previous edge points in the chain. 
when not chained in one of the directions, the corresponding value is set to -1. 
next and prev must be allocated before calling.*/
void chain_edge_points(int * next, int * prev, double * Ex,	double * Ey,double * Gx, double * Gy, int X, int Y)
{
	int x, y, i, j, alt;

	if (next == NULL || prev == NULL || Ex == NULL || Ey == NULL || Gx == NULL || Gy == NULL)
		error("chain_edge_points: invalid input");

	/* initialize next and prev as non linked */
	for (i = 0; i<X*Y; i++) next[i] = prev[i] = -1;

	/* try each point to make local chains */
	for (x = 2; x<(X - 2); x++)   /* 2 pixel margin to include the tested neighbors */
	for (y = 2; y<(Y - 2); y++)
	if (Ex[x + y*X] >= 0.0 && Ey[x + y*X] >= 0.0) /* must be an edge point */
	{
		int from = x + y*X;  /* edge point to be chained			*/
		double fwd_s = 0.0;  /* score of best forward chaining		*/
		double bck_s = 0.0;  /* score of best backward chaining		*/
		int fwd = -1;        /* edge point of best forward chaining */
		int bck = -1;        /* edge point of best backward chaining*/

		/* try all neighbors two pixels apart or less.
		looking for candidates for chaining two pixels apart, in most such cases, 
		is enough to obtain good chains of edge points that	accurately describes the edge.	*/
		for (i = -2; i <= 2; i++)
		for (j = -2; j <= 2; j++)
		{
			int to = x + i + (y + j)*X; /* candidate edge point to be chained */
			double s = chain(from, to, Ex, Ey, Gx, Gy, X, Y);  /* score from-to */

			if (s > fwd_s)	/* a better forward chaining found    */
			{
				fwd_s = s;  /* set the new best forward chaining  */
				fwd = to;
			}
			if (s < bck_s)	/* a better backward chaining found	  */
			{
				bck_s = s;  /* set the new best backward chaining */
				bck = to;
			}
		}

		/* before making the new chain, check whether the target was
		already chained and in that case, whether the alternative
		chaining is better than the proposed one.

		x alt                        x alt
		\                          /
		\                        /
		from x---------x fwd              bck x---------x from

		we know that the best forward chain starting at from is from-fwd.
		but it is possible that there is an alternative chaining arriving
		at fwd that is better, such that alt-fwd is to be preferred to
		from-fwd. an analogous situation is possible in backward chaining,
		where an alternative link bck-alt may be better than bck-from.

		before making the new link, check if fwd/bck are already chained,
		and in such case compare the scores of the proposed chaining to
		the existing one, and keep only the best of the two.

		there is an undesirable aspect of this procedure: the result may
		depend on the order of exploration. consider the following
		configuration:

		a x-------x b
		/
		/
		c x---x d    with score(a-b) < score(c-b) < score(c-d)
		or equivalently ||a-b|| > ||b-c|| > ||c-d||

		let us consider two possible orders of exploration.

		order: a,b,c
		we will first chain a-b when exploring a. when analyzing the
		backward links of b, we will prefer c-b, and a-b will be unlinked.
		finally, when exploring c, c-d will be preferred and c-b will be
		unlinked. the result is just the chaining c-d.

		order: c,b,a
		we will first chain c-d when exploring c. then, when exploring
		the backward connections of b, c-b will be the preferred link;
		but because c-d exists already and has a better score, c-b
		cannot be linked. finally, when exploring a, the link a-b will
		be created because there is no better backward linking of b.
		the result is two chainings: c-d and a-b.

		we did not found yet a simple algorithm to solve this problem. by
		simple, we mean an algorithm without two passes or the need to
		re-evaluate the chaining of points where one link is cut.

		for most edge points, there is only one possible chaining and this
		problem does not arise. but it does happen and a better solution
		is desirable.
		*/
		if (fwd >= 0 && next[from] != fwd &&
			((alt = prev[fwd]) < 0 || chain(alt, fwd, Ex, Ey, Gx, Gy, X, Y) < fwd_s))
		{
			if (next[from] >= 0)     /* remove previous from-x link if one */
				prev[next[from]] = -1;  /* only prev requires explicit reset  */
			next[from] = fwd;         /* set next of from-fwd link          */
			if (alt >= 0)            /* remove alt-fwd link if one         */
				next[alt] = -1;         /* only next requires explicit reset  */
			prev[fwd] = from;         /* set prev of from-fwd link          */
		}
		if (bck >= 0 && prev[from] != bck &&
			((alt = next[bck]) < 0 || chain(alt, bck, Ex, Ey, Gx, Gy, X, Y) > bck_s))
		{
			if (alt >= 0)            /* remove bck-alt link if one         */
				prev[alt] = -1;         /* only prev requires explicit reset  */
			next[bck] = from;         /* set next of bck-from link          */
			if (prev[from] >= 0)     /* remove previous x-from link if one */
				next[prev[from]] = -1;  /* only next requires explicit reset  */
			prev[from] = bck;         /* set prev of bck-from link          */
		}
	}
}

/* apply Canny thresholding with hysteresis

next and prev contain the number of next and previous edge points in the
chain or -1 when not chained. modG is modulus of the image gradient. X,Y is
the image size. th_h and th_l are the high and low thresholds, respectively.

this function modifies next and prev, removing chains not satisfying the
thresholds.
*/
void thresholds_with_hysteresis(int * next, int * prev,	
	double * modG,	int X, int Y,double th_h, double th_l)
{
	int * valid;
	int i, j, k;

	/* check input */
	if (next == NULL || prev == NULL || modG == NULL)
		error("thresholds_with_hysteresis: invalid input");

	/* get memory */
	valid = (int *)xmalloc(X * Y * sizeof(int));
	for (i = 0; i<X*Y; i++) valid[i] = FALSE;

	/* validate all edge points over th_h or connected to them and over th_l */
	for (i = 0; i<X*Y; i++)   /* prev[i]>=0 or next[i]>=0 implies an edge point */
	if ((prev[i] >= 0 || next[i] >= 0) && !valid[i] && modG[i] >= th_h)
	{
		valid[i] = TRUE; /* mark as valid the new point */

		/* follow the chain of edge points forwards */
		for (j = i; j >= 0 && (k = next[j]) >= 0 && !valid[k]; j = next[j])
		if (modG[k] < th_l)
		{
			next[j] = -1;  /* cut the chain when the point is below th_l */
			prev[k] = -1;  /* j must be assigned to next[j] and not k,
						   so the loop is chained in this case */
		}
		else
			valid[k] = TRUE; /* otherwise mark the new point as valid */

		/* follow the chain of edge points backwards */
		for (j = i; j >= 0 && (k = prev[j]) >= 0 && !valid[k]; j = prev[j])
		if (modG[k] < th_l)
		{
			prev[j] = -1;  /* cut the chain when the point is below th_l */
			next[k] = -1;  /* j must be assigned to prev[j] and not k,
						   so the loop is chained in this case */
		}
		else
			valid[k] = TRUE; /* otherwise mark the new point as valid */
	}

	/* remove any remaining non-valid chained point */
	for (i = 0; i<X*Y; i++)   /* prev[i]>=0 or next[i]>=0 implies edge point */
	if ((prev[i] >= 0 || next[i] >= 0) && !valid[i])
		prev[i] = next[i] = -1;

	/* free memory */
	free((void *)valid);
}

/* create a list of chained edge points composed of 3 lists
x, y and curve_limits; it also computes N (the number of edge points) and
M (the number of curves).

x[i] and y[i] (0<=i<N) store the sub-pixel coordinates of the edge points.
curve_limits[j] (0<=j<=M) stores the limits of each chain in lists x and y.

x, y, and curve_limits will be allocated.

example:

curve number k (0<=k<M) consists of the edge points x[i],y[i]
for i determined by curve_limits[k] <= i < curve_limits[k+1].

curve k is closed if x[curve_limits[k]] == x[curve_limits[k+1] - 1] and
y[curve_limits[k]] == y[curve_limits[k+1] - 1].
*/
void list_chained_edge_points(double ** x, double ** y, int * N,int ** curve_limits,
	int * M,int * next, int * prev,	double * Ex, double * Ey, int X, int Y)
{
	int i, k, n;

	/* initialize output: x, y, curve_limits, N, and M

	there cannot be more than X*Y edge points to be put in the output list:
	edge points must be local maxima of gradient modulus, so at most half of
	the pixels could be so. when a closed curve is found, one edge point will
	be put twice to the output. even if all possible edge points (half of the
	pixels in the image) would form one pixel closed curves (which is not
	possible) that would lead to output X*Y edge points.

	for the same reason, there cannot be more than X*Y curves: the worst case
	is when all possible edge points (half of the pixels in the image) would
	form one pixel chains. in that case (which is not possible) one would need
	a size for curve_limits of X*Y/2+1. so X*Y is enough.

	(curve_limits requires one more item than the number of curves.
	a simplest example is when only one chain of length 3 is present:
	curve_limits[0] = 0, curve_limits[1] = 3.)
	*/
	*x = (double *)xmalloc(X * Y * sizeof(double));
	*y = (double *)xmalloc(X * Y * sizeof(double));
	*curve_limits = (int *)xmalloc(X * Y * sizeof(int));
	*N = 0;
	*M = 0;

	/* copy chained edge points to output */
	for (i = 0; i<X*Y; i++)   /* prev[i]>=0 or next[i]>=0 implies an edge point */
	if (prev[i] >= 0 || next[i] >= 0)
	{
		/* a new chain found, set chain starting index to the current point
		and then increase the curve counter */
		(*curve_limits)[*M] = *N;
		++(*M);

		/* set k to the beginning of the chain, or to i if closed curve */
		for (k = i; (n = prev[k]) >= 0 && n != i; k = n);

		/* follow the chain of edge points starting on k */
		do
		{
			/* store the current point coordinates in the output lists */
			(*x)[*N] = Ex[k];
			(*y)[*N] = Ey[k];
			++(*N);

			n = next[k];   /* save the id of the next point in the chain */

			next[k] = -1;  /* unlink chains from k so it is not used again */
			prev[k] = -1;

			/* for closed curves, the initial point is included again as
			the last point of the chain. actually, testing if the first
			and last points are equal is the only way to know that it is
			a closed curve.

			to understand that this code actually repeats the first point,
			consider a closed chain as follows:  a--b
			|  |
			d--c

			let us say that the algorithm starts by point a. it will store
			the coordinates of point a and then unlink a-b. then, will store
			point b and unlink b-c, and so on. but the link d-a is still
			there. (point a is no longer pointing backwards to d, because
			both links are removed at each step. but d is indeed still
			pointing to a.) so it will arrive at point a again and store its
			coordinates again as last point. there, it cannot continue
			because the link a-b was removed, there would be no next point,
			k would be -1 and the curve is finished.
			*/

			k = n;  /* set the current point to the next in the chain */
		} while (k >= 0); /* continue while there is a next point in the chain */
	}
	(*curve_limits)[*M] = *N; /* store end of the last chain */
}

/* chained, sub-pixel edge detector. based on a modified Canny non-maximal
suppression and a modified Devernay sub-pixel correction.

input:

image        : the input image
X,Y          : the size of the input image
sigma        : standard deviation sigma for the Gaussian filtering
(if sigma=0 no filtering is performed)
th_h         : high gradient threshold in Canny's hysteresis
th_l         : low gradient threshold in Canny's hysteresis

output:

x,y          : lists of sub-pixel coordinates of edge points
curve_limits : the limits of each curve in lists x and y
N            : number of edge points
M            : number of curves

the input is a XxY graylevel image given as a pointer to an array of doubles
such that image[x+y*X] is the value at coordinates x,y
(for 0 <= x < X and 0 <= y < Y).

the output are the chained edge points given as 3 allocated lists: x, y and
curve_limits. also the numbers N (size of lists x and y) and M (number of
curves).

x[i] and y[i] (0<=i<N) store the sub-pixel coordinates of the edge points.
curve_limits[j] (0<=j<=M) stores the limits of each chain in lists x and y.

example:

curve number k (0<=k<M) consists of the edge points x[i],y[i]
for i determined by curve_limits[k] <= i < curve_limits[k+1].

curve k is closed if x[curve_limits[k]] == x[curve_limits[k+1] - 1] and
y[curve_limits[k]] == y[curve_limits[k+1] - 1].
*/
void devernay(double ** x, double ** y, int * N, int ** curve_limits,int * M,
	uchar * image, uchar * gauss, int X, int Y, double sigma, double th_h, double th_l)
{
	double * Gx = (double *)xmalloc(X * Y * sizeof(double));     /* grad_x */
	double * Gy = (double *)xmalloc(X * Y * sizeof(double));     /* grad_y */
	double * modG = (double *)xmalloc(X * Y * sizeof(double));   /* |grad| */
	double * Ex = (double *)xmalloc(X * Y * sizeof(double));     /* edge_x */
	double * Ey = (double *)xmalloc(X * Y * sizeof(double));     /* edge_y */
	int * next = (int *)xmalloc(X * Y * sizeof(int));			 /* next point in chain */
	int * prev = (int *)xmalloc(X * Y * sizeof(int));			 /* prev point in chain */

	if (sigma == 0.0) compute_gradient(Gx, Gy, modG, image, X, Y);
	else
	{
		gaussian_filter(image, gauss, X, Y, sigma);
		compute_gradient(Gx, Gy, modG, gauss, X, Y);
	}

	compute_edge_points(Ex, Ey, modG, Gx, Gy, X, Y);

	chain_edge_points(next, prev, Ex, Ey, Gx, Gy, X, Y);

	thresholds_with_hysteresis(next, prev, modG, X, Y, th_h, th_l);

	list_chained_edge_points(x, y, N, curve_limits, M, next, prev, Ex, Ey, X, Y);

	/* free memory */
	free((void *)Gx);
	free((void *)Gy);
	free((void *)modG);
	free((void *)Ex);
	free((void *)Ey);
	free((void *)next);
	free((void *)prev);
}




////---------²Î¿¼----------
//void GetGuassFilter(double* pfGaussFilter, const int iFilterHeight, const int iFilterWidth, const double fSigma)
//{
//	for (int i = 0; i < iFilterHeight; i++)
//	{
//		for (int j = 0; j < iFilterWidth; j++)
//		{
//			pfGaussFilter[i*iFilterHeight + j] = 1.0 / exp(((i - iFilterHeight / 2) ^ 2 + (j - iFilterWidth) ^ 2));
//		}
//		cout << endl;
//	}
//	double sum = 0;
//	for (int j = 0; j < iFilterWidth*iFilterHeight; j++)
//	{
//		sum += pfGaussFilter[j];
//	}
//	for (int j = 0; j < iFilterHeight*iFilterWidth; j++)
//	{
//		pfGaussFilter[j] = pfGaussFilter[j] / sum;
//	}
//}
//void GaussFilt2D(uchar* pSrc, uchar* pDst, const int iHeight, const int iWidth,
//	const int iFilterHeight, const int iFilterWidth, const double fSigma)
//{
//	double* pfGaussFilter = new double[iFilterHeight*iFilterWidth];
//	GetGuassFilter(pfGaussFilter, iFilterHeight, iFilterWidth, fSigma);
//
//	double sum = 0;
//	int index = 0;
//	double current_src = 0;
//	for (int i = -iFilterHeight; i < iHeight - 1; i++)
//	for (int j = -iFilterWidth; j < iWidth - 1; j++)
//	{
//		sum = 0;
//		index = 0;
//
//		for (int m = i - iFilterHeight / 2; m < i + iFilterHeight / 2 + 1; m++)
//		{
//			for (int n = j - iFilterWidth / 2; n<j + iFilterWidth / 2 + 1; n++)
//			{
//				if (m<0 || n<0 || m>iHeight - 1 || n>iWidth - 1)
//					continue;
//				current_src = (double)pSrc[m*iWidth + n];
//				sum += current_src*pfGaussFilter[index++];
//
//			}
//		}
//		if (i + iFilterHeight / 2>0 && i + iFilterHeight / 2<iHeight - 1 && j + iFilterWidth / 2 > 0 && j + iFilterWidth / 2 < iWidth - 1)
//			pDst[(i + iFilterHeight / 2)*iWidth + j + iFilterWidth / 2] = (uchar)sum;
//
//		if (pDst[i*iWidth + j]<0)
//			pDst[i*iWidth + j] = 0;
//		if (pDst[i*iWidth + j]>255)
//			pDst[i*iWidth + j] = 255;
//	}
//}
//void GaussFilt2DFull(uchar* pSrc, uchar* pDst, const int iHeight, const int iWidth,
//	const int iFilterHeight, const int iFilterWidth, const double fSigma)
//{
//	double* pfGaussFilter = new double[iFilterHeight*iFilterWidth];
//	GetGuassFilter(pfGaussFilter, iFilterHeight, iFilterWidth, fSigma);
//	double sum = 0;
//	int index = 0;
//	double current_src = 0;
//	for (int i = 0; i < iHeight - 1; i++)
//	for (int j = 0; j < iWidth - 1; j++)
//	{
//		sum = 0;
//		index = 0;
//
//		for (int m = i - iFilterHeight / 2; m < i + iFilterHeight / 2; m++)
//		for (int n = j - iFilterWidth / 2; n < j + iFilterWidth / 2; n++)
//		{
//			if (m<0 || n<0 || m>iHeight - 1 || n>iWidth - 1)
//				continue;
//			current_src = (double)pSrc[m*iWidth + n];
//			sum += current_src*pfGaussFilter[index++];
//
//		}
//		//if(i+iFilterHeight/2>0&&i+iFilterHeight/2<iHeight-1&&j+iFilterWidth/2>0&&j+iFilterWidth/2<iWidth-1)
//		pDst[i*iWidth + j] = (uchar)sum;
//
//		//pDst[i*iWidth+j]=(uchar)sum;
//		if (pDst[i*iWidth + j]<0)
//			pDst[i*iWidth + j] = 0;
//		if (pDst[i*iWidth + j]>255)
//			pDst[i*iWidth + j] = 255;
//	}
//}