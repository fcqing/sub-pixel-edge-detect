#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <opencv2/opencv.hpp>    
#include <opencv2/imgproc/imgproc.hpp>    
#include <opencv2/core.hpp>  
#include <opencv2/highgui.hpp>  

void error(char * msg);
void * xmalloc(size_t size);
int greater(double a, double b);
double dist(double x1, double y1, double x2, double y2);
void gaussian_kernel(double * kernel, int n, double sigma, double mean);
void gaussian_filter(uchar* image, uchar* out, int iHeight, int iWidth, double sigma);
double chain(int from, int to, double * Ex, double * Ey, double * Gx, double * Gy, int X, int Y);
void compute_gradient(double * Gx, double * Gy, double * modG, uchar * image, int X, int Y);
void compute_edge_points(double * Ex, double * Ey, double * modG,
	double * Gx, double * Gy, int X, int Y);
void chain_edge_points(int * next, int * prev, double * Ex,
	double * Ey, double * Gx, double * Gy, int X, int Y);
void thresholds_with_hysteresis(int * next, int * prev,
	double * modG, int X, int Y, double th_h, double th_l);
void list_chained_edge_points(double ** x, double ** y, int * N,int ** curve_limits, 
	int * M,int * next, int * prev,	double * Ex, double * Ey, int X, int Y);
void devernay(double ** x, double ** y, int * N, int ** curve_limits, int * M,
	uchar * image, uchar * gauss, int X, int Y, double sigma, double th_h, double th_l);

////---------
//void GetGuassFilter(double* pfGaussFilter, const int iFilterHeight, const int iFilterWidth, const double fSigma);
//void GaussFilt2D(uchar* pSrc, uchar* pDst, const int iHeight, const int iWidth,
//	const int iFilterHeight, const int iFilterWidth, const double fSigma);
//void GaussFilt2DFull(uchar* pSrc, uchar* pDst, const int iHeight, const int iWidth,
//	const int iFilterHeight, const int iFilterWidth, const double fSigma);