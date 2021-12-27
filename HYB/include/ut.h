#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define XWIDTH 256
#define YWIDTH 256
#define MAXVAL 65535

#if ((0x100 & 0xf) == 0x0)
#	define I_M_LITTLE_ENDIAN 1
	#define swap(mem) (( (mem) & (short int)0xff00) >> 8) +	( ((mem) & (short int)0x00ff) << 8)
#else
	#define I_M_LITTLE_ENDIAN 0
	#define swap(mem) (mem)
#endif

#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))

#ifndef KTYPE
	#define KTYPE float
#endif

//professors routines for pgm file management 
void write_pgm_image( void *image, int maxval, int xsize, int ysize, const char *image_name);
void read_pgm_image( void **image, int *maxval, int *xsize, int *ysize, const char *image_name);
void swap_image( void *image, int xsize, int ysize, int maxval );
void * generate_gradient( int maxval, int xsize, int ysize );
void * generate_random( int maxval, int xsize, int ysize );


//my routines

//kernels
void uniform_kernel(KTYPE *mat, size_t x, size_t y);
void weighted_kernel(KTYPE *mat, size_t x, size_t y, KTYPE f);
void gaussian_kernel(KTYPE *mat, size_t x, size_t y, KTYPE sigma_sq);
KTYPE *normalize(void *kimage,size_t xkernel,size_t ykernel,int maxval);

//Convolution
void OMP_swap_image( void *image, int xsize, int ysize, int maxval );
void OMP_MPIConvolve(unsigned short int *image, unsigned short int *blurred, int xsize, int ysize,KTYPE *convolution_matrix, int xconv, int yconv, int space_up, int space_down);





