#include "ut.h"
#include <string.h>
#include <time.h>
#include <omp.h>

int main(int args, char** argv)
{

	clock_t t0, tIO, tcalc, twrite;
	t0=clock();

	/*********************************************************************************************************************
	*																																																										 *
	*																			I/0 MANAGEMENT - INITIALIZATION																								 *
	*																																																										 *
	**********************************************************************************************************************/
	char* usage = "Usage: ./blur [kernel-type] {x-kernel-size} {y-kernel-size} {additional-kernel-param} [input-file] {output-file}";
	#define MAX_ARGS 7
	#define MIN_ARGS 4
	
	if(args<MIN_ARGS) 
	{
		printf("Too few arguments in function blur. %s\n",usage);
		return 1;
	}
	else if(args>MAX_ARGS)
	{
		printf("Too many arguments in function blur. %s\n",usage);
		return 1;
	}
	
	int arg_counter = 0;
	int xkernel,ykernel;
	KTYPE f, *kernel;

	int ktype = atoi(argv[++arg_counter]);
	  
	/************************************************************************************	
	 kernel choice: 0->"uniform" 1->"weighted" 2->"gaussian" 3->"imported from pgm file"
	*************************************************************************************/
  
	switch(ktype)
	{	
		case 0: case 1: case 2:

			xkernel = atoi(argv[++arg_counter]);
			ykernel = atoi(argv[++arg_counter]);
			if(!(xkernel%2) || !(ykernel%2))
			{
				printf("Kernel dimensions must be odd integers. Dimensions given were %dx%d\n",xkernel,ykernel);
				return 3;
			}
					
			kernel = (KTYPE *)malloc(xkernel*ykernel*sizeof(KTYPE));

			switch(ktype)
			{
				case 0:
					if(args > MAX_ARGS-1)
					{
						printf("Too many arguments in function blur. %s\n",usage);
						free(kernel);
						return 1;
					}
					uniform_kernel(kernel,xkernel,ykernel);
					printf("Using Mean Kernel of dimension %dx%d\n",xkernel,ykernel);
				break;
					
				case 1:
					if(args < MAX_ARGS-1)
					{
						printf("Too few arguments in function blur, additional parameter for the kernel needed.\n%s\n",usage);
						free(kernel);
						return 1;
					}
									
					f = atof(argv[++arg_counter]);
					weighted_kernel(kernel,xkernel,ykernel,f);
					printf("Using Weighted Kernel of dimension %dx%d and f = %f\n",xkernel,ykernel,f);
				break;
					
				case 2:
					if(args > MAX_ARGS-1)
					{
						printf("Too many arguments in function blur. %s\n",usage);
						free(kernel);				
						return 1;
					}
								
					// sigma for the gaussian function ---> here it has been (arbitrarily) chosen to be the maximum half dimension of the kernel			
					int s = max((xkernel/2),(ykernel/2));
					gaussian_kernel(kernel,xkernel,ykernel,s*s);
					printf("Using Gaussian Kernel of dimension %dx%d (s: %d)\n",xkernel,ykernel,s);
				break;
			}
		break;

		case 3:
			if(args > MAX_ARGS-2)
			{
				printf("Too many arguments in function blur. No dimension needed for pgm kernel.\n %s\n",usage);
				return 1;
			}
			
			int kmaxval;
			void *kimage;
			
			read_pgm_image(&kimage, &kmaxval, &xkernel, &ykernel, argv[++arg_counter]);
			kernel = normalize(kimage,xkernel,ykernel,kmaxval);
			free(kimage);

		break;

		default:
			printf("Kernel id must be either 0,1,2 for automatic generation or 3 for pgm file. Given id was %d",ktype);
			return 2;
	}
	
	/********************
	 input name setting 
	********************/
	
	char* input_name = argv[++arg_counter];
	
	int namelen = strlen(input_name);
	if(namelen > 4 && !strcmp(&input_name[namelen-4],".pgm")) { printf("Trying to blur image \"%s\"\n",input_name); }
	else
	{
		printf("Input file name must end in \".pgm\". Given file name was %s\n",input_name);
		free(kernel);
		return 4; //exit
	}
	
	
	/**********************************************************************************************************************
	*																																																										  *
	*																			  BLURRING - OMP COMMUNICATION																								  *
	*																																																										  *
	**********************************************************************************************************************/

	void *image; //pointer to be used to store the image 
	int maxval, xsize, ysize; //useful info about image
	
	read_pgm_image(&image, &maxval, &xsize, &ysize, input_name);
	
	//the following code is structured to work with 16bit images, so it wont work with 8 bits images (yet)
	if(maxval<255)
	{
		printf("8bit pictures not supported (yet).\n");
		free(kernel);
		return 5;
	}
	
	tIO = clock();

	void *blurred = malloc(xsize*ysize*sizeof(short unsigned int)); 	//piece of memory to memorize the blurred image 
	
	#pragma omp parallel
	{	
		//check endianism - eventually swap
  	if ( I_M_LITTLE_ENDIAN ) OMP_swap_image(image, xsize, ysize, maxval);
	
  	OMP_Convolve((unsigned short int*)image, (unsigned short int*)blurred, xsize, ysize ,kernel, xkernel, ykernel);	//actual convolution
    
  	// swap the endianism again
  	if ( I_M_LITTLE_ENDIAN ) OMP_swap_image(blurred , xsize, ysize, maxval);
	}
	
	//free the matrix resources and image vector
  free(kernel);
  free(image);

	tcalc = clock();

	/********************
	 output name setting 
	********************/

	char *output_name;
	if(args > arg_counter+1) output_name = argv[++arg_counter];		
	else 
	{
		input_name[strlen(input_name)-4]='\0';
		char charf[20];
		char out[50] = "";
			
		sprintf(charf,"%e",f);
		charf[1] = charf[2];
		charf[2] = '\0';
		
		if(ktype-1) sprintf(out,"%s.bb_%d_%dx%d.omp.pgm",input_name,ktype, xkernel, ykernel);
		else 				sprintf(out,"%s.bb_1_%dx%d_%s.omp.pgm",input_name, xkernel, ykernel,charf);
				
		output_name = out;
	}
	
	write_pgm_image(blurred, maxval, xsize, ysize, output_name);
	printf("Blurred image was succesfully stored in the file \"%s\"\n",output_name);
	
	twrite = clock();
	//free other resources
	free(blurred);
	
	printf("Walltime timings. Input: %lfs, Calculation (threads avg): %lfs, Output: %lfs. Total: %lfs\n",(double)(tIO-t0)/CLOCKS_PER_SEC,(double)(tcalc-tIO)/CLOCKS_PER_SEC/omp_get_max_threads(),(double)(twrite-tcalc)/CLOCKS_PER_SEC,(double)(twrite-t0)/CLOCKS_PER_SEC);
	
	return 0;
}
