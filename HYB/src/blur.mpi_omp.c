#include <stdio.h>
#include <mpi.h>
#include <omp.h>
#include <string.h>
#include "ut.h"


int main(int args, char** argv)
{	
	int rank, size;
	int mpi_provided_thread_level;

	MPI_Init_thread( &args, &argv, MPI_THREAD_FUNNELED, &mpi_provided_thread_level);
	if ( mpi_provided_thread_level < MPI_THREAD_FUNNELED ) 
	{
		printf("a problem arise when asking for MPI_THREAD_FUNNELED level\n");
		MPI_Finalize();
		exit(1);
	}
	
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
	//time measurements
	double t0, tIO, tcomm, tcalc, tcomm2;
	t0 = MPI_Wtime();

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
		if(!rank) printf("Too few arguments in function blur. %s\n",usage);
		MPI_Finalize();
		return 1;
	}
	else if(args>MAX_ARGS)
	{
		if(!rank) printf("Too many arguments in function blur. %s\n",usage);
		MPI_Finalize();
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
				if(!rank) printf("Kernel dimensions must be odd integers. Dimensions given were %dx%d\n",xkernel,ykernel);
				MPI_Finalize();
				return 3;
			}
					
			kernel = (KTYPE *)malloc(xkernel*ykernel*sizeof(KTYPE));

			switch(ktype)
			{
				case 0:
					if(args > MAX_ARGS-1)
					{
						if(!rank) printf("Too many arguments in function blur. %s\n",usage);
						free(kernel);
						MPI_Finalize();
						return 1;
					}
					uniform_kernel(kernel,xkernel,ykernel);
					if(!rank) printf("Using Mean Kernel of dimension %dx%d\n",xkernel,ykernel);
				break;
					
				case 1:
					if(args < MAX_ARGS-1)
					{
						if(!rank) printf("Too few arguments in function blur, additional parameter for the kernel needed.\n%s\n",usage);
						free(kernel);
						MPI_Finalize();
						return 1;
					}
									
					f = atof(argv[++arg_counter]);
					weighted_kernel(kernel,xkernel,ykernel,f);
					if(!rank) printf("Using Weighted Kernel of dimension %dx%d and f = %f\n",xkernel,ykernel,f);
				break;
					
				case 2:
					if(args > MAX_ARGS-1)
					{
						if(!rank) printf("Too many arguments in function blur. %s\n",usage);
						free(kernel);				
						MPI_Finalize();
						return 1;
					}
								
					// sigma for the gaussian function ---> here it has been (arbitrarily) chosen to be the maximum half dimension of the kernel			
					int s = max((xkernel/2),(ykernel/2));
					gaussian_kernel(kernel,xkernel,ykernel,s*s);
					if(!rank) printf("Using Gaussian Kernel of dimension %dx%d (s: %d)\n",xkernel,ykernel,s);
				break;
			}
		break;

		case 3:
			if(args > MAX_ARGS-2)
			{
				if(!rank) printf("Too many arguments in function blur. No dimension needed for pgm kernel.\n %s\n",usage);
				MPI_Finalize();
				return 1;
			}
			
			int kmaxval;
			void *kimage;
			
			read_pgm_image(&kimage, &kmaxval, &xkernel, &ykernel, argv[++arg_counter]);
			kernel = normalize(kimage,xkernel,ykernel,kmaxval);
			free(kimage);

		break;

		default:
			if(!rank) printf("Kernel id must be either 0,1,2 for automatic generation or 3 for pgm file. Given id was %d",ktype);
			MPI_Finalize();
			return 2;
	}
	
	/********************
	 output name setting 
	********************/
	
	char* input_name = argv[++arg_counter];
	
	int namelen = strlen(input_name);
	if(namelen > 4 && !strcmp(&input_name[namelen-4],".pgm")) { if(!rank) printf("Trying to blur image \"%s\"\n",input_name); }
	else
	{
		if(!rank) printf("Input file name must end in \".pgm\". Given file name was %s\n",input_name);
		free(kernel);
		MPI_Finalize();
		return 4; //exit
	}
	

	/*********************************************************************************************************************
	*																																																										 *
	*																			  BLURRING - MPI COMMUNICATION																								 *
	*																																																										 *
	**********************************************************************************************************************/
	
	//pointer to be used to store the image to blur
	void *image = NULL;
	
	//useful info about image (unpacked, just declaration)  
	int maxval, xsize, ysize;
	
	//useful info (packed)
	int image_parameters[3];
	
	//master process reads the image from memory and distributes the interesting (packed) information to all via an MPI Broadcast
	if(!rank)	read_pgm_image(&image, image_parameters, image_parameters+1, image_parameters+2, input_name);
	MPI_Bcast(image_parameters,3,MPI_INT,0,MPI_COMM_WORLD);
	
	//unpacking of info
	maxval = image_parameters[0];
	xsize = image_parameters[1];
	ysize = image_parameters[2];
	
	//the following code is structured to work with 16bit images, so it wont work with 8 bits images (yet)
	if(maxval<255)
	{
		if(!rank) printf("8bit pictures not supported (yet).\n");
		free(kernel);
		MPI_Finalize();
		return 5;
	}
	
	//dimesion of the halo layer above and below the chunks (at least in non-pathological case, and before correction)
	//expressed in terms of number of pixels
	int halo_size = (ykernel/2)*xsize; 
	
	/*******************************************************************************************************************************
		WORK SUBDIVISION
	  
	  To ease the job of subdiving the image read from file and the reunification after the blurring, it has been decided to split
	  it in contiguous portions of memory, that is in rows. The master process blurres the first bunch of rows of the image, the first
	  slave the second and so on. Note that this is easier to implement but implies that the border in split unevenly between the 
	  processes. Since working on the border requires more FLOPs (need for renormalization), the workload is also distributed in a 
	  (possibly) uneven manner: if the total number of processes does not divide the the number of rows, the work is distributed 
	  evenly among the "middle processes" and the remaing amount of work is divided between master and last (the least workload is 
	  always given to the master) according to the subsequent formulas: */

		#define INNER_WORKLOAD(x,y,procs)  ( ( (y%procs) ? (y+procs-y%procs) : (y) )/procs * x )
		#define FIRST_WORKLOAD(x,y,procs)  ((y-(procs-2)*(INNER_WORKLOAD(x,y,procs)/x))/2 * x)
		#define LAST_WORKLOAD(x,y,procs)   (FIRST_WORKLOAD(x,y,procs) + ((y-(procs-2)*(INNER_WORKLOAD(x,y,procs)/x))%2 * x))
	
	/*If the number of processes does divide the the number of rows on the other hand, the workload is split evenly
	
	*********************************************************************************************************************************/

	//this represents the workload of each process
	int workload;
	//time took to do I/O stuff
	tIO = MPI_Wtime();
		
	/*
	* MASTER's TASK
	*/
	
	if(!rank)
	{
		//requests for the asynchronous communication
		MPI_Request *requests = (MPI_Request *)malloc((size-1)*sizeof(MPI_Request));
		
		//start and end position wrt the image pointer of the the chunks of memeory to be sent to each slave.
		//chunk is the number of pixels that each must receive, and includes halo layers.			
		int i, start, end, chunk;
				
		workload = FIRST_WORKLOAD(xsize,ysize,size);
		int inner = INNER_WORKLOAD(xsize,ysize,size);
		int last = LAST_WORKLOAD(xsize,ysize,size);
		
		for(i=1;i<size;i++)
		{			
			//this takes into account cases where start might be negative (which will be an error) and prevents it.
			//the same goes for the end 
			start = max(0,workload - halo_size + inner*(i-1));
			end = (rank != size-1)? min(xsize*ysize,workload + halo_size + inner*i) : xsize*ysize;
			
			chunk = end-start;
			
			MPI_Isend((unsigned short int *)image + start, chunk, MPI_UNSIGNED_SHORT, i, 0, MPI_COMM_WORLD, &requests[i-1]);
		}		
			
		/******************************************		
		*Gatherv requirements setup:
		******************************************/
		
		int *recv_size = (int *)malloc(size*sizeof(int));
		int *displs = (int *)malloc(size*sizeof(int));
		
		//this macro calculates the size of workload given the rank
		#define PROPER_SIZE(rank) (!(rank))?(workload):(((rank) == size-1)?(last):(inner))
	
		displs[0] = 0;
		for(i=1;i<size;i++) displs[i] = displs[i-1] + (recv_size[i-1] = PROPER_SIZE(i-1));
		recv_size[size-1] = last;
				
		/*****************************************/
		
		//now chunk refers to the masters own chunk of data, and also here we must prevent illegal indices
		chunk = min(xsize*ysize,workload+halo_size);
		
		//space up and down respresents the quantity of data contained in the upper and lower halo layers: they might differ from halo_size
		//if illegal indices have been found earlier, so we must correct. Of course the master only has a lower halo layer so up must be 0.
		int space_up = 0, space_down = chunk - workload;
		
		//allocatetes the space for the complete blurred image to be stored (also used for the local part of the master to be stored directly)
		void *blurred = malloc(sizeof(unsigned short int)*xsize*ysize);
		
		//Hereafter the buffer image is used in the convolution, so we must wait that all have recevied the correct image before we can modify it.
		MPI_Waitall(size-1, requests, MPI_STATUSES_IGNORE);
		free(requests);
		
		//time took to communicate
		tcomm = MPI_Wtime();
		
		/********************************************************
		* BLURRING - BEGIN
		*********************************************************/
		
		#pragma omp parallel
		{
			if ( I_M_LITTLE_ENDIAN) OMP_swap_image(image, xsize, chunk/xsize, maxval);
			
			//this function does the convolution of image and stores the result in blurred. 
			//It hadles different sizes of the two by means of the space up and down counters.
			OMP_MPIConvolve((unsigned short int *)image, (unsigned short int *)blurred, xsize, workload/xsize, kernel, xkernel, ykernel, space_up/xsize, space_down/xsize);
			
			if ( I_M_LITTLE_ENDIAN) OMP_swap_image(blurred, xsize, workload/xsize, maxval);
		}	
		free(image);
		/********************************************************
		* BLURRING - END
		*********************************************************/
     
     tcalc = MPI_Wtime();
     
    //recombination of the split image is done by means of Gatherv function. Here the sendbuffer is MPI_IN_PLACE since the master works already
    //in the array of the complete image by construction of the algorithm.    	
		MPI_Gatherv(MPI_IN_PLACE, workload, MPI_UNSIGNED_SHORT, blurred, recv_size, displs, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
		
		tcomm2 = MPI_Wtime(); 
		
		free(recv_size);
		free(displs);
		
		
		
		/*******************************************************************************************
		* Here the name of the output file and printing of the image is handled (only by the master) 
		********************************************************************************************/
		
		char *output_name;
		if(args > arg_counter+1) output_name = argv[++arg_counter];
		else 
		{
		input_name[strlen(input_name)-4]='\0';
		char charf[20];
		char out[100] = "";
			
		sprintf(charf,"%e",f);
		charf[1] = charf[2];
		charf[2] = '\0';
		
		if(ktype-1) sprintf(out,"%s.bb_%d_%dx%d.mpi_omp.pgm",input_name,ktype, xkernel, ykernel);
		else 				sprintf(out,"%s.bb_1_%dx%d_%s.mpi_omp.pgm",input_name, xkernel, ykernel,charf);
				
		output_name = out;
		}
	
		write_pgm_image(blurred, maxval, xsize, ysize, output_name);
		printf("Blurred image was succesfully stored in the file \"%s\"\n",output_name);
		
		free(blurred);
	}
	
	/*
	* SLAVES' TASK
	*/
	
	else
	{
		//the workload of each is establishe by means of the above definitions
		workload  = (rank != size-1) ? INNER_WORKLOAD(xsize,ysize,size) : LAST_WORKLOAD(xsize,ysize,size);
		
		int space_up, space_down, start, end;
		
		//this is totally analogous to the earlier discussion
		start = FIRST_WORKLOAD(xsize,ysize,size) - halo_size + INNER_WORKLOAD(xsize,ysize,size)*(rank-1);  //rank specific
		if(start<0)	space_up = halo_size + start, start = 0;
		else space_up = halo_size;
		
		if(rank != size-1)
		{
			end = FIRST_WORKLOAD(xsize,ysize,size) + halo_size + INNER_WORKLOAD(xsize,ysize,size)*rank;
			if(end>xsize*ysize) space_down = halo_size + xsize*ysize - end, end = xsize*ysize;
			else space_down = halo_size;
		}
		else end = xsize*ysize, space_down = 0;

		int chunk = end-start;
					
		//allocates space for the local copy of the image chunk to be stored			
		void *local_image = malloc(sizeof(unsigned short int)*chunk);
		
		//waits until the image is received
		MPI_Recv(local_image, chunk, MPI_UNSIGNED_SHORT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		//time took to communicate
		tcomm = MPI_Wtime();
		
		/********************************************************
		* BLURRING - BEGIN
		*********************************************************/
		
		//allocates space for the blurred local copy to be stored
		void *blurred = malloc(sizeof(unsigned short int)*workload);
		
		#pragma omp parallel
		{
			if ( I_M_LITTLE_ENDIAN) OMP_swap_image(local_image, xsize, chunk/xsize, maxval);
		
			OMP_MPIConvolve((unsigned short int *)local_image, (unsigned short int *)blurred, xsize, workload/xsize, kernel, xkernel, ykernel, space_up/xsize, space_down/xsize);
		
			if ( I_M_LITTLE_ENDIAN) OMP_swap_image(blurred, xsize, workload/xsize, maxval);
		}
		free(local_image);
		/********************************************************
		* BLURRING - END
		*********************************************************/
		
		//time for computation
		tcalc = MPI_Wtime(); 

		//sends back the local copy to master		
		MPI_Gatherv((unsigned short int *)blurred, workload, MPI_UNSIGNED_SHORT, NULL, NULL, NULL, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
		
		//time for the second communication
		tcomm2 = MPI_Wtime(); 
		free(blurred);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	printf("[%d] Walltime timings. I/0: %fs, Scattering: %fs, Calculation: %fs, Gathering: %fs. Total: %fs\n",rank,tIO-t0,tcomm-tIO,tcalc-tcomm,tcomm2-tcalc,tcomm2-t0);
	free(kernel);
	

	MPI_Finalize();
	return 0;
}
