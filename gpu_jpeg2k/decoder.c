/* 
Copyright 2009-2013 Poznan Supercomputing and Networking Center

Authors:
Milosz Ciznicki miloszc@man.poznan.pl

GPU JPEG2K is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GPU JPEG2K is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with GPU JPEG2K. If not, see <http://www.gnu.org/licenses/>.
*/
/**
 * 	@file decoder.c
 *	@brief This the decoder main file. It loads JPEG 2000 file/stream, executes all modules in order to produce output in desired format.
 *
 *  @author Jakub Misiorny <misiorny@man.poznan.pl>
 *  @author Milosz Ciznicki
 */
//#include "threads/ThreadPool.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "config/arguments.h"
#include "config/init_device.h"
#include "config/help.h"
#include "preprocessing/preprocess_gpu.h"
#include "types/image_types.h"
#include "types/buffered_stream.h"
#include "types/image.h"


#include "print_info/print_info.h"

#include "dwt/dwt.h"
#include "dwt/iwt_1d.h"

#include "tier1/dequantizer.h"
#include "tier1/coeff_coder/gpu_coder.h"

#include "tier2/codestream.h"
#include "tier2/buffer.h"
#include "file_format/boxes.h"

#include "klt/klt.h"
#include <pthread.h>
#include <stdint.h>
#include "threads/thpool.h"

/**
 * @brief Main decoder function. It all starts here.
 * Input parameters:
 * 1) input JPEG 2000 image
 *
 * @return 0 on success
 */
int GPU_Decode(void *arg)
{	

//	println_start(INFO);
	type_image *img = (type_image *)malloc(sizeof(type_image));
	memset(img, 0, sizeof(type_image));
	img->in_file="artificial.j2k";
	img->out_file="artificial.bmp";
	//img->out_file=strcat("artificial.bmp",(char*)num);
	//img->conf_file="../config_files/lossy.config";
	img->conf_file=NULL;
	// check_args_dec(img);
	// if((parse_args(argc, argv, img) == ERROR) || (check_args_dec(img) == ERROR))
	// {
	// 	fprintf(stderr, "Error occurred while parsing arguments.\n");
	// 	fprintf(stdout, "%s", help);
	// 	return 1;
	// }

	type_parameters *param = (type_parameters*)malloc(sizeof(type_parameters));
	if((parse_config(img->conf_file, param) == ERROR) || (check_config(param) == ERROR)) {
		fprintf(stderr, "Error occurred while parsing configuration file.\n");
		fprintf(stdout, "%s", help);
		return 1;
	}
	
	//init_device(param);
	free(param);
	FILE *fsrc = fopen(img->in_file, "rb");

	

	if (!fsrc) {
		fprintf(stderr, "Error, failed to open %s for reading\n", img->in_file);
		return 1;
	}

	type_tile *tile;
	int i;

	if(strstr(img->in_file, ".jp2") != NULL) {
		println(INFO, "It's a JP2 file");

		//parse the JP2 boxes
		jp2_parse_boxes(fsrc, img);
		fclose(fsrc);

		println_var(INFO, "Number of tiles : %d", img->num_tiles);

		// Do decoding for all tiles
		for(i = 0; i < img->num_tiles; i++) {
			tile = &(img->tile[i]);
			/* Decode data */
			decode_tile(tile);
			/* Dequantize data */
			dequantize_tile(tile);
			/* Do inverse wavelet transform */
			iwt(tile);
		}

		if(img->use_mct == 1) {
			// lossless decoder
			if(img->wavelet_type == 0) {
				color_decoder_lossless(img);
			}
			else { //lossy decoder
				color_decoder_lossy(img);
			}
		} else if (img->use_part2_mct == 1) {
			decode_klt(img);
		} else {
			if(img->sign == UNSIGNED) {
				idc_level_shifting(img);
			}
		}
	}
	else {//It is not a JP2 file.
		type_buffer *src_buff = (type_buffer *) malloc(sizeof(type_buffer));

		init_dec_buffer(fsrc, src_buff);
		fclose(fsrc);
		long int start_t2;
		start_t2 = start_measure();
		decode_codestream(src_buff, img);
		free(src_buff->data);
		free(src_buff);
		printf("t2 %ld ms\n", stop_measure(start_t2)/1000);


	//	get_next_box(fsrc);
	printf("num_tiles %d \n", img->num_tiles);
		// Do decoding for all tiles
		for(i = 0; i < img->num_tiles; i++)	{
			tile = &(img->tile[i]);
			/* Decode data */
			long int start_t1;
			start_t1 = start_measure();
			decode_tile(tile);
			printf("t1 %ld ms\n", stop_measure(start_t1)/1000);

			/* Dequantize data */
			long int start_dequantion;
			start_dequantion = start_measure();
			dequantize_tile(tile);
			printf("dequantize %ld ms\n", stop_measure(start_dequantion)/1000);

			/* Do inverse wavelet transform */
			long int start_iwt;
			start_iwt = start_measure();
			iwt(tile);
			printf("iwt %ld ms\n", stop_measure(start_iwt)/1000);
		}

		if(img->use_mct == 1) {
			// lossless decoder
			if(img->wavelet_type == 0) {
				color_decoder_lossless(img);
			}
			else {  //lossy decoder
				color_decoder_lossy(img);
			}
		} else if (img->use_part2_mct == 1) {
			decode_klt(img);
		} else {
			if(img->sign == UNSIGNED) {
				idc_level_shifting(img);
			}
		}
	}

	save_image(img);
	free(img);
}
void* thread( void *arg )  
{  
    printf( "This is a thread and arg = %d.\n", *(int*)arg);  
    *(int*)arg = 0;  
    return arg;  
}  
int main(int argc, char **argv){
	
	type_parameters *param = (type_parameters*)malloc(sizeof(type_parameters));
	param->param_device=0;
	init_device(param);
	
	 pthread_t th;
	int ret;  
    int arg = 10;  
    int *thread_ret = NULL;  
	// GPU_Decode((void*)1);
    // ret = pthread_create( &th, NULL, GPU_Decode, &arg );  
    // if( ret != 0 ){  
    //     printf( "Create thread error!\n");  
    //     return -1;  
    // }  
    // printf( "This is the main process.\n" );  
    // pthread_join( th, (void**)&thread_ret );
	// printf( "thread_ret = %d.\n", *thread_ret);  
	int __argc = argc;
	char **__argv=argv;
	puts("Making threadpool with 4 threads");
	threadpool thpool = thpool_init(12);
	puts("Adding 40 tasks to threadpool");
	int i;
	for (i=0; i<100; i++){
	thpool_add_work(thpool, GPU_Decode, (void*)(uintptr_t)i);
	};

	thpool_wait(thpool);
	puts("Killing threadpool");
	thpool_destroy(thpool);
}