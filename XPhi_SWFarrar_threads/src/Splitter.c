/*
 * Splitter.c
 *
 *  Created on: 17/03/2015
 *      Author: galvez
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <string.h>
#include <fcntl.h>
#include "Splitter.h"

static off_t get_filesize(const char * filename){
	FILE* fp;
	off_t file_size = -1;
	const off_t bad_size = -1;

	fp = fopen(filename, "r");
	if (fp == NULL) {
	  /* Handle error */
		return bad_size;
	}
	if (fseek(fp, 0 , SEEK_END) != 0) {
	  /* Handle error */
		return bad_size;
	}

	file_size = ftell(fp);
	if (file_size == -1) {
	  /* Handle error */
		return bad_size;
	}
	fclose(fp);
	return file_size;
}

static off_t look_around(FILE * fp, off_t pos){
	int offset_relative;
	int error;
	for(offset_relative=0;; offset_relative++){
		error =
		fseek(fp, pos + offset_relative, SEEK_SET);
		if (fgetc(fp) == '>') return pos + offset_relative;
		fseek(fp, pos - offset_relative, SEEK_SET);
		if (fgetc(fp) == '>') return pos - offset_relative;
		if (error) return -1;
	}
	return -1;
}

off_t * split(const char * filename, int16_t num_parts){
	int cur_part=0;
	off_t * slices = (off_t *)malloc(num_parts*2*sizeof(off_t));
	off_t file_size = get_filesize(filename);
	off_t slice_size = file_size / num_parts;
	FILE *fp = fopen(filename, "r");
	slices[0] = 0;
	for(cur_part=1; cur_part < num_parts; cur_part++){
		off_t pos_init_seq = look_around(fp, slice_size * cur_part);
		// Set end position of previous slice
		slices[2*cur_part - 1] = pos_init_seq - 1;
		// Set start position of current slice
		slices[2*cur_part] = pos_init_seq;
	}
	slices[2*num_parts - 1] = file_size;
	fclose(fp);
	return slices;
}
