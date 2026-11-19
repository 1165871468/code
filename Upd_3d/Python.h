#ifndef PYTHON_H
#define PYTHON_H

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

void source_time_adj(int Ns, int Nr, double dt, double T);
void filter_bandpass(const char *input_file, const char *output_file);
int wave_window(int Ns, int Nr, double t_start, double t_end);

#endif
