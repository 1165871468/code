#include "Python.h"

void source_time_adj(int Ns, int Nr, double dt, double T)
{
    int ret;
    char command[512];

    ret = snprintf(command, sizeof(command),
                   "cd ./Upd && pythonx source_time_adj.py %d %d %.15g %.15g",
                   Ns, Nr, dt, T);

    if (ret < 0 || ret >= (int)sizeof(command)) {
        fprintf(stderr, "Error: command too long -- source_time_adj\n");
        exit(EXIT_FAILURE);
    }

    ret = system(command);
    if (ret != 0) {
        fprintf(stderr, "Error: source_time_adj failed, ret=%d\n", ret);
        exit(EXIT_FAILURE);
    }
}

void filter_bandpass(const char *input_file, const char *output_file)
{
    int ret;
    char command[512];

    ret = snprintf(command, sizeof(command), "cd ./Upd && pythonx filter_bandpass.py %s %s && cd ../", input_file, output_file);

    if (ret < 0 || ret >= (int)sizeof(command)) {
        fprintf(stderr, "Error: command too long -- filter_bandpass\n");
        exit(EXIT_FAILURE);
    }

    ret = system(command);

    if (ret == -1) {
        fprintf(stderr, "Error: system() failed: %s -- filter_bandpass\n",
                strerror(errno));
        exit(EXIT_FAILURE);
    }

    if (ret != 0) {
        fprintf(stderr,
                "Error: filter_bandpass.py failed, return code = %d\n",
                ret);
        exit(EXIT_FAILURE);
    }
}

int wave_window(int Ns, int Nr, double t_start, double t_end)
{
    char input_file[512], command[1024];
    int status;

    if (t_end <= t_start) {
        fprintf(stderr, "Error: invalid time window %.6f-%.6f\n", t_start, t_end);
        return -1;
    }

    snprintf(input_file, sizeof(input_file), "../../Output/Receivers_syn/wave_s%d_r%d.txt", Ns, Nr);

    snprintf(command, sizeof(command),
             "cd ./Upd && pythonx wave_window.py \"%s\" %.17g %.17g",
             input_file, t_start, t_end);

    status = system(command);

    if (status == -1 || !WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        fprintf(stderr, "Error: wave_window failed (s%d_r%d): %s\n",
                Ns, Nr, status == -1 ? strerror(errno) : "Python execution failed");
        return -1;
    }

    return 0;
}
