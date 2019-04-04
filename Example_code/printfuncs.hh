#pragma once

void Print_Positions(double **position, int dim, int num_p, int current_step, double delta_t, int modulo_const, FILE *output_file, int print_flag);
void Print_Distances(double **position, int dim, int num_p, int current_step, double delta_t, int modulo_const, FILE *distance_file, int print_flag);
void Print_Radii(vector<double> &Ri, vector<double> &divR, vector<double> &vcap_array, int num_p, int current_step, double delta_t, int modulo_const, FILE *output_file, int print_flag);
int Print_Arrays(container_poly &con1, int max_num_p, int loopsims, int snapshot_counter, const vector<vector<int> > &neighbor_array, const vector<vector<double> > &area_array, FILE *area_file, FILE *adj_file);
