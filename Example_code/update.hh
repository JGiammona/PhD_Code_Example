#pragma once
#include "simdefs.hh"

//Do force numerical integration for one timestep
void Update(double **position, int dim, int num_p, double delta_t, double sdt, double alpha, double beta, vector<double> &Ri, vector<double> &r0i, vector<double> &a, double a2, double noise_int, double ellipse_a, double ellipse_b, double ellipse_c, double ellipse_percent, double pot_strength, const vector<vector<int> > &neighbor_array, int pot_toggle, int pot_type, FILE *maxdisp_file, double elapsed_time, double &maxforce_ever);

//Produce envelope force on cells
struct coord Envelope_Potential(struct coord position, double ellipse_a, double ellipse_b, double ellipse_c, double R, double r0, double alpha, double beta, double a, double a2, int num_p, double ellipse_percent, double pot_strength, int pot_type);

//Cause a cell to divide into two cells
int Do_Divisions(double **position, int num_p, int rule_counter, int div_counter, int max_num_p, int dim, vector<double> &Ri, vector<double> &r0i, vector<double> &divR, vector<double> &a, double ellipse_a, double ellipse_b, double ellipse_c, double ellipse_percent, int pot_toggle, vector<double> &time_left, vector<double> &time2div, double div_time, double clock_noise, vector<int> &cell_parent, double div_noise, vector<int> &num_divs_i, double time_offset);

//Functional form of force between cells
double force_ij(double dist, double alpha, double beta, double a, double R, double r0); //Calculate the force between two points

//Use overlap between cells and use that to increase cell volume
void FindVolumeOverlap(double **position, vector<double> &Ri, vector<double> &divR, int max_num_p, int dim, vector<double> &vcap_array, const vector<vector<int> > &neighbor_array);
double CapVolume(double R1, double R2, double dist);

//Use Voro++ to find cell neighbors and areas of contact
void Analyze_Voro_Container(container_poly &con1, vector<vector<int> > &neighbor_array, vector<vector<double> > &area_array);

//Remove pair from neighbor and area arrays if area is too small or neighbor pair goes in only one direction
void CleanUp_Arrays(int max_num_p, vector<vector<int> > &neighbor_array, vector<vector<double> > &area_array);

//Determine if cells have reached the tethrahedron equilibrium
int Array_Is_Equil(double **position, int dim, int max_num_p, vector<vector<int> > &neighbor_array, vector<double> &Ri);

//Checks if time_left has gone to zero or below
int Time_To_Divide(vector<double> &time_left, int num_p);

//Decrements entry in time_left vector by delta_t
void Increment_Cell_Clock(vector<double> &time_left, int num_p, double delta_t);
