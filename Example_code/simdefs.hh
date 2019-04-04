#pragma once
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "voro++.hh"
#include <vector>
#include <iostream>
#include <string>
#include <numeric>
#include <algorithm>

using namespace std;
using namespace voro;

#define PI 3.14159265

struct coord{
    double x;
    double y;
    double z;
};

struct update_params{
    int dim;
    int num_p;
    double delta_t;
    double sdt;
    double alpha;
    double beta;
    double epsilon;
    double R;
    double a;
    double noise_int;
    struct coord box;
    double gamma;
    double shear;
    double elapsed_time;
};

double Ranf(void); //This generates a uniform random number from 0 to 1
double Gauss_Rand(double mean, double stddev); // This generates a normally distrubited random number
struct coord Div_Direction(int i); //Causes ordered or random division
struct coord Div_Direction2(int i, double div_noise); //Causes T division
struct coord Sphere_Rand(void); //Generate vector with three Gaussian numbers to generate points uniformly distributed on sphere

void Check_File(FILE *file1);

double rel_dist(struct coord pt1, struct coord pt2); //Find the distance between two vectors
struct coord rel_dir(struct coord pt1, struct coord pt2); //Find a vector pointing from 1 to 2
double dot_product(struct coord pt1, struct coord pt2); //Calculate the dot product of two vectors
struct coord add_vector(struct coord pt1, struct coord pt2); //Add vectors
struct coord scale_vector(struct coord pt1, double scale); //Scale vector by scalar
struct coord zero_vector(struct coord pt1); //Set vector to zero
double sign(double x);

double** alloc_2d(int m, int n); //These let us create arrays of various sizes
int** alloc_2d_int(int m, int n);

