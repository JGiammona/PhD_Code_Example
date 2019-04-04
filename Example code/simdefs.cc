#include "simdefs.hh"

//Generate uniform rand num from 0 to 1
double Ranf(){
    double num;
    
    num = rand()/((double) RAND_MAX);
    
    return num;
}

//Generate normally distributed random number
double Gauss_Rand(double mean, double stddev){
    double x1, x2, w, y1, y2;
    
    do {
        x1 = 2.0 * Ranf() - 1.0;
        x2 = 2.0 * Ranf() - 1.0;
        w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;
    
    return (stddev*y1)+mean;
}

//Generates a vector to cause ordered division
struct coord Div_Direction(int i){
    
    struct coord dir_vect;
    double x,y,z, len;
    
    x = 0.;
    y = 0.;
    z = 0.;
    
    if(i == 1){
        x = 1.;
    }
    else if(i == 2){
        y = 1.;
    }
    else if(i == 3){
        z = 1.;
    }
    
    dir_vect.x = x;
    dir_vect.y = y;
    dir_vect.z = z;
    
    len = sqrt((x*x)+(y*y)+(z*z));
    
    if(len < 0.1){
        printf("Length is too small");
    }
    
    dir_vect = scale_vector(dir_vect,(1/len));
    
    return dir_vect;
}

//Generates a vector to cause T division (which is wild type), div_noise in radians
struct coord Div_Direction2(int i, double div_noise){
    
    struct coord dir_vect;
    double x,y,z, len;
    double theta, phi;
    
    phi = (2.0*PI)*Ranf(); //in radians
    theta = div_noise*Ranf(); //in radians
    
    x = 0.;
    y = 0.;
    z = 0.;
    
    if(i == 1){
        x = cos(theta);
        
        y = sin(theta)*sin(phi);
        z = sin(theta)*cos(phi);
    }
    else if(i == 2){
        y = cos(theta);
        
        x = sin(theta)*sin(phi);
        z = sin(theta)*cos(phi);
    }
    else if(i == 3){
        z = cos(theta);
        
        x = sin(theta)*sin(phi);
        y = sin(theta)*cos(phi);
    }
    
    dir_vect.x = x;
    dir_vect.y = y;
    dir_vect.z = z;
    
    len = sqrt((x*x)+(y*y)+(z*z));
    
    if(len < 0.1){
        printf("Length is too small");
    }
    
    dir_vect = scale_vector(dir_vect,(1/len));
    
    return dir_vect;
}

//Generates a vector of three random gaussian random numbers to produce a uniformly distributed set of points on the sphere
struct coord Sphere_Rand(){
    
    struct coord rand_vect;
    double x,y,z, len;
    
    x = Gauss_Rand(0,1);
    y = Gauss_Rand(0,1);
    z = Gauss_Rand(0,1);
    
    rand_vect.x = x;
    rand_vect.y = y;
    rand_vect.z = z;
    
    len = sqrt((x*x)+(y*y)+(z*z));
    
    rand_vect = scale_vector(rand_vect,(1/len));
    
    return rand_vect;
}

//Checks if file was made
void Check_File(FILE *file1){
    if (file1 == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
    return;
}


//Find the distance between two points
double rel_dist(struct coord pt1, struct coord pt2){
    double dist, sqr_dist;
    sqr_dist = pow((pt1.x-pt2.x),2)+pow((pt1.y-pt2.y),2)+pow((pt1.z-pt2.z),2);
    dist = sqrt(sqr_dist);
    return dist;
}

//Find the vector pointing in the direction from pt1 to pt2
struct coord rel_dir(struct coord pt1, struct coord pt2){
    struct coord rel_vector;
    double rel_len;
    rel_vector.x = pt2.x-pt1.x;
    rel_vector.y = pt2.y-pt1.y;
    rel_vector.z = pt2.z-pt1.z;
    //Normalize vector
    rel_len = sqrt(dot_product(rel_vector,rel_vector));
    rel_vector = scale_vector(rel_vector,(1/rel_len));
    return rel_vector;
}

//Calculate the dot product
double dot_product(struct coord pt1, struct coord pt2){
    double dot;
    dot = (pt1.x*pt2.x)+(pt1.y*pt2.y)+(pt1.z*pt2.z);
    return dot;
}



//Add two vectors together
struct coord add_vector(struct coord pt1, struct coord pt2){
    struct coord sum_vector;
    sum_vector.x = pt1.x+pt2.x;
    sum_vector.y = pt1.y+pt2.y;
    sum_vector.z = pt1.z+pt2.z;
    return sum_vector;
}

//Multiply a vector by a scalar
struct coord scale_vector(struct coord pt1, double scale){
    pt1.x=pt1.x*scale;
    pt1.y=pt1.y*scale;
    pt1.z=pt1.z*scale;
    return pt1;
}

//Set all coords of a vector to zero
struct coord zero_vector(struct coord pt1){
    pt1.x=0;
    pt1.y=0;
    pt1.z=0;
    return pt1;
}

//Find the sign of a double
double sign(double x){
    
    if (x > 0) return 1.;
    if (x < 0) return -1.;
    return 0.;
    
}

//Create 2D array
double** alloc_2d(int m, int n){
    
    double** position = (double**)malloc(m * sizeof(double*));
    if (position == 0)
    {
        printf("ERROR: Out of memory for position\n");
        return position;
    }
    
    int k;
    for (k=0; k < m; ++k) {
        position[k] = (double*)malloc(n * sizeof(double));
        if (position[k] == 0)
        {
            printf("ERROR: Out of memory for position[%i]\n",k);
            return position;
        }
    }
    
    return position;
}

//Create 2D int array
int** alloc_2d_int(int m, int n){
    
    int** position = (int**)malloc(m * sizeof(int*));
    if (position == 0)
    {
        printf("ERROR: Out of memory for position\n");
        return position;
    }
    
    int k;
    for (k=0; k < m; ++k) {
        position[k] = (int*) malloc(n * sizeof(int));
        if (position[k] == 0)
        {
            printf("ERROR: Out of memory for position[%i]\n",k);
            return position;
        }
    }
    
    return position;
}


