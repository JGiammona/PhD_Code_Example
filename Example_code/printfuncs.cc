#include "simdefs.hh"
#include "printfuncs.hh"

//Print positions of particles
void Print_Positions(double **position, int dim, int num_p, int current_step, double delta_t, int modulo_const, FILE *output_file, int print_flag){
    
    int i,k;
    //int equil_steps = 1./delta_t;
    
    if (current_step%(modulo_const) == 0 || print_flag == 1){
        
        fprintf(output_file, " %.5e", current_step*delta_t );	/* current time */
        
        for (k=0; k<num_p; k++) {
            for (i=0; i<dim; i++) {
                fprintf(output_file, " ,%.5e", position[k][i] );	/* current positions */
            }
        }
        
        fprintf(output_file, "\n");
        
    }
    
    return;
}

//Print distances of particles
void Print_Distances(double **position, int dim, int num_p, int current_step, double delta_t, int modulo_const, FILE *distance_file, int print_flag){
    
    int i,j,k;
    //int equil_steps = 1./delta_t;
    double dist;
    vector<coord> pt(num_p);
    
    if (current_step%(modulo_const) == 0 || print_flag == 1){
        
        //Initialize points
        for (k=0; k<num_p; k++) {
            pt[k] = zero_vector(pt[k]); //Initialize to zero
            for (i=0; i<dim; i++) {
                if (i == 0) {
                    pt[k].x = position[k][i];
                }
                if (i == 1) {
                    pt[k].y = position[k][i];
                }
                if (i == 2) {
                    pt[k].z = position[k][i];
                }
            }
        }
        
        fprintf(distance_file, " %.5e", current_step*delta_t );	/* current time */
        
        //Loop over all particles i
        for (i=0; i<num_p; i++) {
            //Loop over all other particles j
            for (j=i; j<num_p; j++) {
                
                if (j==i) {
                    continue;
                }
                
                dist = rel_dist(pt[i],pt[j]);
                fprintf(distance_file, ", %i ,%i ,%.5e", i, j, dist );
                
            }
        }
        fprintf(distance_file, "\n");
        
    }//End if
    
    return;
}

//Print radii of particles
void Print_Radii(vector<double> &Ri, vector<double> &divR, vector<double> &vcap_array, int num_p, int current_step, double delta_t, int modulo_const, FILE *output_file, int print_flag){
    
    int k;
    //int equil_steps = 1./delta_t;
    
    if (current_step%(modulo_const) == 0 || print_flag == 1){
        
        fprintf(output_file, " %.5e", current_step*delta_t );	/* current time */
        
        for (k=0; k<num_p; k++) {
            fprintf(output_file, " ,%.5e  ,%.5e  ,%.5e", Ri[k], divR[k], vcap_array[k]);	/* current radii */
        }
        
        fprintf(output_file, "\n");
        
    }
    
    return;
}

//Output Voro++ vertices and edges for plotting in Mathematica or gnuplot
int Print_Arrays(container_poly &con1, int max_num_p, int loopsims, int snapshot_counter, const vector<vector<int> > &neighbor_array, const vector<vector<double> > &area_array, FILE *area_file, FILE *adj_file){
    
    for (int outiter = 0; outiter < max_num_p; outiter++) {
        for (int initer = 0; initer < max_num_p; initer++) {
            fprintf(area_file, " %.5e", area_array[outiter][initer]);
            if (initer < max_num_p - 1) {
                fprintf(area_file, ",");
            }
            fprintf(adj_file, " %i", neighbor_array[outiter][initer]);
            if (initer < max_num_p - 1) {
                fprintf(adj_file, ",");
            }
        }
        fprintf(area_file, "\n");
        fprintf(adj_file, "\n");
    }
    
    char voro_buffer[80];
    //con.draw_particles("irregular_p1.dat");
    snprintf(voro_buffer, sizeof(char) * 80, "dodec_sim%i_%i.dat", loopsims, snapshot_counter);
    con1.draw_cells_gnuplot(voro_buffer);
    snapshot_counter++;
    
    return snapshot_counter;
}

