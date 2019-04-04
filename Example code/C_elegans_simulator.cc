
//  C_elegans_simulator.c
//
//
//  Created by James Giammona
//  This program simulates the development of the early C. elegans embryo using the Euler-Maruyama process for a particle in 3d


#include "simdefs.hh"
#include "printfuncs.hh"
#include "update.hh"
#include "vorofuncs.hh"


/*************************** main program ****************************/
int
main (void)
{
    printf("Starting Program\n");
    srand(time(NULL)); //This is necessary to generate random numbers
    
    char directory_buffer[100];
    snprintf(directory_buffer, sizeof(char) * 100, "C_elegans_test"); //Name the directory
    
    //Create directory (if none available)
    struct stat st = {0};
    
    if (stat(directory_buffer, &st) == -1) {
        mkdir(directory_buffer, 0777);
    }
    
    //Change to directory
    int change_dir;
    change_dir = chdir(directory_buffer);
    
    if (change_dir != 0) {
        printf("Error changing directory!\n");
        exit(1);
    }
    
    //We want div time, R, and div rule loops here
    //int potloopmax = 1;
    double pot_vals[] = {/*0.1,1.,*/10.};
    int potloop = 0;
    //int divloopmax = 1;
    double div_vals[] = {/*0.1,1.,*/10.};
    int divloop = 0;
    //int ruleloopmax = 1;
    const char* rules[] = {/*"order",*/ /*"order",*/ "Tdiv"};
    int ruleloop = 0;
    
    //Loops are time offset, time noise, and max div angle
    int offsetloopmax = 3;
    double time_offsets[] = {-0.1,0.,0.1};
    int timenoiseloopmax = 1;
    double time_noises[] = {0.,0.1};
    int divangleloopmax = 1;
    double div_angles[] = {0.,(PI/4.)};
    const char* div_angle_names[] = {"0",".25Pi"};
    
    int totloop = offsetloopmax*timenoiseloopmax*divangleloopmax;
    int totloopcounter = 0;
    
    for (int offsetloop = 0; offsetloop < offsetloopmax; offsetloop++) { //Loop through time offsets
        for (int timenoiseloop = 0; timenoiseloop < timenoiseloopmax; timenoiseloop++) { //Loop through timing noise values
            for (int divangleloop = 0; divangleloop < divangleloopmax; divangleloop++) { //Loop through division angle noise values
                
                cout << totloopcounter+1 << " out of " << totloop << endl;
                totloopcounter++;
                
                double AR = 1.6;
                double current_v = 1.2;

                double current_b = pow(current_v*(1./AR),(1./3.));
                //cout << "V:" << current_v << " b:" << current_b << endl;
                
                snprintf(directory_buffer, sizeof(char) * 100, "Toff%g_ClockN%g_Angle%s_AR%g_V%g_r0.7_N4_div%g_n5e-5_%s_pot%g",time_offsets[offsetloop], time_noises[timenoiseloop], div_angle_names[divangleloop],AR,current_v,div_vals[divloop],rules[ruleloop],pot_vals[potloop]); //Name the directory
                
                if (stat(directory_buffer, &st) == -1) {
                    mkdir(directory_buffer, 0777);
                }
                
                //Change to directory
                change_dir = chdir(directory_buffer);
                
                if (change_dir != 0) {
                    printf("Error changing directory!\n");
                    exit(1);
                }

                //Declare variables
                //double area_a, area_b;
                int snapshot_counter;
                
                double div_time = div_vals[divloop];
                double delta_t = 1e-3; //Timestep size
                double sdt = sqrt(delta_t); //Square root of timestep (used in update calc)
                int div_steps = (div_time)/delta_t;
                int div_modulo_const = div_steps/10;
                int dim = 3; //Spatial Dimension of simulation
                int max_divs = 2;
                int num_p_start = 1;
                int num_p = num_p_start; //How many points there are currently
                int max_num_p = num_p_start*(pow(2,max_divs)); //Max number of point
                
                int loopsimsmax = 2; //Num sims
                
                //Constants in force calculation
                double alpha = 4.0; //Scale of repulsive potential
                double beta = 3.0;
                double R = 1.; //Compared to equiliburim radius r
                //double r0 = r_vals[Rloop]*R;
                double r0 = 0.7*R;
                double a_const = 0.01; // Width of transition region in fermi function
                double a2 = 0.01; //Size of fermi cutoff for external potential
                //noise_int = 0.00005; //Noise intensity (related to kT in equilibrium)
                double noise_int = 5e-5; //Noise intensity (related to kT in equilibrium)
                //double clock_noise = div_time/10.;
                double clock_noise = time_noises[timenoiseloop]; //Want a relative error
                double div_noise = div_angles[divangleloop]; //Division angle noise in radians (uniformly from 0 to max angle)
                double time_offset = time_offsets[offsetloop]*div_time;
                //double noise_int = 0.;
                //Ellipsoid constants
                int pot_toggle = 1;
                int pot_type = 5; // 5 is repulsive confinement, 6 is adhesive confinement
                double ellipse_a = AR*current_b;//0.5 is to go from being normalized to 2r to r
                double ellipse_b = current_b;
                double ellipse_c = ellipse_b;
                double ellipse_percent = 0.9; //Percent of volume of ellipse that point volume should add up to
                double pot_strength = pot_vals[potloop];
                double maxforce_ever = 0.;
                
                char pos_buffer[40], dist_buffer[40], area_buffer[40], adj_buffer[40], force_buffer[40], radii_buffer[40], cell_parent_buffer[40];
                
                double **position;
                //Make max_num_p x max_num_p array
                vector< vector<int> > neighbor_array;
                int neighbor_array_init = 0;
                neighbor_array.resize(max_num_p, vector<int>(max_num_p,neighbor_array_init));
                
                //Make max_num_p x max_num_p array
                vector< vector<double> > area_array;
                double area_array_init = 0.;
                area_array.resize(max_num_p, vector<double>(max_num_p,area_array_init));
                
                //Make vectors for other important quantities
                vector<double> Ri(max_num_p);
                vector<double> divR(max_num_p);
                vector<double> r0i(max_num_p);
                vector<double> vcap_array(max_num_p);
                vector<double> a(max_num_p);
                vector<double> time2div(max_num_p);
                vector<double> time_left(max_num_p);
                vector<int> cell_parent(max_num_p);
                vector<int> num_divs_i(max_num_p);
                
                //File for cell parent/child info
                snprintf(cell_parent_buffer, sizeof(char) * 40, "cell_parent.csv");
                FILE *cell_parent_file = fopen(cell_parent_buffer, "w");
                Check_File(cell_parent_file);
        
                //Run loopsimmax number of simulations for each parameter value
                for (int loopsims = 0; loopsims < loopsimsmax; loopsims++) {
                    
//                    if(loopsims%10 == 0){
//                        cout << "Sim " << loopsims+1 << endl;
//                    }
                    
                    //Define data files
                    snprintf(pos_buffer, sizeof(char) * 40, "pos_%i.csv",loopsims);
                    snprintf(dist_buffer, sizeof(char) * 40, "dist_%i.csv",loopsims);
                    snprintf(area_buffer, sizeof(char) * 40, "area_%i.csv",loopsims);
                    snprintf(adj_buffer, sizeof(char) * 40, "adj_%i.csv",loopsims);
                    snprintf(force_buffer, sizeof(char) * 40, "maxforce_%i.csv",loopsims);
                    snprintf(radii_buffer, sizeof(char) * 40, "radii_%i.csv",loopsims);
                    
                    FILE *pos_file = fopen(pos_buffer, "w");
                    FILE *dist_file = fopen(dist_buffer, "w");
                    FILE *area_file = fopen(area_buffer, "w");
                    FILE *adj_file = fopen(adj_buffer, "w");
                    FILE *force_file = fopen(force_buffer, "w");
                    FILE *radii_file = fopen(radii_buffer, "w");
                    
                    Check_File(pos_file);
                    Check_File(dist_file);
                    Check_File(area_file);
                    Check_File(adj_file);
                    Check_File(force_file);
                    Check_File(radii_file);

                    //Allocate memory for position array
                    position = alloc_2d(max_num_p,dim);
                    
                   
                    //Initialize variables and vectors
                    num_p = num_p_start;
                    for (int init_loop1=0; init_loop1<max_num_p; init_loop1++) {
                        
                        cell_parent[init_loop1] = 0;
                        num_divs_i[init_loop1] = 0;
                        Ri[init_loop1] = R;
                        r0i[init_loop1] = r0;
                        divR[init_loop1] = R;
                        a[init_loop1] = a_const;
                        //Put in r0 initialization
                        vcap_array[init_loop1] = 0.;
                        time2div[init_loop1] = fabs(div_time);
                        if(init_loop1 == 1){ //Add offset to P1 cell relative to AB cell
                            time2div[init_loop1] = fabs(div_time + time_offset);
                        }
                        time2div[init_loop1] = fabs(time2div[init_loop1] + Gauss_Rand(0.,(clock_noise*time2div[init_loop1]))); //Add a relative error to the division time
                        time_left[init_loop1] = time2div[init_loop1];
                        
                        for (int init_loop2=0; init_loop2<dim; init_loop2++) {
                            position[init_loop1][init_loop2] = -999.; //Initialize array
                        }
                    }
                    
                    //Initialize arrays
                    for (int init_loop1=0; init_loop1<max_num_p; init_loop1++) {
                        for (int init_loop2=0; init_loop2<max_num_p; init_loop2++) {
                            area_array[init_loop1][init_loop2] = 0.;
                            neighbor_array[init_loop1][init_loop2] = 0;
                        }
                    }
                    
                    //Initialize positions
                    position[0][0] = 0. + Gauss_Rand(0, (ellipse_b/10.));			// initial x1 value
                    position[0][1] = 0. + Gauss_Rand(0, (ellipse_b/10.));			// initial y1 value
                    position[0][2] = 0. + Gauss_Rand(0, (ellipse_b/10.));			// initial z1 value
                    
            //        position[1][0] = 1.;			// initial x1 value
            //        position[1][1] = 0.;			// initial y1 value
            //        position[1][2] = 0.;			// initial z1 value
            //        
            //        position[2][0] = 2.;			// initial x1 value
            //        position[2][1] = 0.;			// initial y1 value
            //        position[2][2] = 0.;			// initial z1 value
                    
            //        position[3][0] = 0.6;			// initial x1 value
            //        position[3][1] = 0.6;			// initial y1 value
            //        position[3][2] = 0.;			// initial z1 value
                    
                    // Create a container with the geometry given above. This is bigger
                    // than the particle packing itself.
                    container_poly con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                                  false,false,false,8);
                    
                    // Create the "initial shape" wall class and add it to the container
                    wall_initial_shape(wis);
                    con.add_wall(wis);
                    
                    
                    snapshot_counter = 0;
                    maxforce_ever = 0.;
                    
                    //for (int div_counter=0; div_counter<(max_divs+1); div_counter++) { //It's max_divs +1 so you actually get the right number of divisions
                    int step_counter = 0;
                    int stop_token = 0;
                    int div_counter = 0;
                    int print_flag = 0;
                    int print_counter = 10;
                    
                    //Loop over timesteps, main loop for each simulation
                    while(stop_token == 0){
                    
                        if(Time_To_Divide(time_left, num_p) == 1){
                            
                            print_flag = 1;
                            //print_flag = 0; // Turning this flag on records the ten timesteps after division
                            
                            if(num_p == 1){
                                div_counter = 1;
                            }
                            if(num_p >= 2){
                                div_counter = 2;
                            }
                            
                            
                        //Do the divisions
                            num_p = Do_Divisions(position, num_p, 2 /*ruleloop*/, div_counter, max_num_p, dim, Ri, r0i, divR, a, ellipse_a, ellipse_b, ellipse_c, ellipse_percent, pot_toggle, time_left, time2div, div_time, clock_noise, cell_parent, div_noise, num_divs_i, time_offset);
                            //cout << num_p << endl;
//                            cout << Ri[0] << endl;
//                            cout << Ri[1] << endl;
                            //Add number of particles and default R value
                            wis.init_values(num_p, Ri[0]);
                        }
                
                        //for (int step_counter = 0; step_counter < div_steps; step_counter++) {
                            
                        int current_time = step_counter;
                        
                        con.clear();
                        for (int particle_loop = 0; particle_loop < num_p; particle_loop++) {
                            wis.set_R(particle_loop,Ri[particle_loop]);
                            con.put(particle_loop+1,position[particle_loop][0],position[particle_loop][1],position[particle_loop][2],Ri[particle_loop]);
                        }
                        
                        // Loop over all particles in the container and compute each Voronoi
                        // cell
                        c_loop_all cl(con);
                        //use max_num_p instead of 4
                        
                        
                        //Initialize arrays
                        for (int outiter = 0; outiter < max_num_p; outiter++) {
                            for (int initer = 0; initer < max_num_p; initer++) {
                                neighbor_array[outiter][initer] = 0;
                                area_array[outiter][initer] = 0.;
                            }
                        }
                        
                        //Find neighbors and face areas
                        Analyze_Voro_Container(con, neighbor_array, area_array);
                        CleanUp_Arrays(max_num_p,neighbor_array,area_array);
                        
                        con.clear();
                        
                        //Calculate new Ri's
                        FindVolumeOverlap(position, Ri, divR, num_p, dim, vcap_array, neighbor_array);
                        
                        for (int particle_loop = 0; particle_loop < num_p; particle_loop++) {
                            wis.set_R(particle_loop,Ri[particle_loop]);
                            con.put(particle_loop+1,position[particle_loop][0],position[particle_loop][1],position[particle_loop][2],Ri[particle_loop]);
                        }
                        
                        //Find neighbors and face areas after volume correction
                        Analyze_Voro_Container(con, neighbor_array, area_array);
                        CleanUp_Arrays(max_num_p,neighbor_array,area_array);

                        //Print matrices and Voronoi lines
                        if (step_counter % (div_modulo_const) == 0 || print_flag == 1) {
                            
                            snapshot_counter = Print_Arrays(con, max_num_p, loopsims, snapshot_counter, neighbor_array, area_array, area_file, adj_file);
//                                    cout << Ri[0] << endl;
//                                    cout << Ri[1] << endl;
//                                    cout << Ri[2] << endl;
//                                    cout << Ri[3] << endl;
                        }
                        
                        //Print cell info and update to new values
                        Print_Positions(position, dim, max_num_p, current_time, delta_t, div_modulo_const, pos_file, print_flag);
                        Print_Distances(position, dim, max_num_p, current_time, delta_t, div_modulo_const, dist_file, print_flag);
                        Print_Radii(Ri, divR, vcap_array, max_num_p, current_time, delta_t, div_modulo_const, radii_file, print_flag);
                        Update(position, dim, num_p, delta_t, sdt, alpha, beta, Ri, r0i, a, a2, noise_int, ellipse_a, ellipse_b, ellipse_c, ellipse_percent, pot_strength, neighbor_array, pot_toggle, pot_type, force_file, current_time*delta_t, maxforce_ever);
                        
                        step_counter++;
                        if(print_flag == 1){
                            print_counter--;
                        }
                        
//                            if(print_counter < 10){
//                                cout << print_counter << endl;
//                            }
                        
                        if(print_counter == 0){
                            print_counter = 10;
                            print_flag = 0;
                            
                        }
                        
                        Increment_Cell_Clock(time_left, num_p, delta_t);
                        if(num_p == max_num_p && Time_To_Divide(time_left, num_p) == 1){
                            stop_token = 1;
                        }
                    
                    }//While loop
                    
                    //Print cell_parent_file
                    for(int i=0; i<max_num_p; i++){
                        fprintf(cell_parent_file, "%i, %i", i, cell_parent[i] );
                        
                        if(i < max_num_p-1){
                            fprintf(cell_parent_file, ", ");
                        }
                        else{
                            fprintf(cell_parent_file, "\n");
                        }
                    }

                    //Free memory for position array
                    for (int k=0; k<max_num_p; k++) {
                        free(position[k]);
                    }
                    free(position);
                    
                    fclose(pos_file);
                    fclose(dist_file);
                    fclose(area_file);
                    fclose(adj_file);
                    fclose(force_file);
                    fclose(radii_file);
                    
                } // Numsims loop
                
                fclose(cell_parent_file);
                
                //Add change dir here
                change_dir = chdir("..");
                
                if (change_dir != 0) {
                    printf("Error changing directory!\n");
                    exit(1);
                }
                
            }//Div angle
        } //Clock noise loop
    } //Time offset loop
    
    printf("Program Ending\n");
}

/*************************** end main program ****************************/

