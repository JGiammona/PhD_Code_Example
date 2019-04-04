#include "simdefs.hh"
#include "update.hh"

//This function updates the positions of all the particles using the Euler Maruyama method of numerically integrating a stochastic DE
void Update(double **position, int dim, int num_p, double delta_t, double sdt, double alpha, double beta, vector<double> &Ri, vector<double> &r0i, vector<double> &a, double a2, double noise_int, double ellipse_a, double ellipse_b, double ellipse_c, double ellipse_percent, double pot_strength, const vector<vector<int> > &neighbor_array, int pot_toggle, int pot_type, FILE *maxdisp_file, double elapsed_time, double &maxforce_ever){
    
    int i,j,k, print_token;
    double kappa, force, dist;
    kappa = 1.; //Strength of force(?)
    vector<coord> pt(num_p);
    vector<coord> force_vector(num_p);
    struct coord direction;
    
    print_token = 0;
    
    //Initialize points
    for (k=0; k<num_p; k++) {
        pt[k] = zero_vector(pt[k]); //Initialize to zero
        force_vector[k] = zero_vector(force_vector[k]);
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
    
    double maxforce_sofar = 0.;
    double force_mag = 0.;
    
    //Loop over each particle i
    for (i=0; i<num_p; i++) {
        force_vector[i] = zero_vector(force_vector[i]);
        
        //Loop over all other particles j
        for (j=0; j<num_p; j++) {
            
            if (neighbor_array[i][j]!=0 && rel_dist(pt[i],pt[j]) < (Ri[i]+Ri[j])) {
                dist = rel_dist(pt[i],pt[j]);
                
                //Just for compiling
                double avg_r0 = (r0i[i]+r0i[j])/2.;
                
                force = force_ij(dist,alpha,beta,a[i],((Ri[i]+Ri[j])/2.),avg_r0);
                //force = 0.;
                //(Ri[i]+Ri[j])
            }
            else {
                continue;
            }
            
            if (print_token == 1) { //print token
                printf("Force on %i from %i is %.5e\n",i,j, force );
            }
            //Vector from j to i
            direction = rel_dir(pt[j], pt[i]);
            /*
             if (print_token) {
             printf("Force direction x-%.5e y-%.5e z-%.5e\n", direction.x, direction.y, direction.z );
             }*/
            
            force_vector[i] = add_vector(force_vector[i],scale_vector(direction,force));
            
//             if (print_token) {
//             printf("Force x-%.5e y-%.5e z-%.5e\n", force_vector[i].x, force_vector[i].y, force_vector[i].z );
//             }
            
        }
        
        // length of force_vector
        force_mag = sqrt(dot_product(force_vector[i],force_vector[i]));
        
        if(force_mag > maxforce_sofar){
            maxforce_sofar = force_mag;
        }
        
        //printf("Cell %i in contact: %i,%i,%i,%i\n",i,in_contact[i][0],in_contact[i][1],in_contact[i][2],in_contact[i][3]);
    }
    if(maxforce_sofar*delta_t > 1e-3 && maxforce_sofar>maxforce_ever){
        fprintf(maxdisp_file, "%.5e, %.5e \n",elapsed_time, maxforce_sofar*delta_t );
        maxforce_ever = maxforce_sofar;
    }
    
    //Calculate the new positions from the old positions and forces on the particles
    for (k=0; k<num_p; k++) {
        for (i=0; i<dim; i++) {
            if (i == 0) {
                position[k][i] = position[k][i] + (kappa*delta_t*force_vector[k].x) + (sdt*sqrt(2.*noise_int)*Gauss_Rand(0,1));
                if(pot_toggle == 1){
                    position[k][i] = position[k][i] + (delta_t*(Envelope_Potential(pt[k], ellipse_a, ellipse_b, ellipse_c, Ri[k], r0i[k], alpha, beta, a[k], a2, num_p, ellipse_percent, pot_strength, pot_type).x));
                }
                
            }
            if (i == 1) {
                position[k][i] = position[k][i] + (kappa*delta_t*force_vector[k].y) + (sdt*sqrt(2.*noise_int)*Gauss_Rand(0,1));
                if(pot_toggle == 1){
                    position[k][i] = position[k][i] + (delta_t*(Envelope_Potential(pt[k], ellipse_a, ellipse_b, ellipse_c, Ri[k], r0i[k], alpha, beta, a[k], a2, num_p, ellipse_percent, pot_strength, pot_type).y));
                }
            }
            if (i == 2) {
                position[k][i] = position[k][i] + (kappa*delta_t*force_vector[k].z) + (sdt*sqrt(2.*noise_int)*Gauss_Rand(0,1));
                if(pot_toggle == 1){
                    position[k][i] = position[k][i] + (delta_t*(Envelope_Potential(pt[k], ellipse_a, ellipse_b, ellipse_c, Ri[k], r0i[k], alpha, beta, a[k], a2, num_p, ellipse_percent, pot_strength, pot_type).z));
                }
            }
            
        }
    }
    
    return;
}

//This calculates the force from an elliptical external potential and scales up the size to conserve volume after cell divisions (if volume toggle on). It uses the negative gradient of a potential defined elsewhere
struct coord Envelope_Potential(struct coord position, double ellipse_a, double ellipse_b, double ellipse_c, double R, double r0, double alpha, double beta, double a, double a2, int num_p, double ellipse_percent, double pot_strength, int pot_type){
    
    struct coord force;
    double scaling, xsqr, ysqr, zsqr, asqr, bsqr, csqr, abc, force_const, dist_const, force_strength, exp_constant, exp_constant2, sharpness, force_mag;
    int potential_flag;
    
    potential_flag = pot_type; // 1 is old potential based on repulsive term from the cells themselves
    int volume_toggle = 0; // If on, a,b,c are scaled to contsrain the total volume
    // 2 is the new potential based on a fermi support function convolved with a linear force
    //3 is simlpe linear force
    //4 is adhesion potential (depricated, from ellipsoid potential 2.nb)
    //5 is new confining potential
    //6 is new adhesion potential (from adhesion potential new.nb)
    
    force_strength = pot_strength;
    sharpness = 5.;
    
    abc = ellipse_a*ellipse_b*ellipse_c;
    
    //Sphere radius
    
    scaling = (R)*(pow((num_p)/(ellipse_percent*abc),(1./3.)));
//    scaling2 = pow((num_p),(1./3.));
    
    
    //scaled_a = (a*scaling)/(R/2.);
    
    force = zero_vector(force);
    
    xsqr = position.x*position.x;
    ysqr = position.y*position.y;
    zsqr = position.z*position.z;
    
    
    if(volume_toggle == 1){
        asqr = (ellipse_a*ellipse_a*scaling*scaling);
        bsqr = (ellipse_b*ellipse_b*scaling*scaling);
        csqr = (ellipse_c*ellipse_c*scaling*scaling);
    }
    else{
        //arsqr = (ar*ar);
        asqr = (ellipse_a*ellipse_a);
        bsqr = (ellipse_b*ellipse_b);
        csqr = (ellipse_c*ellipse_c);
    }
    
    //dist_const = sqrt((xsqr/arsqr)+(ysqr)+(zsqr));
    dist_const = sqrt((xsqr/(asqr/bsqr))+(ysqr)+(zsqr));
    
    exp_constant = exp((1.-dist_const)/a2);
    exp_constant2 = a2*(1+exp_constant)*(1+exp_constant);
    
    //1 is old potential based on repulsive term from the cells themselves
    if (potential_flag == 1){
        if (dist_const > .6){
            dist_const = .6;
            force_const = ((pow((1.-dist_const),(-1.-alpha)))*alpha*beta)/((alpha-beta)*dist_const);
        }
        else if (dist_const < .001){
            force_const = 0.;
        }
        else{
            force_const = ((pow((1.-dist_const),(-1.-alpha)))*alpha*beta)/((alpha-beta)*dist_const);
        }
        
        
        force.x = -(position.x/asqr)*force_const;
        force.y = -(position.y/bsqr)*force_const;
        force.z = -(position.z/csqr)*force_const;
        
        force_mag = sqrt(dot_product(force,force));
        
        if (force_mag > 1.0E2){
            
            force = scale_vector(force,1.0E2/force_mag);
            
        }
    }
    
    // 2 is the new potential based on a fermi support function convolved with a linear force
    if (potential_flag == 2){
        
        if (dist_const < .001){
            force.x = 0.;
            force.y = 0.;
            force.z = 0.;
        }
        else{
            force.x = -((force_strength*exp_constant*position.x)/(asqr*exp_constant2))-(force_strength*position.x/(asqr*(1+exp_constant)*dist_const));
            force.y = -((force_strength*exp_constant*position.y)/(bsqr*exp_constant2))-(force_strength*position.y/(bsqr*(1+exp_constant)*dist_const));
            force.z = -((force_strength*exp_constant*position.z)/(csqr*exp_constant2))-(force_strength*position.z/(csqr*(1+exp_constant)*dist_const));
        }
        
    }
    
    //3 is simlpe linear force
    if (potential_flag == 3){
        force.x = -((2*force_strength*position.x)/(asqr));
        force.y = -((2*force_strength*position.y)/(bsqr));
        force.z = -((2*force_strength*position.z)/(csqr));
        
        force_mag = sqrt(dot_product(force,force));
        
        if (force_mag > 5.0E2){
            
            force = scale_vector(force,5.0E2/force_mag);
            
        }
    }
    
    //4 is adhesion potential (depricated, from ellipsoid potential 2.nb)
    if (potential_flag == 4){
        
        
        
        double first_term = 2.*force_strength*alpha*beta*(pow((sharpness*(dist_const-(R/2.))),(-alpha-beta)));
        double second_term = (pow((sharpness*(dist_const-(R/2.))),alpha))+(pow((sharpness*(dist_const-(R/2.))),beta));
        double third_term = dist_const*(alpha-beta)*(2.*dist_const-R);
        
        if (dist_const > ((R/2.))){
            force.x = -((2*force_strength*position.x)/(asqr));
            force.y = -((2*force_strength*position.y)/(bsqr));
            force.z = -((2*force_strength*position.z)/(csqr));
        }
        else if (dist_const < .001){
            force_const = 0.;
            force.x = (position.x/asqr)*force_const;
            force.y = (position.y/bsqr)*force_const;
            force.z = (position.z/csqr)*force_const;
        }
        else{
            force_const = (first_term*second_term)/third_term;
            force.x = (position.x/asqr)*force_const;
            force.y = (position.y/bsqr)*force_const;
            force.z = (position.z/csqr)*force_const;
        }
        
        force_mag = sqrt(dot_product(force,force));
        
        if (force_mag > 5.0E2){
            
            force = scale_vector(force,5.0E2/force_mag);
            
        }
        
    }
    
    //5 is new confining potential
    if (potential_flag == 5){
        
        if (dist_const > .99*ellipse_b){
            dist_const = .99*ellipse_b;
            //force_const = 2./(dist_const*pow((1.-dist_const),3.));
            force_const = 1./(dist_const*pow((ellipse_b-dist_const),2.));
        }
        else if(dist_const < 0.01){
            force_const = 0.;
        }
        else if (dist_const < ((ellipse_b - R)/1.)){
            force_const = 0.;
        }
        else{
            //force_const = 2./(dist_const*pow((1.-dist_const),3.));
            force_const = 1./(dist_const*pow((ellipse_b-dist_const),2.));
        }
        
        force.x = -(position.x/(asqr/bsqr))*force_const*pot_strength;
        force.y = -(position.y/1.)*force_const*pot_strength;
        force.z = -(position.z/1.)*force_const*pot_strength;
        
        force_mag = sqrt(dot_product(force,force));
        
        if (force_mag > 2.0E2){
            force = scale_vector(force,2.0E2/force_mag);
        }
    }
    
    //6 is new adhesion potential (from adhesion potential new.nb)
    if (potential_flag == 6){
        
        if (dist_const > .99*ellipse_b){
            dist_const = .99*ellipse_b;
            force_const = (1./(alpha-beta))*((alpha*beta)/(dist_const))*((pow(r0,alpha)/pow((ellipse_b-dist_const),(alpha+1.)))-(pow(r0,beta)/pow((ellipse_b-dist_const),(beta+1.))));
        }
        else if(dist_const < 0.01){
            force_const = 0.;
        }
        else if (dist_const < ((ellipse_b - R)/1.)){
            force_const = 0.;
        }
        else{
            force_const = (1./(alpha-beta))*((alpha*beta)/(dist_const))*((pow(r0,alpha)/pow((ellipse_b-dist_const),(alpha+1.)))-(pow(r0,beta)/pow((ellipse_b-dist_const),(beta+1.))));
        }
        
        force.x = -(position.x/(asqr/bsqr))*force_const*pot_strength;
        force.y = -(position.y/1.)*force_const*pot_strength;
        force.z = -(position.z/1.)*force_const*pot_strength;
        
        force_mag = sqrt(dot_product(force,force));
        
        if (force_mag > 2.0E2){
            force = scale_vector(force,2.0E2/force_mag);
        }
    }
    
    return force;
}

//This function calculates the positions of the newly divided particles
int Do_Divisions(double **position, int num_p, int rule_counter, int div_counter, int max_num_p, int dim, vector<double> &Ri, vector<double> &r0i, vector<double> &divR, vector<double> &a, double ellipse_a, double ellipse_b, double ellipse_c, double ellipse_percent, int pot_toggle, vector<double> &time_left, vector<double> &time2div, double div_time, double clock_noise, vector<int> &cell_parent, double div_noise, vector<int> &num_divs_i, double time_offset){
    
    int k,l,m, num_sofar, point_check, pass;
    int div_counter1 = 0;
    struct coord div_dir;
    vector<coord> pos(max_num_p);
    double check_dist, dist_const, abc, xsqr, ysqr, zsqr, asqr, bsqr, csqr;
    
    //Shrink cells, they have divided, we are now just placing the new cells
    //Want new volume to be half the old volume so R needs to decrease by 2^(1/3)
    
    num_sofar = num_p;
    
    //Loop through all current cells
    for (k=0; k<num_p; k++) {
        
        //Check if it's time for each cell to divide
        if (time_left[k] > 0.){
            continue;
        }
        
        //Make dividing cell smaller and reset clock
        if (time_left[k] <= 0.){
            time_left[k] = fabs(div_time + Gauss_Rand(0.,(clock_noise*div_time)));
            Ri[k] = Ri[k]/(pow(2,(1./3.)));
            r0i[k] = r0i[k]/(pow(2,(1./3.)));
            a[k] = a[k]/(pow(2,(1./3.)));
            divR[k] = divR[k]/(pow(2,(1./3.)));
        }
        
        point_check = 0; // Reset check for each currrent point
        
        //Place new cells perpendicular from old cells in x, y or z positions based on div_counter
        if(rule_counter == 0){
            
            div_dir = Div_Direction(div_counter);
            
            position[num_sofar][0] = position[k][0] + 1.*Ri[k]*div_dir.x;
            position[num_sofar][1] = position[k][1] + 1.*Ri[k]*div_dir.y;
            position[num_sofar][2] = position[k][2] + 1.*Ri[k]*div_dir.z;
            
            Ri[num_sofar] = Ri[k];
            r0i[num_sofar] = r0i[k];
            a[num_sofar] = a[k];
            divR[num_sofar] = divR[k];
            cell_parent[num_sofar] = k;
            num_divs_i[k] = num_divs_i[k] + 1;
            
        }
        
        //Place new cell randomly around old cell
        if(rule_counter == 1){
            
            int check_counter = 0;
            double min_dist = 1.2*Ri[k];
            
            while (point_check == 0) {
                
                //Initialize positions of all particles so far
                
                for (l=0; l<num_sofar; l++) {
                    pos[l] = zero_vector(pos[l]); //Initialize to zero
                    for (m=0; m<dim; m++) {
                        if (m == 0) {
                            pos[l].x = position[l][m];
                        }
                        if (m == 1) {
                            pos[l].y = position[l][m];
                        }
                        if (m == 2) {
                            pos[l].z = position[l][m];
                        }
                    }
                }
                
                //Place new particle at 0.5 distance from center of old particle
                div_dir = Sphere_Rand();
                
                div_dir = scale_vector(div_dir,min_dist);
                
                pos[num_sofar] = add_vector(pos[k],div_dir);
                
                pass = 1;
                
                abc = ellipse_a*ellipse_b*ellipse_c;
                
                //2*num_p because we already want the ellipse to be big
                //scaling = (R/2.)*(pow((2*num_p)/(ellipse_percent*abc),(1./3.)));
                //scaling2 = pow((num_p),(1./3.));
                
                xsqr = pos[k].x*pos[k].x;
                ysqr = pos[k].y*pos[k].y;
                zsqr = pos[k].z*pos[k].z;
                
                asqr = (ellipse_a*ellipse_a);
                bsqr = (ellipse_b*ellipse_b);
                csqr = (ellipse_c*ellipse_c);
                
                //asqr = (ellipse_a*ellipse_a*scaling*scaling);
                //bsqr = (ellipse_b*ellipse_b*scaling*scaling);
                //csqr = (ellipse_c*ellipse_c*scaling*scaling);
                
                dist_const = sqrt((xsqr/asqr)+(ysqr/bsqr)+(zsqr/csqr));
                
                
                //Check if too close
                for (l=0; l<num_sofar; l++) {
                    
                    check_dist = rel_dist(pos[num_sofar],pos[l]);
                    
                    if (check_dist < min_dist){
                        pass = 0;
                    }
                    if(pot_toggle == 1){
                        if (dist_const > 0.99){
                            pass = 0;
                        }
                    }
                }
                
                if (pass == 1){
                    point_check = 1;
                }
                
                check_counter++;
                if (check_counter >= 1000){
                    check_counter = 0;
                    if(min_dist > 0.1){
                        min_dist = min_dist - (0.1*Ri[k]);
                        cout << "Shrinking Check Dist" << endl;
                    }
                }
                
            } //While loop
            
            //Use div_dir that was checked
            position[num_sofar][0] = position[k][0] + div_dir.x;
            position[num_sofar][1] = position[k][1] + div_dir.y;
            position[num_sofar][2] = position[k][2] + div_dir.z;
            
            Ri[num_sofar] = Ri[k];
            r0i[num_sofar] = r0i[k];
            a[num_sofar] = a[k];
            divR[num_sofar] = divR[k];
            cell_parent[num_sofar] = k;
            num_divs_i[k] = num_divs_i[k] + 1;
            
        } //If rule_counter == 1
        
        //Divide using T-division instead of perpendicular rules
        if(rule_counter == 2){
            
            if (k == 0 && num_divs_i[k] == 0){
                div_counter1 = 1;
            }
            else if (k == 0 && num_divs_i[k] == 1){
                div_counter1 = 2;
            }
            else if (k == 1 && num_divs_i[k] == 0){
                div_counter1 = 1;
            }
            else{
                div_counter1 = 3;
            }
            
            div_dir = Div_Direction2(div_counter1, div_noise);
            
            position[num_sofar][0] = position[k][0] + 1.*Ri[k]*div_dir.x;
            position[num_sofar][1] = position[k][1] + 1.*Ri[k]*div_dir.y;
            position[num_sofar][2] = position[k][2] + 1.*Ri[k]*div_dir.z;
            
            Ri[num_sofar] = Ri[k];
            r0i[num_sofar] = r0i[k];
            a[num_sofar] = a[k];
            divR[num_sofar] = divR[k];
            cell_parent[num_sofar] = k;
            num_divs_i[k] = num_divs_i[k] + 1;
            
        }
        
        num_sofar++;
        
        
    }//Loop over points before division
    
    // update particle number counter
    num_p = num_sofar;
    
    return num_p;
}

//The force funcition (derivative of the potential function)
double force_ij(double dist, double alpha, double beta, double a, double R, double r0){
    double force, part1, part1_prime, part2, part2_prime, part3, d1, d2;
    
    
    //Distance should always be positive
    double dist_sign;
    dist_sign = 1.0;
    if (dist < 0) {
        dist = -dist;
        dist_sign = -1.0;
    }
    
    //This should limit numbers from overflowing
    if (fabs(dist) < 0.2) {
        dist = 0.2;
    }
//    if (fabs(dist) > (2.*R)+0.2) {
//        dist = (2.*R)+0.2;
//    }
    
    //The force calculaton
    part1 = -(pow(2.*r0,alpha)*(beta+(exp((2.*(r0-R))/a)*(beta+(2.*r0/a)))))/(pow(dist,alpha));
    part1_prime = pow(2.*r0,alpha)*(-alpha)*((beta+(exp((2.*(r0-R))/a)*(beta+(2.*r0/a))))/(pow(dist,alpha+1.)));
    part2 = (pow(2.*r0,beta)*(alpha+(exp((2.*(r0-R))/a)*(alpha+(2.*r0/a)))))/(pow(dist,beta));
    part2_prime = pow(2.*r0,beta)*(beta)*((alpha+(exp((2.*(r0-R))/a)*(alpha+(2.*r0/a))))/(pow(dist,beta+1.)));

    part3 = 1./(1.+(exp((dist-(2.*R))/a)));
    
    d1 = (1./a)*(exp((dist-(2.*R))/a))*(1./(alpha-beta))*1./(pow(1.+(exp((dist-(2.*R))/a)),2.))*(part1+part2);
    d2 = (1./(alpha-beta))*part3*(part1_prime+part2_prime);
    
    force = -(d1+d2);
    
    //Limiting the absolute magnitude of the force
//    if (fabs(force) < 1.0E-7) {
//        force = 0.0;
//    }
    if (fabs(force) > 2.0E2) {
        force = 2.0E2;
    }
    
    return force;
}

//Use overlap between cells and use that to increase cell volume
void FindVolumeOverlap(double **position, vector<double> &Ri, vector<double> &divR, int max_num_p, int dim, vector<double> &vcap_array, const vector<vector<int> > &neighbor_array){
    
    int i,j,k;
    vector<double> cap_volumes(max_num_p);
    double dist;
    vector<coord> pt(max_num_p);
    
    for (k=0; k<max_num_p; k++) {
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
    
    //Initialize
    for (i = 0; i < max_num_p; i++) {
        cap_volumes[i] = 0.;
    }
    
    //Find cap volumes
    
    for (i = 0; i < max_num_p; i++) {
        for (j = 0; j < max_num_p; j++) {
            
            if(i == j){
                continue;
            }
            
            dist = rel_dist(pt[i],pt[j]);
            //cout << "dist " << dist << endl;
            
            if (dist > (divR[i]+divR[j])) {
                continue;
            }
            
            if (dist < 0.01) {
                continue;
            }
            
            if (neighbor_array[i][j] == 0) {
                continue;
            }
            //cout << "Cap volume from particle " << j << " is "<< CapVolume(R,dist) << endl;
            cap_volumes[i] = cap_volumes[i] + CapVolume(divR[i],divR[j],dist);
            //printf("Point %i, Cap%i: %.5e \n",i,j,CapVolume(divR[i],divR[j],dist) );
            
        }
        
        //Find rescaled R and save in Ri
        Ri[i] = pow((divR[i]*divR[i]*divR[i])+(3./(4.*PI))*cap_volumes[i],(1./3.));
        vcap_array[i] = cap_volumes[i];
        
//                printf("Point %i: New radius %.5e \n",i, Ri[i] );
//                printf("Point %i: Total Cap volume %.5e \n",i, vcap_array[i] );
//                printf("Point %i: Original volume %.5e \n",i, ((4.*PI)/3.)*pow(divR[i],3.) );
//                printf("Point %i: Cap volume percentage %.5e \n",i, vcap_array[i]*100./(((4.*PI)/3.)*pow(divR[i],3.)) );
    }
    
    
    return;
}

//Determine volume of overlap between two spheres of radius R1 and R2 at distance dist
double CapVolume(double R1, double R2, double dist){
    
    double cap_volume, a, h1, x1, s, area;
    int R1_obtuse = 0;
    int R2_obtuse = 0;
    x1 = 0;
    a = 0;
    s = 0;
    area = 0;
    
    if( ((R1*R1) + (dist*dist)) < (R2*R2) ){ //Check if R1-dist angle is obtuse
        R1_obtuse = 1;
        //Cell 1 cap should be hemisphere
        h1 = R1;
    }
    
    else if( ((R2*R2) + (dist*dist)) < (R1*R1) ){ //Check if R2-dist angle is obtuse
        R2_obtuse = 1;
        //Cell 1 cap should be cut at cell 2 position, h = r1 - dist
        h1 = R1 - dist;
    }
    else{
        //Use Heron's formula to determine the area of the triangle made by R1, R2 and dist as base
        s = (R1 + R2 + dist)/2.; //semi-perimeter
        area = sqrt(s*(s-R1)*(s-R2)*(s-dist));
        
//        cout << "R1 " << R1 << endl;
//        cout << "R2 " << R2 << endl;
//        cout << "s " << s << endl;
//        cout << "area " << area << endl;
        
        //Altitude which we call a can be found by using area = 1/2 bh with base = dist
        a = (2.*area)/dist;
        
        //We use the Pythagorean theorem to find x1 and x2 which add to dist
        if( ((R1*R1)-(a*a)) < 0.){
            cout << "Distance error in cap volume" << endl;
        }
        
        x1 = sqrt( ((R1*R1)-(a*a)) );
        h1 = R1 - x1;
    }
//    printf("x %.5e \n", x1 );
//    printf("h %.5e \n",h1 );
//    printf("a %.5e \n",a );
    
    cap_volume = (PI/3.)*h1*h1*((3.*R1)-h1);
    
    return cap_volume;
}

//Use Voro++ to find cell neighbors and areas of contact
void Analyze_Voro_Container(container_poly &con1, vector<vector<int> > &neighbor_array, vector<vector<double> > &area_array){
    
    voronoicell_neighbor d;
    c_loop_all c2(con1);
    int id, neighbor_loop;
    vector<int> neigh;
    vector<double> face_areas;
    
    if(c2.start()) do if(con1.compute_cell(d,c2)) {
        id=c2.pid();
        // Gather information about the computed Voronoi cell
        d.neighbors(neigh);
        d.face_areas(face_areas);
        //volume = c.volume();
        //cout << volume << endl;
        //cout << "Particle: " << id << endl;
        
        if (neigh.size()!=face_areas.size()){
            cout << "# of Neighbors != Areas" << endl;
        }
        
        for (neighbor_loop=0; neighbor_loop<neigh.size();neighbor_loop++){
            //neighbor_array[particle_counter][neighbor_loop] = 0;
            if (neigh[neighbor_loop]!=0) {
                //cout << "Neighbor " << neigh[neighbor_loop] << " Side Area ";
                //cout << face_areas[neighbor_loop] << " Face index " << neighbor_loop << endl;
                
                if (id != (neigh[neighbor_loop])) {
                    //cout << "Id " << id << " Neigh " << neigh[neighbor_loop] << endl;
                    
                    neighbor_array[id-1][neigh[neighbor_loop]-1] = 1;
                    area_array[id-1][neigh[neighbor_loop]-1] = face_areas[neighbor_loop];
                }
            }
        }
        //cout << endl;
        
    } while (c2.inc());
    
    return;
}

//Remove pair from neighbor and area arrays if area is too small or neighbor pair goes in only one direction
void CleanUp_Arrays(int max_num_p, vector<vector<int> > &neighbor_array, vector<vector<double> > &area_array){
    
    int neigh_a, neigh_b;
    double area_a, area_b;
    
    //Clean up the area and adjacency matrices
    for (int outiter = 0; outiter < max_num_p; outiter++) {
        for (int initer = 0; initer < max_num_p; initer++) {
            
            area_a = area_array[outiter][initer];
            area_b = area_array[initer][outiter];
            
            neigh_a = neighbor_array[outiter][initer];
            neigh_b = neighbor_array[initer][outiter];
            
            if (neigh_a != neigh_b) {
                neighbor_array[outiter][initer] = 0;
                neighbor_array[initer][outiter] = 0;
            }
            
            if (area_a < 0.01 || area_b < 0.01) {
                neighbor_array[outiter][initer] = 0;
                neighbor_array[initer][outiter] = 0;
            }
            
            
        }
    }

    
    return;
}

//Determine if cells have reached the tethrahedron equilibrium
int Array_Is_Equil(double **position, int dim, int max_num_p, vector<vector<int> > &neighbor_array, vector<double> &Ri){
    int is_equil = 0;
    int edge_counter = 0;
    int dist_counter = 0;
    
    vector<coord> pt(max_num_p);
    
    //Initialize points
    for (int k=0; k<max_num_p; k++) {
        pt[k] = zero_vector(pt[k]); //Initialize to zero
        for (int i=0; i<dim; i++) {
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

    
    for (int i = 0; i < max_num_p; i++) {
        for (int j = 0; j < max_num_p; j++) {
            
            if(rel_dist(pt[i],pt[j]) < (Ri[i]+Ri[j])){
                dist_counter += 1;
            }
            edge_counter += neighbor_array[i][j];
            
        }
    }
    
    if(edge_counter == 12 && dist_counter == 16){
        is_equil = 1;
    }
    if(edge_counter > 12){
        cout << "Too many 1's in your array!" << endl;
    }
    
    return is_equil;
}

//Checks if time_left has gone to zero or below
int Time_To_Divide(vector<double> &time_left, int num_p){
    
    int divide_indicator = 0;
    
    for(int i = 0; i < num_p; i++){
        if(time_left[i] <= 0.){
            divide_indicator = 1;
        }
    }
    
    return divide_indicator;
}

//Decrements entry in time_left vector by delta_t
void Increment_Cell_Clock(vector<double> &time_left, int num_p, double delta_t){
    
    for(int i = 0; i < num_p; i++){
        time_left[i] = time_left[i] - delta_t;
    }
    
}
