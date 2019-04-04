#pragma once
# include "simdefs.hh"

// Set up constants for the Voro++ container geometry
const double x_min=-20,x_max=20;
const double y_min=-20,y_max=20;
const double z_min=-20,z_max=20;

// Golden ratio constants for making dodecahedron
const double Phi=0.5*(1+sqrt(5.0));
const double phi=0.5*(1-sqrt(5.0));

// Set up the number of blocks that the container is divided
// into.
const int n_x=1,n_y=1,n_z=1;

// Create a wall class that, whenever called, will replace the Voronoi cell
// with a prescribed shape, in this case a dodecahedron
class wall_initial_shape : public wall {
public:
    wall_initial_shape() {
        
        counter = 0;
        num_p = 1;
        R.reserve(10);
        
    };
    bool point_inside(double x,double y,double z) {return true;}
    bool cut_cell(voronoicell &c,double x,double y,double z) {
        
        // Set the cell to be equal to the dodecahedron
        double scaling;
        //cout <<"Counter "<< counter << endl;
        //cout <<"R "<< R[counter] << endl;
        scaling = R[counter]/0.4755; //0.4755 is Ri for unscaled dodecaheron
        // Create a dodecahedron
        v.init(-3,3,-3,3,-3,3);
        //Should this be bigger for bigger R?
        v.plane(0,scaling*(Phi/2.),scaling*(1/2.));v.plane(0,scaling*(-Phi/2.),scaling*(1/2.));v.plane(0,scaling*(Phi/2.),scaling*(-1/2.));
        v.plane(0,scaling*(-Phi/2.),scaling*(-1/2.));v.plane(scaling*(1/2.),0,scaling*(Phi/2.));v.plane(scaling*(-1/2.),0,scaling*(Phi/2.));
        v.plane(scaling*(1/2.),0,scaling*(-Phi/2.));v.plane(scaling*(-1/2.),0,scaling*(-Phi/2.));v.plane(scaling*(Phi/2.),scaling*(1/2.),0);
        v.plane(scaling*(-Phi/2.),scaling*(1/2.),0);v.plane(scaling*(Phi/2.),scaling*(-1/2.),0);v.plane(scaling*(-Phi/2.),scaling*(-1/2.),0);
        c=v;
        counter++;
        if (counter >= (num_p)){
            counter = 0;
        }
        return true;
    }
    bool cut_cell(voronoicell_neighbor &c,double x,double y,double z) {
        
        // Set the cell to be equal to the dodecahedron
        double scaling;
        //cout <<"Counter "<< counter << endl;
        //cout <<"R "<< R[counter] << endl;
        scaling = R[counter]/0.4755; //0.4755 is Ri for unscaled dodecaheron
        // Create a dodecahedron
        v.init(-3,3,-3,3,-3,3);
        //Should this be bigger for bigger R?
        v.plane(0,scaling*(Phi/2.),scaling*(1/2.));v.plane(0,scaling*(-Phi/2.),scaling*(1/2.));v.plane(0,scaling*(Phi/2.),scaling*(-1/2.));
        v.plane(0,scaling*(-Phi/2.),scaling*(-1/2.));v.plane(scaling*(1/2.),0,scaling*(Phi/2.));v.plane(scaling*(-1/2.),0,scaling*(Phi/2.));
        v.plane(scaling*(1/2.),0,scaling*(-Phi/2.));v.plane(scaling*(-1/2.),0,scaling*(-Phi/2.));v.plane(scaling*(Phi/2.),scaling*(1/2.),0);
        v.plane(scaling*(-Phi/2.),scaling*(1/2.),0);v.plane(scaling*(Phi/2.),scaling*(-1/2.),0);v.plane(scaling*(-Phi/2.),scaling*(-1/2.),0);
        c=v;
        counter++;
        if (counter >= (num_p)){
            counter = 0;
        }
        return true;
    }
    
    //Allow input radius to scale dodecahedron
    bool set_R(int index, double val){
        
        R[index] = val;
        
        return true;
    }
    
    bool init_values(int num_val, double R_init){
        num_p = num_val;
        
        //Reset counter
        counter = 0;
        //Set all the radii to 1
        for (int i=0; i<num_p; i++) {
            R.push_back(R_init);
        }
        
        return true;
    }
    
private:
    voronoicell v;
    int counter;
    int num_p;
    vector<double> R;
};

// Want to add some of the container stuff I do in main in here
