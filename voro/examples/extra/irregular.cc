// Irregular packing example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#include "voro++.hh"
#include <vector>
#include <iostream>
#include <string>
using namespace std;
using namespace voro;

// Set up constants for the container geometry
const double x_min=-6,x_max=6;
const double y_min=-6,y_max=6;
const double z_min=-3,z_max=9;

// Golden ratio constants
const double Phi=0.5*(1+sqrt(5.0));
const double phi=0.5*(1-sqrt(5.0));

// Set up the number of blocks that the container is divided
// into.
const int n_x=5,n_y=5,n_z=5;


// Create a wall class that, whenever called, will replace the Voronoi cell
// with a prescribed shape, in this case a dodecahedron
class wall_initial_shape : public wall {
	public:
		wall_initial_shape() {

			// Create a dodecahedron
			v.init(-2,2,-2,2,-2,2);
			v.plane(0,Phi,1);v.plane(0,-Phi,1);v.plane(0,Phi,-1);
			v.plane(0,-Phi,-1);v.plane(1,0,Phi);v.plane(-1,0,Phi);
			v.plane(1,0,-Phi);v.plane(-1,0,-Phi);v.plane(Phi,1,0);
			v.plane(-Phi,1,0);v.plane(Phi,-1,0);v.plane(-Phi,-1,0);
		};
		bool point_inside(double x,double y,double z) {return true;}
		bool cut_cell(voronoicell &c,double x,double y,double z) {

			// Set the cell to be equal to the dodecahedron
			c=v;
			return true;
		}
		bool cut_cell(voronoicell_neighbor &c,double x,double y,double z) {

			// Set the cell to be equal to the dodecahedron
			c=v;
			return true;
		}
	private:
		voronoicell v;
};

int main() {
    
    voronoicell_neighbor c;
    vector<int> neigh;
    double x,y,z;
    int id, particle_counter;
    
    char p_buffer[40], v_buffer[40];

	// Create a container with the geometry given above. This is bigger
	// than the particle packing itself.
	container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);

	// Create the "initial shape" wall class and add it to the container
	wall_initial_shape(wis);
	con.add_wall(wis);

	// Import the irregular particle packing
	//con.import("pack_irregular");
    for (int j = 0; j < 2; j++) {
        
        particle_counter = 1;
        
        if(j == 0){
            snprintf(p_buffer, sizeof(char) * 40, "irregular_p1.gnu");
            snprintf(v_buffer, sizeof(char) * 40, "irregular_v1.gnu");
        }
        if(j == 1){
            snprintf(p_buffer, sizeof(char) * 40, "irregular_p2.gnu");
            snprintf(v_buffer, sizeof(char) * 40, "irregular_v2.gnu");
        }
        
        
        con.put(1,1.,0.,0.);
        con.put(2,-1.,0.,0.);
        con.put(3,0.,((double)j)+1.,0.);
        //con.put(4,0.,0.1,0.);
        
    //    con.put(1,-1,0.,0.);
    //    con.put(2,0.,0.5,0.);
    //    con.put(3,0.,-0.5,0.);
    //    con.put(4,1.,0.,0.);

        
        // Save the particles and Voronoi cells in POV-Ray format (couldn't get POV to work so using gnu)
        con.draw_particles(p_buffer);
        con.draw_cells_gnuplot(v_buffer);
        
        // Loop over all particles in the container and compute each Voronoi
        // cell
        c_loop_all cl(con);
        if(cl.start()) do if(con.compute_cell(c,cl)) {
            cl.pos(x,y,z);id=cl.pid();
            // Gather information about the computed Voronoi cell
            c.neighbors(neigh);
            
            cout << "Particle: " << particle_counter << " Neighbor";
            for (int i=0; i<neigh.size();i++){
                if (neigh[i]!=0) {
                    cout << " " << neigh[i];
                }
            }
            cout << endl;
            particle_counter++;
        } while (cl.inc());
        
        con.clear();
        
    }
}
