// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#include "voro++.hh"
using namespace voro;

// Set up constants for the container geometry
const double a_const=1,b_const=1,c_const=1;

// Set up constants for the container geometry
const double x_min=-a_const-0.5,x_max=a_const+0.5;
const double y_min=-b_const-0.5,y_max=b_const+0.5;
const double z_min=-c_const-0.5,z_max=c_const+0.5;
const double cvol=3.1415;

// Set up the number of blocks that the container is divided into
const int n_x=6,n_y=6,n_z=6;

// Set the number of particles that are going to be randomly introduced
const int particles=10;

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

//// This class creates a custom toroidal wall object that is centered on the
//// origin and is aligned with the xy plane. It is derived from the pure virtual
//// "wall" class. The "wall" class contains virtual functions for cutting the
//// Voronoi cell in response to a wall, and for telling whether a given point is
//// inside the wall or not. In this derived class, specific implementations of
//// these functions are given.
//class wall_ellipsoid : public wall {
//public:
//    
//    // The wall constructor initializes constants for the major and
//    // minor axes of the torus. It also initializes the wall ID
//    // number that is used when the plane cuts are made. This is
//    // only tracked with the voronoicell_neighbor class and is
//    // ignored otherwise. It can be omitted, and then an arbitrary
//    // value of -99 is used.
//    wall_ellipsoid(double ia,double ib, double ic,int iw_id=-99)
//    : w_id(iw_id), a(ia), b(ib), c(ic) {};
//    
//    // This returns true if a given vector is inside the torus, and
//    // false if it is outside. For the current example, this
//    // routine is not needed, but in general it would be, for use
//    // with the point_inside() routine in the container class.
//    bool point_inside(double x,double y,double z) {
//        double temp=((x*x)/(a*a))+((y*y)/(b*b))+((z*z)/(c*c));
//        return temp<1.;
//    }
//    
//    // This template takes a reference to a voronoicell or
//    // voronoicell_neighbor object for a particle at a vector
//    // (x,y,z), and makes a plane cut to to the object to account
//    // for the toroidal wall
//    template<class vc_class>
//    inline bool cut_cell_base(vc_class &c,double x,double y,double z) {
//        double orad=sqrt(x*x+y*y);
//        double odis=orad-mjr;
//        double ot=odis*odis+z*z;
//        
//        // Unless the particle is within 1% of the major
//        // radius, then a plane cut is made
//        if(ot>0.01*mnr) {
//            ot=2*mnr/sqrt(ot)-2;
//            z*=ot;
//            odis*=ot/orad;
//            x*=odis;
//            y*=odis;
//            return c.nplane(x,y,z,w_id);
//        }
//        return true;
//    }
//    
//    // These virtual functions are called during the cell
//    // computation in the container class. They call instances of
//    // the template given above.
//    bool cut_cell(voronoicell &c,double x,
//                  double y,double z) {return cut_cell_base(c,x,y,z);}
//    bool cut_cell(voronoicell_neighbor &c,double x,
//                  double y,double z) {return cut_cell_base(c,x,y,z);}
//private:
//    // The ID number associated with the wall
//    const int w_id;
//    // The major radius of the torus
//    const double a;
//    // The minor radius of the torus
//    const double b;
//    // The minor radius of the torus
//    const double c;
//};

int main() {
	int i;
	double x,y,z;

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block
	container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);
    // Add a cylindrical wall to the container
    wall_sphere sph(0,0,0,1,10);
    con.add_wall(sph);

	// Randomly add particles into the container
	for(i=0;i<particles;i++) {
        x = 10;y = 10;z=10;
        while(x*x+y*y+z*z>1){
            x=x_min+rnd()*(x_max-x_min);
            y=y_min+rnd()*(y_max-y_min);
            z=z_min+rnd()*(z_max-z_min);
        }
		con.put(i,x,y,z);
	}

	// Sum up the volumes, and check that this matches the container volume
	double vvol=con.sum_cell_volumes();
	printf("Container volume : %g\n"
	       "Voronoi volume   : %g\n"
	       "Difference       : %g\n",cvol,vvol,vvol-cvol);

	// Output the particle positions in gnuplot format
	con.draw_particles("test_p.gnu");

	// Output the Voronoi cells in gnuplot format
	con.draw_cells_gnuplot("test_v.gnu");
}
