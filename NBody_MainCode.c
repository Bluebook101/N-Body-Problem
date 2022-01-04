/*
-----------------------------------------------------------------------------------------------------------------------------------
PHYM004
HW2 - The N-Body Problem
Date: 11/12/20
Author: Joe Salkeld 670008743
-----------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------WELCOME-------------------------------------------------------------

This piece of code is written to read in data from a text file of N bodies and evolve the system using two different methods of
numerical analysis. The format of the text file should tab seaprated and on each line will contain the name of the body followed
by its mass, its 3 position coordinates and its 3 velocity coordinates all tab separated. An example of the file format for the 
Sun-Earth-Moon system is below:

# Name	Mass		x			y		z		Vx		Vy			Vz	
Earth	5.9726e24	1.496e11	0.0		0.0		0.0		2.98e4		0.0
Sun		1.9885e30	0.0			0.0		0.0		0.0		0.0			0.0
Moon	7.342e22	1.503844e11	0.0		0.0		0.0		3.0312e4	0.0

To compile the code, the following command was used:

$ gcc NBody.c -Wall -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas -lm -o nbody.exe

After compiling the code, the user then has two different options of numerical solvers. These are called as follows:

-v x		using the Velocity-Verlet solver with a timestep of x days
-r x		using the RK4 solver with a maximum timestep of x days

And would be run like this:

$ ./nbody.exe -v 1

The operations in this code are called using a getoptlong() function in the main, adapted from Dr C.D.H Williams' code.

The code then outputs a text file with updated positions for the bodies for every timestep for N_POINTS, as defined at the top of 
the code. This is then plotted by Gnuplot and the results are displayed.

The code also calculates the orbital period and the error in the total energy of the system, and outputs these values to the 
console.

The code comes with three example text files along with 3 script files to plot them:

Sun_Earth_Moon.txt and SEM.script for the Sun-Earth-Moon system
Solar_System.txt and SolarSystem.script for the Solar System
Figure_of_8.txt and 3Body.script for the stable figure of eight loop

N.b. In order for the figure of eight to plot, G_CONST must be set to 1 and the timestep must be 0.0000001 or smaller

-----------------------------------------------------------------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <getopt.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>


#define MAX_LINE_SIZE 	250 //maximum length of a line
#define ITEMS_LINE    	8 //maximum number of values per line
#define MAX_NAME_SIZE   32 //maximum size of body name
#define DAYS_SECS 86400 //constant of number of seconds in a day
#define G_CONST 6.67e-11 //gravitational constant


#define INPUT_FILE "Sun_Earth_Moon.txt" //input file which contains the initial conditions
#define GNUPLOT_SCRIPT "./SEM.script" //script file which contains plotting choices
#define DATAFILE   "Evolution.dat" //name of the output datafile 
#define N_POINTS 10000 //number of points the code runs for


static int maxbodies = 0; //global variable used for the number of bodies found in the code

typedef enum coords { X, Y, Z, N_COORDS } Coords;

/*
This structure contains variables for storing all the initial information for each body 
*/
typedef struct body {
  char name[MAX_NAME_SIZE];
  double mass;
  double r[N_COORDS]; // displacement 
  double v[N_COORDS]; // velocity 
  double a[N_COORDS]; // acceleration 
} Body;

/*
This structure contains variables for storing all the information for orbital period calculations 
*/
typedef struct period {
	double yprev; //previous y-displacement value
	double ynew; //new y-displacement value
	double av_period; //stores the average orbital period
	double laps; //keeps track of how many times the body has completed an orbit
	double orbit; //stores how many timesteps the most recent orbit took
} Period;

/*
This structure contains variables for storing all the information for energy calculations 
*/
typedef struct sys_energy{
	double Etotal;
	double Ek;
	double Ep;
	double Eprev; //previous timestep's total energy
	double av_error; //stores the average error for the total energy
	double percerr; //stores the current error for the total energy
	double day; //keeps track of how many timesteps into the code the function is 
} Energy;

/*
-----------------------------------------------------------------------------------------------------------------------------------
Read File Function - Adapted from CDH Williams ReadOrbits.c file, this function reads the input file to calculate the number
of bodies in the file, and then allocated enough memory accordingly. It then loops through these bodies storing all the data
in the structure Body. Error checks are in place to make sure the code reads in the data correctly. 
-----------------------------------------------------------------------------------------------------------------------------------
*/
static struct body * readfile(){
	
	char line[MAX_LINE_SIZE];
	char nameBuf[MAX_LINE_SIZE];
	FILE *input = fopen( INPUT_FILE, "r" );

	if(!input) {
		fprintf(stderr, "Error: Could not open file '%s'.\n",INPUT_FILE);
		exit(1);
	}
	
	//looping through the file to count the number of bodies
	while(fgets(line, MAX_LINE_SIZE, input)){
		if(line[0] != '#'){
			maxbodies += 1;
		}
	}
	
	rewind(input);
  
	Body *bodies = malloc(maxbodies * sizeof(Body));

	int bodyN = 0; // Number of bodies successfully read
	
	while (bodyN < maxbodies && fgets(line, MAX_LINE_SIZE, input)){
		if(line[0] != '#'){
			int nFound = sscanf(line,"%s %lg %lg %lg %lg %lg %lg %lg", nameBuf, &bodies[bodyN].mass, &bodies[bodyN].r[X], 
			&bodies[bodyN].r[Y], &bodies[bodyN].r[Z], &bodies[bodyN].v[X], &bodies[bodyN].v[Y], &bodies[bodyN].v[Z]);
			
			if(nFound == ITEMS_LINE){
				strncpy(bodies[bodyN].name,nameBuf,MAX_NAME_SIZE);
				bodyN++;
			} 
			//Error check if the number of values on a line does not equal the expected amount
			else{
				fprintf(stderr, "Unknown format: %s\n",line);
				exit(1);
			}
		}
	}
	
	fclose(input);
	return bodies;
}

/*
-----------------------------------------------------------------------------------------------------------------------------------
Acceleration Calculation Function - This function takes in the values stored in the structure Body and uses them to calculate the
acceleration of the system. It then updates the values of the acceleration stored in the structure and passes the pointer to the 
structure back to the code.
-----------------------------------------------------------------------------------------------------------------------------------
*/
static struct body * acc(struct body *bodies){
	
	for (int x=0; x<maxbodies; x++){
		//all acceleration is initialised to zero
		bodies[x].a[X] = 0;
		bodies[x].a[Y] = 0;
		bodies[x].a[Z] = 0;
		for (int y=0; y<maxbodies; y++){
			if(y != x){
				
				//value of the modulus of the distance between the two bodies is stored
				double modr = sqrt(pow((bodies[x].r[X] - bodies[y].r[X]),2) + pow((bodies[x].r[Y] - bodies[y].r[Y]),2) + pow((bodies[x].r[Z] - bodies[y].r[Z]),2));
				
				double ax = -((bodies[x].r[X] - bodies[y].r[X]) * G_CONST * bodies[y].mass) / (pow(modr,3)); 
				
				double ay = -((bodies[x].r[Y] - bodies[y].r[Y]) * G_CONST * bodies[y].mass) / (pow(modr,3));
				
				double az = -((bodies[x].r[Z] - bodies[y].r[Z]) * G_CONST * bodies[y].mass) / (pow(modr,3));
				
				//total acceleration is updated with each body until all bodies have been looped through
				bodies[x].a[X] += ax;
				bodies[x].a[Y] += ay;
				bodies[x].a[Z] += az;
			}
		}	
	}
	
	return bodies;
}

/*
-----------------------------------------------------------------------------------------------------------------------------------
Total System Energy Calculation Function - This function takes in pointers to the structures of Body and Energy as well as the
current point number and uses them to calculate the total energy of the system, before returning a pointer to the structure Energy. 
It calculates the Kinetic and the Gravitational Potential Energy separately, then calculates the Total Energy before working out 
the average error between the Total Energy and the previous timestep's Total Energy. A running average of the percentage error in 
the Total Energy is kept track of, which is output at the end of the program to the console.
-----------------------------------------------------------------------------------------------------------------------------------
*/
static struct sys_energy * energy_func(struct body *bodies, struct sys_energy *energy, double point){
	
	energy->Etotal = 0; //initialising the total energy of the system
	
	for(int x=0; x<maxbodies; x++){
		
		energy[x].Ek = 0.5 * bodies[x].mass * (pow(bodies[x].v[X],2) + pow(bodies[x].v[Y],2) + pow(bodies[x].v[Z],2));
		
		energy[x].Ep = 0; //intialising the potential energy
		
		for (int y=0; y<maxbodies; y++){
			if(y != x){
				//value of the modulus of the distance between the two bodies is stored
				double modr = sqrt( pow((bodies[x].r[X] - bodies[y].r[X]),2) + pow((bodies[x].r[Y] - bodies[y].r[Y]),2) + pow((bodies[x].r[Z] - bodies[y].r[Z]),2) );
				
				//running total of the potential energy for body x
				energy[x].Ep = energy[x].Ep - (( G_CONST * bodies[x].mass * bodies[y].mass) / modr );
			}	
		}
		//Total energy. The potential energy is multiplied by 0.5 to account for counting each body twice in the loop
		energy->Etotal = energy->Etotal + 0.5 * energy[x].Ep + energy[x].Ek;
	}

	//initialising values if this is the first time the fucntion has been called
	if(point == 0){
		energy->Eprev = energy->Etotal;
		energy->av_error = 0;
	}
	
	//percentage error and average percentage error totals being calculated
	else{
		energy->percerr = sqrt(pow(((energy->Etotal - energy->Eprev) / (energy->Eprev)),2)) * 100;
		energy->Eprev = energy->Etotal;
		energy->av_error = ((energy->av_error *(point -1)) + energy->percerr) / point;
	}
	
	return energy;
}

/*
-----------------------------------------------------------------------------------------------------------------------------------
Orbital Period Function - This function takes in pointers to the structures Body and Period, the timestep input by the 
user at runtime and the current point number. The code checks if the y value of the displacement has passed its initial position, 
and increments a value if it has, storing the days it took to reach that point. The average orbital period is stored for bodies 
which pass their initial position more than once in the N_POINTS. It returns a pointer to the structure Period.
-----------------------------------------------------------------------------------------------------------------------------------
*/
static struct period * orbital(struct body *bodies, struct period *orbit_period, double timestep, double point){
	
	double time = timestep / DAYS_SECS; //converts the timestep into days.
	
	//initialising values if this is the first time the function has run
	if(point == 0){
		for(int x=0; x<maxbodies;x++){
			orbit_period[x].yprev = bodies[x].r[Y]; 
			orbit_period[x].av_period = 0;
			orbit_period[x].laps = 0;
		}
		return orbit_period;
	}
	
	//runs through each body and calculates whether it has completed and orbit and stores values if it has
	for(int x=0; x<maxbodies;x++){
		orbit_period[x].ynew = bodies[x].r[Y];
	
		if(orbit_period[x].ynew >= 0 && orbit_period[x].yprev < 0){
			orbit_period[x].laps += 1; //increments lap counter
			orbit_period[x].orbit = (time * point) / orbit_period[x].laps; //orbital period in days
			//updating the average orbital period
			orbit_period[x].av_period = ((orbit_period[x].av_period *(orbit_period[x].laps -1)) + orbit_period[x].orbit) / orbit_period[x].laps ;
			orbit_period[x].yprev = orbit_period[x].ynew;
		}
		else{
			orbit_period[x].yprev = orbit_period[x].ynew;
		}
	}
	
	return orbit_period;
}

/*
-----------------------------------------------------------------------------------------------------------------------------------
Velocity Verlet Method - This function takes in the value of the timestep as specified by the user and uses the Velocity-Verlet 
method of numerical analysis to evolve the system. It first reads the input file before beginning the analysis. The function 
calculates the acceleration of the bodies, before updating the velocity of each body a half timestep. It then updates the position 
of each body and recalculates the acceleration, before updating the velocity another half timestep. This process is looped for 
N_POINTS and at each point the position values are output to a separate text file to be plotted at the end. The values for the 
orbital period and energy of the system are output to the console at the end too.
-----------------------------------------------------------------------------------------------------------------------------------
*/
static void verlet(double timestep){
	
	//Reading in of input file and allocation of memory
	struct body *bodies = readfile(); 
	struct period *orbit_period = malloc(maxbodies * sizeof(Period)); 
	struct sys_energy *energy = malloc(maxbodies * sizeof(Energy));
  
	FILE *fdat = fopen (DATAFILE, "w");
	
	//Exit if the file cannot be opened
	if (!fdat) {
		exit(1);
	}
	
	fprintf(fdat, "# X\tY\tZ\n" );
	for(int x=0; x<maxbodies; x++){
			fprintf(fdat, "%.8e\t%.8e\t%.8e\t", bodies[x].r[X], bodies[x].r[Y], bodies[x].r[Z] ); //printing the initial positions to file
		}
	fprintf(fdat, "\n");
	
	//Initialising values for the orbital period and system energy
	orbit_period = orbital(bodies, orbit_period, timestep, 0);	
	energy = energy_func(bodies, energy, 0);
	
	//Initial acceleration calculation
	bodies = acc(bodies);	
	
	for (int d = 1; d <= N_POINTS; d++ ) { 
   
		//First step of velocity verlet: new half time step velocity

		for (int x=0 ; x<maxbodies; x++){
			for (int y=0; y<N_COORDS; y++){
				bodies[x].v[y] = bodies[x].v[y] + (0.5 * timestep * bodies[x].a[y]);
			}
		}		
   
		//Second step of velocity verlet: new positions
	
		for (int x=0 ; x<maxbodies; x++){
			for (int y=0; y<N_COORDS; y++){
				bodies[x].r[y] = bodies[x].r[y] + (bodies[x].v[y] * timestep);		
			}
		}	
   
		//Third step of velocity verlet: rederive aceleration
		
		bodies = acc(bodies);

		//Fourth step of velocity verlet: correct velocities accordingly
	
		for (int x=0 ; x<maxbodies; x++){
			for (int y=0; y<N_COORDS; y++){
				bodies[x].v[y] = bodies[x].v[y] + (0.5 * timestep * bodies[x].a[y]);
			}
		}
		
		//Energy of system 
		energy = energy_func(bodies, energy, d);
	
		//Orbital Period
		
		orbit_period = orbital(bodies, orbit_period, timestep, d);		
	
		//Printing new positions to file
		for(int x=0; x<maxbodies; x++){
			fprintf(fdat, "%.8e\t%.8e\t%.8e\t", bodies[x].r[X], bodies[x].r[Y], bodies[x].r[Z] );
		}
		fprintf(fdat,"\n");
	}
	
	
	//Outputting the orbital period of all bodies to the console
	for(int x=0; x<maxbodies ; x++){
		if(orbit_period[x].av_period == 0){
			printf("%s did not complete a full orbit\n", bodies[x].name);
		}
		else{
		printf("The average Orbital Period of %s is %lg\n", bodies[x].name, orbit_period[x].av_period);
		}
	}
	
	//Outputting the percentage error in Total energy to the console
	printf("The average percentage error in the Total Energy is %.3lg %%", energy->av_error);
	

	fclose(fdat);  
	free(bodies);   
	free(orbit_period);
	free(energy);
}

/*
-----------------------------------------------------------------------------------------------------------------------------------
ODE Solver Function - This function is part of the RK4 driver, and reads in the values for position and velocity
in a 1D array and then calculates the velocity and acceleration. The function returns a value GSL_SUCCESS if it solves the 
ODE correctly. It is based on the example code given in the GSL Library for ODE solvers.
-----------------------------------------------------------------------------------------------------------------------------------
*/
static int ode (double t, const double y[], double f[], void* params){
	
	(void)(t); //Avoids unused parameter warning 
	struct body *values = params;
  
	for(int i = 0; i <maxbodies ; i++){
	  
		f[i*6] = y[(i*6)+3]; //Vx
		f[(i*6)+1] = y[(i*6)+4]; //Vy
		f[(i*6)+2] = y[(i*6)+5]; //Vz
		f[(i*6)+3] = 0; //Ax
		f[(i*6)+4] = 0; //Ay
		f[(i*6)+5] = 0; //Az
	  
	//Loop to calculate total acceleration due to all bodies
		for (int j=0; j<maxbodies; j++){
			if(j != i){
				//value of the modulus of the distance between the two bodies is stored
				double modr = sqrt( pow((y[i*6] - y[j*6]),2) + pow((y[(i*6)+1] - y[(j*6)+1]),2) + pow((y[(i*6)+2] - y[(j*6)+2]),2) );
				
				//equation for velocity derivative (acceleration) is updated.
				f[(i*6)+3] += -( G_CONST * (values+j)->mass * (y[i*6] - y[j*6])) / (pow(modr,3));
				f[(i*6)+4] += -( G_CONST * (values+j)->mass * (y[(i*6)+1] - y[(j*6)+1])) / (pow(modr,3));
				f[(i*6)+5] += -( G_CONST * (values+j)->mass * (y[(i*6)+2] - y[(j*6)+2])) / (pow(modr,3));
			}
		}
	}
	
	return GSL_SUCCESS;
}

/*
-----------------------------------------------------------------------------------------------------------------------------------
Runge Kutta 4 Function - This function is part of the GSL standard library, using the ODE solver, specifically the RK4 method, to
create a series of ODE for all the bodies and solving them to evolve the system. It is called from the getoptlong() function in the
main and reads in the textfile of bodies before rewriting them as a one dimensional array which contains [x,y,z,vx,vy,vz] for all 
bodies one after the other. It then uses the GSL library to solve this equation with the function above, and outputs the updated
positions to the text file which is used by Gnuplot to plot the orbit. The energy of the system and the orbital period are also 
calculated and output to the console. 
-----------------------------------------------------------------------------------------------------------------------------------
*/
static void runge (double usertime){
	
	void *nul = NULL;
	double t=0;  
	double hstart = usertime * 0.1;

	//Reading in the input file and allocating memory accordingly
	struct body *bodies = readfile();
	struct period *orbit_period = malloc(maxbodies * sizeof(Period));
	struct sys_energy *energy = malloc(maxbodies * sizeof(Energy));

	FILE *fdat = fopen (DATAFILE, "w");
	
	//Exits if the file cannot be opened
	if (!fdat) {
		exit(1);
	}
	
	
	fprintf(fdat, "# X\tY\tZ\n" );
	for(int x=0; x<maxbodies; x++){
		fprintf(fdat, "%.8e\t%.8e\t%.8e\t", bodies[x].r[X], bodies[x].r[Y], bodies[x].r[Z]); //printing the initial positions to file
	}
	fprintf(fdat, "\n");
	
	//Initialising values for the orbital period and system energy
	orbit_period = orbital(bodies, orbit_period, usertime, 0);	
	energy = energy_func(bodies, energy, 0);
	
	for (int i = 1; i <= N_POINTS; i++){
		
		int body_values = 2*N_COORDS;//defining the number of values for each body
		double bodydata[body_values*maxbodies]; //defining a 1D array for the data to be rewritten into
		
		/*
		This loops through the two separate arrays of position and velocity for each body and stores them in the one dimensional 
		array which has been defined, so that the array contains [Rx1,Ry1,Rz1,Vx1,Vy1,Vz1,Rx2,Ry2,Rz2,Vx2,Vy2,Vz2,...] for a 
		system in three coordinates.
		*/
		for (int x = 0; x<maxbodies ; x++){
			for(int z = 0; z<=N_COORDS; z+=N_COORDS){
				for (int y = 0 ; y<N_COORDS ; y++){
					if (z==0){
						bodydata[(body_values*x)+y+z] = bodies[x].r[y];
					}
					else{
						bodydata[(body_values*x)+y+z] = bodies[x].v[y];
					}
				}
			}
		}
  
		gsl_odeiv2_system sys = {ode, nul, body_values*maxbodies, bodies}; //ODE system is setup 
		
		//RK4 method is selected with max errors of 1e-8 and an hstart value of 0.1*usertime
		gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4, hstart, 1e-8, 1e-8); 
  
		double ti = i * usertime; //user defined timestep to solve between
	
		int status = gsl_odeiv2_driver_apply (d, &t, ti, bodydata); //ODE system is solved
	 
		//Error check if the solver failed
		if (status != GSL_SUCCESS){
			printf ("error, return value=%d\n", status);
			break;
        } 
		
		/*
		Same loop as before the solver but reversed to allow the 1D array to be stored into the two separate arrays in
		the Body structure, which allows for orbital period and energy calculation, as well as printing the updated 
		positions and updating values to be read in when the loop repeats.
		*/
		for (int x = 0; x<maxbodies ; x++){
			for(int z = 0; z<=N_COORDS; z+=N_COORDS){
				for (int y = 0 ; y<N_COORDS ; y++){
					if (z==0){
						bodies[x].r[y] = bodydata[(body_values*x)+y+z];
					}
					else{
						bodies[x].v[y] = bodydata[(body_values*x)+y+z] ;
					}
				}
			}
		}		
		
		//driver is freed before the loop begins again
		gsl_odeiv2_driver_free (d);
		
		//Energy of system 
		energy = energy_func(bodies, energy, i);
		
		//Orbital Period
		orbit_period = orbital(bodies, orbit_period, usertime, i);	
		
		
		//Printing new positions to file
		for(int x=0;x<maxbodies;x++){
			fprintf(fdat, "%.8e\t%.8e\t%.8e\t", bodies[x].r[X], bodies[x].r[Y], bodies[x].r[Z] );
		}
		fprintf(fdat,"\n");
	}
	
	//Outputting the orbital period of all bodies to the console
	for(int x=0; x<maxbodies ; x++){
		if(orbit_period[x].av_period == 0){
			printf("%s did not complete a full orbit\n", bodies[x].name);
		}
		else{
		printf("The average Orbital Period of %s is %lg\n", bodies[x].name, orbit_period[x].av_period);
		}
	}
	
	//Outputting the percentage error in Total energy to the console
	printf("The average percentage error in the Total Energy is %.3lg %%", energy->av_error);

	fclose(fdat);  
	free(bodies); 
	free(orbit_period);
	free(energy);
}


/*
-----------------------------------------------------------------------------------------------------------------------------------
Main Section - This is the main function which contains the getoptlong() function allowing the choice of solver to be input by the 
user at runtime. It is adapted from Dr C.D.H Williams' code used in mat_gen.c, as well as the online GNU C library. The timestep 
the user inputs is in days, and the code converts this to seconds to be used by the solvers.
-----------------------------------------------------------------------------------------------------------------------------------
*/
int main(int argc,char ** argv) {
	
    static int normal_flg;
	char command[PATH_MAX];
	double usertime = 0; //User inputs timestep in days//
	double time = 0; 
	
//getoptlong() function and options are adapted from Dr. C.D.H. Williams code "mat_gen.c"

    while (1) {
        static struct option long_options[] = {
            // These options set flags. 
            {"verbose", no_argument,      &normal_flg, 1},
            {"normal", no_argument,       &normal_flg, 0},
            // These options don’t set a flag, they are distinguished by their indices. 
            {"verlet", required_argument,  0, 'v'},
			{"runge-kutta", required_argument,  0, 'r'},
            {0, 0, 0, 0}
        };
        
        // getopt_long needs somewhere to store its option index. 
        int option_index = 0;
        
        int c = getopt_long( argc, argv, "r:v:", long_options, &option_index );
        
        // End of options is signalled with '-1' 
        if (c == -1)
            break;
        
        switch (c) {
            case 0:
                // If this option sets a flag we have nothing to do. 
                if (long_options[option_index].flag != 0) {
                    break;
                }
                printf ("option %s", long_options[option_index].name);
                if (optarg)
                    printf (" with arg %s", optarg);
                printf ("\n");
                break;
				
            case 'v':
                usertime = atof(optarg); //value for timestep in days read in
				time = usertime * DAYS_SECS; //converts from days to seconds
				verlet(time);
				break;	
			
			case 'r':
				usertime = atof(optarg); //value for the maximum timestep
				time = usertime * DAYS_SECS; //converts from days to seconds
				runge(time);
				break;
        }
    }
    
    // Prints any remaining command line arguments (not options).
    if (optind < argc) {
        fprintf (stderr, "Error: Unrecognised arguments: ");
        while (optind < argc) {
            fprintf (stderr, "%s ", argv[optind++]);
        }
        fprintf (stderr, "\n");
    }
	
	//Command to invoke Gnuplot and plot the output data.
	snprintf(command, sizeof(command), "%s", GNUPLOT_SCRIPT );
	system(command);
	
	return(0);
}