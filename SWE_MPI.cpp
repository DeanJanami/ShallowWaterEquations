///////////////////////////////////////////////
//
// Dean Pakravan: 757389
// Qijie Li:      927249
//
// Assignment 1: AHPC
//
// Method:		Finite Difference Method with sixth order central differences
//				and Fourth Order Runge-Kutta Method, parallelized with MPI
//
// Compilation:	mpicxx SWE_MPI.cpp -o SWE_MPI
//
// Execution:	mpirun  -n 4 SWE_MPI
//              OR
//              srun --ntasks=4 --nodes=1 --ntasks-per-node=4 ./SWE_MPI
//              OR
//              srun -n 4 -N 1 --ntasks-per-node=4 ./SWE_MPI
//
/////////////////////////////////////////////////////////////////////////////

#include <cstring>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <math.h>

using namespace std;

// Simulation Parameters
const double x_min   =  0.00;
const double x_max   = 100.00;
const double y_min   =  0.00;
const double y_max   = 100.00;
const double t_min   =  0.00;
const double t_max   = 100.00;
const double Delta_x =  1.00;
const double Delta_y =  1.00;
const double Delta_t =  0.10;
const int	 N_x     =  (x_max-x_min)/Delta_x;
const int	 N_y     =  (y_max-y_min)/Delta_y;
const int	 N_t     =  (t_max-t_min)/Delta_t;
const double g       =  9.81;
// MPI param
const int    N_D     = 2;
const int    X       = 0;
const int    Y       = 1;
const int    numElementsPerBlock = 1; 

// Function declarations
void f(double** kx, double** ky, double** kh, double** Vx, double** Vy, double** h, int myN_x, int myN_y,
         int iStart, int iEnd, int jStart, int jEnd);
void write(MPI_File& file, double** h, int l, int myN_x, int myN_y);
void exchange(double** phi, int myN_x, int myN_y, MPI_Status status, MPI_Comm Comm2D, MPI_Datatype stridetype, 
                int leftNeighbor, int rightNeighbor, int bottomNeighbor, int topNeighbor);

int main(int argc, char** argv)
{
    // ID for each processor
    int          myID            = 0;
    int          N_Procs         = 0;      // Number of processors
    double       wtime           = 0.0;	  // Record the starting time

	MPI_Init(&argc, &argv);								
	MPI_Comm_rank(MPI_COMM_WORLD, &myID);				// Find out this process' rank
	MPI_Comm_size(MPI_COMM_WORLD, &N_Procs);            // Find out how many processes in total

    // Pre-allocate datatypes
    fstream      file;
    MPI_Status   status;
	MPI_Datatype stridetype;
	MPI_Comm     Comm2D;

	int          N               = sqrt(N_Procs);
	int          dimensions[N_D] = {N, N};
	int          isPeriodic[N_D] = {0, 0};
	int          myCoords[N_D]   = {0, 0};

    // Compute the number of grid points
	int          myN_x           = N_x/(N);
	int          myN_y           = N_y/(N);

    int          myiStart        = 0;
	int          myiEnd          = 0;
	int          myjStart        = 0;
	int          myjEnd          = 0;

    // Indexes
    double       t               = 0;
	int          k               = 0;
	int          j               = 0;
	int          i               = 0;
	int          l               = 0;

    int          reorder         = 1;
	int          leftNeighbor    = 0;
	int          rightNeighbor   = 0;
	int          bottomNeighbor  = 0;
	int          topNeighbor     = 0;
    char         myFileName[64];


    if(myID==0)
	{
		wtime = MPI_Wtime();
	}

	// Allocate arrays - 2D arrays are much easier to deal with in MPI
    double** Vx         = new double* [myN_x+6];
    double** Vy         = new double* [myN_x+6];
    double** h          = new double* [myN_x+6];
    double** tempVx     = new double* [myN_x+6];
    double** tempVy     = new double* [myN_x+6];
    double** temph      = new double* [myN_x+6];
    Vx[0]               = new double  [(myN_x+6)*(myN_y+6)];
    Vy[0]               = new double  [(myN_x+6)*(myN_y+6)];
    h[0]                = new double  [(myN_x+6)*(myN_y+6)];
    tempVx[0]           = new double  [(myN_x+6)*(myN_y+6)];
    tempVy[0]           = new double  [(myN_x+6)*(myN_y+6)];
    temph[0]            = new double  [(myN_x+6)*(myN_y+6)];

    double** k1x        = new double* [myN_x+6];
    double** k2x        = new double* [myN_x+6];
    double** k3x        = new double* [myN_x+6];
    double** k4x        = new double* [myN_x+6];
    k1x[0]              = new double  [(myN_x+6)*(myN_y+6)];
    k2x[0]              = new double  [(myN_x+6)*(myN_y+6)];
    k3x[0]              = new double  [(myN_x+6)*(myN_y+6)];
    k4x[0]              = new double  [(myN_x+6)*(myN_y+6)];

    double** k1y        = new double* [myN_x+6];
    double** k2y        = new double* [myN_x+6];
    double** k3y        = new double* [myN_x+6];
    double** k4y        = new double* [myN_x+6];
    k1y[0]              = new double  [(myN_x+6)*(myN_y+6)];
    k2y[0]              = new double  [(myN_x+6)*(myN_y+6)];
    k3y[0]              = new double  [(myN_x+6)*(myN_y+6)];
    k4y[0]              = new double  [(myN_x+6)*(myN_y+6)];
        
    double** k1h        = new double* [myN_x+6];
    double** k2h        = new double* [myN_x+6];
    double** k3h        = new double* [myN_x+6];
    double** k4h        = new double* [myN_x+6];
    k1h[0]              = new double  [(myN_x+6)*(myN_y+6)];
    k2h[0]              = new double  [(myN_x+6)*(myN_y+6)];
    k3h[0]              = new double  [(myN_x+6)*(myN_y+6)];
    k4h[0]              = new double  [(myN_x+6)*(myN_y+6)];

    for(int i=1, ii=myN_y+6; i<myN_x+6; i++, ii+=(myN_y+6))
	{
		Vx[i]           = &Vx[0][ii];
        Vy[i]           = &Vy[0][ii];
		h[i]            = &h[0][ii];
		tempVx[i]       = &tempVx[0][ii];
		tempVy[i]       = &tempVy[0][ii];
		temph[i]        = &temph[0][ii];

        k1x[i]          = &k1x[0][ii];
        k2x[i]          = &k2x[0][ii];
        k3x[i]          = &k3x[0][ii];
        k4x[i]          = &k4x[0][ii];

        k1y[i]          = &k1y[0][ii];
        k2y[i]          = &k2y[0][ii];
        k3y[i]          = &k3y[0][ii];
        k4y[i]          = &k4y[0][ii];

        k1h[i]          = &k1h[0][ii];
        k2h[i]          = &k2h[0][ii];
        k3h[i]          = &k3h[0][ii];
        k4h[i]          = &k4h[0][ii];
	}

    MPI_Cart_create(MPI_COMM_WORLD, N_D, dimensions, isPeriodic, reorder, &Comm2D);
	MPI_Comm_rank(Comm2D, &myID);
	MPI_Cart_coords(Comm2D, myID, N_D, myCoords);
	MPI_Cart_shift(Comm2D, X, 1, &leftNeighbor,   &rightNeighbor);
	MPI_Cart_shift(Comm2D, Y, 1, &bottomNeighbor, &topNeighbor);

    // Adjust if change to dirichlet BC
	myiStart = 3;
	myiEnd   = myN_x+3;
	myjStart = 3;
	myjEnd   = myN_y +3;
    cout  << "iStart = " << myiStart <<  ", iEnd = " << myiEnd << ", jStart = " << myjStart <<  ", jEnd = " << myjEnd << endl;

    // Create a new datatype to store values on an x boundary
    // x is not contigous in memory, hence we create this new data type
    MPI_Type_vector(myN_x, numElementsPerBlock, myN_y+6, MPI_DOUBLE, &stridetype);
    MPI_Type_commit(&stridetype);

    // Set initial condition
	for(i=3; i<myN_x+3; i++)
	{
        for(j=3; j<myN_x+3; j++) {
            Vx[i][j] = 0.0;
            Vy[i][j] = 0.0;
            h[i][j]  = 1 + 0.5*exp((-1.0/25.0)*(pow(((i+myN_x*myCoords[X])-33),2) +pow(((j+myN_y*myCoords[Y])-33),2)))
                         + 0.5*exp((-1.0/25.0)*(pow(((i+myN_x*myCoords[X])-83),2) +pow(((j+myN_y*myCoords[Y])-23),2)))
                         + 0.5*exp((-1.0/25.0)*(pow(((i+myN_x*myCoords[X])-53),2) +pow(((j+myN_y*myCoords[Y])-73),2)));
        }

	}

    // Time marching loop
	for(l=0; l<N_t; l++)
	{
        
        t += Delta_t;

        // Re-adjust each of the processors' ghost points
        exchange(Vx, myN_x, myN_y, status, Comm2D, stridetype, leftNeighbor, rightNeighbor, bottomNeighbor, topNeighbor);
        exchange(Vy, myN_x, myN_y, status, Comm2D, stridetype, leftNeighbor, rightNeighbor, bottomNeighbor, topNeighbor);
        exchange(h, myN_x, myN_y, status, Comm2D, stridetype, leftNeighbor, rightNeighbor, bottomNeighbor, topNeighbor);     
     

        // First step
		f(k1x, k1y, k1h, Vx, Vy, h, myN_x, myN_y, myiStart, myiEnd, myjStart, myjEnd);

        for(j = myjStart; j<myjEnd; j++) 
        {
            for(i = myiStart; i<myiEnd; i++)
            {
                tempVx[i][j]	= Vx[i][j] + Delta_t/2*k1x[i][j];
                tempVy[i][j]	= Vy[i][j] + Delta_t/2*k1y[i][j];
                temph[i][j]	    = h[i][j]  + Delta_t/2*k1h[i][j];
            }
        }
        exchange(tempVx, myN_x, myN_y, status, Comm2D, stridetype, leftNeighbor, rightNeighbor, bottomNeighbor, topNeighbor);
        exchange(tempVy, myN_x, myN_y, status, Comm2D, stridetype, leftNeighbor, rightNeighbor, bottomNeighbor, topNeighbor);
        exchange(temph, myN_x, myN_y, status, Comm2D, stridetype, leftNeighbor, rightNeighbor, bottomNeighbor, topNeighbor);

   
		// Second step
        f(k2x, k2y, k2h, tempVx, tempVy, temph, myN_x, myN_y, myiStart, myiEnd, myjStart, myjEnd);

        for(j = myjStart; j<myjEnd; j++) 
        {
            for(i = myiStart; i<myiEnd; i++)
            {
                tempVx[i][j]	= Vx[i][j] + Delta_t/2*k2x[i][j];
                tempVy[i][j]	= Vy[i][j] + Delta_t/2*k2y[i][j];
                temph[i][j]	    = h[i][j]  + Delta_t/2*k2h[i][j];
            }
        }
        exchange(tempVx, myN_x, myN_y, status, Comm2D, stridetype, leftNeighbor, rightNeighbor, bottomNeighbor, topNeighbor);
        exchange(tempVy, myN_x, myN_y, status, Comm2D, stridetype, leftNeighbor, rightNeighbor, bottomNeighbor, topNeighbor);
        exchange(temph, myN_x, myN_y, status, Comm2D, stridetype, leftNeighbor, rightNeighbor, bottomNeighbor, topNeighbor);
    
        // Third step
        f(k3x, k3y, k3h, tempVx, tempVy, temph, myN_x, myN_y, myiStart, myiEnd, myjStart, myjEnd);

        for(j = myjStart; j<myjEnd; j++) 
        {
            for(i = myiStart; i<myiEnd; i++)
            {
                tempVx[i][j]	= Vx[i][j] + Delta_t*k3x[i][j];
                tempVy[i][j]	= Vy[i][j] + Delta_t*k3y[i][j];
                temph[i][j]	    = h[i][j]  + Delta_t*k3h[i][j];
            }
        }
        exchange(tempVx, myN_x, myN_y, status, Comm2D, stridetype, leftNeighbor, rightNeighbor, bottomNeighbor, topNeighbor);
        exchange(tempVy, myN_x, myN_y, status, Comm2D, stridetype, leftNeighbor, rightNeighbor, bottomNeighbor, topNeighbor);
        exchange(temph, myN_x, myN_y, status, Comm2D, stridetype, leftNeighbor, rightNeighbor, bottomNeighbor, topNeighbor);

    
		// Fourth step
        f(k4x, k4y, k4h, tempVx, tempVy, temph, myN_x, myN_y, myiStart, myiEnd, myjStart, myjEnd);

        for(j = myjStart; j<myjEnd; j++) 
        {
            for(i = myiStart; i<myiEnd; i++)
            {
                Vx[i][j]	=  Vx[i][j] + Delta_t*(k1x[i][j]/6 + k2x[i][j]/3 
                    + k3x[i][j]/3 + k4x[i][j]/6);
                Vy[i][j]	=  Vy[i][j] + Delta_t*(k1y[i][j]/6 + k2y[i][j]/3 
                    + k3y[i][j]/3 + k4y[i][j]/6);
                h[i][j]	    =  h[i][j] + Delta_t*(k1h[i][j]/6 + k2h[i][j]/3 
                    + k3h[i][j]/3 + k4h[i][j]/6);

            }
        }

        if (l==5) {
            sprintf(myFileName, "part_Process_%d_%d.csv", myCoords[X], myCoords[Y]);
            file.open(myFileName, ios::out);
            for(i=3; i<myN_x+3; i++)
            {
                for(j=3; j<myN_y+3; j++)
                {
                    file << (i-3)+myN_x*myCoords[X] << ", " << (j-3)+myN_y*myCoords[Y] << ", " << h[i][j] << "\n";
                }
            }
        }
    
        // Write the solution - change l if desire different time-step
        if (l==100) {
	    	write(file, h, l, myN_x, myN_y);
        }
        if(myID==0)
        {
	    	cout << "t = " << t << endl;
        }
	}
          file.close();


    // Deallocate arrays
    delete [] Vx;
	delete [] Vy;
	delete [] h;
    delete [] tempVx;
    delete [] tempVy;
    delete [] temph;
	delete [] k1x;
	delete [] k2x;
	delete [] k3x;
	delete [] k4x;

    delete [] k1y;
	delete [] k2y;
	delete [] k3y;
	delete [] k4y;

    delete [] k1h;
	delete [] k2h;
	delete [] k3h;
	delete [] k4h;

	if(myID==0)
	{
		wtime	= MPI_Wtime() - wtime;	// Record the end time and calculate elapsed time
		cout << "Simulation took " << wtime/N_t << " seconds per time step with " << N_Procs << " processes" << endl;
	}

	MPI_Finalize();

	return 0;
}

void    f(double** kx, double** ky, double** kh, double** Vx, double** Vy, double** h, int myN_x, int myN_y,
            int iStart, int iEnd, int jStart, int jEnd) {

	for(int i=iStart; i<iEnd; i++)
	{
        for (int j = jStart; j< jEnd; j++) 
        {

            int n_1 = i, n_2 = i, n_3 = i, n_4 = i, n_5 = i, n_6 = i;

            int m_1 = j, m_2 = j, m_3 = j, m_4 = j, m_5 = j, m_6 = j;


            kx[i][j]    = (-g/Delta_x)*((-1.0/60.0)*h[n_1-3][j] + (3.0/20.0)*h[n_2-2][j] + (-3.0/4.0)*h[n_3-1][j] + (3.0/4.0)*h[n_4+1][j] + (-3.0/20.0)*h[n_5+2][j] + (1.0/60.0)*h[n_6+3][j])
                                + (-Vx[i][j]/Delta_x)*((-1.0/60.0)*Vx[n_1-3][j] + (3.0/20.0)*Vx[n_2-2][j] + (-3.0/4.0)*Vx[n_3-1][j] + (3.0/4.0)*Vx[n_4+1][j] + (-3.0/20.0)*Vx[n_5+2][j] + (1.0/60.0)*Vx[n_6+3][j])
                                + (-Vy[i][j]/Delta_y)*((-1.0/60.0)*Vx[i][m_1-3] + (3.0/20.0)*Vx[i][m_2-2] + (-3.0/4.0)*Vx[i][m_3-1] + (3.0/4.0)*Vx[i][m_4+1] + (-3.0/20.0)*Vx[i][m_5+2] + (1.0/60.0)*Vx[i][m_6+3]);


            ky[i][j]	= (-g/Delta_y)*((-1.0/60.0)*h[i][m_1-3] + (3.0/20.0)*h[i][m_2-2] + (-3.0/4.0)*h[i][m_3-1] + (3.0/4.0)*h[i][m_4+1] + (-3.0/20.0)*h[i][m_5+2] + (1.0/60.0)*h[i][m_6+3])
                                + (-Vx[i][j]/Delta_x)*((-1.0/60.0)*Vy[n_1-3][j] + (3.0/20.0)*Vy[n_2-2][j] + (-3.0/4.0)*Vy[n_3-1][j] + (3.0/4.0)*Vy[n_4+1][j] + (-3.0/20.0)*Vy[n_5+2][j] + (1.0/60.0)*Vy[n_6+3][j])
                                + (-Vy[i][j]/Delta_y)*((-1.0/60.0)*Vy[i][m_1-3] + (3.0/20.0)*Vy[i][m_2-2] + (-3.0/4.0)*Vy[i][m_3-1] + (3.0/4.0)*Vy[i][m_4+1] + (-3.0/20.0)*Vy[i][m_5+2] + (1.0/60.0)*Vy[i][m_6+3]);


            kh[i][j]    =  (-1/Delta_x)*((-1.0/60.0)*Vx[n_1-3][j]*h[n_1-3][j] + (3.0/20.0)*Vx[n_2-2][j]*h[n_2-2][j] + (-3.0/4.0)*Vx[n_3-1][j]*h[n_3-1][j] + (3.0/4.0)*Vx[n_4+1][j]*h[n_4+1][j] + (-3.0/20.0)*Vx[n_5+2][j]*h[n_5+2][j] + (1.0/60.0)*Vx[n_6+3][j]*h[n_6+3][j])
                                + (-1/Delta_y)*((-1.0/60.0)*Vy[i][m_1-3]*h[i][m_1-3] + (3.0/20.0)*Vy[i][m_2-2]*h[i][m_2-2] + (-3.0/4.0)*Vy[i][m_3-1]*h[i][m_3-1] + (3.0/4.0)*Vy[i][m_4+1]*h[i][m_4+1] + (-3.0/20.0)*Vy[i][m_5+2]*h[i][m_5+2] + (1.0/60.0)*Vy[i][m_6+3]*h[i][m_6+3]);

        }
    }
}

void exchange(double** phi, int myN_x, int myN_y, MPI_Status status, MPI_Comm Comm2D, MPI_Datatype stridetype,
             int leftNeighbor, int rightNeighbor, int bottomNeighbor, int topNeighbor) {

    for (int i = 0; i < 3; i++) {
        // Connect the inner boundaries
        MPI_Sendrecv(&(phi[i+3][3]),      myN_y, MPI_DOUBLE,  leftNeighbor,   0, &(phi[myN_x+i+3][3]), myN_y, MPI_DOUBLE, rightNeighbor,  0, Comm2D, &status);
        MPI_Sendrecv(&(phi[myN_x+i][3]),  myN_y, MPI_DOUBLE,  rightNeighbor,  0, &(phi[i][3]),         myN_y, MPI_DOUBLE, leftNeighbor,   0, Comm2D, &status);
        MPI_Sendrecv(&(phi[3][3+i]),      1,     stridetype,  bottomNeighbor, 0, &(phi[3][myN_y+i+3]), 1,     stridetype, topNeighbor,    0, Comm2D, &status);
        MPI_Sendrecv(&(phi[3][myN_y+i]),  1,     stridetype,  topNeighbor,    0, &(phi[3][i]),         1,     stridetype, bottomNeighbor, 0, Comm2D, &status);

        // // Connect the outer boundaries
        MPI_Sendrecv(&(phi[myN_x+i][3]),  myN_y, MPI_DOUBLE,  leftNeighbor,   0, &(phi[i][3]),         myN_y, MPI_DOUBLE, rightNeighbor,  0, Comm2D, &status);		
        MPI_Sendrecv(&(phi[i+3][3]),      myN_y, MPI_DOUBLE,  rightNeighbor,  0, &(phi[myN_x+3+i][3]), myN_y, MPI_DOUBLE, leftNeighbor,   0, Comm2D, &status);
        MPI_Sendrecv(&(phi[3][myN_y+i]),  1,     stridetype,  bottomNeighbor, 0, &(phi[3][i]),         1,     stridetype, topNeighbor,    0, Comm2D, &status);
        MPI_Sendrecv(&(phi[3][i+3]),      1,     stridetype,  topNeighbor,    0, &(phi[3][myN_y+3+i]), 1,     stridetype, bottomNeighbor, 0, Comm2D, &status);
    }
}

