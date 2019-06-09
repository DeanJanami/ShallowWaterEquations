/////////////////////////////////////////////////////////////////////////////
//
// Dean Pakravan: 757389
// Qijie Li:      927249
//
// Applied High Performance Computing: Assignment 1
//
// The Shallow Water Equations
//
// Method:		Finite Difference Method with sixth order central differences
//				and Fourth Order Runge-Kutta Method, parallelized with OpenMP
//
// Compilation:	module load gcc
//              g++ SWE_op.cpp -o SWE_op -fopenmp
//              OR
//              bgxlc++ SWE_op.cpp -o SWE_op -qsmp=omp
//
// Execution:	export OMP_NUM_THREADS=1, 2, 4, 8, 16, etc
//				./SWE_op
//              OR
//              srun --ntasks=1 --nodes=1 --ntasks-per-node=1 ./SWE_op
//              OR
//              srun -n 1 -N 1 --ntasks-per-node=1 ./SWE_op
//
/////////////////////////////////////////////////////////////////////////////

#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h>
#include <omp.h>

using namespace std;

// Simulation Parameters
const double x_min   =  0.00;
const double x_max   = 100.00;
const double y_min   =  0.00;
const double y_max   = 100.00;
const double t_min   =  0.00;
const double t_max   = 100.00;
const double Delta_x =  0.100;
const double Delta_y =  0.100;
const double Delta_t =  0.20;
const int	 N_x     =  (x_max-x_min)/Delta_x+1;
const int	 N_y     =  (y_max-y_min)/Delta_y+1;
const int	 N_t     =  (t_max-t_min)/Delta_t+1;
const double g       =  9.81;

// Function declarations
void f(double*** k, double*** phi);
void write(fstream& file, double*** phi, int N_x, int N_y);

int main(int argc, char** argv)
{
	// Simulation parameters
	double  t            = 0;
	// Indexes
	int     k            = 0;
	int     j            = 0;
	int     i            = 0;
	int     l            = 0;
    fstream		file;
	double  wtime     = omp_get_wtime();			// Record the starting time
	int     N_Threads = omp_get_max_threads();


    // Allocate arrays	
	double*** phi     = new double**[N_x];
    double*** tempPhi = new double**[N_x];
    double*** k1      = new double**[N_x];
    double*** k2      = new double**[N_x];
    double*** k3      = new double**[N_x];
    double*** k4      = new double**[N_x];
	for(i =0; i<N_x; i++){
	    phi[i]     = new double*[N_y];
        tempPhi[i] = new double*[N_y];
        k1[i]      = new double*[N_y];
        k2[i]      = new double*[N_y];
        k3[i]      = new double*[N_y];
	    k4[i]      = new double*[N_y];
		for(j =0; j<N_y; j++){
			phi[i][j]     = new double[3];
            tempPhi[i][j] = new double[3];
            k1[i][j]      = new double[3];
            k2[i][j]      = new double[3];
            k3[i][j]      = new double[3];
	        k4[i][j]      = new double[3];
			for(k = 0; k<3;k++){
				phi[i][j][k]     = 0;
                tempPhi[i][j][k] = 0;
                k1[i][j][k]      = 0;
                k2[i][j][k]      = 0;
                k3[i][j][k]      = 0;
                k4[i][j][k]      = 0;
			}
		}
	}
	// Vx[i][j] == phi[i][j][0]
    // Vy[i][j] == phi[i][j][1]
    //  h[i][j] == phi[i][j][2]

    // Set initial condition
	for(i=0; i<N_x; i++)
	{
        for(j=0; j<N_x; j++) {
            phi[i][j][0] = 0.0;
            phi[i][j][1] = 0.0;
            phi[i][j][2]  = 1 + 0.5*exp((-1.0/25.0)*(pow((i-30),2) +pow((j-30),2)))
			                  + 0.5*exp((-1.0/25.0)*(pow((i-80),2) +pow((j-20),2)))
                              + 0.5*exp((-1.0/25.0)*(pow((i-50),2) +pow((j-70),2)));
        }

	}

	file.open("SWE_op.csv", ios::out);
    // We write down the initial conditions
	write(file, phi, N_x,N_y);

    // Time marching loop
	#pragma omp parallel default(shared) private(i,j,k,l)
	for(l=0; l<N_t; l++)
	{
		#pragma omp single
		{
		t += Delta_t;
		}

        // First step
		f(k1, phi);
		
		#pragma omp for schedule(static)
        for (k = 0; k < 3; k++) 
        {
            for(j = 0; j<N_y; j++) 
            {
                for(i = 0; i<N_x; i++)
                {
                    tempPhi[i][j][k]	= phi[i][j][k] + Delta_t/2*k1[i][j][k];
                }
            }
        }
		// Second step
        f(k2, tempPhi);
		
		#pragma omp for schedule(static)
        for (k = 0; k < 3; k++) 
        {
            for(j = 0; j<N_y; j++) 
            {
                for(i = 0; i<N_x; i++)
                {
                    tempPhi[i][j][k]	= phi[i][j][k] + Delta_t/2*k2[i][j][k];
                }
            }
        }
        // Third step
		f(k3, tempPhi);

		#pragma omp for schedule(static)
        for (k = 0; k < 3; k++) 
        {
            for(j = 0; j<N_y; j++) 
            {
                for(i = 0; i<N_x; i++)
                {
                    tempPhi[i][j][k]	= phi[i][j][k] + Delta_t*k3[i][j][k];
                }
            }
        }
		// Fourth step
		f(k4, tempPhi);

		#pragma omp for schedule(static)
        for (k = 0; k < 3; k++) 
        {
            for(j = 0; j<N_y; j++) 
            {
                for(i = 0; i<N_x; i++)
                {
                    phi[i][j][k]	=  phi[i][j][k] + Delta_t*(k1[i][j][k]/6 + k2[i][j][k]/3 
                                        + k3[i][j][k]/3 + k4[i][j][k]/6);
                }
            }
        }
		#pragma omp single
		{

		// Write down the solution
		write(file, phi, N_x, N_y);


		cout << "t = " << t << endl;
		}
	}

	file.close();

    // Deallocate arrays
    delete [] tempPhi;
	delete [] k1;
	delete [] k2;
	delete [] k3;
	delete [] k4;
	delete [] phi;

	wtime	= omp_get_wtime() - wtime;	// Record the end time and calculate elapsed time
	cout << "Simulation took " << wtime/100 << " seconds per time step with " << N_Threads << " threads" << endl;

	return 0;
}

void	f(double*** k, double*** phi)
{

	#pragma omp for nowait
	for(int i=0; i<N_x; i++)
	{
        for (int j = 0; j< N_y; j++) 
        {
            // Periodic boundary conditions on i
            int n_1 = i, n_2 = i, n_3 = i, n_4 = i, n_5 = i, n_6 = i;
            if (i==2) {
                n_1 = N_x+2;
            }
            else if (i==1) {
                n_1 = N_x+1;
                n_2 = N_x+1;
            }
            else if (i==0) {
                n_1 = N_x;
                n_2 = N_x;
                n_3 = N_x;
            }
            else if (i==N_x-3) {
                n_6 = -3;
            }
            else if (i==N_x-2) {
                n_5 = -2;
                n_6 = -2;
            }
            else if (i==N_x-1) {
                n_4 = -1;
                n_5 = -1;
                n_6 = -1;
            }
            // Periodic boundary conditions on j
            int m_1 = j, m_2 = j, m_3 = j, m_4 = j, m_5 = j, m_6 = j;
            if (j==2) {
                m_1 = N_y+2;
            }
            else if (j==1) {
                m_1 = N_y+1;
                m_2 = N_y+1;
            }
            else if (j==0) {
                m_1 = N_y;
                m_2 = N_y;
                m_3 = N_y;
            }
            else if (j==N_y-3) {
                m_6 = -3;
            }
            else if (j==N_y-2) {
                m_5 = -2;
                m_6 = -2;
            }
            else if (j==N_y-1) {
                m_4 = -1;
                m_5 = -1;
                m_6 = -1;
            }
            // NewKx
            k[i][j][0] = (-g/Delta_x)*((-1.0/60.0)*phi[n_1-3][j][2] + (3.0/20.0)*phi[n_2-2][j][2] + (-3.0/4.0)*phi[n_3-1][j][2] + (3.0/4.0)*phi[n_4+1][j][2] + (-3.0/20.0)*phi[n_5+2][j][2] + (1.0/60.0)*phi[n_6+3][j][2])
                                + (-phi[i][j][0]/Delta_x)*((-1.0/60.0)*phi[n_1-3][j][0] + (3.0/20.0)*phi[n_2-2][j][0] + (-3.0/4.0)*phi[n_3-1][j][0] + (3.0/4.0)*phi[n_4+1][j][0] + (-3.0/20.0)*phi[n_5+2][j][0] + (1.0/60.0)*phi[n_6+3][j][0])
                                + (-phi[i][j][1]/Delta_y)*((-1.0/60.0)*phi[i][m_1-3][0] + (3.0/20.0)*phi[i][m_2-2][0] + (-3.0/4.0)*phi[i][m_3-1][0] + (3.0/4.0)*phi[i][m_4+1][0] + (-3.0/20.0)*phi[i][m_5+2][0] + (1.0/60.0)*phi[i][m_6+3][0]);

 			// NewKy
            k[i][j][1]	= (-g/Delta_y)*((-1.0/60.0)*phi[i][m_1-3][2] + (3.0/20.0)*phi[i][m_2-2][2] + (-3.0/4.0)*phi[i][m_3-1][2] + (3.0/4.0)*phi[i][m_4+1][2] + (-3.0/20.0)*phi[i][m_5+2][2] + (1.0/60.0)*phi[i][m_6+3][2])
                                + (-phi[i][j][0]/Delta_x)*((-1.0/60.0)*phi[n_1-3][j][1] + (3.0/20.0)*phi[n_2-2][j][1] + (-3.0/4.0)*phi[n_3-1][j][1] + (3.0/4.0)*phi[n_4+1][j][1] + (-3.0/20.0)*phi[n_5+2][j][1] + (1.0/60.0)*phi[n_6+3][j][1])
                                + (-phi[i][j][1]/Delta_y)*((-1.0/60.0)*phi[i][m_1-3][1] + (3.0/20.0)*phi[i][m_2-2][1] + (-3.0/4.0)*phi[i][m_3-1][1] + (3.0/4.0)*phi[i][m_4+1][1] + (-3.0/20.0)*phi[i][m_5+2][1] + (1.0/60.0)*phi[i][m_6+3][1]);

            // NewKh
            k[i][j][2] = (-1/Delta_x)*((-1.0/60.0)*phi[n_1-3][j][0]*phi[n_1-3][j][2] + (3.0/20.0)*phi[n_2-2][j][0]*phi[n_2-2][j][2] + (-3.0/4.0)*phi[n_3-1][j][0]*phi[n_3-1][j][2] + (3.0/4.0)*phi[n_4+1][j][0]*phi[n_4+1][j][2] + (-3.0/20.0)*phi[n_5+2][j][0]*phi[n_5+2][j][2] + (1.0/60.0)*phi[n_6+3][j][0]*phi[n_6+3][j][2])
                                + (-1/Delta_y)*((-1.0/60.0)*phi[i][m_1-3][1]*phi[i][m_1-3][2] + (3.0/20.0)*phi[i][m_2-2][1]*phi[i][m_2-2][2] + (-3.0/4.0)*phi[i][m_3-1][1]*phi[i][m_3-1][2] + (3.0/4.0)*phi[i][m_4+1][1]*phi[i][m_4+1][2] + (-3.0/20.0)*phi[i][m_5+2][1]*phi[i][m_5+2][2] + (1.0/60.0)*phi[i][m_6+3][1]*phi[i][m_6+3][2]);

        }
    }
}

void    write(fstream& file, double*** phi, int N_x, int N_y)
{
	//phi[i][j][2] are the height values
	for(int i=0; i<N_x; i++)
	{
        for (int j = 0; j < N_y; j++)
        {
		    file << i << ", " << j << ", " << phi[i][j][2] << endl;
        }
	}
	return;
}