// This is the code for Multiple Allocation Hub and Spoke Network Design for the model described in:
// Hamachar 1996. Efficient algorithms for the uncapacitated multiple allocation p-hub median problem. Location Science, 4(3), 139-154.
#pragma warning(disable : 4996)//For Visual Studio 2012
#include<stdio.h>
#include<conio.h>
#include<iostream>
#include<fstream>
#include<iosfwd>
#include<string>
#include <deque>
#include <sstream>
#include <time.h>
#include <stdlib.h>

#include <ilcplex/ilocplex.h>
#include <ilconcert/ilosys.h>

ILOSTLBEGIN

typedef IloArray<IloBoolVarArray> array2d;//Creating a 2d array of x variables
typedef IloArray<array2d> array3d;//Creating a 3d array of x variables
typedef IloArray<array3d> array4d;//Creating a 4d array of x variables

//***************************For LP Relaxation ********************************
typedef IloArray<IloNumVarArray> array2dNumVar;//Creating a 2d array of x variables
typedef IloArray<array2dNumVar> array3dNumVar;//Creating a 3d array of x variables
typedef IloArray<array3dNumVar> array4dNumVar;//Creating a 4d array of x variables
//*****************************************************************************
typedef IloArray<IloNumArray> TwoDMatrix;
typedef IloArray<TwoDMatrix> ThreeDMatrix;
typedef IloArray<ThreeDMatrix> FourDMatrix;

int main(int argc, char **argv)
{
	IloEnv env;
	int N=12;//N=no. of nodes
	

	TwoDMatrix w(env);//flow data between pairs of nodes
		
	IloNum eps;//later assigned as eps = cplex.getParam(IloCplex::EpInt); EpInt is the tolerance gap for integer variables

	try
	{
		///////// DATA FILE READING //////////

		//const char* data_filename  = "Data_HS_CAB_5node.dat";
		//const char* data_filename  = "Data_HS_CAB_10node.dat";
		//const char* data_filename = "Data_HS_CAB_15node.dat";
		const char* data_filename = "Data.dat";
		if (argc > 1)
		{
			data_filename = argv[1];
		}
		fstream datafile;
		datafile.open(data_filename, ios::in);

		if (!datafile)
		{
			cerr << "ERROR: could not open file " << data_filename << " for reading" << endl;
			cerr << "usage:   " << argv[0] << " <datafile>" << endl;
			throw(-1);
		}

		datafile >> w;



		datafile.close();

	

		//=================Decleare Variables===============================
		IloModel model(env);
		//array4dNumVar y(env, N);//xijkm = amount of the flow originated at i passing via hubs k and m in that order
		//IloBoolVarArray z(env, N); //zk = 1 
		// for lp relaxation
		array2d y(env, N);
		//=======================================================================
		TwoDMatrix y_val(env, N);//to hold values of x
				//=======================================================================

		for (int i = 0; i < N; i++)
		{
			y[i] = IloBoolVarArray(env, N);
			y_val[i] = IloNumArray(env, N);
			
		}
		
		// Objective Function: Minimize: sum {i in 1..N}{j in i..N}{k in 1..N}{m in 1..N}
		//lambda[i][j][k][m]*(c_collect*dist[i][k]+c_tranship*dist[k][m]+c_distribute*dist[m][j])*x[i][j][k][m]

		IloExpr Obj(env); // Creates an expression with the name Obj (Objective)		


		for (int i = 0; i < N-1; i++)
		{
			for (int j = i+1; j < N; j++)
			{
				
				Obj += w[i][j] * y[i][j];
								
			}
			
		}
		
		
		// model.add is a function to add constraints and objectives in the CPLEX environment
		model.add(IloMaximize(env, Obj)); // IloMinimize is used for minimization problems
		Obj.end(); // Clear out the memory allocated for the expression 


		//Constraint 2: for {i in 1..N}{j in 1..N}{k in 1..N}: x[i][j][k][m] <= z[k] 
		for (int i=0; i < N-2; i++)
		{
			for (int j = i+1; j < N-1; j++)
			{
				for (int k = j+1; k < N; k++)
				{
					model.add(y[i][j] + y[i][k] + y[k][j] <= 2);
					
				}
			}
		}
		
		//Constarint 4: sum{ k in 1..N }(z[k][k]) = p;
		for (int i = 0; i < N - 1; i++)
		{
			for (int j = i + 1; j < N ; j++)
			{
				for (int k=0; k<N;k++)
										
				{
					if (k != i && k != j)
					{
						model.add(y[i][j] - y[i][k] - y[k][j] <= 0);
					}
				}
			}
		}

		for (int i = 0; i < N - 1; i++)
		{
			for (int j = i + 1; j < N; j++)
			{
				model.add(y[i][j] = y[j][i]);
			}
		}


		//==============================================================================

		//Optimize
		IloCplex cplex(model);
		int N_Variables = cplex.getNcols();
		int N_Constraints = cplex.getNrows();
		cout << "No. of Original Variables = " << N_Variables << endl;
		cout << "No. of Original Constraints = " << N_Constraints << endl;
		cplex.exportModel("Breaks.lp");
		//system("pause");
		//cplex.setOut(env.getNullStream()); //This is to supress the output of Branch & Bound Tree on screen
		//cplex.setWarning(env.getNullStream()); //This is to supress warning messages on screen
		eps = cplex.getParam(IloCplex::EpInt);

		clock_t t1, t2;
		t1 = clock();
		cplex.solve();//solving the MODEL
		t2 = clock();
		float diff((float)t2 - (float)t1);
		float seconds = diff / CLOCKS_PER_SEC;
		cout << "TOTAL CPUTIME = " << seconds << " secs" << endl;

		if (cplex.getStatus() == IloAlgorithm::Infeasible) // if the problem is infeasible
		{
			cout << "Problem is Infeasible" << endl;
		}

		// Print results
		cout << "Maximum Cost = " << cplex.getObjValue() << endl;

		cout << "Nodes \t Hubs: " << endl;
		// for lp relaxation
		
			
		cout << "i j y[i][j]" << endl;
		for (int i = 0; i < N-1; i++)
		{
			for (int j = i+1; j < N; j++)
			{
				if (cplex.getValue(y[i][j]) > 0)
						{
							cout << i << " " << j <<  " " << cplex.getValue(y[i][j]) << endl;
						}
						else
						{
							y[i][j] = 0;
						}
					}
				}
			}
		

		

	//end of of try block
	catch (IloException& ex)
	{
		cerr << "Error: " << ex << endl;
	}
	catch (...)
	{
		cerr << "Error" << endl;
	}
	return 0;
	env.end();
}