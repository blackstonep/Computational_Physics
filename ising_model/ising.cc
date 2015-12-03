#include "nr3.h"
#include "ludcmp.h"
#include "ran.h"

using namespace std; 

time_t thyme = time(0); 
double masterSeed = double(thyme);



double expo(int energy, double temp) {
	return exp(- double(energy) / temp );
};

struct Cake {
	int m; 
	double T; 
	MatInt *state;

	VecDoub *magstats;
	VecDoub *engstats; 


	//Constructor
	Cake(int m_in, double T_in) {
		m = m_in;
		T = T_in;


		Ran ran(masterSeed+1.0);

		//Clear stat vectors
		magstats = new VecDoub(2);
		engstats = new VecDoub(2);
		(*engstats)[0] = 0.0; (*engstats)[1] = 0.0; 
		(*magstats)[0] = 0.0; (*magstats)[1] = 0.0;

		//Initialize state matrix
		state = new MatInt(m,m);
		for (int i=0; i<m; i++) {
			for (int j=0; j<m; j++) {
				(*state)[i][j] = 1;

//				double random = ran.doub(); 
//				if (random < 0.5) {
//					(*state)[i][j] = -1;
//				}
//				else {
//					(*state)[i][j] = 1;
//				}

			};
		};
	};

	//Destructor
	~Cake() {
		cout << "destructing cake...\n";
		delete state;
		delete magstats;
		delete engstats;
	};

	//Print state of cake
	void printstate() {

		//Cycle through rows
		for (int i=0; i<m; i++) {

			//Cycle through columns
			for (int j=0; j<m; j++) {
				if ((*state)[i][j]==1){
					cout << (*state)[i][j]<< "     ";
				}
				else {
					cout << (*state)[i][j]<<"    ";
				};
			};
			cout<<endl;
		};
	};

	void flip(int row, int col) {
		(*state)[row][col] *= -1;
	}


	//Should we flip it?!
	bool decide(int row, int col, double random) {
		

		//Energy of current spin state
		double E;

		//Account for periodic boundary cdts
		if (row == 0) {
			if (col == 0) {
				E = - (*state)[0][0] * ( (*state)[m-1][0] + (*state)[0][m-1]
							+ (*state)[0][1] + (*state)[1][0] );
			}
			else if (col == m-1) {
				E = - (*state)[0][col] * ( (*state)[0][m-2] + (*state)[0][0] 
							+ (*state)[1][m-1] + (*state)[m-1][m-1] );
			}
			else {
				E = - (*state)[0][col] * ( (*state)[0][col+1] + (*state)[0][col-1] 
							+ (*state)[1][col] + (*state)[m-1][col] );
			}
		}
		else if (row == m-1) {
			if (col == 0) {
				E = - (*state)[m-1][0] * ( (*state)[m-1][1] + (*state)[m-1][m-1] 
							+ (*state)[m-2][0] + (*state)[0][0] );
			}
			else if (col == m-1) {
				E = - (*state)[m-1][m-1] * ( (*state)[m-1][0] + (*state)[m-1][m-2] 
							+ (*state)[m-2][m-1] + (*state)[0][m-1] ); 
			}
			else {
				E = - (*state)[m-1][col] * ( (*state)[m-1][col-1] + (*state)[m-1][col+1] 
							+ (*state)[m-2][col] + (*state)[0][col] );
			}
		}
		else if (col == 0) {
			E = - (*state)[row][0] * ( (*state)[row][1] + (*state)[row][m-1] 
						+ (*state)[row-1][0] + (*state)[row+1][0]);
		}
		else if (col == m-1) {
			E = - (*state)[row][m-1] * ( (*state)[row][0] + (*state)[row][m-2] 
						+ (*state)[row-1][m-1] + (*state)[row+1][m-1]);
		}
		else {
			E = - (*state)[row][col] * ( (*state)[row-1][col] + (*state)[row+1][col]
						+ (*state)[row][col-1] + (*state)[row][col+1] );
		};

		//dE = Efin - Einit = (-E) - E = - 2 * E
		double dE = - 2 * E;

		if(  expo(dE, T) > random) {return true;}
		else {return false;}
	}

	//Calculate energy of single state
	double energy(Cake &object, int dim) {	
		int intergy = 0;
		double energy;
		for (int i=0; i<dim; i++) {
			for (int j=0; j<dim; j++) {
				if (i==dim-1) {
					if (j==dim-1) {
						intergy += - (*state)[i][j] * ((*state)[i][0] + (*state)[0][j] ); 
					}
					else {
						intergy += - (*state)[i][j] * ((*state)[i][j+1] + (*state)[0][j] );
					}
				}
				else if (j==0) {
					intergy += - (*state)[i][j] * ((*state)[i][0] + (*state)[i+1][j] );
				}
				else {
					intergy += - (*state)[i][j] * ((*state)[i+1][j] + (*state)[i][j+1] );
				}
				//Doublize energy
				energy = double(intergy) / double(dim*dim);
			}
		}
		return energy;
	}

	//Calculate magn. of single state
	double magnet(Cake &object, int dim) {
		int mag = 0; 
		for (int i=0; i<dim; i++) {
			for (int j=0; j<dim; j++) {
				mag += (*state)[i][j];
			}
		}
		return ( double(mag) / double(dim * dim) );
	}

	//Record <E>, <M>, and mean-squares
	void stats(double engstat, double magstat, double total) {
		(*engstats)[0] += engstat / total ;
		(*engstats)[1] += engstat * engstat / total ; 
		(*magstats)[0] += magstat / total; 
		(*magstats)[1] += magstat * magstat / total;
	}

	void sweep(Cake &object, int sweep_num, int dim, int stat_cutoff) {

		//Set random number seqs.
		double col_seed = masterSeed + 5.0; 
		double row_seed = masterSeed + 2.0; 
		Ran ranrow(row_seed);
		Ran rancol(col_seed); 
		double decide_seed = masterSeed;
		Ran ran_decide(decide_seed);

		//Make energy and magn. data files
		ofstream energyfile("energy.dat");
		ofstream magnetfile("magnet.dat"); 			
		
		int step = dim*dim;

		for (int i=0; i<sweep_num; i++) {
			for (int spin_num=0; spin_num < step; spin_num++ ) {
		
				//Choose random element
				int ran_row = int(floor(dim * ranrow.doub() ) );
				int ran_col = int(floor(dim * rancol.doub() ) ); 
				double ran_dec = ran_decide.doub();

				//Flip if it needs flippin'
				if (object.decide(ran_row, ran_col, ran_dec ) ) {
					object.flip(ran_row, ran_col);
				}
			};

			//Output to data files
			double magnet = object.magnet(object, dim);
			magnetfile << i << ", " << magnet << endl;

			double energy = object.energy(object, dim);
			energyfile << i << ", " << energy << endl;

			//Update stat vectors
			if (i>=stat_cutoff) {
				object.stats(energy, magnet, double(sweep_num - stat_cutoff ));
			};
		};
	}

	//Output <E>, <M>, and std deviations
	VecDoub printstats(Cake &object) {
		VecDoub stat(4);

		stat[0] = (*engstats)[0];
		stat[2] = (*magstats)[0];

		stat[1] = sqrt( (*engstats)[1] 
			- (*engstats)[0] * (*engstats)[0] );
		stat[3] = sqrt( (*magstats)[1] 
			- (*magstats)[0] * (*magstats)[0] ); 

		return stat;
	}
};



int main() {

	int dim; 
	cout << "Enter dimension: ";
	cin >> dim; 
	VecDoub stats(4);
	double test_temp;
	cout<<"Set plotting temp: ";
	cin>>test_temp;
	Cake cake1(dim, test_temp); 
	cake1.printstate();

/*	int sweep_num = 2500;
	cake1.sweep(cake1, sweep_num, dim, sweep_num-1000);

	//Testing stat vector system
	stats = cake1.printstats(cake1);

	cout<<"Average Energy = " << stats[0] <<endl;
	cout<<"Standard Energy Deviation: " << 
		stats[1] <<endl;
	cout << "Average Magnetization = " << stats[2] << endl;
	cout << "Standard Mag. Deviation: " << 
		stats[3]  << endl;

	ofstream tempdat("tempdata.dat"); 
*/
//	for (double s=5.0; s>0.0; s-=.1) {
//		VecDoub tempstats(4); 
//
//		//Coordinate transformation
//			//to focus on temp = 2.37
//		double A = .135; 
//		double B = 13.0/625.0 - 15.0*A/2.0;
//		double C = 25.0*A/2.0 + 112.0/125.0;
//		double temp = s * (A*s*s + B*s + C);
//
//		Cake cake(dim, temp);
//		cake.sweep(cake, sweep_num, dim, sweep_num-100);
//
//		tempstats = cake.printstats(cake);
//
//		tempdat << temp << ", " << tempstats[0] << ", " 
//			<< tempstats[1] << ", " << tempstats[2] 
//			<< ", " << tempstats[3] <<endl;
//			
//	};


	return 0;
};

