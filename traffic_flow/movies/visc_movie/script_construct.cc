#include "nr3.h"

using namespace std;

int main() {

	int cycles = 3;

	ofstream outfile;
	outfile.open("movie-nonviscous");

	outfile << "set xrange[0:1]\n";
	outfile << "set yrange[0.6:1.1]\n";
	outfile << "set grid\n";

	for (int i=0; i<cycles; i++) {
		for (int t=0; t<50; t+= 1) {
			outfile << "t = " << t << endl;
			outfile << "plot '0.7_data_t='.t.'.dat' w l" << endl;
			outfile << "pause 0.3" << endl;
		}
	}

	outfile.close();

	outfile.open("movie-viscous"); 

	outfile << "set xrange[0:1]\n";
	outfile << "set yrange[0.6:1.1]\n";
	outfile << "set grid\n";

	for (int i=0; i<cycles; i++) {
		for (int t=0; t<50; t+= 1) {
			outfile << "t = " << t << endl;
			outfile << "plot 'visc=0.17_rho_bar=0.7_t='.t.'.dat' w l" << endl;
			outfile << "pause 0.3" << endl;
		}
	}

	outfile.close();
	return 0;
}