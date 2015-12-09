#include "nr3.h"

using namespace std;

int main() {

	int cycles = 3;

	ofstream outfile;
	outfile.open("movie-traffic");

	outfile << "set xrange[0:1]\n";
	outfile << "set yrange[0.5998:0.601]\n";
	outfile << "set grid\n";

	for (int i=0; i<cycles; i++) {
		for (int t=0; t<50; t+= 1) {
			outfile << "t = " << t << endl;
			outfile << "plot 'trash_data_t='.t.'.dat' w l" << endl;
			outfile << "pause 0.05" << endl;
		}
	}

	outfile.close();

	return 0;
}