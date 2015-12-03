#include <iostream>
#include <fstream>

using namespace std;

int main() {
	int a, b;
	ofstream file("output.dat");
	for (a=1; a<6; a++) {
		for (b=1; b<6; b++) {
			file << a << " " << b << endl;
		}
	}
	file.close();

	return 0;
}
