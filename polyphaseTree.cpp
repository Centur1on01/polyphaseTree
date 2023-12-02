#define _SILENCE_NONFLOATING_COMPLEX_DEPRECATION_WARNING
#define _USE_MATH_DEFINES

#include <iostream>
#include <complex>
#include <string>
#include <fstream>
#include <vector>
#include <deque>


using namespace std;

class polyphaseTree {

	int numOutput, log2num, order;
	string fileName;
	double coeff;
	vector<double> coeffTable;
	deque<complex<double>> output;

public:

	polyphaseTree(const string& fileName = "coeffTable.txt", const int& numOutput = 2) : fileName(fileName), numOutput(numOutput), order(0) {
		ifstream file(fileName, ios::in);

		if (file.is_open())
			while (!file.eof()) {
				file >> coeff;
				coeffTable.push_back(coeff);
				++order;
			}
		file.close();

		log2num = log2(numOutput);

		output.resize(numOutput, { 0.0, 0.0 });
	}



	void next(const deque<complex<double>>& signal) {
		//double divider = numOutput / 2;
		//int jj = 0;

		for (int i = 0; i < numOutput; ++i) { // iteration for outputs

			for (int k = 0; k < log2num; ++k) { // iteration through sections

				complex<double> outputTemp = output[i];
				output[i] = 0;

				// check if needed multiplication by -1 on subfilter output
				// correct iteration through coeffTable
				int evenSubfilter = 0;
				if (i % (int)pow(2, k + 1) >= (int)pow(2, k + 1) - (int)pow(2, k + 1) / 2) {
					evenSubfilter = 1;
				}

				// defining which section output is being calculated
				if (k == 0) {
					for (int l = 0; l < order / 2; ++l) {
						output[i] += signal[i] * coeffTable[order - 1 - (2 * l) - evenSubfilter];
					}
				}
				else {
					for (int l = 0; l < order / 2; ++l) {
						output[i] += outputTemp * coeffTable[order - 1 - (2 * l) - evenSubfilter];
					}
				}

				// multiplying by -1
				if (evenSubfilter) {
					output[i] *= -1;
					if (output[i].imag() == -0) {
						output[i].imag(output[i].imag() * -1.);
					}
				}
			}
		}

	}



	complex<double> getResult() {
		complex<double> temp = output.front();
		output.pop_front();
		output.push_back(0);
		return temp;
	}

};

int main()
{
	ofstream sineout("sine.txt");
	ofstream resout("result.txt");
	ofstream reschk("resultcheck.txt");
	ofstream respcm("res.pcm", ios::binary), sinepcm("sinepcm.pcm", ios::binary);
	complex<double> sine;
	deque<complex<double>> sinedeque;
	const complex<double> jj(0.0, 1.0);

	int numOutput = 2;
	polyphaseTree filter("coeffTable.txt", numOutput);
	polyphaseTree filtercheck("coeffTable.txt", numOutput);


	ifstream file;
	file.open("out_mod_x8_up12000.pcm", std::ios::binary);
	//file.open("sine_tone.dat", std::ios::binary);

	vector<complex<int16_t>> signalsine;

	int16_t realPart, imagPart;
	while (file.read(reinterpret_cast<char*>(&realPart), sizeof(realPart)) &&
		file.read(reinterpret_cast<char*>(&imagPart), sizeof(imagPart))) {
		signalsine.emplace_back(realPart, imagPart);
		//file.write(reinterpret_cast<const char*>(&realPart), sizeof(float));
		//file.write(reinterpret_cast<const char*>(&imagPart), sizeof(float));
		//file.write((const char*)&realPart, sizeof(int16_t));
		//file.write((const char*)&imagPart, sizeof(int16_t));
	}


	vector<complex<double>> signalsine_double;
	for (int i = 0; i < signalsine.size(); ++i) {
		signalsine_double.emplace_back(0., 0.);
		signalsine_double[i].real((double)signalsine[i].real());
		signalsine_double[i].imag((double)signalsine[i].imag());
	}


	//ofstream testInput("test.txt");




	for (double i = 0; i < signalsine.size(); ++i) {




		if ((int)i % numOutput == 0 && i != 0) {
			filtercheck.next(sinedeque);
			sinedeque.clear();
			for (int j = 0; j < numOutput; ++j) {
				complex<double> result = filter.getResult();
				complex<double> resultcheck = filtercheck.getResult();
				if (j == numOutput - 1) {
					reschk << resultcheck.real() << (resultcheck.imag() < 0 ? "" : "+") << resultcheck.imag() << "i" << endl;
				}
				else {
					reschk << resultcheck.real() << (resultcheck.imag() < 0 ? "" : "+") << resultcheck.imag() << "i" << "\t\t";
				}
			}


		}

		//sine = 64. * exp(10000. / 2000000. * i * 2. * M_PI * jj) + 64. * exp(-10000. / 2000000. * i * 2. * M_PI * jj);
		sine = signalsine_double[i];
		if (i < 1000)
			cout << signalsine[i] << endl;
		sinedeque.push_back(sine);

		sineout << sine.real() << (sine.imag() < 0 ? "" : "+") << sine.imag() << "i" << endl;

		//sineout << sine.real() << (sine.imag() < 0 ? "" : "+") << sine.imag() << "i" << endl;

	}

	system("pause");
}
