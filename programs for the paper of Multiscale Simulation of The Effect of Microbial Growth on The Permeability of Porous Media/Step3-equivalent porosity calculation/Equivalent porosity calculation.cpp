# include<iostream> 
# include<cmath>
#include<math.h>
# include<cstdlib>
# include<iomanip>
# include<fstream>
# include<sstream>
# include<string>
# include <stdio.h>
# include <stdlib.h>


using namespace std;
const int N = 196;//Look at the "Permeability from LBM.dat" how many data there are in this file , put two more

double f(double x);
double k[200];
double kLBM[N];
double epsilon[N];
int i;

int main()
{


	for (int i = 0; i <= 199; i++)//199
		k[i] = f(i * 1.0 / 200);

	int p2;

	ifstream fin1("Permeability from LBM.dat");
	for (p2 = 0; p2 < N; p2++)
		fin1 >> kLBM[p2];

	epsilon[0] = 0;

	for (int j = 0; j < N; j++)
	{
		int go = 1;
		i = 0;
		do
		{
			i = i + 1;
			if (k[i] > kLBM[j]) {
				go = 0;
				epsilon[j] = ((k[i] - kLBM[j]) * (i - 1) * 1.0 / 200 + (kLBM[j] - k[i - 1]) * (i) * 1.0 / 200) / (k[i] - k[i - 1]);//Linear interpolation
			}

		} while (go == 1 && i <= 199);//199
	}
	ofstream   out1("The calculated equivalent porosity.dat ", ios::app);
	for (i = 0; i < N; i++)
		out1 << epsilon[i] << endl;

	return 0;

}

double f(double x)
{
	double d = 1.95E-4, c = 1.0 / 150;
	double k = pow(x, 3) / pow(1 - x, 2) * d * d * c;
	return k;
}