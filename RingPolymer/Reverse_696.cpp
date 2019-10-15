#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main()
{
  ifstream fin;
  fin.open("coord_2.xyz");
  string junk;
  getline(fin,junk);
  getline(fin,junk);
  double c,x[346],y[346],z[346];
  for (int i=0;i<346;i++)
    {
      fin >> c >> x[i] >> y[i] >> z[i];
    }

  ofstream fout;
  fout.open("coord_2_reversed.xyz");

  fout << "346\nAtoms\n";
  for (int i=345;i>=0;i--)
    {
      fout << "1 " << x[i] << "\t" << y[i] << "\t" << z[i] << "\n";
    }
  fin.close();
  fout.close();
  return 0;
}
