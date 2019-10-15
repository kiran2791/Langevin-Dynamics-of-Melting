#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main()
{
  ifstream fin;
  fin.open("coord_6_2.xyz");
  string junk;
  getline(fin,junk);
  getline(fin,junk);
  double c,x[6],y[6],z[6];
  for (int i=0;i<6;i++)
    {
      fin >> c >> x[i] >> y[i] >> z[i];
    }

  ofstream fout;
  fout.open("coord_6_reversed.xyz");

  fout << "6\nAtoms\n";
  for (int i=5;i>=0;i--)
    {
      fout << "1 " << x[i] << "\t" << y[i] << "\t" << z[i] << "\n";
    }
  fin.close();
  fout.close();
  return 0;
}
