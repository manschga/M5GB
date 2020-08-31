#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <ctime>
#include <utility>
#include <cstdlib>
#include <fstream>
#include "other_fct.h"
#include "lab_poly.h"
#include "M5GB.h"
#include <unordered_map>
//#include <bits/stdc++.h>
//#include <boost/functional/hash.hpp>

using namespace std;

//#define PRINT

//need to check for monomials with p != 2
int main(int argc, char *argv[])
{
  int n = 20;
  int d = 10;
  binom_table = Create_binom_table(n+d,n+d);
  cout << "binom table created!" << endl;
  table = Create_table(n,d); //check, maybe neednot d = n!
  cout << "size: " << table.size() << endl;
  int p = 2;
  //const string input_file("Examples/medium_ideal.txt");
  //const string input_file("Examples/test_ideal.txt");
  //const string input_file("Examples/small_test.txt");
  //const string input_file("Examples/test_ideal_orig.txt");
  //const string input_file("Examples/example_2.txt");
  //const string input_file("Examples/big_example.txt");
  const string input_file("Examples/n_20.txt");
  //const string input_file("Examples/ToyExample-type4-n10-seed3");
  //const string input_file("Examples/ToyExample-type4-n15-seed0");
  //const string input_file("katsura_4.txt");
  //const string input_file("Examples/example_31.txt");

  vector<polynomial> system_ = Create_poly_system(input_file,n,p);
  //#define PRINT
  
  #ifdef PRINT
  cout << "Polynomials before reduce:" << endl;
  for (int i = 0; i < system_.size(); ++i)
    {
      vector<vector<int>> poly_i = vec_from_poly(system_[i],n);
	for (int j = 0; j < poly_i.size(); ++j)
	  {
	    bool is_zero = true;
	    for (int k = 0; k < poly_i[j].size(); ++k)
	      {
		if ((poly_i[j])[k]>0 )
		{
		  is_zero = false;
		  cout << char(97+k) << (poly_i[j])[k];
		//		    cout << "x_" << k+1 << "^" << ((poly_i[j])[k] );
		}
	      }
	    if (is_zero)
	      cout << "1";
	    if (j < poly_i.size()-1)
	    cout << "+";
	  }
      //      cout << basis[i] << endl;
	cout << "," << endl;
    }
  #endif

  cout << "M5GB start!!!!!" << endl << endl;

  
  //  Basis solution = M5GB(system_,n);
    
  vector<polynomial> basis = M5GB(system_,n);//Poly(solution,n,p); //remove signatures

  
  /*
  //cout << "not reduced basis:" << endl;
  for (int i = 0; i < basis.size(); ++i)
    {
      basis[i].update();
      //cout << basis[i] << endl;
    }

  int size = 0;
  int size_before = 1;
  size = 0;
  size_before = 1;
  while(size != size_before)
    {
      size_before = basis.size();
      reduce(basis); //get the (unique) reduced gröbner basis
      //cout << "in here! " << endl;
      size = basis.size();//maybe the size drops
    }
  sort(basis);
  */

  sort(basis);
  ///*
  cout << "Gröbner Basis polynomials:" << endl;
  for (int i = 0; i < basis.size(); ++i)
    {
      vector<vector<int>> poly_i = vec_from_poly(basis[i],n);
      for (int j = 0; j < poly_i.size(); ++j)
	  {
	    bool is_zero = true;
	    cout << basis[i].Get_monomials()[j].Get_coeff();
	    for (int k = 0; k < poly_i[j].size(); ++k)
	      {
		if ((poly_i[j])[k]>0 )
		{
		  is_zero = false;
		  cout << char(97+k) << (poly_i[j])[k];
		//		    cout << "x_" << k+1 << "^" << ((poly_i[j])[k] );
		}
	      }
	    if (is_zero)
	      cout << "1";
	    if (j < poly_i.size()-1)
	    cout << "+";
	  }
	cout << "," << endl;
    }
  //*/
 

  //cout << "basis size:" << basis.size() << endl;

     return 0;
}

