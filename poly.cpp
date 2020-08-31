#include "poly.h"
#include "other_fct.h"
#include <algorithm>
#include <sstream>
#include<vector>
#include <fstream>
#include <string>


using namespace std;

monomial::monomial()
{
  _term = {};
  _coeff = 0;
  _p = 0;
}

monomial::monomial(const term_int term, const int coeff, const int p)
{
  _term = term;
  _coeff = coeff;
  _p = p;
}

//for singular
ostream& operator<<(ostream& s, const monomial& rhs)
{
  //if (rhs._coeff > 1)
  cout << rhs._coeff << "*";
  cout << rhs._term;
  /*
  for (int i = 0; i < rhs._term.size()-1; ++i)
    {
      if (rhs._term[i] > 1)
	cout << char(97+i) << rhs._term[i];
      if (rhs._term[i] == 1)
	cout << char(97+i);
    }
  if (rhs._term[rhs._term.size()-1] > 1)
    cout << char(97+rhs._term.size()-1) << rhs._term[rhs._term.size()-1];
    if (rhs._term[rhs._term.size()-1] == 1)
      cout << char(97+rhs._term.size()-1);
  */
  return s;
}


/* for sage
ostream& operator<<(ostream& s, const monomial& rhs)
{
  //if (rhs._coeff > 1)
  cout << rhs._coeff;
  for (int i = 0; i < rhs._term.size()-1; ++i)
    {
      if (rhs._term[i] > 1)
	cout << "*" << "x_" << i+1 << "^" << rhs._term[i];
      if (rhs._term[i] == 1)
	cout << "*" << "x_" << i+1;
    }
  if (rhs._term[rhs._term.size()-1] > 1)
    cout << "*" << "x_" << rhs._term.size() << "^" << rhs._term[rhs._term.size()-1];
    if (rhs._term[rhs._term.size()-1] == 1)
      cout << "*" << "x_" << rhs._term.size();
  return s;
}
*/

polynomial::polynomial()
{
  _monomials = {};
  _p = 0;
}

polynomial::polynomial(int n, int p)//construct zero polynomial
{
  //  vector<int> zeros(n,0);
  _monomials = {monomial(0,0,p)};
  _p = p;
}


polynomial::polynomial(const vector<monomial>& mon)
{
  _monomials = mon;
  _p = mon[0].Get_p();
}

ostream& operator<<(ostream& s, const polynomial& rhs)
{
  for (int i = 0; i < rhs.size()-1; ++i)
    {
      cout << rhs._monomials[i] << "+";
    }
  cout << rhs._monomials[rhs.size()-1] <<",";
  return s;
}

void lcm(pair<int,vector<int>>& u, const monomial& LM_f, pair<int,vector<int>>& v, const monomial& LM_g)
{
  vector<int> term_f = table[LM_f.Get_term()].first;
  vector<int> term_g = table[LM_g.Get_term()].first;
  vector<int> term_u,term_v; //terms for u,v
  lcm(term_u,term_f,term_v,term_g);//get lcm on term level
  int c_u = Inverse(LM_f.Get_coeff(),LM_f.Get_p());//construct coeff of u,v
  int c_v = Inverse(LM_g.Get_coeff(),LM_g.Get_p());
  u = make_pair(c_u,term_u); //monomials in vector form
  v = make_pair(c_v,term_v);
  return;
}

/*
polynomial spol(const polynomial& f, const polynomial& g)//needed?
{
  monomial LM_f = f.Get_LM();
  monomial LM_g = g.Get_LM();
  pair<int,vector<int>> u,v;
    lcm(u,LM_f,v,LM_g);
    polynomial(f*u-g*v);
}
*/

term_vec get_remainder(term_vec term, term_vec divisor)//if divisor divides term, then term = divisor*remainder, is only on term level!
{
  //cout << "in remainder " << endl;
  //  int rem_coeff = mon.Get_coeff()*Inverse(divisor.Get_coeff(),divisor.Get_p());
  /*
  cout << "term_size: " << term.size() << endl;
  cout << "divisor size: " << divisor.size() << endl;
  cout << "in remainder " << endl;
  */
  term_vec rem_vec = {};
  if (term.size() != divisor.size())
    {
      cout << "Error in Remainder computation!" << endl;
    }
  for (int i = 0; i < term.size(); ++i)
    {
      //cout << "before emplace " << endl;
      rem_vec.emplace_back(term[i]-divisor[i]);
      //cout << "after emplace " << endl;
    }
  return rem_vec;
}

vector<polynomial> Create_poly_system(const string& input_file, int n, int p)
{
  vector<polynomial> system = {};
  vector<int> input_vec;
  string line;
  ifstream myfile (input_file);
  if (myfile.is_open())
  {
    while (getline (myfile,line))
    {
      input_vec = {};
      int i = 0;
      while (i < line.size())
	{
	  input_vec.emplace_back(line[i]-48); //from char to int
	  i = i+2;
	}
      input_vec.pop_back();
      system.emplace_back(poly_from_vec(input_vec,n,p));
      //cout << "test function:" << poly_from_vec(input_vec,n,p);
    }
    myfile.close();
  }
  else cout << "Unable to open file";
  return system;
}


polynomial poly_from_vec(const vector<int>& vec, int n, int p)
{
  //  monomial constant(0,vec[vec.size()-1],p);
  vector<monomial> monomials = {};
  for (int i = 0; i < vec.size(); ++i)
    {
      if ((vec[i])%p > 0)
	{
	  monomials.emplace_back(monomial(vec.size()-i-1,(vec[i]),p));//new term with coefficient
	}
    }
  return polynomial(monomials);//.update();
}

//case p = 2
vector<vector<int>> vec_from_poly(const polynomial& poly, int n)
{
  vector<vector<int>> poly_vec = {};
  //assume poly is correctly ordered, else call poly.update() -> no more const
  //first extra, get size
  if (poly.size() == 0)
    {
      cout << "empty polynomial!" << endl;
      return {};
    }
  //  vector<int> vec(max_size,0);
  for (int i = 0; i < poly.Get_monomials().size(); ++i)
    {
      if (poly.Get_monomials()[i].Get_coeff() > 0)
	{
	  vector<int> term_i = table[poly.Get_monomials()[i].Get_term()].first;
	  poly_vec.emplace_back(term_i);
	}
    }
  return poly_vec;
}
//*/

  
  /*if (term[0] == 0)
    {
  //else term = x_l^i*....x_n^r
      term[non_zero_idx] --;
      term[non_zero_idx-1] ++;//new term = x_1^k*x_(l-1)x_l^(i-1)*....x_n^r;//term[non_zero_idx-1] - 1; //new term = x_1^k*x_(l-1)x_l^(i-1)*....x_n^r
      return;
      }*/
  //else
  //{
  


  
  /*deglex!
  //get last non-zero entry
  int non_zero_idx = 0;
  for (int i = term.size()-1; i >= 0; --i)
    {
      if (term[i] > 0)
	{
	  non_zero_idx = i;
	  break;
	}
    }
  if (non_zero_idx == 0)//term = x_1^degree
    {
      term[0] = 0;
      term[term.size()-1] = degree+1; //new term = x_n^(degree+1)
      return;
    }
  //else term = x_1^.....x_(l-1)^(k)*x_l^(m)
  term[non_zero_idx-1] ++;
  int temp = term[non_zero_idx];
  term[non_zero_idx] = 0;
  term[term.size()-1] = temp-1; //new term = x_1^.....x_(l-1)^(k+1)*x_n^(m-1)
  */

int reduce(vector<polynomial>& G)
{
  int reductions = 0;
  int size = 0;
  int size_before = 1;
  while(size != size_before)
    {
      size_before = G.size();
      
      bool changed = true;
      while (changed)
	{
	  changed = false;
	  int size = G.size();
	  int i = 0;
	  while (i < G.size())
	    {
	      int j = 0;
	      while (j < G.size())
		{
		  if (j != i)
		    {
		      reduce(G[i],G[j],changed,reductions);//reduce G[i] by all other in G
		      int LT_G = G[i].Get_LT();
		      int LC_G = G[i].Get_LC();
		      if (LT_G == 0 && LC_G == 0) //is zero-polynomial
			{
			  G.erase (G.begin()+i);
			  size--;
			  j = G.size();
			  i = G.size();
			}
		    }
		  ++j;
		}
	      ++i;
	    }
	}
      //lc 1
      /* //works for p = 2
	 for (int i = 0; i < G.size(); ++i)
	 {
	 G[i] = G[i]*Inverse(G[i].Get_LC(),G[i].Get_p());
	 }
      */
      size = G.size();//maybe the size drops
    }
  return reductions;
  //sort(G);
}

void reduce(polynomial &G_i, const polynomial &G_j, bool &changed, int& reductions)
{
  term_vec LT_j_vec = table[G_j.Get_LT()].first;
  for (int terms = 0; terms < G_i.size(); ++terms)
    {
      monomial monom = G_i.Get_monomials()[terms];
      int coeff = monom.Get_coeff();
      vector<int> monom_vec = table[monom.Get_term()].first;
      if (isDivisible(monom_vec,LT_j_vec))
	{
	  changed = true;
	  reductions ++;
	  /*
	  for (int i = 0; i < monom_vec.size(); ++i)
	    {
	      cout << monom_vec[i] << ",";
	    }
	  cout << "LT_G_j " << endl;
	  for (int i = 0; i < LT_j_vec.size(); ++i)
	    {
	      cout << LT_j_vec[i] << ",";
	    }
	  cout << endl;
	  changed = true;
	  /*
	  cout << "G_i before: " << G_i << endl << endl;
	  cout << "G_j: " << polynomial({monomial(get_remainder(monom.Get_term(),G_j.Get_LT()),monom.Get_coeff()*Inverse(G_j.Get_LC(),monom.Get_p()),monom.Get_p())}) << endl << endl;
	  */
	  polynomial vG_j = G_j.Multiply(get_remainder(monom_vec,LT_j_vec)); //term multi
	  //vG_j = G_j.Multiply(monom.Get_coeff()*Inverse(G_j.Get_LC(),monom.Get_p()),true,table);//constant multi
	  int c = coeff*Inverse(vG_j.Get_LC(),vG_j.Get_p());
	  /*


	  cout << "G_i: " << G_i << endl;
	  */
	  G_i = G_i - vG_j.Multiply(c,true);
	  //cout << "G_i after : " << G_i<< endl << endl;
	}
    }
}

polynomial Tail(const polynomial& f) //returns all except leading monomial
{
  
  //cout << "Tail computation begin!" << endl;
  vector<monomial> tails = f.Get_monomials();
  //int n = tails[0].Get_term().size();
  if (tails.size() == 0)
    {
      cout << "Error in Tail! Empty polynomial!" << endl;
    }
  int p = tails[0].Get_p();
  //cout << "size of tails before removing: " << tails.size() << endl;
  tails.erase(tails.begin()); //remove LM
  //cout << "size of tails after removing: " << tails.size() << endl;
  if (tails.size() > 0)
    return polynomial(tails);
  else
    {// cout << "possible seg fault in Tail!" << endl;
      return polynomial({monomial(0, 0, p)});//construct zero polynomial
    }
}

void sort(vector<polynomial>& basis)
{
  sort(basis.begin(), basis.end(),
       [](const polynomial &one , const polynomial &two)
       {
	 return is_Smaller(one,two);
       });
}

bool is_Smaller(const polynomial &one , const polynomial &two)
{
  if (one.Get_LT() < two.Get_LT())
    {
      return true;
    }
  if (two.Get_LT() < one.Get_LT())
    {
      return false;
    }
  //else
  return is_Smaller(Tail(one),Tail(two));
}
