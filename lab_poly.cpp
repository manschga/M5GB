#include "lab_poly.h"
//#include "assignments.h"
#include <algorithm>
#include <sstream>

using namespace std;
//default constructor
Lab_poly::Lab_poly()
{
  _poly = polynomial();
  _signature = signature();
  _p = 0;
}

//constructor
Lab_poly::Lab_poly(polynomial poly, signature sign)
{
  _poly = poly;
  _signature = sign;
  _p = _poly.Get_p();
}


ostream& operator<<(ostream& s, const Lab_poly& rhs)
{
  cout << "(" << rhs.Get_poly() << "," << rhs.Get_sign()<< ")";
  return s;
}

///*
Lab_poly spol(const Lab_poly& F, const Lab_poly& G)
{
  pair<int,vector<int>> u,v;
  //cout << "in spol-comp!" << endl;
  lcm(u,F.Get_LM(),v,G.Get_LM());
  return F.Multiply(u.second).Multiply(u.first,true) - G.Multiply(v.second).Multiply(v.first,true);
}
//*/

pair<Lab_poly,Lab_poly> spair(const Lab_poly& F, const Lab_poly& G)
{
  #ifdef Null
  cout << "Spair construction: " << endl;
  cout << "F: " << F << endl;
  cout << "G: " << G << endl;
  #endif
  pair<int,vector<int>> u,v;
  //  monomial u,v;
  lcm(u,F.Get_LM(),v,G.Get_LM());
  //for case p = 2 works, else need coeff!
  #ifdef Null
  cout << "uF: " << F.Multiply(u.second).Multiply(u.first,true);
  cout << "vG:"  << G.Multiply(v.second).Multiply(v.first,true);
  #endif
  return make_pair(F.Multiply(u.second).Multiply(u.first,true),G.Multiply(v.second).Multiply(v.first,true));
}

vector<polynomial> poly(const vector<Lab_poly>& G)
{
  vector<polynomial> basis = {};
  for (int g = 0; g < G.size(); ++g)
    {
      basis.emplace_back(G[g].Get_poly());
    }
  return basis;
}

void sign_sort(vector<pair<Lab_poly,Lab_poly>>& P)
{
  sort(P.begin(), P.end(),
       [](const pair<Lab_poly,Lab_poly> &one , const pair<Lab_poly,Lab_poly> &two)
       {
	 return !Smaller_pair(one,two);//compare signatures
       });
  
  return;
}

bool Smaller_pair(const pair<Lab_poly,Lab_poly>& one, const pair<Lab_poly,Lab_poly>& two)
{
  if (one.first.Get_sign() == two.first.Get_sign())
    {
      if (one.second.Get_sign() == two.second.Get_sign())
	{
	  return (one.second.Get_sign() < two.second.Get_sign());
	}
    }
  //else
  return (one.first.Get_sign() < two.first.Get_sign());
}

//test
/*
Flag_poly::Flag_poly()
{
  const polynomial poly;
  const signature sign;
  Flag_poly FP(poly,sign);
  _mac_poly = FP._mac_poly;
  _Flag = FP._Flag;
}

Flag_poly::Flag_poly(const polynomial& poly, const signature& flag)
{
  //determine degree of polynomial
  vector<int> LT = poly.Get_LT();
  int d = 0;
  for (int i = 0; i < LT.size(); ++i)
    {
      d += LT[i];
    }
  int max_size = lower_degrees(d+1,LT.size())+1; //of degree d
  
  if (poly.Get_monomials().size() == 0)
    cout << "Error in Flag_poly construction! Empty poly!" << endl;
  int n = poly.Get_monomials()[0].Get_term().size();
  
  _mac_poly = vec_from_poly(poly, n, max_size);
  _Flag = term_to_int(flag.Get_term(),n);
}
*/




//old one
///*
Flag_poly::Flag_poly()
{
  _Lab_poly = Lab_poly();
  _Flag = signature();
}

//constructor
Flag_poly::Flag_poly(Lab_poly poly, signature flag)
{
  _Lab_poly = poly;
  _Flag = flag;
}
//*/






