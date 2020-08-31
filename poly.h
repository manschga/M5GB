#ifndef POLY_H_INCLUDED
#define POLY_H_INCLUDED

#include<vector>
#include <algorithm>
#include "other_fct.h"

using namespace std;

extern Table table;



class monomial{
 public:
  //const vector<pair<vector<int>,vector<int>>>& table
  //= Create_table(n,d);
  //default constructor
  monomial();
  //constructor
  //monomial(const vector<int>& term, int coeff, int p);
  
  monomial(const term_int term, const int coeff, const int p);
  
  term_int Get_term() const
  {
    return this->_term;
  }
  int Get_coeff() const
  {
    return this->_coeff;
  }
  int Get_p() const
  {
    return this->_p;
  }
  
  /*
    int Get_degree() const
    {
    int d = 0;
    for (int i = 0; i < this->_term.size(); ++i)
    {
    d += this->_term[i];
    }
    }
  */
  
  bool operator < (const monomial &rhs) const
  {
    if (this->_term < rhs._term)
    {
      return true;
      }
    if (this->_term == rhs._term)
      {
	//same terms
	return (this->_coeff < rhs._coeff);
      }
    return false;
}

bool operator > (const monomial &rhs) const
{
  if (this->_term > rhs._term)
    {
      return true;
      }
  if (this->_term == rhs._term)
    {
      //same terms
	return (this->_coeff > rhs._coeff);
      }
  return false;
}

 bool operator == (const monomial &rhs) const
 {
   bool one = (*this > rhs);
   bool two = (*this<rhs);
   bool result = (!(*this>rhs)) && (!(*this<rhs));
   return result;
 }

 

 monomial Multiply(const term_vec& t) const//multiply with term t in vector formulation
{
  /*
  cout << "new multiplication!" << endl;
  cout << "t in multply: " << endl;
  for (int i = 0; i < t.size(); ++i)
    {
      cout << t[i] << ",";
    }
  cout << endl;
  */
  int result = this->_term; //monomial in integer format
  for (int i = 0; i < t.size(); ++i)
	{
	  for (int j = 0; j < t[i]; ++j) //read up from mult table
	    {
	      //cout << "result before multiplication: " << result <<  endl;
	      result = table[result].second[i]; //multiply with x_i+1
	      //cout << "result after multiplication: " << result << endl;
	  // cout << "before term" << new_term[i] << endl;
	  //cout << "t[i]" << t[i] << endl;
	      //	  new_term[i] += t_vec[i];
	  //cout << "after term" << new_term[i] << endl;
	    }
	}
  return monomial(result,this->_coeff,this->_p);
}

 monomial Multiply(const term_int t, bool is_constant) const//multiply with constant / term, can remove table for constant case
{
  if (is_constant)//multiply with constant
    {
      return monomial(this->_term,(this->_coeff*t) % (this->_p),this->_p);
    }
  else {
    term_vec t_vec = table[t].first; //get vector repr of t
    term_int result = this->_term; //monomial in integer format
    for (int i = 0; i < t_vec.size(); ++i)
      {
	for (int j = 0; j < t_vec[i]; ++j) //read up from mult table
	  {
	    
	    result = table[result].second[i]; //multiply with x_i+1
	  }
      }
    return monomial(result,this->_coeff,this->_p);
  }
}

/*
monomial operator*(const monomial m) const//multiply with monomial
{
  
  vector<int> new_term = this->_term;
  //add term powers!
  if (m._term.size() != new_term.size())
    {
      cout << "Multiplying monomials of different size!" << endl;
    }
  for (int i = 0; i < this->_term.size(); ++i)
    {
      new_term[i] += m._term[i];
    }
  return monomial(new_term,(this->_coeff * m._coeff) % (this->_p),this->_p);
}
*/
   
   friend ostream& operator<< (ostream& s, const monomial& rhs);

   //friend ostream& operator<<(ostream& s, const monomial& rhs);

 private:
  term_int _term; //maybe long int
  int _coeff;
  int _p;
};

ostream& operator<<(ostream& s, const monomial& rhs);

//typedef set<monomial> Poly;

class polynomial{
 public:
  //default constructor
  polynomial();

  //constructor
  polynomial(const vector<monomial>& mon);

  //zero poly constructor
  polynomial(int n, int p);

  /*
  int Get_degree() const
  {
    if ((this->_monomials).size() == 0)
      {
	  vector<int> empty = {};em
	  cout << "error, zero polynomial!" << endl;
      }
    return (this->_monomials)[0].Get_degree();
  }
  */
  
  const monomial& Get_LM() const
    {
      if ((this->_monomials).size() == 0)
	{
	  vector<int> empty = {};
	  //std::cout << "error, zero polynomial!" << endl;
	  //	  return empty;
	}
      //    return (this->_coeffs).size()-1;
      return (this->_monomials)[0];
    }

  

  int Get_LT() const
    {
      return (this->_monomials)[0].Get_term();
    }

  int Get_LC() const
    {
      return ((this->_monomials)[0]).Get_coeff();
    }

  int Get_p() const
  {
    return (this->_monomials)[0].Get_p();
  }

  const vector<monomial>& Get_monomials() const
    {
      return (this->_monomials);
    }

  int size() const
  {
    return (this->_monomials.size());
  }

  polynomial Multiply(const term_vec& t) const//multiply with term t in vector formulation
  {
    vector<monomial> mons = this->_monomials;
    for (int i = 0; i < mons.size(); ++i)
      {
	//cout << "multiplication takes place!" << endl;
	mons[i] = mons[i].Multiply(t); //is_constant
	//determines whether term or constant multiplication
      }
    return polynomial(mons);
  }
    
  polynomial Multiply(const int t, bool is_constant) const//multiply with constant / term
      {
	//if (is_constant)
	//{
	    vector<monomial> mons = this->_monomials;
	    for (int i = 0; i < mons.size(); ++i)
	      {
		mons[i] = mons[i].Multiply(t,is_constant); //is_constant
		//determines whether term or constant multiplication
	      }
	    return polynomial(mons);
	    }
	    // }
	    /*
	else
	  {
	    vector<monomial> mons = this->_monomials;
	    vector<monomial> product = {};
	    for (int i = 0; i < mons.size(); ++i)
	      {
		//cout << "before: " << endl;
		//mons[i] = mons[i]*t;
		product.emplace_back((mons[i]*t));
		//cout << "after" << endl;
		//cout << (mons[i])*t << endl;
	      }
	    //cout << "got here" << endl;
	    polynomial test(product);
	    //cout << "after test" << endl;
	    return test;
	    //      return polynomial(mons);
	  }
      }
	    */


/*
    polynomial operator*(const vector<int> t) const//multiply with term t
    {
      //cout << "polynomial level" << endl;
    }
*/
/*
    polynomial operator*(const monomial m) const//multiply with monomial m
    {
      vector<monomial> mons = this->_monomials;
      for (int i = 0; i < mons.size(); ++i)
	{
	  mons[i] = mons[i]*m;
	}
      return polynomial(mons);
    }
*/

  polynomial operator+(const polynomial& rhs) const
  {
    //polynomial test_1 = *this;
    //polynomial test_2 = rhs;
    //cout << "Add " << test_1 << "and" << test_2 << endl;

    //compute *this+rhs = f+g (sum of polynomials)
    vector<monomial> f_mon = this->_monomials;

    if (this->Get_LC() == 0)
      {
	return rhs;
      }
    if (rhs.Get_LC() == 0)
      {
	return *this;
      }
    
    vector<monomial> g_mon = rhs._monomials;
    reverse(f_mon.begin(),f_mon.end());
    reverse(g_mon.begin(),g_mon.end());
    vector<monomial> h_mon = {};

    int f_i = f_mon.size()-1;
    int g_i = g_mon.size()-1;
    h_mon.reserve(f_i+g_i+2);
    vector<monomial> merged(f_i+g_i+2);
    std::merge(f_mon.begin(),f_mon.end(),g_mon.begin(),g_mon.end(),merged.begin());

    if (merged.size() > 1)
      {
	for (int i = 0; i < merged.size()-1; ++i)
	  {
	    //if (count(merged[i], merged[i+1], merged[i+1].Get_term()) == 1) 
	    if (merged[i].Get_term() != merged[i+1].Get_term())
	      {
		h_mon.push_back(merged[i]);
	      }
	    else
	      {
		++i;// do sum of both in one step
		//push_back sum of coeffs, case p != 2
	      }
	  }
	{
	  if (merged[merged.size()-2].Get_term() < merged[merged.size()-1].Get_term())
	    {
	      h_mon.push_back(merged[merged.size()-1]);
	    }
	}
      }
    else
      {
	h_mon = merged;
      }
    reverse(h_mon.begin(),h_mon.end());
    /*
    
   
    while (f_i >= 0 || g_i >= 0)//go through all terms of f and g
      {
	if (f_i < 0) //all terms of f done, avoid index crash!
	{
	  h_mon.emplace_back(g_mon[g_i]);
	  g_i--;
	}
      else if (g_i < 0)
	{
	  h_mon.emplace_back(f_mon[f_i]);
	  f_i--;
	}
	else
	{
	  if (f_mon[f_i].Get_term() == g_mon[g_i].Get_term())
	    {//same term appears in both, add up coefficients
	      h_mon.emplace_back(monomial(f_mon[f_i].Get_term(), (f_mon[f_i].Get_coeff() + g_mon[g_i].Get_coeff())%rhs._p, f_mon[f_i].Get_p()));
	      f_i--;
	      g_i--;
	    }
	  else
	    {
	      if (f_mon[f_i].Get_term()< g_mon[g_i].Get_term())
		{
		  h_mon.emplace_back(f_mon[f_i]);
		  f_i--;
		}
	      else
		{
		  h_mon.emplace_back(g_mon[g_i]);
		  g_i--;
		}
	    }
      }
      }
    */
    if (h_mon.size() == 0)
      {
	int p = this->Get_p();
	polynomial h({monomial(0,0,p)});
	//cout << "equals" << h << endl;
	return h;
      }
    polynomial h(h_mon);
    //cout << "in addtion, p:" << h.Get_p() << endl;
    //h = h.update();
    //reorder, remove zeros, do a function!
    return h;
  }

  polynomial operator-(const polynomial& rhs) const
  {
    
    //construct -rhs and call f + -rhs
    vector<monomial> minus_mon_rhs = {};
    for (int i = 0; i < rhs.Get_monomials().size(); ++i)
      {
	minus_mon_rhs.emplace_back(monomial(rhs.Get_monomials()[i].Get_term(),(rhs.Get_p()-rhs.Get_monomials()[i].Get_coeff()) % rhs.Get_p(),rhs.Get_monomials()[i].Get_p()));
      }
    polynomial minus_rhs(minus_mon_rhs);
    
    //cout << "succesful until subtraction!" << endl;
    return (*this + minus_rhs);
  }

  bool operator == (const polynomial &rhs) const
  {
    if (this->Get_monomials().size() > 0 && rhs.Get_monomials().size()>0)
      {
	if (! (this->Get_LM() == rhs.Get_LM()))
	  {
	    return false;
	  }
	if (this->Get_LC() == 0 && rhs.Get_LC() == 0)
	  return true;
	polynomial lhs_tail = *this - polynomial({this->Get_LM()});
	polynomial rhs_tail = rhs - polynomial({rhs.Get_LM()});
	return (lhs_tail == rhs_tail);
      }
    else
      {
	if (this->Get_monomials().size() == 0 && rhs.Get_monomials().size() == 0)
	  {
	    return true;
	  }
      }
    return false;
  }

   bool operator < (const polynomial &rhs) const
  {
    if (this->Get_monomials().size() > 0 && rhs.Get_monomials().size()>0)
      {
	if (this->Get_LM() == rhs.Get_LM())
	  {
	    polynomial lhs_tail = *this - polynomial({this->Get_LM()});
	    polynomial rhs_tail = rhs - polynomial({rhs.Get_LM()});
	    return (lhs_tail < rhs_tail);
	  }
	else
	  {
	    return this->Get_LM() < rhs.Get_LM();
	  }
      }
    return false;
  }

  const polynomial update()//removes monomials with 0 coeff and reorders in decreasing terms
  {
    polynomial temp = *this;
    vector<monomial> temp_monom = temp._monomials;
    vector<monomial> new_vec = {};
    for (int i = 0; i < temp_monom.size(); ++i)
      {
	if (temp_monom[i].Get_coeff() != 0)
	  {
	    new_vec.emplace_back(temp_monom[i]);
	    //	    temp_monom.erase(temp_monom.begin() + i);
	  }
      }
    ///*
    //cout << "update called!" << endl;
    sort(new_vec.begin(), new_vec.end(),
       [](const monomial &one , const monomial &two)
       {
	 return (one.Get_term() > two.Get_term());
       });
    //cout << "end of update!" << endl;
    //*/

    /*
    //cout << "got also here!" << endl;
    if (new_vec.size() == 0)
      {
	//cout << "attention! zero polynomial!" << endl;
	int n = temp_monom[0].Get_term().size();
	vector<int> zero_pol(n,0);
	int p = temp_monom[0].Get_p();
	new_vec.emplace_back(monomial(zero_pol,0,p));
      }
    */
    if (new_vec.size() == 0)
      {
	int p = temp_monom[0].Get_p();
	new_vec.emplace_back(monomial(0,0,p));
      }
    //cout << "end of update!" << endl;
    return polynomial(new_vec);
  }
    /*
    cout << "update called! " << endl;
    for (int i = 0; i < this->_monomials.size(); ++i)
      {
	if ((((this->_monomials)[i]).Get_coeff()) == 0)
	  {
	    this->_monomials.erase(this->_monomials.begin() + i);
	  }
      }
  
  /*
~polynomial() 
    { 
        cout << "Destructor for poly called " << endl;  
    }
  */ 
  
  friend ostream& operator<<(ostream& s, const polynomial& rhs);
  
 private:
  vector<monomial>_monomials; //all terms with non-zero coefficient
  int _p;
};

  ostream& operator<<(ostream& s, const polynomial& rhs);

void lcm(pair<int,vector<int>>& u, const monomial& LM_f, pair<int,vector<int>>& v, const monomial& LM_g);

polynomial spol(const polynomial& f, const polynomial& g);

vector<int> get_remainder(vector<int> term, vector<int> divisor);//if divisor divides term, then term = divisor*remainder

vector<polynomial> Create_poly_system(const string& input_file, int n, int p);

polynomial poly_from_vec(const vector<int>& vec, int n, int p);

void reduce(polynomial &G_i, const polynomial &G_j, bool &changed, int &reduction);

//void reduce(polynomial &G_i, const polynomial &G_j, bool &changed);

int reduce(vector<polynomial>& G);

polynomial Tail(const polynomial& f); //returns all except leading monomial


vector<vector<int>> vec_from_poly(const polynomial& poly, int n);

void sort(vector<polynomial>& basis);

bool is_Smaller(const polynomial &one , const polynomial &two);


#endif
