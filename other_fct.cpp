#include "other_fct.h"
//#include "assignments.h"
#include <algorithm>
#include <sstream>

Table table;
Binom_table binom_table;


using namespace std;

vector<int> Subtract(const vector<int>& f, const vector<int>& g,int p)//computes f-g over F_p
{
  int max_size = max(f.size(),g.size());
  int min_size = max(f.size(),g.size());
  vector<int> result(max_size,0);
  for (int i = 0; i < min_size; ++i)
    {
      result[i] = (p+f[i] - g[i]) % p;
    }
  if (f.size() > g.size())
    {
      for (int i = g.size(); i < max_size; ++i)
	{
	  result[i] = f[i];
	} 
    }
  else
    {
       for (int i = g.size(); i < max_size; ++i)
	{
	  result[i] = (p-g[i]) % p;
	} 
    }
  return result;
}

vector<int> Add(const vector<int>& f, const vector<int>& g,int p)//computes f-g over F_p
{
  int max_size = max(f.size(),g.size());
  int min_size = max(f.size(),g.size());
  vector<int> result(max_size,0);
  for (int i = 0; i < min_size; ++i)
    {
      result[i] = (f[i] + g[i]) % p;
    }
  if (f.size() > g.size())
    {
      for (int i = g.size(); i < max_size; ++i)
	{
	  result[i] = f[i];
	} 
    }
  else
    {
       for (int i = g.size(); i < max_size; ++i)
	{
	  result[i] = g[i];
	} 
    }
  return result;
}

bool term_is_Zero(vector<int> term)
{
  for (int i = 0; i < term.size(); ++i)
    {
      if (term[i] > 0)
	{
	  return false;
	}
    }
  return true;
}

int LT(const vector<int>& poly)
{
  for (int i = poly.size()-1; i >= 0; --i)
    {
      if (poly[i] >0)
	return i;
    }
  cout << "Zero polynomial!" << endl;
  return 0;
}

vector<int> Tail(const vector<int>& poly)
{
  vector<int> tail = poly;
  //improve this by LT extra!
  int LT_ = LT(poly);
  tail[LT_] = 0;
  return tail;
}

int Inverse(int c, int p)//compute c^-1 in F_p, can do a table for large fields
{
  int c_copy = c;
  if (c == 0)
    {
      cout << "Next polynomial added to input set!" << endl;
      //cout << "Error! No inverse to 0 exists!" << endl;
      return 0;
    }
  for (int i = 1; i < p; ++i)
    {
      if (c == 1)
	{
	  return i;
	}
      else
	{
	  c = (c + c_copy)%p;
	}
    }
  cout << "Error, no inverse found! Sure that p is prime?" << endl;
  return 0;
}
int binom(int n, int k)   //from geeksforgeeks.org
{
  if (n < k)
    return 0;
  if (k < 0 || n < 0)
    return 0;
    int res = 1;  
  
    // Since C(n, k) = C(n, n-k)  
    if ( k > n - k )  
        k = n - k;  
  
    // Calculate value of  
    // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]  
    for (int i = 0; i < k; ++i)  
    {  
        res *= (n - i);  
        res /= (i + 1);  
    }  
  
    return res;  
}  

/*
int binom(int n,int k)
{
  int numerator = 1;
  int denominator = 1;
  //cout << "n: " << n << endl;
  if (n < k)
    return 0;
  if (k < 0 || n < 0)
    return 0;
  for (int i = n-k+1; i <= n; ++i)
    {
      numerator *= i;
    }
  for (int i = 2; i <= k; ++i)
    {
      denominator *= i;
      //cout << "denominator in loop: " << denominator << endl;
    }
  //cout << "numerator: " << numerator << endl;
  //cout << "denominator: " << denominator << endl;
  /*
  if (numerator % denominator != 0)
    {
      cout << "Error in binom computation!!!" << endl;
    }

  //cout << "nearly end of binom!" << endl;
  return (numerator/denominator);
}
*/

vector<vector<term_int>> Create_binom_table(int n, int d)
{
  vector<vector<term_int>> Binom_table = {};
  for (int i = 0; i <= n; ++i)
    {
      vector<term_int> table_n = {};
      for (int j = 0; j <= d; ++j)
	{
	  table_n.push_back(binom(i,j));
	}
      Binom_table.push_back(table_n);
    }
  return Binom_table;
}


term_int term_to_int(term_vec term, int n) //computes entry in macaulay matrix
{
  //cout << "into term_to_int" << endl;
  //get degree first
  int d = 0;
  for (int i = 0; i < term.size(); ++i)
    {
      d += term[i];
    }
  //cout << "degree: " << d << endl;
  //int mac = binom(n+d-1,n);//get all polynomials with degrees up to d-1
  int mac = 0;
  if (n+d > 0)
    {
      mac += (binom_table[n+d-1])[n];//binom(n+d-1,n);//get all polynomials with degrees up to d-1
    }
  //cout << "after lower degrees: " << mac << endl;
  mac ++; //go to x_n^d

  for (int i = n-1; i >= 0; i--)
    {
      d -= term[i]; //get degree in x_1,...,x_i
      if (i+d > 0 && d > 0)
	{
	  mac += (binom_table[i+d-1])[d-1];//binom(i+d-1,d-1); //compute # of polynomials of degree d-1 in (x_1,...,x_i,x_i+1)
	  //mac += binom(i+d-1,d-1);
	}
    }
  //left here!
  return mac;
}

vector<int> create_degrees(int d, int n)
{
  vector<int> degrees = {};
  for (int i = 0; i < d; ++i)
    {
      degrees.push_back(lower_degrees(i,n));
    }
return degrees;
}

int lower_degrees(int d, int n) //gets the number of x_1^(d-1) (largest term of degree d-1), is highly inefficient!
{
  //cout << "into lower degrees!" << endl;
  int number = 0;
  for (int i = 0; i < d; ++i)
    {
      number += binom(n+i-1,i);
    }
  return number;
}


bool is_Smaller(const vector<int>& t, const vector<int>& u) //checks if t == u for t,u terms, order is degrevlex
{
  //cout << "is smaller called! " << endl;
  /*for (int i = 0; i < u.size(); ++i)
    {
      cout << u[i] << " ";
    }
  cout << endl;
  cout << "t:";
  for (int i = 0; i < t.size(); ++i)
    {
      cout << t[i] << " ";
    }
  cout << endl;
  */

  if (t.size() != u.size())
    {
      cout << "Error! wrong sizes in is_smaller" << endl;
      cout << "t.size: " << t.size() << endl;
      cout << "u.size: " << u.size() << endl;
    }
  int deg_t = 0;
  int deg_u = 0;
  for (int i = 0; i < u.size(); ++i)
    {
      deg_t += t[i];
      deg_u += u[i];
    }
   
  if (deg_t < deg_u)
    return true;
  if (deg_t > deg_u)
    return false;

  //case deg(t) = deg(u)
  int n = t.size();
  for (int i = n-1; i >= 0; --i)
    {
      if (t[i] > u[i])
	return true;
      if (t[i] < u[i])
	return false;
    }
  //case equal
  return true;
}


bool is_Equal(const vector<int>& t, const vector<int>& u) //checks if t < u for t,u terms, order is deglex
{
  if (t.size() != u.size())
    {
      cout << "Error! wrong sizes in is_equal!" << endl;
    }

  int n = t.size();
  for (int i = 0; i < t.size(); ++i)
    {
      if (t[i] != u[i])
	return false;
    }
  return true;//they are equal
}



//need to change!
void lcm(vector<int>& u, const vector<int>& LT_f, vector<int>& v, const vector<int> LT_g)
{
  int n = LT_f.size();
  u.assign(n, 0);
  v.assign(n, 0);
  for (int i = 0; i < n; ++i)
    {
      if (LT_f[i] > LT_g[i])
	{
	  v[i] = LT_f[i] - LT_g[i];
	}
      else
	{
	  u[i] = LT_g[i] - LT_f[i];
	}
    }
  return;
}

/*
void reorder_poly(vector<vector<int>>& f)
{
  //reorder such that terms are ordered in decreasing way!
  sort(f.begin(), f.end(),
       [](const vector<int> &one , const vector<int> &two)
       {
	 return !is_Smaller(one,two);
       });
}
*/

bool isDivisible(const vector<int>& term, const vector<int>& divisor)
{
  /*
  if (term.size() != divisor.size())
    {
      cout << "term size: " << term.size() << endl;
      cout << "divisor size: " << divisor.size() << endl;
      cout << "Error!" << endl;
    }
  */
  int n = term.size();
  for (int i = 0; i < n; ++i)
    {
      if(term[i] < divisor[i])
	{
	  return false;
	}
    }
  return true;
}


  
vector<pair<term_vec,vector<term_int>>> Create_table(int n, int d)//multiplication table
{
  vector<pair<term_vec,vector<term_int>>> table = {};
  term_int max_size = binom(n+d-1,n)+1;
  table.reserve(max_size);
  term_vec term(n,0);
  vector<term_int> multis(n,0);
  for (int j = 0; j < n; ++j)
    {
      multis[j] = n-j;
    }
  table.emplace_back(make_pair(term,multis));
  for (term_int i = 1; i < max_size; ++i)
    {
      next_term(term);
      for (int j = 1; j <= n; ++j)
	{
	  term_vec term_xj = vector_add(term,j); //new term = term*x_j
	  multis[j-1] = term_to_int(term_xj,n)-1; //might be too slow!,
	  //check for over size!
	}
      table.emplace_back(make_pair(term,multis));
    }
  return table;
}

term_vec vector_add(const term_vec& term, int j) //new term = term*x_j
{
  vector<int> new_term = term;
  new_term[j-1] ++;
  return new_term;
}

void next_term(vector<int>& term) //degree revers lex, change if needed!
{
  //get degree
  int degree = 0;
  for (int i = 0; i < term.size(); ++i)
    {
      degree += term[i];
    }
  if (degree == 0) //term = 1
    {
      term[term.size()-1] = 1; //new term = x_n
      return;
    }
  //get frist non-zero entry which is not x_1
  int non_zero_idx = 0;
  for (int i = 1; i < term.size(); ++i)
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
  int temp = term[0];
  term[0] = 0;
  term[non_zero_idx] --;
  term[non_zero_idx-1] ++;
  term[non_zero_idx-1] += temp;

      //}

  return;
}



  

