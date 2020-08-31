#include "signature.h"
#include <algorithm>
#include <sstream>
#include<vector>

signature::signature()
{
  _term = 0;
  _index = 0;
}

signature::signature(int term, int index)
{
  _term = term;
  _index = index;
}


//construct a signature which is larger than everything else
signature::signature(int m)
{
  _term = 1;
  _index = m+1; //makes sure this signature is greater than all possibly signatures
}

const signature& max(const signature& one, const signature& two)
{
  if (Smaller_sign (one,two))
    {
      return two;
    }
  //else
  if (Smaller_sign (two, one))
    {
      return one;
    }
  //else same signature!
  //cout << "Possible signature drop!" << endl;
  return one;  
}

const signature& min(const signature& one, const signature& two)
{
  if (Smaller_sign (one,two))
    {
      return one;
    }
  //else
  if (Smaller_sign (two, one))
    {
      return two;
    }
  //else same signature!
  //cout << "Possible signature drop!" << endl;
  return one;  
}



ostream& operator<<(ostream& s, const signature& rhs)
{
  cout << rhs.Get_term() << "*e_" << rhs.Get_index();
  /*
  bool is_zero = true;
  int j = 0;
  while (is_zero && j < rhs._term.size())
    {
      if (rhs._term[j] != 0)
	{
	  is_zero = false;
	  break;
	}
      ++j;
    }

  if (is_zero)
    {
      cout << "1*e_" << rhs._index;
    }
  else
    {
  for (int i = 0; i < rhs._term.size(); ++i)
    {
      if (rhs._term[i] > 0)
      cout << "x_" << i+1 << "^" <<rhs._term[i];
    }
  cout << "*e_" << rhs._index;
  }*/
  return s;
}


bool Smaller_sign (const signature& S, const signature& T)
{
  return (S<T);
  #ifdef LTPOT
  if (S.Get_term() < T.Get_term())
    {
      return true;
    }
  if (S.Get_term() > T.Get_term())
    {
      return false;
    }
  //else same term
  if (S.Get_index() < T.Get_index())
    return true;
  #endif
  
  //if (S.Get_term().size() != T.Get_term().size())
    //cout << "problem in smaller sign!" << endl;
    //cout << "after smaller sign!" << endl;
  return false; //same signature as false 
    }
