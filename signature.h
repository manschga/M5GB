#ifndef SIGNATURE_H_INCLUDED
#define SIGNATURE_H_INCLUDED

#include<vector>
#include <algorithm>
#include "other_fct.h"


//const vector<pair<vector<int>,vector<int>>>& table

#define LTPOT
using namespace std;

class signature{
 public:
  signature();

  #ifdef DPOT
  int n = 10; //change always -> do better!
  vector<int>  degrees = create_degrees(n,20); //upper found for degree, have to change for larger examples
#endif
  
  signature(int term, int index);
  
  //constructor for "infinity" signature
  signature(int m);
  
  int Get_index() const
  {
    return this->_index;
  }
  term_int Get_term() const
    {
      return this->_term;
    }


#ifdef DPOT
  int degree(term_int term) const
  {
    for (int i = 1; i < degrees.size(); ++i)
      {
	if (degrees[i] > term)
	  {
	    return i-1;
	  }
      }
    //cout << term << endl;
    //cout << "Degree bound overflow!" << endl;
    return 0;
  }
#endif

  bool operator < (const signature &rhs) const
  {
    //cout << "< called with " << this->_term << " " << rhs._term << endl;
#ifdef LTPOT
    if (this->_term == rhs._term)
      {
    	return (this->_index < rhs._index);
      }
    else
      {
	return (this->_term < rhs._term);
      }
    #endif
    #ifdef DPOT
    //cout << "< degrees: " << degree(this->_term) << " " << degree(rhs._term) << endl;
    if (degree(this->_term) == degree(rhs._term))
	  {
	    if (this->_index == rhs._index)
	      {
		return (this->_term < rhs._term);
	      }
	    else
	      {
		return this->_index < rhs._index;
	      }	    
	  }
    else
      {
	return degree(this->_term) < degree(rhs._term);
      }
    #endif
    
    return false;
    //same signature here!
  }

  bool operator == (const signature &rhs) const
  {
    if (this->_index == rhs._index && this->_term == rhs._term)
      {
	return true;
      }
    return false;
  }

  bool operator > (const signature &rhs)
  {
    return ((rhs < *this) && !(*this == rhs));
  }

  signature Multiply(const term_vec& t) const//multiply with term t in vector formulation
{
  term_int result = this->_term; //monomial in integer format
  for (int i = 0; i < t.size(); ++i)
	{
	  for (int j = 0; j < t[i]; ++j) //read up from mult table
	    {
	      result = table[result].second[i]; //multiply with x_i+1
	      // cout << "before term" << new_term[i] << endl;
	  //cout << "t[i]" << t[i] << endl;
	  //cout << "after term" << new_term[i] << endl;
	    }
	}
  return signature(result,this->_index);
}

  signature Multiply(const term_int t) const//multiply with term t in integer formulation
{
  term_vec t_vec = table[t].first; //get vector repr of t
  return this->Multiply(t_vec);
}

  /*
  signature operator*(const vector<int> t) const//multiply with term
    {
      vector<int> current_term = this->_term;
      if (t.size() != current_term.size())
	{
	  cout << "Error in signature multiplication!" << endl;
	}
      for (int i = 0; i < t.size(); ++i)
	{
	  current_term[i] += t[i];
	}
      return signature(current_term,this->_index);
    }
  */

  friend ostream& operator<<(ostream& s, const signature& rhs);
  
 private:
  int _index;
  term_int _term;
};


bool Smaller_sign (const signature& S, const signature& T);

const signature& max(const signature& one, const signature& two);

const signature& min(const signature& one, const signature& two);

ostream& operator<<(ostream& s, const signature& rhs);

#endif
