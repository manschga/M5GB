#ifndef LAB_POLY_H_INCLUDED
#define LAB_POLY_H_INCLUDED

#include<array>
#include <iostream>
#include<vector>
#include "other_fct.h"
#include "poly.h"
#include "signature.h"
#include <algorithm>

//extern Table table;

//const vector<pair<vector<int>,vector<int>>>& tableconst vector<pair<vector<int>,vector<int>>>& table
using namespace std;

class Lab_poly{
  public:
  //default constructor
  Lab_poly();
  //real constructor!
  Lab_poly(polynomial poly, signature sign);
  
  int Get_p() const
  {
    return this->_p;
  }
  
  const polynomial& Get_poly() const
  {
    return this->_poly;
  }

  const monomial& Get_LM() const
  {
    #ifdef DEBUG
    if ((this->_poly).Get_monomials().size() == 0)
      {
	cout << "Error! 0-polynomial!" << endl;
      }
    #endif
    return ((this->_poly).Get_monomials())[0];
  }

  int Get_LT() const
  {
    //cout << "LT called!" << endl;
    //not necessary!
    if (this->_poly.size() == 0)
      {
	cout << "term empty!" << endl;
	return 0;
      }
    return (this->_poly).Get_monomials()[0].Get_term();
  }
    
  int Get_LC() const
  {
    return (this->_poly).Get_monomials()[0].Get_coeff();
    //    return (this->_coeffs)[(this->_coeffs).size()-1];
  }

  const signature& Get_sign() const
  {
    return this->_signature;
  }

  int Get_signature_term() const
    {
      return (this->_signature).Get_term();
    }
  
  int Get_index() const
  {
    return (this->_signature).Get_index();
  }

  /*
  bool operator < (const Lab_poly& G) const //change that!
  {
    if (this->Get_index() < G.Get_index())
      {
	return true;
      }
    if (this->Get_index() > G.Get_index())
      {
	return false;
      }
    return (this->Get_sign().Get_term() < G.Get_sign().Get_term());
  }
  */
    

  bool is_Smaller(const Lab_poly& G) const//ratio-order, written with monomials, check that! Unprobable Error source?
  {
    if (this->Get_index() < G.Get_index())
      {
	return true;
      }
    if (this->Get_index() > G.Get_index())
      {
	return false;
      }
    //same index here
    monomial lhs = G.Get_LM().Multiply(this->_signature.Get_term(), false); //= LM(g)S(F)
    monomial rhs = this->Get_LM().Multiply(G._signature.Get_term(), false); //= LM(f)S(G)
    if (lhs < rhs)
      {
	return true;
      }
    if (rhs > lhs)
      {
	return false;
      }
    return (this->_signature.Get_term() < G._signature.Get_term()); //S(F) < S(G))
  }


  /*  
  bool operator < (const Lab_poly& G) const //ratio-order, written with monomials, check that! Unprobable Error source?
  {
    if (this->Get_index() < G.Get_index())
      {
	return true;
      }
    if (this->Get_index() > G.Get_index())
      {
	return false;
      }
    //same index here
    if (G.Get_LM()*this->_signature.Get_term() < this->Get_LM()*G._signature.Get_term()) //LM(g)S(F) < LM(f)S(G)
      {
	return true;
      }
    if (G.Get_LM()*this->_signature.Get_term() > this->Get_LM()*G._signature.Get_term()) //LM(g)S(F) > LM(f)S(G)
      {
	return false;
      }
    return (this->_signature.Get_term() < G._signature.Get_term()); //S(F) < S(G))
  }
  */

  Lab_poly operator+(const Lab_poly& G) const
 {
   //create new polynomial poly(F) + poly(G)
   polynomial h = this->_poly + G._poly;
  
   //create new signature max(S(F),S(G))
   signature new_signature = max(this->_signature,G._signature);

   return Lab_poly(h,new_signature);
      }

  Lab_poly operator+(const polynomial& g) const //same signature, just add poly
  {
    //create new polynomial poly(F) + g with S(F)
    
    return Lab_poly(this->_poly + g,this->_signature);
  }

    Lab_poly operator-(const polynomial& g) const //same signature, just add poly
  {
    //create new polynomial poly(F) + g with S(F)
    
    return Lab_poly(this->_poly - g,this->_signature);
  }
  

  Lab_poly operator-(const Lab_poly& G) const
  {
     
    //create new polynomial poly(F) - poly(G)
    polynomial new_poly = this->_poly - G._poly;
    //create new signature max(S(F),S(G))
    signature new_signature = max(this->_signature,G._signature);
    
    return Lab_poly(new_poly,new_signature);
  }

  Lab_poly Multiply(const int t, bool is_constant) const //multiply with constant / term
  {
    if (is_constant) //mutliply with constant
      {
	return Lab_poly(this->_poly.Multiply(t,is_constant),this->_signature);
      }
    else //multiply with term
      {
	return Lab_poly(this->_poly.Multiply(t,is_constant),this->_signature.Multiply(t));
      }
  }

  Lab_poly Multiply(const vector<int>& t) const//multiply with monomial in vector vector
  {
      /*
      cout << "multiply lab poly with term" << endl;
      cout << "polynomial:" << (this->_poly)*t;
      cout << "signature: " << (this->_signature)*t;
      */
    return Lab_poly((this->_poly).Multiply(t),(this->_signature).Multiply(t));
    }

  /*
  Lab_poly operator*(const monomial m) const//multiply with monomial
    {
      return Lab_poly(this->_poly*m,this->_signature*m.Get_term());
    }
  */

  void update()//remove zeros and reorder terms in decreasing way in polynomial
  {
    _poly =  this->_poly.update();
  }

  /*
  ~Lab_poly() 
    { 
      //cout << "Destructor for lab_poly called " << endl;  
    }
  */ 
    
  friend ostream& operator<<(ostream& s, const Lab_poly& rhs);

  private:
  //  vector<int> _coeffs; //coefficient in macaulay matrix/consider each given term with coeff + array of size n for sparse systems
  polynomial _poly; //polynomial with vector of vectors a_0,t,...,a_n,t with a_0,t = c_t, a_i,t = degx_i(t)
  signature _signature; //todo: first integer: corrsponds to term, translation via macaulay matrix, 2nd integer to index i of e_i
  int _p; //prime for F_p
};

bool Smaller_pair(const pair<Lab_poly,Lab_poly>& one, const pair<Lab_poly,Lab_poly>& two);

ostream& operator<<(ostream& s, const Lab_poly& rhs);

pair<Lab_poly,Lab_poly> spair(const Lab_poly& F, const Lab_poly& G);

Lab_poly spol(const Lab_poly& F, const Lab_poly& G);

vector<polynomial> poly(const vector<Lab_poly>& G);

void sign_sort(vector<pair<Lab_poly,Lab_poly>>& P);

/* test
class Flag_poly{
 public:
  //default constructor
  Flag_poly();
  //real constructor!
  Flag_poly(const polynomial& poly, const signature& flag);
  
   vector<int> Get_mac_poly() const
  {
    return (this->_mac_poly);
  }

    int Get_Flag() const
  {
    return (this->_Flag);
  }
  
 private:
    vector<int> _mac_poly;
    int _Flag; //signature number!
};
*/

//old one
///*
class Flag_poly{
 public:
  //default constructor
  Flag_poly();
  //real constructor!
  Flag_poly(Lab_poly poly, signature flag);
  
   Lab_poly Get_Lab_poly() const
  {
    return (this->_Lab_poly);
  }

    signature Get_Flag() const
  {
    return (this->_Flag);
  }
  
 private:
  Lab_poly _Lab_poly;
  signature _Flag;
};
//*/


#endif // LAB_POLY_H_INCLUDED
