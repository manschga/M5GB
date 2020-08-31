#include "lab_poly.h"
#include <algorithm>
#include <chrono>
#include <vector>
#include "M5GB.h"
#include "other_fct.h"
#include <set>
#include <unordered_map> 

using namespace std;

//todo: normalize G, case p=2 should work 
//some printing messages on
//define DEBUG

//#define SHOW_PAIRS

//all printing messages on
//#define PRECISE_DEBUG

//time measurement on
//#define TIME

//for counting the values
#define statistics

//define rewrite order, exactly one of those should be defined
//#define NUMBER
#define RATIO

//define order extension
//#define DPOT
#define LTPOT

//test to reduce a better polynomial with same signature
#define LOWER_POL

//#define SMALLEST_LT //-> 
//#define LARGE_GEN -> nice
//#define LARGEST_LT //-> nice
#define NUMBER_OF_TERMS //-> very nice for fast red, not for total red.!
//#define SMALLEST_REMAINDER //-> very nice
//#define LARGEST_FLAG

//install different criteria

int p,m;

#ifdef statistics
int full_red_steps = 0;
int number_of_reductions = 0;
int straight_m_reduces = 0;
int m_re_reduced = 0;
int g_reduced = 0;
int rewrite_counter = 0;
int spol_is_zero = 0;
int syz_counter = 0;
int actual_red = 0;
int zero_reductions = 0;
int BB_counter = 0;
#endif

Fast_conversion_map G_conv;
Fast_conversion_map M_conv;
LT_sig_map LT_G;
LT_sig_map LT_M;
Basis G_;
Multiples M;
Flags Flags_T = {};//for irreducible terms
Flags Flags_M = {};
Generations Generations_G = {};
Generations Generations_M = {};
Generations Generations_T = {};
Syzygies H;
Signatures Spol_sign;
Spair_map Spairs_map;
Spairs P;
vector<int> degrees;
int d = 5;

//////////////////////////////////////////////////////  main functions /////////////////////////////////////////////////


vector<polynomial> M5GB(vector<polynomial>& system, int n) //table might overflow -> need to redefine inside
{
  auto start = std::chrono::steady_clock::now();
  degrees = create_degrees(n,10); //upper found for degree, have to change for larger examples
  Initialize(system,n);  

#ifdef DEBUG
  cout << "Size of P: " << P.size() << endl;
#endif
  int polys_left = m-2;
  int debug = 1000000;
  int maxiter = debug;
  int last_non_zero = maxiter-debug;
  //signature last_syz(table.size(),m);

  while (Spol_sign.size() > 0 && debug > 0)
    {
#ifdef DEBUG
      //cout << "iteration " << maxiter-debug << endl;
#endif
      if (Spol_sign.size() > 0)
	{
	  signature smallest_sign  = *Spol_sign.begin();
	  /*
	  */
	  Spol_sign.erase(Spol_sign.begin());//remove first item
	  if (Syzygy(smallest_sign,H))
	     {
	       syz_counter ++; //to be precise, take away all spairs of that syzygy
	       debug --;
	       continue;
	     }
	  if (smallest_sign.Get_term() != 0)
	    {
	      auto same_sign_pairs = Spairs_map.equal_range(smallest_sign);
	      for (auto pair = same_sign_pairs.first; pair!= same_sign_pairs.second; ++pair)
		{
		  Lab_poly uF = pair->second.first;
		  //Lab_poly vG = pair->second.second.first;
		  /*
		   if(Syzygy(vG.Get_sign(),H))//meets syzygy criterion
		    {
#ifdef DEBUG
		      cout << "syzygy criterion works for " << vG.Get_sign() << endl;
#endif
		      syz_counter ++;
		      debug --;
		      continue;
		    }
		  */
		   if (isRewritable(uF,pair->second.second)) //changed order, check!, could change first taking maximum w.r.t. rewrite order!
		    {
		      rewrite_counter++;
#ifdef DEBUG
		      cout << "Rewritable criterion works for " << uF << endl;
#endif
		      debug --;
		      continue;
		    }
		   /*
		    if (isRewritable(vG,pair->second.second.second))//Generations_G.find(G.Get_LT())->second)) //changed order, check!, could change first taking maximum w.r.t. rewrite order!)) //changed order, check!
		    {
		      rewrite_counter++;
#ifdef DEBUG
		      cout << "Rewritable criterion works for " << vG << endl;
#endif
		      debug --;
		      continue;
		    }
		   */
		    //Lab_poly Spol = uF - vG;
		   term_vec out = table[smallest_sign.Get_term()].first;
		   for (auto &it : out)
		     {
		       //cout << it;
		     }
		   //cout << endl;
#ifdef LOWER_POL
		   Lab_poly Spol = Low_pol(smallest_sign,uF,pair->second.second,n);
#else
		   Lab_poly Spol = uF;//variant
#endif
		   
		  /*
		  if (Spol.Get_LC() == 0)//is 0-polynomial
		    {
		      H[Spol.Get_sign().Get_index()].emplace_back(Spol.Get_sign().Get_term());//add syzygy signature
		      spol_is_zero ++;
		      debug --;
		      break;
		    }
		  */
#ifdef PRECISE_DEBUG
		  cout << "uF:" << uF<< endl;
		  //cout << "vG:" << vG << endl;
		  cout << "Next signature: " << Spol.Get_sign() << endl;
	  if (((maxiter - debug) % 1000) == 1)
	    {
	      cout << "size of G:" << G_.size() << endl;
	      cout << "iteration " << maxiter - debug << endl;
	      cout << "First Buchberger criterion works: " << BB_counter << endl;
	      cout << "spol is zero: " << spol_is_zero << endl;
	      cout << "rewritten: " << rewrite_counter << endl;
	      cout << "syzygy crit: " << syz_counter << endl;
	      cout << "actually reduced: " << actual_red << endl;
	      cout << "signature: " << Spol.Get_sign() << endl;
	    }
#endif
	  if (((maxiter - debug) % 100) == 33)
	    cout << "signature: " << Spol.Get_sign() << endl;
	  //if not rewriteable
	    
	  actual_red ++;
	  //cout << "signature: " << Spol.Get_sign() << endl;
	  //cout << Spol << endl << endl << endl;
	  if (actual_red == 1500)
	    {
	      debug = 0;
	      continue;
	    }
	  Full_S_Reduce(Spol,n);
	  //cout << " --> " << Spol << endl << endl << endl;
#ifdef DEBUG
	  cout << " --> " << Spol << endl << endl << endl;
#endif
	  if((Spol.Get_poly()).Get_LC() == 0)// 0-polynomial
	    {
	      //cout << maxiter-debug << endl;
	      debug = 0;
	      continue;
	      zero_reductions ++;
	      debug --;
	      H[Spol.Get_sign().Get_index()].emplace_back(Spol.Get_sign().Get_term());//add syzygy signature
	      break; //this signature already reduces to 0
	    }
	  else
	    {
	      last_non_zero = maxiter-debug;
	      Update_syz(Spol);//add new koszul syzygies
	      Update_pairs(Spol,system,n,BB_counter);//system needed? check!
#ifdef PRECISE_DEBUG
	      if (G_.find(Spol.Get_LT()) != G_.end())
		cout << "Same leading term! Error!" << endl;
#endif
	      Add_to_M(Spol);
	      Add_to_G(Spol);
	      /*
	      G_[Spol.Get_LT()] = Spol;
	      M[Spol.Get_LT()] = Spol;
	      Generations_M[Spol.Get_LT()] = G_.size(); //generation for M
	      Generations_G[Spol.Get_LT()] = G_.size(); //generation for G
	      Flags_M[Spol.Get_LT()] = Create_Flag(Spol.Get_poly()); //create flag for M
	      */
	      //Update_step(Spol); -> not needed since done in lazy manner!
	    }
	  debug --;
	  break; //one element of this signature already examined
		}
	    }
	}
    }
  vector<polynomial> basis = Poly(G_,n,p); //remove signatures

  sort(basis);
  Interreduce(basis);
  
#ifdef statistics
  cout << "Number of examined S-polynomials: " << maxiter - debug << endl;
  cout << "discarded since zero:" << spol_is_zero << endl;
  cout << "Discarded by Rewriteable criterion: " << rewrite_counter << endl;
  cout << "discarded by syzygy criterion: " << syz_counter << endl;
  cout << "actually reduced: " << actual_red << endl;
  cout << "Number of full reduction steps:" << full_red_steps << endl;
  cout << "Number of total reductions:" << number_of_reductions << endl;
  cout<< "Number of fast reductions: " <<  straight_m_reduces << endl;
  cout<< "Number of re-reductions: " <<  m_re_reduced << endl;
  cout<< "Number of g-reductions: " <<  g_reduced << endl;
  cout<< "Size of M: " << M.size() << endl;
  cout<< "Size of G: " << G_.size() << endl;
  cout << "Number of Zero reductions:" << zero_reductions << endl;
  cout << "First Buchberger criterion works: " << BB_counter << endl;
  cout << "Last non-zero: " << last_non_zero << endl;
#endif
  if (debug == 0)
    cout << "maximal number of iterations reached!!!"<< endl << endl;
  auto end = chrono::steady_clock::now();
  cout << "Elapsed time in Seconds : "
	 << chrono::duration_cast<chrono::milliseconds>(end - start).count()
	 << " milliseconds " << endl;
  return basis;
  //return G_;
}

void Add_to_M(const Lab_poly& F)
{
  term_int LT = F.Get_LT();
  M_conv[LT] = make_pair(table[LT].first,F.Get_sign());
  M.erase(LT);
  M[LT] = F.Get_poly();
  Generations_M.erase(LT);
  Generations_M[LT] = G_.size();
  Flags_M.erase(LT);
  Flags_M[LT] = Create_Flag(F.Get_poly());
  LT_M.erase(LT);
  LT_M[LT] = F.Get_sign();
  return;
}

//#define SMALLEST_LT// -> 
//#define LARGE_GEN// -> nice
//#define LARGEST_LT //-> nice
//#define NUMBER_OF_TERMS// -> very nice
//#define SMALLEST_REMAINDER// -> very nice
//#define LARGEST_FLAG

//construct a labelled polynomial which is better reduced than uF with same signature
Lab_poly Low_pol(const signature& smallest_sign, const Lab_poly& uF, const int generation, int n)
{
  term_int LT_low = 0;
  term_int remainder_low = 0;
  Lab_poly Low = uF;
  bool not_yet = true;
#ifdef LARGEST_LT //no sense to do largest lt heuristic here
  return uF;
#endif
  for (auto& it : LT_M)
    {
      if (smallest_sign.Get_index() == it.second.Get_index())
	{
	  if (isDivisible(table[smallest_sign.Get_term()].first,table[it.second.Get_term()].first))
	    {
	      polynomial test = M.find(it.first)->second.Multiply(get_remainder(table[smallest_sign.Get_term()].first,table[it.second.Get_term()].first));
	      
#ifdef LARGE_GEN
	      if (not_yet) //still uF is best
		{

	       if (Generations_G.find(it.first)->second > generation) //different heuristics to check
		    {
		      not_yet = false;
		      Low = Lab_poly(test,smallest_sign);
		      LT_low = it.first;
		    }
		}
	      else
		{
		  if (Generations_G.find(it.first)->second > Generations_G.find(LT_low)->second) //different heuristics to check
		     {
		       Low = Lab_poly(test,smallest_sign);
		       LT_low = it.first;
		     }
		}
#endif
	      //*
#ifdef NUMBER_OF_TERMS
	      if (not_yet) //still uF is best
		{
		  if (test.size() <= uF.Get_poly().size()) //different heuristics to check
		    {
		      not_yet = false;
		      Low = Lab_poly(test,smallest_sign);
		      LT_low = it.first;
		    }
		}
	      else
		{
		  if (test.size() < Low.Get_poly().size()) //different heuristics to check
		     {
		       Low = Lab_poly(test,smallest_sign);
		       LT_low = it.first;
		     }
		}
#endif
	      //*/
#ifdef SMALLEST_REMAINDER
	      if (not_yet) //still uF is best
		{
		  if (test.size() < uF.Get_poly().size()) //different heuristics to check
		    {
		      not_yet = false;
		      Low = Lab_poly(test,smallest_sign);
		      LT_low = it.first;
		      remainder_low = term_to_int(get_remainder(table[smallest_sign.Get_term()].first,table[it.second.Get_term()].first),n);
		    }
		}
	      else
		{
		  if(term_to_int(get_remainder(table[smallest_sign.Get_term()].first,table[it.second.Get_term()].first),n) > remainder_low)
		    {
		       Low = Lab_poly(test,smallest_sign);
		       LT_low = it.first;
		       remainder_low = term_to_int(get_remainder(table[smallest_sign.Get_term()].first,table[it.second.Get_term()].first),n);
		     }
		}
#endif
#ifdef SMALLEST_LT
	      if (not_yet) //still uF is best
		{
		  if (test.Get_LT() <= uF.Get_LT()) //different heuristics to check
		    {
		      not_yet = false;
		      Low = Lab_poly(test,smallest_sign);
		      LT_low = it.first;
		    }
		}
	      else
		{
		  if (test.Get_LT() <= Low.Get_LT()) //different heuristics to check
		     {
		       Low = Lab_poly(test,smallest_sign);
		       LT_low = it.first;
		     }
		}
#endif
	    }
	}
    }
	  return Low;
	  /*

  
  for (auto& it : LT_M) 
    {
      if (Lab_pol_M.find(it.first)->second == true) //check if correct label!
	{
      if (smallest_sign.Get_index() == it.second.Get_index())
	{
	  if (isDivisible(table[smallest_sign.Get_term()].first,table[it.second.Get_term()].first))
	    {
	      polynomial test = M.find(it.first)->second.Multiply(get_remainder(table[smallest_sign.Get_term()].first,table[it.second.Get_term()].first));
#ifdef high_gen
	      if (not_yet) //still uF is best
		{

	       if (Generations_M.find(it.first)->second > generation) //different heuristics to check
		    {
		      not_yet = false;
		      Low = Lab_poly(test,smallest_sign);
		      LT_low = it.first;
		    }
		}
	      else
		{
		  if (Generations_M.find(it.first)->second > Generations_M.find(LT_low)->second) //different heuristics to check
		     {
		       Low = Lab_poly(test,smallest_sign);
		       LT_low = it.first;
		     }
		}
#endif

#ifdef smallest_size
	      if (not_yet) //still uF is best
		{
		  if (test.size() <= uF.Get_poly().size()) //different heuristics to check
		    {
		      not_yet = false;
		      Low = Lab_poly(test,smallest_sign);
		      LT_low = it.first;
		    }
		}
	      else
		{
		  if (test.size() < Low.Get_poly().size()) //different heuristics to check
		     {
		       Low = Lab_poly(test,smallest_sign);
		       LT_low = it.first;
		     }
		}
#endif
#ifdef smallest_remainder
	      if (not_yet) //still uF is best
		{
		  if (test.size() <= uF.Get_poly().size()) //different heuristics to check
		    {
		      not_yet = false;
		      Low = Lab_poly(test,smallest_sign);
		      LT_low = it.first;
		      remainder_low = term_to_int(get_remainder(table[smallest_sign.Get_term()].first,table[it.second.Get_term()].first),n);
		    }
		}
	      else
		{
		  if(term_to_int(get_remainder(table[smallest_sign.Get_term()].first,table[it.second.Get_term()].first),n) > remainder_low)
		    {
		       Low = Lab_poly(test,smallest_sign);
		       LT_low = it.first;
		       remainder_low = term_to_int(get_remainder(table[smallest_sign.Get_term()].first,table[it.second.Get_term()].first),n);
		     }
		}
#endif
	    }
	    }
	}
    }
  return Low;
	  */
}

void Add_to_G(const Lab_poly& F)
{
  term_int LT = F.Get_LT();
  G_conv[LT] = make_pair(table[LT].first,F.Get_sign());
  G_[LT] = F.Get_poly();
  Generations_G[LT] = G_.size();
  LT_G[LT] = F.Get_sign();
  return;
}
  
void Full_S_Reduce(Lab_poly &F, const int n)
{
  //number_of_reductions++;
#ifdef PRECISE_DEBUG
  cout << "doing full-reduce with " << F << endl;
#endif
  signature T = F.Get_sign();
  polynomial f = F.Get_poly();
  polynomial f_ = f;
  for (auto& monom : f.Get_monomials())
    {
#ifdef PRECISE_DEBUG
      cout << "t = " << monom.Get_term() << endl;
#endif
      int t = monom.Get_term();
       int c = monom.Get_coeff();
       if (Flags_T.find(t) != Flags_T.end() &&  (Flags_T.find(t)->second > T || Flags_T.find(t)->second == T))
	 {
         if (Generations_T.find(t)->second == G_.size())//t is known to be irreducible -> nothing needs to be done!
	     {
	       continue;
	     }
	   else if (!Check_TGFT(t,T))//if not reducible by G, possible improvement: check only for newer elements, nothing for counting argument!
	     {
	       Generations_T.erase(t);
	       Generations_T[t] = G_.size();
	       Flags_T.erase(t);//update signature flag since new generation, improvement: take only newer ons and do min with old flag
	       Flags_T[t] = Create_term_flag(t);
	     }
	 }
       //bool compatible = false; -> improvement to save some div by G operations
      auto itr = LT_M.find(t);
      if (itr != LT_M.end()) //element in M found
	{
	  Lab_poly m_ = Update_m(Lab_poly(M.find(itr->first)->second,itr->second),T,n); //bring m in right form such that it is tail-irreducible again
	  f_ = f_ - m_.Get_poly().Multiply(c,true); //f - c_t(f)*m_t
	  number_of_reductions++;
	  full_red_steps ++;
	}
      else if (Check_TGFT(t,T))//if reducible by G
	{
	  Flags_T.erase(t);//t is no longer irreducible (only first one needed, rest for less storage) 
	  Generations_T.erase(t);
	  Lab_poly vM = Div_by_M(t,T);//find best divisor in M, improvment: map G -> M, only need to check that m's on divisibility, if not by M -> reduce by G, inside implemented!
	  Lab_poly m_ = Tail_S_Reduce(vM,T,n); //create new element in M such that it is tail-irreducible
	  f_ = f_ - m_.Get_poly().Multiply(c,true); //f - c_t(f)*m_t
	  number_of_reductions++;
	  full_red_steps ++;
	}
      else
	{
	  Flags_T.erase(t);//update signature flag since new generation, improvement: take only newer ons and do min with old flag
	  Flags_T[t] = Create_term_flag(t);
	  Generations_T[t] = G_.size();
      }
    }
  F = Lab_poly(f_,T);
#ifdef DEBUG
  //cout << "Full_reduce element: " << F << endl;
#endif
  return;
}

signature Create_term_flag(const term_int t)
{
  vector<int> t_ = table[t].first;
  signature flag(table.size(),m); //cheap workaround, infinite signature
  for (auto& it : G_conv)
    {
      if (isDivisible(t_,it.second.first))//compute function to get remainder straight
	{
	  signature candidate = it.second.second.Multiply(get_remainder(t_, it.second.first));//if divisor divides term, then term = divisor*remainder;
	  flag = min(flag,candidate);
	}
    }
  return flag;
}

signature Create_Flag(const polynomial& f)
{
  signature flag(table.size(),m); //cheap workaround, infinite signature
  for (auto monom = ++f.Get_monomials().begin(); monom != f.Get_monomials().end(); ++monom) //leave out LT, need to stay at that datastructure, otherwise for-loop! (poss error source)
    {
      term_int t = monom->Get_term();
      if (Generations_T.find(t)->second < G_.size()) //only do something then
	flag = min(flag,Create_term_flag(t));
      else
	flag = min(flag,Flags_T.find(t)->second);
    }
  return flag;
}

Lab_poly Tail_S_Reduce(const Lab_poly& F_,const signature& T, const int n)
{
#ifdef PRECISE_DEBUG
  cout << "doing Tail-reduce with polynomial " << F_ << endl;
#endif
  Lab_poly F = F_;
  signature S_F = F.Get_sign();
  //signature flag(table.size(),m);
  int rec_depth = 1; //cheap workaround until now
  if (rec_depth < 10000)
    {
      rec_depth++;
      polynomial f = F_.Get_poly();
      polynomial f_ = f;
      for (auto monom = ++f.Get_monomials().begin(); monom != f.Get_monomials().end(); ++monom) //leave out LT, need to stay at that datastructure, otherwise for-loop! (poss error source)
     {
       int t = monom->Get_term();
       int c = monom->Get_coeff();
#ifdef PRECISE_DEBUG
       cout << "t = " << monom->Get_term() << endl;
#endif
       if (Flags_T.find(t) != Flags_T.end() &&  (Flags_T.find(t)->second > T || Flags_T.find(t)->second == T )) 
	 {
	   if (Generations_T.find(t)->second == G_.size())//t is known to be irreducible -> nothing needs to be done!
	     {
	       continue;
	     }
	   else if (!Check_TGFT(t,T))//if not reducible by G, possible improvement: check only for newer elements, nothing for counting argument!
	     {
	       Generations_T.erase(t);
	       Generations_T[t] = G_.size();
	       Flags_T.erase(t);//update signature flag since new generation, improvement: take only newer ons and do min with old flag
	       Flags_T[t] = Create_term_flag(t);
	       continue;
	     }
	 }
       auto itr = LT_M.find(t);
      if (itr != LT_M.end()) //element found, check on order if checking by div by G is faster?
	{
	  Lab_poly m_ = Update_m(Lab_poly(M.find(itr->first)->second,itr->second),T,n); //bring m in right form such that it is tail-irreducible again
	  f_ = f_ - m_.Get_poly().Multiply(c,true); //f - c_t(f)*m_t
	  S_F = max(S_F,m_.Get_sign());
	  number_of_reductions++;
	}
      else if (Check_TGFT(t,T))//if reducible by G
	{
	  Flags_T.erase(t);//t is no longer irreducible (only first one needed, rest for less storage) 
	  Generations_T.erase(t);
	  
	  Lab_poly vM = Div_by_M(t,T);//find best divisor in M, improvment: map G -> M, only need to check that m's on divisibility
	  Lab_poly m_ = Tail_S_Reduce(vM,T,n); //create new element in M such that it is tail-irreducible
	  f_ = f_ - m_.Get_poly().Multiply(c,true); //f - c_t(f)*m_t
	  S_F = max(S_F,m_.Get_sign());
	  number_of_reductions++;	
	}
       else
	 {
	   //not yet reducible -> keep t as irreducible term!
	   Flags_T.erase(t);//update signature flag since new generation, improvement: take only newer ons and do min with old flag
	   Flags_T[t] = Create_term_flag(t);
	   Generations_T.erase(t);
	   Generations_T[t] = G_.size();//set current generation when t was updated the last time
	 }
     }
      F = Lab_poly(f_,S_F);
      Add_to_M(F);
      return F;
    }
  
  cout << "Too small recursion depth" << endl;
   
  return Lab_poly();
}

//////////////////////////////////////////////////////////////////////////// Criteria ///////////////////////////////////////////////////////////////////////////////

//check if LT(F), LT(G) are disjoint, if so return true -> discard spair F,G
bool First_Buchberger(const vector<int>& LT_F, const vector<int>& LT_G)
{
  int n = LT_F.size();
  for (int i = 0; i < n; ++i)
    {
      if (LT_F[i] > 0 && LT_G[i] > 0)
	{
	  return false;
	}
    }
  return true;
}

//checks if signature S is a known syzygy signature by checking if
//there exists some h in H such that h divides S
bool Syzygy(const signature& sign, const Syzygies& H)
{
  int i = sign.Get_index();
  //  int sign_term = S.Get_term();
  //checks if S is known syzygy, order H by term, then faster!
  for (auto it = H[i].begin(); it != H[i].end(); ++it)
    {
      if (*it == sign.Get_term()) //fast check, if H > S then H does not divide S
	{
	  return true;
	}
    }
  for (auto it = H[i].begin(); it != H[i].end(); ++it)
    {
      if (*it < sign.Get_term()) //fast check, if H > S then H does not divide S
	{
	  if (isDivisible(table[sign.Get_term()].first,table[*it].first)) //look it up or faster by memory?
	    {
	      //cout << "signature H_i: " << H[i] << endl;
	      //cout << "signature from Spol: " << S << endl;
	      return true;
	    }
	}
    }
  return false;
}

bool isRewritable(const Lab_poly& F, const int generation)
{
  signature T = F.Get_sign();
  for (auto it : LT_G)
    {
      if (it.second.Get_index() == T.Get_index())
	{
	  if (isDivisible(table[T.Get_term()].first,table[it.second.Get_term()].first))
	    {
#ifdef NUMBER
	      //rewrite order num: F < G if F created before G, still some error in here!
	      if (generation < Generations_G.find(it.first)->second)
		{
		  /*
		  cout << "Gen(F) = " << generation << endl;
		  cout << "Gen(G) = " << Generations_G.find(it.first)->second << endl;
		  */
		  return true;
		}
	      
#endif
#ifdef RATIO
	      //define vG = v*G_i if S(vG_i) = S(F);
	      vector<int> v = get_remainder(table[T.Get_term()].first,table[it.second.Get_term()].first);
	      //rewrite order, for ratio, only need: LT(g) < LT(f)
	      if (F.Get_LT() > G_.find(it.first)->second.Multiply(v).Get_LT())
		{
#ifdef DEBUG
		  cout << F << "is rewritable by " << it.second << endl;
		  cout << F << "is rewritable by " << it.second.Multiply(v) << endl;
#endif
		  return true;
		}
#endif
	    }
	}
    }
  return false;
}

/////////////////////////////////////////////////////////////// Update functions ////////////////////////////////////////////////////////////////////////////////

Lab_poly Update_m(const Lab_poly& m, const signature& T, const int n)//needs huge improvment -> ineffective!
{
  term_int LT = m.Get_LT();
  int generation = Generations_M[LT];
  signature Flag = Flags_M[LT];
  //4 cases; check if current generation and flag >= T
  if (generation == G_.size() && (Flag > T || Flag == T))
    {
      straight_m_reduces++;
      return m;
    }
  //m_re_reduced++;
  
  if (generation < G_.size())
    {
      /*
      Lab_poly F = m;
      signature S_F = m.Get_sign();
      polynomial f = m.Get_poly();
      polynomial f_ = f;
      for (auto monom = ++f.Get_monomials().begin(); monom != f.Get_monomials().end(); ++monom) //leave out LT, need to stay at that datastructure, otherwise for-loop! (poss error source)
	{
	  int t = monom->Get_term();
	  int c = monom->Get_coeff();
	  if (Flags_T.find(t) != Flags_T.end() &&  (Flags_T.find(t)->second > T || Flags_T.find(t)->second == T )) 
	    {
	      if (Generations_T.find(t)->second == G_.size())//t is known to be irreducible -> nothing needs to be done!
		{
		  continue;
		}
	      else if (!Check_TGFT(t,T,generation))//if not reducible by G, possible improvement: check only for newer elements, nothing for counting argument!
		{
		  Generations_T.erase(t);
		  Generations_T[t] = G_.size();
		  Flags_T.erase(t);//update signature flag since new generation, improvement: take only newer ons and do min with old flag
		  Flags_T[t] = Create_term_flag(t);
		  continue;
		}
	    }
	  auto itr = LT_M.find(t);
	  if (itr != LT_M.end()) //element found, check on order if checking by div by G is faster?
	    {
	      Lab_poly m_ = Update_m(Lab_poly(M.find(itr->first)->second,itr->second),T,n); //bring m in right form such that it is tail-irreducible again
	      f_ = f_ - m_.Get_poly().Multiply(c,true); //f - c_t(f)*m_t
	      S_F = max(S_F,m_.Get_sign());
	      number_of_reductions++;
	    }
	  else if (Check_TGFT(t,T,generation))//if reducible by G
	{
	  Flags_T.erase(t);//t is no longer irreducible (only first one needed, rest for less storage) 
	  Generations_T.erase(t);
	  
	  Lab_poly vM = Div_by_M(t,T);//find best divisor in M, improvment: map G -> M, only need to check that m's on divisibility
	  Lab_poly m_ = Tail_S_Reduce(vM,T,n); //create new element in M such that it is tail-irreducible
	  f_ = f_ - m_.Get_poly().Multiply(c,true); //f - c_t(f)*m_t
	  S_F = max(S_F,m_.Get_sign());
	  number_of_reductions++;	
	}
	  else
	    {
	      //not yet reducible -> keep t as irreducible term!
	      Flags_T.erase(t);//update signature flag since new generation, improvement: take only newer ons and do min with old flag
	      Flags_T[t] = Create_term_flag(t);
	      Generations_T.erase(t);
	      Generations_T[t] = G_.size();//set current generation when t was updated the last time
	    }
	}
      //add_to?
      F = Lab_poly(f_,S_F);
      Add_to_M(F);
      return F;
      */
      return Tail_S_Reduce(m,T,n); //need only to check for newer generations -> improvement, but not for counting argument!
    }
  if (generation == G_.size() && Flag < T)
    {
        m_re_reduced++;
	return Tail_S_Reduce(m,T,n); //need only to check for newer generations -> improvement, but not for counting argument!); //need only to check for those terms -> improvement, but not for counting argument!
    }
  cout << "Error occured! element newer than current generation!" << endl;
  return Lab_poly();
}

//updates the flags of M_ by checking whether t = vLT(f)  for some term t in T(g) where g in M,
//if S(vF) < S-flag(g), set S-flag(g) = S(vF)
//after this function, all g in M_ are S-tail-reduced w.r.t. G and their own S-flag, maybe not checking & always calling tail-reduce on M is faster?
/*
void Update_step(const Lab_poly& F)
{
  for (auto it : M) 
    {
      polynomial g_i = it.second.Get_poly();
      for (int terms = 0; terms < g_i.size(); ++terms) //for all t \in T(f)
	{
	  vector<int> t = table[g_i.Get_monomials()[terms].Get_term()].first;
	  if (isDivisible(t,table[F.Get_LT()].first))
	    {
	      vector<int> v = get_remainder(t, table[F.Get_LT()].first);//if divisor divides term, then term = divisor*remainder
	      Lab_poly vF = F.Multiply(v);
	      if (Smaller_sign(vF.Get_sign(),it.second.Get_sign())) //check if S(vF) < S-Flag(M_i)
		{
		  M[it.first] = Lab_poly(it.second.Get_poly(),vF.Get_sign());//do a setter to increase speed!
		}
	    }
	}
    }
  for (auto it = Irreducible.begin(); it != Irreducible.end(); ++it)
    {
      if (isDivisible(table[*it].first,table[F.Get_LT()].first))
	{
	  Irreducible.erase(*it);
	}
    }
    return;
}
*/

//adds {(uF,vg): g in G} to P, where LM(uF) = LM(vg), change pairs to (vg,uF) if S(vg) > S(uf), remove pair if same signature!
//improvement later:
void Update_pairs(const Lab_poly& F, const vector<polynomial>& system, const int n, int& BB_counter)
{
  //cout << "pair update: " << LT_G.size() << endl;
  for (auto& it : LT_G) 
    {
      //cout << Generations_G.find(it.first)->second << endl;
      if (Generations_G.find(it.first)->second <= n) //do not add spols coming from field equations!
	{
	  //cout << "Field equation detected!" << endl;
	  continue;
	}      
      
      if (First_Buchberger(table[F.Get_LT()].first,table[it.first].first))
	{
	  //inefficient!
	  /*
	  pair<Lab_poly,Lab_poly> Spair = spair(F,Lab_poly(G_.find(it.first)->second,it.second));
	  Lab_poly Spol = Spair.first-Spair.second;
	  H[Spol.Get_sign().Get_index()].emplace_back(Spol.Get_sign().Get_term());
	  */
	  BB_counter ++;
	  continue;
	}
     
      //*/
      pair<Lab_poly,Lab_poly> Spair = spair(F,Lab_poly(G_.find(it.first)->second,it.second));
      pair<pair<Lab_poly,int>,pair<Lab_poly,int>> spair_gen = make_pair(make_pair(Spair.first,G_.size()+1),make_pair(Spair.second,Generations_G.find(it.first)->second));
      if (Spair.first.Get_sign() == Spair.second.Get_sign()) // remove pair if same signature! Get counter for statistic
	{
	  continue;
	  //cout << "singular pair!" << endl;
	}
      else
	{
	  if (Spair.first.Get_sign() < Spair.second.Get_sign()) 
	    {
	      //change order to ensure first > second
	      Spair = make_pair(Spair.second,Spair.first);
	      spair_gen = make_pair(spair_gen.second,spair_gen.first);
	    }
	  //insert correctly

	  Spol_sign.insert(Spair.first.Get_sign()); //not efficient, only inserted when signature not yet in
	  Spairs_map.insert(make_pair(Spair.first.Get_sign(),spair_gen.first)); //not efficient, only inserted when signature not yet in
	  //auto correct_insert = P.insert(Spair);
	  /*
	  if (!correct_insert.second)
	    {
	      std::cout << "Failed to insert element!" << std::endl;
	      cout << Spair.first.Get_sign() << endl;
	      for (auto it = P.begin(); it != P.end(); ++it)
		{
		  cout << it->first.Get_sign() << endl;
		}
	      
	    }
	  */
	}
    }

      
#ifdef DEBUG
    for (auto it = P.begin(); it != P.end(); ++it)
    {
      cout << it->first.Get_sign() << ","  << it->second.Get_sign() << endl;
    }
#endif

  //cout << "size of P after update: " << P.size() << endl;
  return;
}

/*
//adds {spol(F,g): g in G} to P
//do an improvement later by using criteria to discard some S-polynomials already here (non-critical pairs,...)
void Update_spol(const Lab_poly& F, vector<Lab_poly>& P, const vector<Lab_poly>& G)
{
  for (int g = 0; g < G.size(); ++g)
    {
      Lab_poly Spol = spol(F,G[g]);
      Spol.update();
      P.emplace_back(Spol);
    }
  return;
}
*/


//adds the signatures {S(LT(g)*f-LT(f)*g, g in G} to H
void Update_syz(const Lab_poly& Spol)
{
  for (auto& it : LT_G) 
    {
      //cout << "update syz" << Spol.Get_sign().Multiply(it.first) << it.second.Multiply(Spol.Get_LT())<< endl;
      signature temp = max(Spol.Get_sign().Multiply(it.first),it.second.Multiply(Spol.Get_LT()));
      H[temp.Get_index()].emplace_back(temp.Get_term());
    }
  return;
}



/////////////////////////////////////////////////////////////////////////////////////////// Init functions //////////////////////////////////////////////////////////////

void Initialize(vector<polynomial>& system, int n) //table might overflow -> need to redefine inside
{
  //sort(system); //test if this speeds up
  p = system[0].Get_p();
  Add_field_eq(system,n);
  //m = system.size();
  Interreduce(system);
  sort(system.begin()+n, system.end());
  m = system.size();
  cout << "n: " << n << endl;
  cout << "p: " << p << endl;
  cout << "m: " << m << endl;

  vector<Lab_poly> F = Basic_poly(system,n);
  Init_LT_G(F);
  G_ = Init_G(F);
  M = G_;
  LT_M = LT_G;
  P = Init_Spairs(F,Spol_sign,Spairs_map,n);//check that, error source? comparison?
  H = Init_Syzygies(F); //create trivial syzygy signatures
  Init_generations(F);
  
  #ifdef NUMBER
  //Numbers = Init_numbers(F);
  #endif
  
  #ifdef PRECISE_DEBUG
  cout << "Trivial syzygies: " << endl;
  for (auto it = H.begin(); it != H.end(); ++it)
    {
      for (auto it2 = it->begin(); it2 != it->end(); ++it2)
	cout << *it2 << endl;
    }
#endif
    cout << "Initialization done!" << endl;
    return;
}

void Interreduce(vector<polynomial>& system)
{
  int size = 0;
  int size_before = 1;
  while(size != size_before)
    {
      cout << "reduction step: " << system.size() << endl;
      size_before = system.size();
      int reductions = reduce(system); //get the (unique) reduced gröbner basis
      cout << "after reduce: " << system.size() << endl;
      size = system.size();//maybe the size drops
      number_of_reductions += reductions;
    }
  cout << "end of reduction!" << endl;
  //sort(system); //test if this speeds up
  reverse(system.begin(),system.end());
}


//adds the field equations x_i^p - x_i to the polynomial system for all i = 1,...,n
void Add_field_eq(vector<polynomial>& system,int n)
{
  //assume field_eq are not in system yet!
  for (int i = 0; i < n; ++i)
    {
      //add x_i^p - x_i
      vector<int> x_p(n,0);
      x_p[i] = p;
      monomial one(0,1,p); //one-polynomial
      monomial mon_xp = one.Multiply(x_p);//multiply with x^p
      monomial minus_one(0,p-1,p); //-1-polynomial
      vector<int> x_1(n,0);//or n+1?
      x_1[i] = 1;
      monomial mon_x1 = minus_one.Multiply(x_1);

      //cout << "field_eq_i: "<< polynomial({mon_xp,mon_x1}) << endl;
      system.emplace_back(polynomial({mon_xp,mon_x1}));
    }
  //reverse(system.begin(),system.end());
    return;
}

//turns the polynomial system into a labelled polynomial system by f_i -> F_i
//where F_i = (f_i,e_i) for i = 1,...,m


vector<Lab_poly> Basic_poly(const vector<polynomial>& system_, int n)
{
  vector<Lab_poly> F = {};
  F.reserve(m);//(m);
  for (int i = 0; i < m; i++)
    {
      #ifdef LTPOT
      Lab_poly F_i(system_[i],signature(system_[i].Get_LT(),i));
      #endif
      #ifdef DPOT
      Lab_poly F_i(system_[i],signature(0,i));
      #endif
      //Lab_poly F_i(system_[i],signature(0,i)); -> is worse by a test!
      F.emplace_back(F_i);
    }
  return F;
}

Basis Create_basis(const vector<Lab_poly>& F)
{
  Basis basis;
  for (auto it = F.begin(); it != F.end(); ++it)
    {
      basis[it->Get_LT()] = it->Get_poly();
    }
  return basis;
}

//returns the set {F_2,...F_m} where F = {F_1,...,F_m}
Spairs F_2(const vector<Lab_poly>& F, int n, int p)
{
  Spairs P = {};

  vector<int> zero(n,0);
  Lab_poly Zero(polynomial({monomial(0,0,p)}),signature(0,0));
  //cout << "Zero: " << Zero << endl;
  for (int i = 1; i < F.size(); ++i)
    {
      P.insert(make_pair(F[i],Zero));
      //    P.emplace_back(
    }
  //reverse(P.begin(),P.end()); 
  return P;
}

//return poly(G)
vector<polynomial> Poly(const Basis& G, int n, int p)
{
  vector<polynomial> G_poly = {};
  G_poly.reserve(G_.size());// = {};
  for (auto it : G_) 
    {
      G_poly.emplace_back(it.second);
    }
  return G_poly;
}

//initializes hash map with LT(F) -> F
Basis Init_G(const vector<Lab_poly>& F)
{
  Basis G_ = {};
  for (auto it = F.begin(); it != F.end(); ++it)
    {
      G_[it->Get_LT()] = it->Get_poly(); //check: what if same LT in G_?
    }
  return G_;
}

//initializes vector with all elading terms of G
void Init_LT_G(const vector<Lab_poly> & F)
{
  for (auto it = F.begin(); it != F.end(); ++it)
    {
      LT_G[it->Get_LT()] = it->Get_sign();
      G_conv[it->Get_LT()] = make_pair(table[it->Get_LT()].first,it->Get_sign());
      M_conv[it->Get_LT()] = make_pair(table[it->Get_LT()].first,it->Get_sign());
    }
  //return LT_G;
}

//init spol(f_i,f_j) for all i,j
Spairs Init_Spairs(const vector<Lab_poly>& F, Signatures& Spol_sign, Spair_map& Spairs_map, int n)
{
  Spairs P = {};
  for (int i = n-1; i < F.size(); ++i)
    {
      for (int j = i+1; j < F.size(); ++j)
	{
	  pair<Lab_poly,Lab_poly> Spair = spair(F[j],F[i]);
	  Spol_sign.insert(spair(F[j],F[i]).first.Get_sign()); //not efficient, only inserted when signature not yet in
	  pair<pair<Lab_poly,int>,pair<Lab_poly,int>> spair_gen = make_pair(make_pair(Spair.first,j+1),make_pair(Spair.second,i+1));//check that, error source?
	  Spairs_map.insert(make_pair(Spair.first.Get_sign(),spair_gen.first));//
	  //Spairs_map.insert(make_pair(spair(F[j],F[i]).first.Get_sign(),spair(F[j],F[i]))); //not efficient, only inserted when signature not yet in
	  //P.insert(spair(F[j],F[i]));
	}
    }

#ifdef PRECISE_DEBUG
  for (auto it = P.begin(); it != P.end(); ++it)
    {
      cout << it->first.Get_sign() << "," <<  it->second.Get_sign() << endl;
    }
#endif
  return P;
}

//add koszul syzygy signatures: LT(f_i)S(f_j) for j < i if ascending sorted! For descending, error source!
Syzygies Init_Syzygies(const vector<Lab_poly> &F)
{
  int k = F.size();
  Syzygies H(k);
  for (int i = 0; i < k; ++i)
    {
      H[i].reserve(i);
    }
  
  for (int i = 0; i < F.size(); ++i)
    {
      for (int j = i+1; j < F.size(); ++j)
	{
	  H[j].emplace_back(F[j].Multiply(F[i].Get_LT(),false).Get_sign().Get_term()); //should i do without LT? check after!
	}
    }
  return H;
}

/*
Number_lookup Init_numbers(const vector<Lab_poly> &F)
{
  Number_lookup Numbers = {};
  for (auto it = F.begin(); it != F.end(); ++it)
    {
      Numbers[it->Get_LT()] = Numbers.size(); //check: what if same LT in G_?
    }
  return Numbers;
}
*/

void Init_generations(const vector<Lab_poly> &F)
{
  for (auto it = F.begin(); it != F.end(); ++it)
    {
      Generations_G[it->Get_LT()] = Generations_G.size()+1;
      Generations_M[it->Get_LT()] = Generations_M.size()+1;
    }
}




///////////////////////////////////////////////////////////////////////////////// Finding functions /////////////////////////////////////////////////////////////////////////////////


 #if 0
//if there exists some g in G such that t = LT(vg), S(vg) < T, return vg
//else return G[0], should never come to this since by Check_TGFT not possible!
auto Find_in_G(const vector<int>& t, const Basis& G_, const signature& T)
{
  for (auto& it : LT_G) 
    {
      if (isDivisible(t,table[it->first].first))//compute function to get remainder straight
	{
	  //vector<int> v = get_remainder(t, table[it->first].first));//if divisor divides term, then term = divisor*remainder;
	  Lab_poly vG = it->second.Multiply(get_remainder(t, table[it->first].first));//if divisor divides term, then term = divisor*remainder;
	  if (Smaller_sign(vG.Get_sign(),T)) //check if S(vG) < T
	    {
	      return vG;
	    }
	}
    }
 
int g = 0;
  //  vector<int> temp = G[g].Get_LT();
 int size_G = G.size();
 while (g < G.size())
    {
      /*
      cout << "t: " << endl;
      for (int i = 0; i < t.size(); ++i)
	{
	  cout << t[i] << ",";
	}
      */
	 
      if (isDivisible(t,table[G[g].Get_LT()].first))
	{
	  //cout << "G in found G: " << G[g] << endl;
	  vector<int> v = get_remainder(t, table[G[g].Get_LT()].first);//if divisor divides term, then term = divisor*remainder
	  Lab_poly vG = G[g].Multiply(v);
	  //cout << "vG after multi: " << G[g] << endl;
	  if (Smaller_sign(vG.Get_sign(),T)) //check if S(vG) < T
	    {
	      
	      return vG;
	    }
	}
      ++g;
    }
 //debug check: since called after TGFT, we know that some g in G must exist with that property!
  cout << "Error in function Find_in_G! No element in G found!" << endl;
  return G[0];
}
#endif


//if there exists some g in G with t = LT(vg), S(vg) < T, then return the index of g and set compatible = true
//else if there exists some g in G with t = LT(vg) S(vg) >= T, then return the index of g and set compatible = false
//else g = G.size(), compatible = false
bool Check_TGFT(const term_int t_, const signature& T)
{
  #ifdef PRECISE_DEBUG
  cout << "TGFT called!" << endl;
  #endif
  for (auto& it : G_conv) 
    {
      //*
      if (t_ == it.first)
	{
	  signature test = it.second.second;//if divisor divides term, then term = divisor*remainderv);//check: needed whole multiple?
	  if (Smaller_sign(test,T)) //check if S(vG) < T
	    {
	      return true;
	    }
	}
    }
  for (auto& it : G_conv)
    {
      if (t_ > it.first)
	{
	  if (isDivisible(table[t_].first,it.second.first))
	    {
	      signature test = it.second.second.Multiply(get_remainder(table[t_].first, it.second.first));//if divisor divides term, then term = divisor*remainderv);//check: needed whole multiple?
	      if (Smaller_sign(test,T)) //check if S(vG) < T
		{
		  return true;
		}
	    }
	}
    }
  return false;
}

bool Check_TGFT(const term_int t_, const signature& T, const int generation)//check only from generation to end
{
  for (auto& it : G_conv) 
    {
      //*
      if (t_ == it.first)
	{
	  signature test = it.second.second;//if divisor divides term, then term = divisor*remainderv);//check: needed whole multiple?
	  if (Smaller_sign(test,T)) //check if S(vG) < T
	    {
	      return true;
	    }
	}
    }
  for (auto& it : G_conv)
    {
      if (t_ > it.first && Generations_G.find(t_)->second > generation)
	{
	  if (isDivisible(table[t_].first,it.second.first))
	    {
	      signature test = it.second.second.Multiply(get_remainder(table[t_].first, it.second.first));//if divisor divides term, then term = divisor*remainderv);//check: needed whole multiple?
	      if (Smaller_sign(test,T)) //check if S(vG) < T
		{
		  return true;
		}
	    }
	}
    }
  return false;
}


Lab_poly Div_by_M(const term_int t, const signature& T)
{
  g_reduced++;
  Lab_poly current = Lab_poly();
  polynomial Old_M = polynomial();
  bool m_found = false;
  vector<int> t_ = table[t].first;
    if (!m_found)
    {
      for (auto& it : LT_M) 
	{
	  if (isDivisible(table[t].first,table[it.first].first))
	    {
	      term_vec v = get_remainder(table[t].first, table[it.first].first);//if divisor divides term, then term = divisor*remainder);
	      signature test = it.second.Multiply(v);
	      if (Smaller_sign(test,T)) //check if S(vG) < T
		{
		  m_found = true;
		  polynomial New_M = M.find(it.first)->second;
		  //Lab_poly candidate(G_.find(it.first)->second.Multiply(v),test);
		  if (is_more_reduced(New_M,Old_M)) //check different criteria
		    {
		      current = Lab_poly(New_M.Multiply(v),test);
		      Old_M = New_M;
		      //return candidate;
		    }
		}
	    }
	}
    }
    if (!m_found)
    {
      for (auto& it : LT_G) 
	{
	  if (isDivisible(table[t].first,table[it.first].first))
	    {
	      term_vec v = get_remainder(table[t].first, table[it.first].first);//if divisor divides term, then term = divisor*remainder);
	      signature test = it.second.Multiply(v);
	      if (Smaller_sign(test,T)) //check if S(vG) < T
		{
		  m_found = true;
		  polynomial New_M = G_.find(it.first)->second;
		  //Lab_poly candidate(G_.find(it.first)->second.Multiply(v),test);
		  if (is_more_reduced(New_M,Old_M)) //check different criteria
		    {
		      current = Lab_poly(New_M.Multiply(v),test);
		      Old_M = New_M;
		      //return candidate;
		    }
		}
	    }
	}
      }
  return current;
}

//#define SMALLEST_LT -> 
//#define LARGE_GEN -> nice
//#define LARGEST_LT //-> nice
//#define LARGEST_FLAG
//#define NUMBER_OF_TERMS -> very nice
//#define SMALLEST_REMAINDER -> coincides with largest_LT


bool is_more_reduced(const polynomial& one, const polynomial& two)
{
  if (two.Get_monomials().size() == 0) //starting off
    return true;
  term_int LT_1 = one.Get_LT();
  term_int LT_2 = two.Get_LT();
#ifdef LARGE_GEN
  return (Generations_M.find(LT_1)->second > Generations_M.find(LT_2)->second); //break ties arbitrarily
#endif
#ifdef LARGEST_LT
  return (LT_1 > LT_2);
#endif
#ifdef SMALLEST_REMAINDER
  return (LT_1 > LT_2);
#endif
#ifdef SMALLEST_LT
  return (LT_1 < LT_2);
#endif
#ifdef LARGEST_FLAG
  return (Flags_M.find(LT_1)->second > Flags_M.find(LT_2)->second);
#endif
#ifdef NUMBER_OF_TERMS
  return (LT_1 < LT_2);
  //return (M.find(LT_1)->second.Get_monomials().size() < M.find(LT_2)->second.Get_monomials().size());
#endif  
  
  return false;
}


//if there exists an g in M such that t = LT(g), which is tail-reduced w.r.t. G and signature T, return the index of g and flag_compatible = true
//else if there exists an g in M such that t = LT(g), but g is not tail-reduced w.r.t. G and signature T, return the index of g and flag_compatible = false
//else g = M_.size(), flag_compatible = false

#if 0
auto Check_in_M(const int t, const signature& T, const Multiples& M, const vector<Flag_poly>& M_, bool& flag_compatible)
{
  flag_compatible = false;
  int g = 0;
  int not_compat_idx = -1;
  auto itr = F_M.find(t);
  if (itr != F_M.end()) //element found
    {
      if (Smaller_sign(T,itr->second.Get_sign()))//if T < signature flag -> compatible
	{
	  flag_compatible = true;
	}
    }
  return *itr;

  
  while (g < M_.size())
    {
      if (t == M_[g].Get_Lab_poly().Get_LT())
	{
	  /*
	  if (T.Get_term().size() != M_[g].Get_Flag().Get_term().size())
	    {
	      cout << "Error in check in M!" << endl;
	      cout << "T: " << T.Get_term().size() << endl;
	      cout << "M_[g]: " << M_[g].Get_Flag() << endl;
	    }
	  */
	  if (Smaller_sign(T,M_[g].Get_Flag()))//compatible
	    {
	      //found a real reducer!
	      flag_compatible = true;
	      return g;	      
	    }
	  else
	    {
	      //found an element, but we would need to reduce it, maybe find a compatible one! -> not possible since unique? -> need to do exchange instaed of add!
	      not_compat_idx = g;
	    }
	}
      ++g;
    }
  if ((!flag_compatible) && not_compat_idx >= 0)//if no compatible one found, but at least one non-compatible found
    {
      g = not_compat_idx;
    }
  return g;

}
#endif

