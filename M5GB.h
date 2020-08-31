#ifndef M5GB_H_INCLUDED
#define M5GB_H_INCLUDED
#include <vector>
#include "lab_poly.h"
#include <set>
#include <unordered_map>
#include <map>

using namespace std;

struct Pair_compare {
  bool operator()(const pair<Lab_poly,Lab_poly>& one, const pair<Lab_poly,Lab_poly>& two) const{
  if (one.first.Get_sign() == two.first.Get_sign())
    {
      if (one.second.Get_sign() == two.second.Get_sign())
	{
	  /*
	  if (one.second.Get_poly() ==  two.second.Get_poly()); //break tie arbitrarily
	  {
	    cout << "Problem encountered!" << endl;
	  }
	  */
	  return (one.first.Get_poly() < two.first.Get_poly()); //break tie arbitrarily
	}
    }
  //else
  return (one.first.Get_sign() < two.first.Get_sign());
}
};

typedef set<signature> Signatures;

typedef vector<vector<term_int>> Syzygies; //H[i] contains all syzygies with index i

typedef multimap<signature,pair<Lab_poly,int>> Spair_map;

typedef set<pair<Lab_poly,Lab_poly>, Pair_compare> Spairs; //corresponds to P


//M
typedef unordered_map<term_int, polynomial> Multiples; //corresponds to M -> as a map from LT_G to increase find?
//G
typedef unordered_map<term_int, polynomial> Basis; //corresponds to G, same structure than Multiples

typedef unordered_map<term_int, bool> correct_label; //indicates if signature of m is still correctly labelled



/////////////////// The following structures are all for G and M maps for faster access //////////////////////////

//G/M:
typedef unordered_map<term_int, signature> LT_sig_map; 

//G/M -> for later improvement
typedef unordered_map<term_int, pair<term_vec,signature>> Fast_conversion_map;//for LT(g) in integer -> LT(g) in vector form, same for M

//G/M/T
typedef unordered_map<term_int, int> Generations;
//first part: LT(M), second contains generation

//M/T
typedef unordered_map<term_int, signature> Flags;
//first part: term t, second contains flag(t)

extern int p;

///////////////////////// Main functions ////////////////
//extern Table table;

vector<polynomial> M5GB(vector<polynomial>& system, int n); //table might overflow -> need to redefine inside
//Basis M5GB(vector<polynomial>& system, int n); //table might overflow -> need to redefine inside

//void Full_S_Reduce(Lab_poly &F, Multiples &M, const Basis& G_, const int n, const int m, Signed_terms& Irreducible);
void Full_S_Reduce(Lab_poly &F, const int n);


Lab_poly Tail_S_Reduce(const Lab_poly& F_,const signature& T, const int n);
//polynomial Tail_S_Reduce(const Lab_poly& F_,const signature& T, Multiples& M, const Basis& G_, int &rec_depth, const int n, const int m, Signed_terms& Irreducible);

/////////////////////////// Criteria ////////////////////

bool First_Buchberger(const vector<int>& LT_F, const vector<int>& LT_G);

bool Syzygy(const signature& sign, const Syzygies& H);
//bool Syzygy(const signature& S, const vector<signature>& H);
//bool Syzygy(const signature& S, const vector<int>& H);
//bool Syzygy(const int sign_term, const vector<int>& H_);

bool isRewritable(const Lab_poly& F, const int generation);
//bool isRewritable(const Lab_poly& F, const Basis& G_);



////////////////////////// Update functions //////////////

Lab_poly Update_m(const Lab_poly& m, const signature& T, const int n);//needs huge improvment -> ineffective!

//void Update_pairs(const Lab_poly& F, Spairs& P, const Basis& G_, const vector<polynomial>& system, int n, int& BB_counter, Signatures& Spol_sign, Spair_map& Spairs_map);
void Update_pairs(const Lab_poly& F, const vector<polynomial>& system, const int n, int& BB_counter);

//void Update_step(const Lab_poly& F, Multiples& M, Signed_terms& Irreducible);
//void Update_step(const Lab_poly& F, Multiples& M);

//void Update_spol(const Lab_poly& F, vector<Lab_poly>& P, const vector<Lab_poly>& G);

//void Update_syz(int i, const Basis& G_, vector<int>& H_);
void Update_syz(const Lab_poly& Spol);
//void Update_syz(const Lab_poly& Spol, const Basis& G_, vector<signature>& H);

signature Create_Flag(const polynomial& f);

signature Create_term_flag(const term_int t);

/////////////////////////// Init functions ///////////////

void Initialize(vector<polynomial>& system, int n); //table might overflow -> need to redefine inside

void Interreduce(vector<polynomial>& system);

void Add_field_eq(vector<polynomial>& system,int n);

vector<Lab_poly> Basic_poly(const vector<polynomial>& system_, int n);

Spairs F_2(const vector<Lab_poly>& F, int n, int p);

vector<Lab_poly> F_2(const vector<Lab_poly>& F);

vector<polynomial> Poly(const Basis& G_,int n,int p);
//vector<polynomial> Poly(const vector<Lab_poly>& G, int n, int p);

Basis Create_basis(const vector<Lab_poly>& F);

Basis Init_G(const vector<Lab_poly>& F);

void Init_generations(const vector<Lab_poly> &F);

void Init_LT_G(const vector<Lab_poly> & F);

Spairs Init_Spairs(const vector<Lab_poly>& F, Signatures& Spol_sign, Spair_map& Spairs_map, int n);
//Spairs Init_Spairs(const vector<Lab_poly>& F);

Syzygies Init_Syzygies(const vector<Lab_poly> &F);
//vector<signature> Init_Syzygies(const vector<Lab_poly> &F);


///////////////////////// Finding functions /////////////

Lab_poly Low_pol(const signature& smallest_sign, const Lab_poly& uF, const int generation, int n);

//Lab_poly Find_in_G(const vector<int>& t, const vector<Lab_poly>& G, const signature& T);//Finds vG in G, with S(vG) < T, t = LT(vG)

//Lab_poly Check_TGFT(const vector<int>& t, const Basis& G_, const signature& T, bool &compatible);
//int Check_TGFT(const vector<int>& t, const vector<Lab_poly>& G, const signature& T, bool &compatible);//checks if term t = LT(vg) for some g in G, if some with S(vg) < T, then compatible = true
bool Check_TGFT(const term_int t, const signature& T);

bool Check_TGFT(const term_int t_, const signature& T, const int generation);//check only from generation to end

Lab_poly Div_by_M(const term_int t, const signature& T);


void Add_to_M(const Lab_poly& F);

void Add_to_G(const Lab_poly& F);

bool is_more_reduced(const polynomial& one, const polynomial& two);

#endif // M5GB_H_INCLUDED
