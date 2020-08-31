#ifndef OTHER_FCT_H_INCLUDED
#define OTHER_FCT_H_INCLUDED
#include<array>
#include <iostream>
#include <vector>
#include <unordered_map>

using namespace std;

typedef long int term_int; //term in integer representation
typedef vector<int> term_vec; //term in vector representation

typedef vector<pair<term_vec,vector<term_int>>> Table;
typedef vector<vector<term_int>> Binom_table;
typedef unordered_map<term_int, term_vec> Monom_table; //term as integer, term as vector
typedef unordered_map<term_int, vector<term_int>> Multiplic_table; //term as integer, vector of all possible multiplications as vector

//typedef vector<pair<vector<int>,vector<int>>> Table;

extern Monom_table monom_table;
extern Multiplic_table multiplic_table;
extern Table table;
extern Binom_table binom_table;

vector<int> Subtract(const vector<int>& f, const vector<int>& g,int p);//computes f-g over F_p

vector<int> Add(const vector<int>& f, const vector<int>& g,int p);//computes f-g over F_p

vector<int> create_degrees(int d, int n);

bool term_is_Zero(vector<int> term);

int Inverse(int c, int p);//compute c^-1 in F_p, can do a table for large fields

bool is_Equal(const vector<int>& t, const vector<int>& u); //checks if t < u for t,u terms, order is deglex

bool is_Smaller(const vector<int>& t, const vector<int>& u); //checks if t < u for t,u terms, order is deglex

void lcm(vector<int>& u, const vector<int>& LT_f, vector<int>& v, const vector<int> LT_g);

bool isDivisible(const vector<int>& term, const vector<int>& divisor);

int binom(int n,int k);

vector<vector<term_int>> Create_binom_table(int n, int d);

term_int term_to_int(term_vec term, int n); //computes entry in macaulay matrix

int lower_degrees(int d, int n); //gets the number of x_1^(d-1) (largest term of degree d-1), is highly inefficient!

int LT(const vector<int>& poly);

vector<int> Tail(const vector<int>& poly);

Table Create_table(int n, int d);//multiplication table

//Monom_table Create_table(int n, int d);//look up for int -> vector representation

//Multiplic_table Create_table(int n, int d);//look up for multiplication of 



term_vec vector_add(const term_vec& term, int j); //new term = term*x_j

void next_term(vector<int>& term); //degree revers lex, change if needed!

#endif
