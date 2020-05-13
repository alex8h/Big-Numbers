#include <iostream>
#include <stdlib.h>
typedef unsigned short BASE;
typedef unsigned  int DBASE;
typedef int BNError;
using namespace std;

namespace BNErrors
{

  static const BNError OK = 0;
  static const BNError MEMORY_ERROR = 1;
  static const BNError INCORRECT_LENGTH = 2;
  static const BNError INCORRECT_POSITION = 3;
  static const BNError NULL_OPERAND = 4;
  static const BNError INCORRECT_SHIFT_VALUE = 5;
  static const BNError OUTPUT_ABROAD = 6;
  static const BNError ERROR = 7;
};

class BigNum{
					BASE *mass;
					int len;
					int H;


			public:
					BigNum();
					BigNum(BASE);
					~BigNum();
					BigNum(BigNum&);

					void creation(int);
					void deallocation();

					int GetLen() const;
					void putv(BASE,int);
					BASE getv(int);
					int sym_to_int(char);
					BigNum correctlen();
					char int_to_sym(int) const;
					BigNum gen_obj();
					BigNum shift_left(BASE);
					BigNum shift_right(BASE);
                    void BNtoINT(int,int,int **,int);
                    BigNum INTtoBN(int,int,int **,int);

					BigNum operator*(BigNum&);
					BigNum operator/(BigNum&);
					BigNum operator%(BigNum&);
					BigNum & operator=(const BigNum&);
					BigNum operator+(BigNum&);
					BigNum operator-(BigNum&);
					BigNum operator^(BigNum&);
					BigNum operator++();
                    BigNum operator--();
					bool operator==(const BigNum&);
					bool operator>(const BigNum&);
					bool operator>=(const BigNum&);
					bool operator<(const BigNum&);
					bool operator<=(const BigNum&);
					bool operator!=(const BigNum&);
					friend std::ostream & operator <<(std::ostream &r,const BigNum &s);
					friend std::istream & operator >>(std::istream &r, BigNum &s);
					void cout_bin();

					int compare(const BigNum&);
//					void GetElement(int, BigNum&);
					BigNum divbase(BigNum&);
					BigNum mul_base(BASE);
					BigNum diff(BigNum&);
					struct DIV div(BigNum&);
					BigNum degree_mod(BigNum&,BigNum&);
					BigNum divmod(BigNum&);
					BigNum mul_karatsuba(BigNum&);
					BigNum barret(BigNum&,BigNum&);
					bool check_prime();
					//BigNum gen_prime(BASE);
					BigNum bit_gen_prime(BASE);
					BigNum inverse_mod(BigNum&);
					int degree_of_two();
					BigNum NOD(BigNum&);
					BigNum NOD2(BigNum&,BigNum&,BigNum&);
					BigNum montg_mod(BigNum&,BASE&);
                    BigNum sqrt();

		};

struct DIV{
		BigNum q;
		BigNum r;
		};
