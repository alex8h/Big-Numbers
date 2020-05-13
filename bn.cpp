#include <string.h>
#include "bn.h"
#define M 1
#define gr 100
typedef unsigned short BASE;
typedef unsigned int DBASE;
using namespace BNErrors;
#define N (sizeof(BASE)*2)   // = 4
#define Base (sizeof(BASE)*8) // = 16

BigNum::BigNum()
			    {
					creation(M);
				}

BigNum::BigNum(BASE k)
				{
					creation(k);
				}


void BigNum::creation(int n)
				{
					len = n;
					H = len + 1;
					mass = new BASE [H];
					for(BASE i = 0; i < H; i++)
							mass[i] = 0;
				}


BigNum::~BigNum()
				{
					deallocation();
				}


void BigNum::deallocation()
				{
					if (mass) delete [] mass;
					mass = NULL;
					H = 0;
					len = 0;
				}


int BigNum::GetLen() const
			{
				return len;
			}


void BigNum::putv(BASE v, int position)
			{
				mass[position] = v;
			}


BASE BigNum::getv(int i)
			{
				return mass[i];
			}


BigNum::BigNum(BigNum &a)
				{
					creation(a.len);
					for(int i = 0; i < len; i++)
					{
						mass[i] = a.mass[i];
					}
				}


BigNum BigNum::sqrt()
			{
				BigNum x,x0,two;
				two.mass[0] = 2;
				x = x0 = *this;
				while(true)
				{
					x0 = x;
					x = *this / x + x;
					x = x / two;
					if(x >= x0)
						return x0;
				}
			}

int BigNum::sym_to_int(char c)
            {
                if(c >= 48 && c <= 57)
                    return c-'0';
                else if (c >= 65 && c <= 70)
                    return c-55;
                else if (c >= 97 && c<= 102)
                    return c-87;
				else throw(OUTPUT_ABROAD);
            }


char BigNum::int_to_sym(int a) const
			{
				if ( a >= 0 && a <= 9)
					return a+48;
				else if ( a >= 10 && a<=15)
					return a+55;
				else
					throw(OUTPUT_ABROAD);
			}


int BigNum::compare(const BigNum &s)
			{
				if(len > s.len) return 1;
				if(len < s.len) return -1;
				for(int i = len-1; i >= 0; i--)
					{
						if(mass[i] > s.mass[i]) return 1;
						if(mass[i] < s.mass[i]) return -1;

					}
				return 0;
			}


bool BigNum::operator==(const BigNum &s)
			{
				if(compare(s) == 0)
					return true;
				return false;

			}


bool BigNum::operator>(const BigNum &s)
			{
				if( compare(s) == 1)
					return true;
				return false;

			}


bool BigNum::operator>=(const BigNum &s)
			{
				if(compare(s) == 0 || compare(s) == 1)
					return true;
				return false;
			}


bool BigNum::operator<(const BigNum &s)
			{
				if(compare(s) == -1)
					return true;
				return false;
			}


bool BigNum::operator<=(const BigNum &s)
			{
				if(compare(s) == -1 || compare(s) == 0)
					return true;
				return false;
			}


bool BigNum::operator!=(const BigNum &s)
			{
				if(compare(s) != 0)
					return true;
				return false;
			}



BigNum BigNum::operator+(BigNum &s)
				{
					DBASE a;
					BASE k = 0;
					BigNum u,v,w;
				    int n = N*4;
					DBASE base = 1;
					base = base<<(n);
					int l = max(len,s.len);
					w.creation(l+1);
					if(*this >= s)
					{
						u = *this;
						v = s;
					}
					else
					{
						u = s;
						v = *this;
					}

					int i,j;
					for(i = 0; i < v.len; i++)
					{
							a = u.mass[i] + v.mass[i] + k;
							w.mass[i] = a%base;
							k = a/base;

					}
					for(j = i; j < u.len; j++)
					{
							a = u.mass[j] + k;
							w.mass[j] = a%base;
							k = a/base;
					}
					for(int y = j; y < w.len; y++)
					{
							a = w.mass[y] + k;
							w.mass[y] = a%base;
							k = a/base;
					}
					w = w.correctlen();
					return w;

				}



void BigNum::BNtoINT(int n,int m,int ** s,int num_str) // кол-во строк,кол-во столб.,матр,номер стр
{

	for(BASE j = 0; j < len; j++)
		s[num_str][j] = mass[j];
}


BigNum BigNum::INTtoBN(int n,int m,int ** s,int num_str)
{
	creation(m);
	for(BASE j = 0; j < m; j++)
	{
		mass[j] = s[num_str][j];
	}
	*this = correctlen();
	return *this;
}


BigNum BigNum::degree_mod(BigNum &s,BigNum &mod)
				{
					BigNum z,q;
					z.creation(mod.len);
					q = *this;
					if(s.mass[0] & 1)
						z = *this;
					else
						z.mass[0] = 1;
					for(int i = 1; i < 16*s.len ; i++)
					{
						//cout<<"  i  "<<i<<endl;

						q = (q*q)%mod;

						//cout<<"q  "<<q<<endl;

						BASE c = i/(N*4);
						BASE r = i%(N*4);
						BASE b = 1<<r;
						//int h = s.mass[c] & b;
						bool e = (s.mass[c] & b) >>r ;
						//cout<<" e "<<e<<endl;
						if(e == true)
						{
							z = (z*q)%mod;
						}
					}
					z = z.correctlen();
					return z;
				}


BigNum BigNum::mul_karatsuba(BigNum &s)
				{
					BigNum w,u0,u1,v0,v1,a,b,c,k;
					BASE l,n;
					l = (max(len,s.len) + 1)>>1;
					u0.creation(l);
					v0.creation(l);

					if(len < l)
						n = len;
					else
						n = l;

					for(BASE i = 0; i < n; i++)
							u0.mass[i] = mass[i];

					if(s.len < l)
						n = s.len;
					else
						n = l;

					for(BASE i = 0; i < n; i++)
						v0.mass[i] = s.mass[i];

					//u1 = *this;
					//v1 = s;
					//u1 = u1.shift_right(l);

					u1.creation(1);
					int m = len - l;
					if(m > 0)
					{
						u1.creation(m);
						for(BASE i = 0; i < u1.len; i++)
							u1.mass[i] = mass[i + l];
					}

					//u1.correctlen();

					v1.creation(1);
					m = s.len - l;
					if(m > 0)
					{
						v1.creation(m);
						for(BASE i = 0; i < v1.len; i++)
							v1.mass[i] = s.mass[i + l];
					}

					//v1 = v1.shift_right(l);

					if(l > gr )
					{
						a = u1.mul_karatsuba(v1);
						b = u0.mul_karatsuba(v0);
						k = u0 + u1;
						c = v0 + v1;
						c = k.mul_karatsuba(c);

					}
					else
					{
						a = u1 * v1;
						b = u0 * v0;
						k = u0 + u1;
						c = v0 + v1;
						c = k * c;
					}
					c = c - a - b;
					a = a.shift_left(l*2);

					for(BASE i = 0; i < b.len; i++)
						a.mass[i] = b.mass[i];

					c = c.shift_left(l);
					w = a + c;
					return w;
				}


BigNum BigNum::shift_left(BASE k)
				{
					if(k == 0)
						return *this;
					BigNum s;
					s.creation(len + k);
					for(BASE i = 0; i < len ; i++)
						s.mass[i + k] = mass[i];
					//s = s.correctlen();
					//deallocation();
					return s;
				}


BigNum BigNum::shift_right(BASE k)
				{
					if(k == 0)
						return *this;
					BigNum s(1);
					int n = len - k;
					if(n <= 0)
						return s;
					s.creation(n);
					for(BASE i = 0; i < s.len; i++)
						s.mass[i] = mass[i + k];
					//deallocation();
					//cout<<"shift  "<<s<<endl;
					s = s.correctlen();
					return s;

				}


BigNum BigNum::barret(BigNum &mod, BigNum &z)
				{
// перед функ. посч. b = 1; z = (b.shift_left(2*len.mod))/m;
					if(*this < mod)
						return *this;
					// 	ЧТО ТУТ ВООБЩЕ ПРОИСХОДИТ???
					BigNum q,r1,r2,r,b,x;
					x = *this;
					cout<<"x = "<<x<<endl;
					cout<<"m = "<<mod<<endl;
					q = this->shift_right(mod.len - 1);
					cout<<"q ="<<q<<endl;
					cout<<"z = "<<z<<endl;
					b = z.shift_right(mod.len + 1);
					cout<<b<<endl;
					q = q * b;
					cout<<q<<endl;
					b.deallocation();
					b.creation(1);
					b.mass[0] = 1;
					b.shift_left(len + 1);
					r1 = *this%b;
					r2 = (q * mod) % b;
					//r2 = r2%b;

					if(r1 >= r2)
						r = r1 - r2;
					else
						r = b + r1 - r2;
					while(r >= mod)
						r = r - mod;
					//r = r.correctlen();
					return r;

				}



bool BigNum::check_prime()
				{
					BASE t = rand()%50 + 100;
					//cout<<"T "<<t<<endl;
					int s = 0; // степень двойки
					BigNum one,two(1),three,p,w,c,r,k,f;
					one.mass[0] = 1;
					two.mass[0] = 2;
					three.mass[0] = 3;
					if(*this == one || *this == two || *this == three)
						return 1;
					p = *this;
					if(p <= three)
						return 0;
					if(!(p.mass[0] & 1))
						return 0;
					w = p - one;	// -1(mod p)
					r = w;
					c = p - two;
					s = r.degree_of_two();
					f.mass[0] = s;
					k = two ^ f;
					r = r / k;

					while(t)
					{
						BigNum b(p.len);
						b = b.gen_obj();
						b = b % w;
						if(b >= two && b <= c)
						{
							BigNum y;
							y = b.degree_mod(r,p);
							if(y != one && y != w)
							{
								for(int j = 1; j < s && y != w; j++)
								{
									y = y.degree_mod(two,p);
									if(y == one)
										return 0;
								}
								if(y != w)
									return 0;
							}

						}
						else
							continue;
						t--;
					}
					return 1;
				}


BigNum BigNum::bit_gen_prime(BASE bits)
				{
					BASE l = bits/Base;
					BASE b = bits - l*Base;
					//cout<<"l "<<l<<endl;
					//cout<<"b "<<b<<endl;
					bool e = false;
					//cout<<"453"<<endl;
					while(!e)
					{

						if(bits%Base == 0)
						{
							creation(l);
							this->gen_obj();
							mass[0] |= 1;
							mass[len - 1] |= 1<<(Base - 1);

						}
						else
						{
							creation(l + 1);
							this->gen_obj();
							mass[0] |= 1;
							mass[len - 1] |= 1<<(b - 1);
							BASE inv = ~0;
							inv = inv>>(Base - b);
							mass[len - 1 ] &= inv;

						}

						//cout<<*this<<endl;
						e = this->check_prime();
					}
					return *this;
				}


BigNum BigNum::montg_mod(BigNum &m, BASE &z)
				{
					BigNum y,one,q;
					one.mass[0] = 1;
					BASE base = ~0;
					y = *this;
					for(BASE i = 0; i < m.len; i++)
					{
						BigNum k;
						BASE u;
						u = (z * y.mass[i]) % base;
						k =  m.mul_base(u);
						cout<<"k  "<<k<<endl;
						k = k.shift_left(i);
						cout<<"k  "<<k<<endl;
						y = y + k;

					}
					y = y.shift_right(m.len);
					if(y >= m)
						y = y - m;
					//y = y.correctlen();
					return y;
				}


int BigNum::degree_of_two()
				{
					int flag = 0;
					BigNum null;
					null.mass[0] = 0;
					if(*this == null)
						return -1;
					for(BASE i = 0; i < len; i++)
						for(BASE j = 0; j < Base; j++)
						{
							BASE mask = 1<<j;
							BASE e = (mass[i] & mask) >> j;
							if(e == 0)
								flag++;
							else
								return flag;
						}
                    return flag;

				}


BigNum BigNum::NOD2(BigNum &b, BigNum &x, BigNum &y)
				{
					BigNum a,null,one;
					null.mass[0] = 0;
					one.mass[0] = 1;
					if(*this == null)
					{
						x = null;
						y = one;
						return b;
					}
					a = *this;
					BigNum x1,y1,d,a1;
					a1 = b % a;
					d = a1.NOD2(a,x1,y1);
					x = b / a;
					x = x * x1;
					if(y1 >= x)
						x = y1 - x;
					else
						x = x - y1;
					y = x1;
					return d;
				}

BigNum BigNum::NOD(BigNum &r)
				{
					BigNum a,a1,b,b1,k,c,c1;
					BigNum p,f1,f2,f3,two,one;
					a = *this;
					b = r;
					BASE flag1 = 0,flag2 = 0,flag3 = 1;
					one.mass[0] = 1;
					two.mass[0] = 2;

					if(!(a.mass[0] & 1))
						flag1 = a.degree_of_two();
					f1.mass[0] = flag1;
					f1 = two ^ f1;
					if(f1 != one)
						a1 = a / f1;
					else
						a1 = a;

					if(!(b.mass[0] & 1))
						flag2 = b.degree_of_two();
					f2.mass[0] = flag2;
					f2 = two ^ f2;
					if(f2 != one)
						b1 = b / f2;
					else
						b1 = b;

					BASE K = min(flag1,flag2);
					k.mass[0] = K;
					while(a1 != b1)
					{
						if(a1 < b1)
						{
							p = a1;
							a1 = b1;
							b1 = p;
						}
						c = a1 - b1;
						//if(!(c.mass[0] & 1))
						flag3 = c.degree_of_two();
						f3.mass[0] = flag3;
						f3 = two ^ f3;
						c1 = c / f3;
						a1 = c1;
					}
					k = two ^ k;
					a1 = k * a1;
					return a1;

				}


BigNum BigNum::inverse_mod(BigNum &m)
				{
					BigNum a0,a1,y0,y1,b,null(1),one(1),q;
					BASE i;
					null.mass[0] = 0;
					one.mass[0] = 1;
					a0 = m;
					a1 = *this;
					y0 = null;
					y1 = one;
					for(i = 0; a0 % a1 != null; i++)
					{
						BigNum u;
						u = a0 % a1;
						q = a0/a1;
						b = a1;
						a1 = q * a1;
						a1 = a0 - a1;
						a0 = b;
						b = y1;
						y1 = y1 * q;
						y1 = y0 + y1;
						y0 = b;
					}
					i--;
					if(i % 2 == 0)
						y1 = m - y1;
					return y1;
				}

BigNum BigNum::gen_obj()
				{
				for(int i = len - 1; i >= 0; i--)
					mass[i] = rand()%65536;
				return *this;
				}


BigNum BigNum::operator^(BigNum &s)
				{
					BigNum w,null,one,two;
					null.mass[0] = 0;
					one.mass[0] = 1;
					two.mass[0] = 2;
					if(s == null)
						return one;
					if(s == one)
						return *this;
					w = *this;
					for(BigNum i = two; i <= s; ++i)
						w = w * *this;
					return w;
				}


BigNum BigNum::operator++()
				{
					BigNum one;
					one.mass[0] = 1;
					*this = *this + one;
					return *this;
				}

BigNum BigNum::operator--()
				{
					BigNum one;
					one.mass[0] = 1;
					*this = *this - one;
					return *this;
				}

BigNum BigNum::operator-(BigNum &s)
				{
					BigNum o;
					o = diff(s);
					o = o.correctlen();
					return o;
				}


BigNum BigNum::diff(BigNum &s)
				{

					if(*this >= s)
					{
						BigNum b,d;
						b.creation(len);
						DBASE a,t;
						DBASE base = 1;
						int n = N*4;
						base = base<<(n);
						BASE k = 0;
						for(int i = 0; i < s.len ;i++)
						{

							a = mass[i] - s.mass[i] - k;
							b.mass[i] = a%base;
							t = base-1-a;
							k = t/base;
						}
						for(int i = s.len; i < len; i++)
						{
							a = mass[i] - k;
							b.mass[i] = a%base;
							t = base-1-a;
							k = t/base;
						}
						return b;
					}
					else
							throw(INCORRECT_LENGTH);


				}


BigNum BigNum::mul_base(BASE q)
				{
					BigNum b,d,qq;
					b.creation(len + 2);
					DBASE a;
					DBASE base = 1;
					base = base<<(N*4);
					BASE k;
					int i,j;
					for(i = 0,k = 0; i < len; i++)
					{
						a = mass[i] * q + b.mass[i] + k;
						b.mass[i] = a % base;
						k = a / base;
					}
					b.mass[i] = k;
					b = b.correctlen();
					return b;

				}


BigNum BigNum::operator*(BigNum &s)
				{
					if(s.len == 1)
					{
						BigNum x = *this;
						x = x.mul_base(s.mass[0]);
						return x;
					}
					BigNum b,d,u,v,w;
					int m = max(len,s.len);
					w.creation(len + s.len);
					DBASE a;
					DBASE base = 1;
					int n = N*4;
					base = base<<(n);
					int i,j;
					BASE k = 0;
					if(*this >= s)
					{
						u = *this;
						v = s;
					}
					else
					{
						u = s;
						v = *this;
					}
					for( j = 0; j < v.len; j++)
					{
							for(i = 0,k = 0; i < u.len; i++)
							{
								a = u.mass[i] * v.mass[j] + w.mass[i+j] + k;
								w.mass[i+j] = a%base;
								k = a/base;

							}

							w.mass[j+i] = k;
					}
					w = w.correctlen();
				    return w;

				}


BigNum & BigNum::operator=(const BigNum &s)
				{
					deallocation();
					creation(s.len);
					for(int i = 0; i < s.len; i++)
						{
							mass[i] = s.mass[i];

						}
					return *this;
				}


BigNum BigNum::divbase(BigNum &v)
				{
					if(v.mass[0] == 0)
						throw(NULL_OPERAND);
					if(v.mass[0] == 1)
						return *this;
					BigNum w;
					w.creation(len);
					DBASE k,base = 1;
					base = base<<(N*4);
					BASE q = 0;
					for(int i = len-1; i >= 0; i--)
					{
						k = base*q + mass[i];
						w.mass[i] = k / v.mass[0];
						q = k % v.mass[0];

					}
					return w;

				}


BigNum BigNum::divmod(BigNum &v)
				{
					BigNum w;
					w.creation(len);

					if(v.mass[0] == 0)
						throw(NULL_OPERAND);
					if(v.mass[0] == 1)
					{
						w = w.correctlen();
						return w;
					}
					DBASE k,base = 1;
					base = base<<(N*4);
					DBASE q = 0;
					for(int i = len-1; i >= 0; i--)
					{

						k = base*q + mass[i];
						q = k % v.mass[0];

					}
					w.mass[0] = q;
					w = w.correctlen();
					return w;
				}


BigNum BigNum::operator/(BigNum &s)
				{
					DIV Z;
					Z = this->div(s);
					return Z.q;
				}


BigNum BigNum::operator%(BigNum &s)
				{
					DIV Z;
					Z = this->div(s);
					return Z.r;
				}


struct DIV BigNum::div(BigNum &s)
				{
					BigNum null;
					null.creation(1);
					*this = this->correctlen();
					s = s.correctlen();
					if(s.len == 1)
					{
							BASE b = s.mass[0];
							DIV Z;
							BigNum x;
							x = *this;
							Z.q = x.divbase(s);
							Z.r = x.divmod(s);
							return Z;
					}
					//cout<<"711  "<<endl;
					DIV Z;
					BigNum w,d;
					DBASE base = 1;
					BASE n = 4*N;
					base = base<<(n);
					d.creation(1);
					d.mass[0] = base / (s.mass[s.len-1] + 1);
					//cout<<"719  "<<endl;
					BigNum u,v;
					u.creation(len);
					v.creation(s.len);
					u = *this*d;
					v = s*d;
					w.creation(u.len);
					if(u.len - len == 0)
					{
						u.len++;
						u.mass[u.len-1] = 0;

					}
					//cout<<"732  "<<endl;
					for(int i = len - s.len; i >= 0; i--)
					{

						DBASE q = (u.mass[i+s.len]*base + u.mass[i+s.len-1]) / v.mass[s.len-1];
						DBASE r = (u.mass[i+s.len]*base + u.mass[i+s.len-1]) % v.mass[s.len-1];

						while(r < base  &&  (q == base || q*v.mass[s.len-2] > base*r + u.mass[i+s.len-2]))
							{
									q--;
									r = r + v.mass[s.len-1];
							}
					//	cout<<"744  "<<endl;
						BigNum k;
						BigNum qq;
						qq.creation(1);
						qq.mass[0] = q;
						k = v.mul_base(q);
						BigNum y;
						y.creation(v.len+1);
						for(int j = i, l = 0; l < y.len; j++,l++)
								y.mass[l] = u.mass[j];
						w.mass[i] = q;

						if(y >= k)
								y = y.diff(k);
						else
						{
								y = k.diff(y);
								v.len++;
								v.mass[v.len-1] = 0;
								w.mass[i] -= 1;
								y = v.diff(y);
								v.len--;

						}
						//cout<<"768  "<<endl;
						for(int j = i,l = 0; l < y.len; j++,l++)
								u.mass[j] = y.mass[l];

					}
					//cout<<"772  "<<endl;
					BigNum R;
					R.creation(s.len);
					for(int i = 0; i < s.len; i++)
						R.mass[i] = u.mass[i];

					if(*this < s)
						Z.r = *this;
					else
						Z.r = R.divbase(d);
					Z.r = Z.r.correctlen();
					Z.q = w;
					Z.q = Z.q.correctlen();
					return Z;
				}


std::ostream & operator <<(std::ostream &r,const BigNum &s)
				{

					int len = s.len;
					BASE q;
					BASE buf;
					int flag = 1;
				    for(int i = len-1; i >= 0; i--)
				    {
						 q = s.mass[i];

						for (int j = N-1; j >= 0 ; j--)
						{
					   		buf = ((q >> (4*j))&15);
							if( flag == 0)
							{
						   		if (buf < 10 )
						    		  r << static_cast<char>(buf+'0');
						   		else
						    	  	r << static_cast<char>(buf+'a'-10);
							}
							else if (buf == 0 && flag == 1)
								continue;
							else
							{	flag = 0;
								j++;
							}
						}

					}
					if (flag == 1)
						r<<static_cast<char>(flag-1+'0');

					return r;

				}


std::istream & operator >>(std::istream &r, BigNum &s)
				{

					char q[1000000];  // создание символьного массива для записи символов
					//с клавиатуры
					//cout<<"Введите число"<<endl;
					cin>>q;     //  ввод строки
					int len = strlen(q);
					int l = (len-1)/N+1;
					//cout<<"Длина: "<<len<<endl;
					//   N = sizeof(BASE)*2 = 8
					s.deallocation(); // очистка объекта
					s.creation(l);  // создание нового объекта нужного размера
					//cout<<s.H<<endl;
					for(int i = 0; i < len; i++)
					{

						int b = s.sym_to_int(q[i]); //запись значений символов в массив
						int f = ((s.mass[(len-i-1)/N])<<4)+b;
						s.mass[(len-i-1)/N] = f;
						//cout<<b<<".";
					}

					return r;
				}


void BigNum::cout_bin()
				{
					for(int i = len - 1; i >= 0 ; i--)
						for(int j = Base - 1; j >= 0; j--)
						{
							BASE mask = 1 << j;
							BASE e = (mass[i] & mask) >> j;
							cout<<e;

						}
					cout<<endl;
				}


BigNum BigNum::correctlen()
                 {
    				bool flag = true;
					for(int i = len-1; i >= 0; i--)
					{
						if(mass[i] == 0)
							len--;
						else
						{
							flag = false;
							break;
						}
					}
					if(flag == true)
					{
						deallocation();
						creation(1);
					}
					else
						H = len + 1;
					return *this;

                 }
