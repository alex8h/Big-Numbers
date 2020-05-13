#include "bn.h"
#include <fstream>
#include <time.h>
using namespace BNErrors;

int main()
{ try {
	bool flag = 1;
	srand(time(NULL));
	BigNum a,c,b,d,f,z;

	for(int i = 0; i < 1001; i++)	
	{
		BASE M = 1000; 
		BASE n = 1 + rand()%M;
		BASE m = 1 + rand()%M;
		a.creation(n);
		b.creation(m);
		c.creation(1);
		d.creation(1);
		f.creation(1);
		DIV Z;
		a = a.gen_obj();
		b = b.gen_obj();
		Z.q = a/b;
		Z.r = a%b;
		c = Z.q*b;
		d = c + Z.r;
		f = a - Z.r;

		if(a == d && f == c &&  Z.r < b)
		{
			flag &= 1;
			//cout<<i<<endl;
		}
		else
		{
			flag &= 0;
			cout<<"a "<<endl<<a<<endl<<endl;
			cout<<"b "<<endl<<b<<endl<<endl;
			cout<<"c "<<endl<<c<<endl<<endl;
			cout<<"f "<<endl<<f<<endl<<endl;
			cout<<"i = "<<i<<endl;
			break;
		}
		a.deallocation();
		b.deallocation();
		c.deallocation();
		d.deallocation();
		Z.q.deallocation();
		Z.r.deallocation();
	}
	if(flag)
		cout<<"Passed"<<endl;
	else
		cout<<"Failed"<<endl;
	  }

catch(BNError error)
       {switch(error)
	{	
		case(MEMORY_ERROR):
		{
			cout<<"MEMORY_ERROR"<<endl;
			break;
		}
		case(INCORRECT_LENGTH):
                {
                        cout<<"INCORRECT_LENGTH"<<endl;
                        break;
                }
	        case(INCORRECT_POSITION):
                {
                        cout<<"INCORRECT_POSITION"<<endl;
                        break;
                }
		case(NULL_OPERAND):
                {
                        cout<<"NULL_OPERAND"<<endl;
                        break;
                }
		case(OUTPUT_ABROAD):
                {
                        cout<<"OUTPUT_ABROAD"<<endl;
                        break;
                }
		case(INCORRECT_SHIFT_VALUE):
                {
                        cout<<"INCORRECT_SHIFT_VALUE"<<endl;
                        break;
                }
	default : 
		{}
        }
	}

}
