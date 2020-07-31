#include "iostream"
#include "stdlib.h"
#include "math.h"
#include "malloc.h"
using namespace std;

typedef struct list
{
	long long *element;
	int step;
}list;//structure for every steps and corresponding elements
int bina[20];
long dec2bin(long);//convert decimal exponent to binary
long extEculidean(long,long,long *,long *,long);
long gcd(long,long);
long cal_inverse(long,long);//a*x+b*y=k,call extEculidean function
long long giant_calcul(long,long,long);//a square modular function
long long squaremodular(long,long,long);//call giant_calcul to calculate giant elements
long long squaremodular_baby(long,long,long,long);//calculate pow of inverse element of baby by squaremodular
long long baby_calcul(long,long,long,long);//modular of fraction
int main()
{
	long a,x,n,b,i,j,m;//a^xmodn=b
    list baby={NULL,0};
	list giant={NULL,0};
	cout<<"please input base, modular, consequence respectively:\n";
	cin>>a;
	cin>>n;
	cin>>b;
	if(gcd(a,n)!=1)//if a and n are relatively prime, there are no solution
	{
		cout<<"NO SOLUTION. base and modular are not relatively prime";
		return 0;
	}
	m=ceil(sqrt(n));//m is the number of times of iterations
	baby.element=(long long*)malloc(m*sizeof(long long));
    giant.element=(long long*)malloc(m*sizeof(long long));//allocate memory for every elements
    long long *baby0=baby.element;
    long long *giant0=giant.element;//record the head node
	for(i=0;i<=m-1;i++)
	{
		*(baby.element++)=baby_calcul(b,a,-i,n);//calculate every element of baby
		baby.step++;//keep pace with corresponding elements
		*(giant.element++)=squaremodular(a,m*i,n);//calculate every element of giant
		giant.step++;//keep pace with corresponding elements
    	//cout<<*(giant.element-1)<<"   giantstep"<<giant.step-1<<"   "<<*(baby.element-1)<<"   babystep"<<baby.step-1<<endl;
	}
    baby.element=baby0;
    giant.element=giant0;//initial pointers
	for(baby.step=0;baby.step<=m-1;baby.step++)//go through all elements of BSGS
	{
		for(giant.step=0;giant.step<=m-1;giant.step++)
		{
			if(*baby.element==*giant.element)//if a baby element equals a giant element, record the solution and return 0
			{
				x=baby.step+m*giant.step;//the solution should be this
				cout<<"\nthe exponent is"<<x<<" which is "<<a<<"^"<<x<<"mod"<<n<<"="<<b;
				return 0;
			}
			if(giant.step==m-1)
			giant.element=giant0;//initial pointers if no solution found in this iteration
			else
			giant.element++;
		    //cout<<*(giant.element)<<"   giantstep"<<giant.step<<"   "<<*(baby.element)<<"   babystep"<<baby.step<<endl;
		}
		if(baby.step==m-1)
		baby.element=baby0;
		else
		baby.element++;
	}
	free(baby.element);
	free(giant.element);//free the memory allocated
}


long dec2bin(long Deci)//convert decimal exponent to binary
{
	int i;	
	for(i=0;i<=19;i++)
	{
		if(Deci)
		{
		bina[i]=Deci%2;
		Deci/=2;//division algorithm
		}
		else
		return 0;
	}
}
 
 
long extEculidean(long a, long b, long *x, long *y, long c)
//when calculating baby steps, using extended eculidean algorithm to solve the inverse of the base
{
    if (b == 1)//when (a*x)modn=1, a and n are relatively prime and there are x and y satisfying ax+py=1
    {
        *x = 1;
        *y = -c/a;
        return a;
    }  
    long result = extEculidean(b, a%b, x, y, a);
    if (c != 0)
    {
        long temp = *x;
        *x = *y;
        *y = temp - c / a*(*y);
    }
    return result;
}
 
 //notice that gcd(a,b)=gcd(b,a%b)
long gcd(long a, long b)
{
    long temp;
    if (a < b)//make sure that a>b
    {
        long temp = b;
        b = a;
        a = temp;
    }
	return b == 0 ? a : gcd(b, a % b); //recursion for gcd
}
 
long cal_inverse(long a, long b)//a*x+b*y=k,call ExtendEculid function
{
    long x;
    long y;

    if (b > a)
    {
        a += b;
        long result = extEculidean(a, b, &x, &y, 0);
        y += x;
        a -= b;
    }
    else
    {
        long result = extEculidean(a, b, &x, &y, 0);
    }
 
    y %= a;
    if (y < 0)
    return y + a;
    else
    return y;
}



long long giant_calcul(long base,long exponent,long modular)//a square modular function
{
	int i;
	long long answer=1;
	long long factors[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};//build a LUT

	for(i=0;i<=19;i++)
	{
		if(i==0)
		factors[i]=base%modular;
		else		
		factors[i]=factors[i-1]*factors[i-1]%modular;//calculate values of every steps
        if(bina[i])//if this digit is one, multiply the corresponding number of factors
        {
        	answer*=factors[i];
		}
	}
	return answer%modular;//modular for the last result
}

long long squaremodular(long base, long exponent,long modular)//call giant_calcul to calculate giant elements
{
	int i;
	for(i=0;i<20;i++)
	bina[i]=0;
	dec2bin(exponent);//divert the exponent to binary for square modular algorithm
	return giant_calcul(base,exponent,modular);//call giant_cacul to calculate giant step
}


long long squaremodular_baby(long b,long base, long exponent,long modular)//calculate pow of inverse element of baby by squaremodular
{
	int i;
	for(i=0;i<20;i++)
	bina[i]=0;//initiate global array bina[i] first when called
	dec2bin(exponent);
	long long answer=1;
	long long factors[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};//initiate look up table

	for(i=0;i<=19;i++)
	{
		if(i==0)
		factors[i]=base%modular;
		else		
		factors[i]=factors[i-1]*factors[i-1]%modular;//calculate values of every steps
        if(bina[i])//if this digit is one, multiply the corresponding number of factors
        {
        	answer*=factors[i];
		}
	}
	return (answer*b)%modular;//modular for the last result and multiply b for baby steps
}

long long baby_calcul(long b, long base,long exponent,long modular)//modular of fraction
{
     long Y;
     for(Y=1;Y<=modular;Y++)//go through for a solution of b*base^(-exponent)mod(modular)=Y
     {
     	if(squaremodular_baby(b,cal_inverse(modular,base),-exponent,modular)==Y)
	    return Y;
	 }
}
