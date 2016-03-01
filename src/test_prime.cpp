#include <stdio.h>
#include <vector>
#include <algorithm>
#include <math.h>

/* Sieve of Eratosthenes */
/*
 * 方法如下：
 * 首先取10以内的质数，{2,3,5,7}
 * 然后计算100以内的质数，从11到100，将每个数和10以内的质数做除余运算
 * 如果余数等于零，表示不是质数。如果都不为零表示是一个质数。
 * 然后在找10000以内的质数，从101到10000，将每个数和100以内的质数做除余运算
 * 一次类推。
 * 为什么这个方法有效，因为N以内的合数（质数以外的数）。可以表示为Q=p*q
 * 首先p,q < N，sqrt(N)以内的质数被我们找出来并且用来除余运算排除掉了大量的数。
 * 这些数都可以写成Q=m*n,其中m和n肯定可以在sqrt(N)的质数中找到。那么反过来说
 * p和q不在sqrt(N)以内的质数中。并且p和q都大于sqrt(N)。这样导致Q>N,因此经过除余
 * 排除过的小于N的合数都被排除了。
 */
typedef std::vector<int> PrimeArray;
PrimeArray p10={2,3,5,7};

bool isPrime(PrimeArray& pa,int n)
{
	for(auto i=pa.begin();i!=pa.end();i++){
		if(n%(*i)==0)
			return false;
	}
	return true;
}

void nextPrime(PrimeArray& pa,int N)
{
	PrimeArray pp;
	int n,pn;
	pn=n=10;
	if( N==1 ){
		pa.resize(p10.size());
		std::copy(p10.begin(),p10.end(),pa.begin());
		return;
	}
	for(int i=2;i<=N;i++){
		pn = n;
		n = n*n;
	}	
	for(int i=pn;i<=n;i++){
		if(isPrime(pa,i)){
			pp.push_back(i);
		}
	}
	auto size = pa.size();
	pa.resize(size+pp.size());
	std::copy(pp.begin(),pp.end(),pa.begin()+size);
}

/*
 * 筛选素数，N=1 10以内的素数，N=2 100以内，N=3 10000
 */
void searchPrime(int N)
{
	PrimeArray pa;
	for(int i=1;i<=N;i++){
		nextPrime(pa,i);
	}
	printf("N=%d prime %d\n",N,pa.size());
	for(auto i=pa.begin();i!=pa.end();i++){
		printf("%d\n",*i);
	}
	printf("N=%d prime %d\n",N,pa.size());
}

int main(int argn,char *argv[])
{
	//searchPrime(2); //10的2次方,打印100以内的全部素数
	//searchPrime(3); //10的4次方,打印1万以内的全部素数
	searchPrime(4); //10的8次方,打印1亿以内的全部素数
}