#include <stdio.h>
#include <vector>
#include <algorithm>
#include <math.h>

/* Sieve of Eratosthenes */
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