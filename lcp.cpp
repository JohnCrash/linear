
int doSolveN(real *A,real *b,int *N,real *x,int n)
{
	
}

int lcp(real *A,real *b,std::vector<real *>& xs,int n)
{
	int i;
	int * N = (int *)malloc(n*sizeof(int));
	int * x = NULL;
	memset(N,0,n);
	do{
		if(!x)
			x = (real*)malloc(2*n*sizeof(real));
		if( doSolveN(A,b,N,x,n) ){
			xs.push_back(x);
			x = NULL;
		}
		//模拟二进制进位
		for(i=0;i<n;i++){
			if( N[i]==0 ){
				N[i] = 1;
				break;
			}else
				N[i] = 0;
		}		
	}while(i<n);
	free(N);
	return 0;
}

int freeLcpSolve(std::vector<real *>& xs)
{
	for(auto i=xs.begin();i!=xs.end();++i)
		free( *i );
	xs.clear();
}