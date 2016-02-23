##目的
用来熟悉一些数学算法，所有的代码都没有进行优化，以便于更好的理解算法本身。
##实现一些代数计算
##编译方法
windows下需要cygwin
```
autoreconf --install
./configure
make
```
###LU分解
考虑到数值精确性，每次将列中绝对值最大的作为主元。
###crout LU方法
crout方法比标准的LU分解要快一倍左右。
###Cholesky分解
[算法描述](http://mathfaculty.fullerton.edu/mathews/n2003/CholeskyMod.html)</br>
[算法](http://www.netlib.org/utk/papers/factor/node9.html)
###解逆矩阵
###解线性方程Ax=b
##微分方程的数值解
###欧拉法数值解
###中点法数值解
###Runge-Kutta法数据值解
###质数筛选算法
[Sieve of Eratosthenes](https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes)
###Gauss-Seidel迭代法解Ax=b
[Gauss–Seidel method](https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method)
###Jacobi迭代法解Ax=b
[Jacobi method](https://en.wikipedia.org/wiki/Jacobi_method#Description)
###线性互补问题(Linear Complementarity Probliem,LCP)
通过穷举得到全部解，因为程序很慢只能用来作为验证其他算法的佐证。</br>
The Projected Gauss-Seidel Method 该方法只在正定阵和主元占优矩阵收敛。

