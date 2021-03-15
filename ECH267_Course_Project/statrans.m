function Phi=statrans(A)

syms s;
[M,N]=size(A); 
I=eye(M); 
B=(s*I-A);
B_1=inv(B);
Phi=ilaplace(B_1); 