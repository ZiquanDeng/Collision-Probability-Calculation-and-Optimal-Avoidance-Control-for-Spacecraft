function y=qzero2(fun,x0)
syms t
format long
f = fun;
df = diff(f,t);
k = 0;
R = vpa(subs(f,'t',x0));
y = x0;
while (abs(R)>1e-16)
    x1 = vpa(x0-subs(f,'t',x0)/subs(df,'t',x0));
    R = x1-x0;
    x0 = x1;
    k = k+1;
end
y = x0;
end
