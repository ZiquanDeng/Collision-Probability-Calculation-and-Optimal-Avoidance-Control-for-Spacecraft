function Q = jifen3(fun,xmin,xmax,ymin,ymax,zmin,zmax)
%UNTITLED4 Summary of this function goes here
%%Detailed explanation goes here
syms a b c
f=fun;
n = 10;
sum = 0;
da = (xmax-xmin)/n;
db = (ymax-ymin)/n;
dc = (zmax-zmin)/n;
x = xmin:da:xmax;
y = ymin:db:ymax;
z = zmin:dc:zmax;
for i=1:n+1
    for j=1:n+1
        for k=1:n+1
            sum = sum+(subs(f,{a,b,c},{x(i),y(j),z(k)}))*(da*db*dc); % sum=sum+f(x(i),y(j),z(k))*(dx*dy*dz);
        end
    end
end
Q = vpa(sum)
end



