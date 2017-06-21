function [x,residue,iteration] = conjgrad(A,b,tolarence)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

x=b;
r=b-A*x;
residue=[];
iteration=[];
if norm(r)<tolarence
    return
end

y=-r;
z=A*y;
s=y'*z;
t=(r'*y)/s;
x=x+t*y;

for k = 1:numel(b);
    r=r-t*z;
    if(norm(r)<tolarence)
        return
    end
    B = (r'*z)/s;
       y = -r + B*y;
       z = A*y;
       s = y'*z;
       t = (r'*y)/s;
       x = x + t*y;
       residue=[residue;r(size(r,1))];
       iteration=[iteration;k];
end
end
    
       

