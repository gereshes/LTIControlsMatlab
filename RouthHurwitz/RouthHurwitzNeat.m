function [ F ] = RouthHurwitzNeat( charEqnVec )
%Calculates the Routh Hurwitz table given an equation vector.
%
%Inputs:
%   charEqnVec - 1 by n vector - Characteristic equation 
%Outputs:
%   F - ceil(n/2) by 2 matrix - Routh Hurwitz matrix
%
% Ari Rubinsztejn
% a.rubin1225@gmail.com
% www.gereshes.com
roots(charEqnVec)
variables=length(charEqnVec);
Z=zeros(variables-2,round(variables/2));
if (mod(variables,2)~=0)
    charEqnVec=[charEqnVec,0];
    variables=length(charEqnVec);
end
R=NaN(2,variables/2);
for c=1:variables
    temp=charEqnVec(c);
    R(c)=temp;
end

F=[R;Z];
for row=3:variables-1
    for col = 1:(variables/2)-1
        F(row,col)=((F(row-1,1)*F(row-2,col+1)-(F(row-2,1)*F(row-1,col+1))))/F(row-1,1);
    end
end

end

