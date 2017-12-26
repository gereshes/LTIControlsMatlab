function [ F,Fdisp ] = RouthHurwitz( charEqnVec)
%Calculates the Routh Hurwitz table given an equation vector. This form can
%handle symbolic variables in the charachteristic equation.
%
%Inputs:
%   charEqnVec - 1 by n vector - Characteristic equation 
%Outputs:
%   F - ceil(n/2) by 2 matrix - Routh Hurwitz matrix
%   Fdisp - ceil(n/2) by 2 Cell - Routh Hurwitz matrix in cell format so it's
%                               easier to read
%
% Note: If there weill be symbolic variables define them using
%   syms before the function and then just enter those variables
%   in the charEqnVec Ex:[1,2,b,1,c]
%
% Ari Rubinsztejn
% a.rubin1225@gmail.com
% www.gereshes.com
variables=length(charEqnVec);
Z=cell(variables-2,round(variables/2));
for c=1:variables-2
    for d=1:round(variables/2)
        Z{c,d}=0;
    end
end
if (mod(variables,2)~=0)
    charEqnVec=[charEqnVec,0];
    variables=length(charEqnVec);
end
R=cell(2,variables/2);
for c=1:variables
    temp=charEqnVec(c);
    R{c}=temp;
end

F=[R;Z];
for row=3:variables-1
    for col = 1:(variables/2)-1
        F{row,col}=(((F{row-1,1}*F{row-2,col+1}-(F{row-2,1}*F{row-1,col+1})))/F{row-1,1});
    end
end
Fdisp=cell(size(F));
for c=1:variables-1
    for d=1:round(variables/2)
%          F{row,col}=simplify(F{row,col});
        Fdisp{c,d}=char(F{c,d});
    end
end
format long
disp(charEqnVec)
disp(Fdisp)
format short
end

