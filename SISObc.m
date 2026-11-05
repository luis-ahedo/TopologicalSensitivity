function [B,C] = SISObc(caso,n)


B =zeros(n,1);
C =zeros(1,n);
[row,col] =ind2sub([n n],caso);
C(row)=1;
B(col)=1;


end