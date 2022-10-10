function[L, R] = Permutations(m)
v = 1:m;
v1 = [v v]; 
A = diag(v1);
v2 = repelem(v,2);
B = diag(v2);
n = 2*m;
s = 1:n;
I = eye(n);
[~,a1] = sortrows(sort(A,2));
[~,b1] = sortrows(sort(B,2));
b1(b1) = s;
[~,a2] = sortrows(sort(A,1).');
[~,b2] = sortrows(sort(B,1).');
b2(b2) = s;
L = I(a1(b1),:);
R = I(:,a2(b2));
end