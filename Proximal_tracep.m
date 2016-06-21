% min_X  0.5*||X-A||_F^2 + r*||X||_*p
% 0< p <= 2
function X = Proximal_tracep(A, r, p)

[n,c] = size(A);
if n<c
    A = A';
end;

[U S V] = svd(A,0);
s = diag(S);

if p == 1
    s1 = max(0,s-r);
else
    k = length(s);
    s1 = zeros(k,1);
    for i=1:k
        s1(i) = findrootp0(s(i), 2*r, p);
    end;
end;

S1 = diag(s1);
X = U*S1*V';

if n<c
    X= X';
end;

end