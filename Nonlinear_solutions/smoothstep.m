function f = smoothstep(x,N)
%     f = -2*x.^3+3*x.^2;
%     f = -2*x.^3+3*x.^2;
    f = zeros(size(x));
    for n=0:N
        f = f + (-1)^n*nchoosek(N+n,n)*nchoosek(2*N+1,N-n)*x.^(N+n+1);
    end
    f(x<=0) = 0;
    f(x>=1) = 1;
end