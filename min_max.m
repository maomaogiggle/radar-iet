function [x_opt,beta] = min_max(N,Vc,s,Ks)
x0 = ones(N,1);
es = ones(N,1);
mu = 1;
for i = 1:4
cvx_begin
    variable x(N,1); variable t;
    minimize (t + mu*(es-2*x0)'*x);
    subject to
                es'*x == Ks;
                0 <= x <= 1;
                for k = 1:size(Vc,2)
                    vsc = Vc(:,k).*conj(s);
                    Vsc = real(vsc*vsc')+1e-5*eye(N);
                    x'*Vsc*x <= t;
                end
                -t <= 0;
cvx_end
x0 = x;
end
[x0_s,I_s] = sort(x0,'descend');
x_opt = zeros(N,1);
x_opt(I_s(1:Ks)) = 1;
beta = (x_opt'*Vsc*x_opt)/(Ks*Ks);
end
