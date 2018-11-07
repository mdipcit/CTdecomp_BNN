% Author: Shunsuke Ono (ono@sp.ce.titech.ac.jp)

function[u] = PeriodicExpansionTrans(Pu, shiftstep)

n = size(Pu);
u = zeros(n(1:2));
for i = 1:n(3)
    for j = 1:n(4)
        u = u + circshift(Pu(:,:,i,j), [(-i+1)*shiftstep, (-j+1)*shiftstep]);
    end
end
u = u/sqrt(n(3)*n(4));
