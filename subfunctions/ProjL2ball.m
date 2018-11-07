% Author: Shunsuke Ono (ono@sp.ce.titech.ac.jp)

function[u] = ProjL2ball(u, f, epsilon)

radius = sqrt(sum(sum(sum((u - f).^2))));
if radius > epsilon
    u = f + (epsilon/radius)*(u - f);
end