% Author: Shunsuke Ono (ono@sp.ce.titech.ac.jp)

function[Pu] = PeriodicExpansion(u, blocksize, shiftstep)

n = size(u);
Pu = zeros(n(1), n(2), blocksize/shiftstep, blocksize/shiftstep);

for i = 1:blocksize/shiftstep
    for j = 1:blocksize/shiftstep
        Pu(:,:,i,j) = circshift(u, [(i-1)*shiftstep, (j-1)*shiftstep]);
    end
end
Pu = Pu/(blocksize/shiftstep);
