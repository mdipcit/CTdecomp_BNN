% Author: Shunsuke Ono (ono@sp.ce.titech.ac.jp)

function[y] = ProxPreBNN(x, gamma, blocksize)

n = size(x);
y = zeros(n);
bnum = [fix(n(1)/blocksize),fix(n(2)/blocksize)];
for i = 1:n(3)
    for j = 1:n(4)
        for k = 1:bnum(1)
            for l = 1:bnum(2)
                [U, S, V] = svd(x(1+blocksize*(k-1):blocksize*k,1+blocksize*(l-1):blocksize*l,i,j));
                y(1+blocksize*(k-1):blocksize*k,1+blocksize*(l-1):blocksize*l,i,j)= U*max(S-gamma,0)*V';
            end
        end
    end
end