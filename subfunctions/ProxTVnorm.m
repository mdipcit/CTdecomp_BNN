% Author: Shunsuke Ono (ono@sp.ce.titech.ac.jp)

function[Du] = ProxTVnorm(Du, gamma)

n = size(Du);
onemat = ones(n(1:2));
thresh = ((sqrt(sum(Du.^2, 3))).^(-1))*gamma;
thresh(thresh > 1) = 1;
coef = (onemat - thresh);

for k = 1:n(3)
        Du(:,:,k) = coef.*Du(:,:,k);
end

