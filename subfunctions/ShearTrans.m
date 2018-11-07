% Author: Shunsuke Ono (ono@sp.ce.titech.ac.jp)

function[x] = ShearTrans(x, theta, direction)

n = size(x);
switch direction
    case 'r'
        for i = 1:n(1)
            x(i,:) = circshift(x(i,:), [0 fix(((i - 1)*theta)/45)]);
        end    
    case 't'
        for j = 1:n(2)
            x(:,j) = circshift(x(:,j), [fix(((j - 1)*theta)/45) 0]);
        end        
    case 'l'
        for i = 1:n(1)
            x(i,:) = circshift(x(i,:), [0 fix(((i - 1)*(-theta))/45)]);
        end
    case 'b'
        for j = 1:n(2)
            x(:,j) = circshift(x(:,j), [fix(((j - 1)*(-theta))/45) 0]);
        end        
end