% Author: Shunsuke Ono (ono@sp.ce.titech.ac.jp)

function[x] = ProjAverageConst(x, bias)
x = x - sum(sum(sum(x)))/(numel(x)) + bias;