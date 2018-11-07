% Author: Shunsuke Ono (ono@sp.ce.titech.ac.jp)

function[u] = ProjDynamicRangeConstraint(u, range)

u(u < range(1)) = range(1);
u(u > range(2)) = range(2);