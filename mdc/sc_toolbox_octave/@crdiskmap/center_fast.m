function [Mnew,A_D] = center_fast(M,wc,zd)
%CENTER_FAST Recentre a crossratio disk map when the preimage is known.
%   [MNEW,A_D] = CENTER_FAST(M,WC,ZD) returns a map conformally
%   equivalent to M whose conformal center is WC, given ZD with
%   EVAL(M,ZD) == WC. CENTER does the same thing but recovers the
%   preimage of WC with CRIMAP0, an ODE solve followed by Newton
%   iterations; here it is already known, so the Moebius transform that
%   CRFIXWC builds can be composed in closed form instead.
%
%   The second output is the preimage under MNEW of the conformal center
%   of M, i.e. what EVALINV(MNEW,CENTER(M)) would return.
%
%   See also CENTER, CRFIXWC, CRSPREAD.

p = polygon(M);
w = vertex(p);
cr = M.crossratio;
Q = M.qlgraph;
n = length(w);

mt_old = M.center{2}(2:5);
u = (mt_old(1)*[zd,0] + mt_old(2)) ./ (mt_old(3)*[zd,0] + mt_old(4));
ul = crspread(u,M.center{2}(1),cr,Q);

for quadnum = 1:n-3
  if isinpoly(wc,w(Q.qlvert(:,quadnum)),1e-4)
    break
  end
end

zc = ul(quadnum,1);
zold = ul(quadnum,2);

z = crembed(cr,Q,quadnum);
a = sign((1-zc'*z(n))/(z(n)-zc));

M.center = {wc,[quadnum 1 zc*a zc' a]};
Mnew = M;

A_D = a*(zold-zc)/(1-zold*zc');
