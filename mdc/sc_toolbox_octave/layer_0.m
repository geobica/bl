function result = layer_0(A,A_D,outfile)
warning('off', 'Octave:possible-matlab-short-circuit-operator');
warning('off', 'Octave:legacy-function');
p = polygon(A);
f = crdiskmap(p);
f = center(f,0);
p = polygon(f);
w = vertex(p);
beta = angle(p) - 1;
z = get(f,'pre');
R_D = A_D*exp(pi*2/3*1i);
L_D = A_D*exp(-pi*2/3*1i);
A_W = eval(f,A_D);
R_W = eval(f,R_D);
L_W = eval(f,L_D);
save("-v7", outfile);
%outfile
result = 0;
end