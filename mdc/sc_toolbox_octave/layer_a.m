function result = layer_a(infile_sc,infile_level,outfile_level_computed)
warning('off', 'Octave:possible-matlab-short-circuit-operator');
warning('off', 'Octave:legacy-function');
S = load(infile_sc);
centers_load = load(infile_level);
centers = centers_load(1).centers;
A_D_load = centers(1,5);
R_D_load = centers(1,6);
L_D_load = centers(1,7);
true_W_center_load = centers(1,2);
anchor_str_load = centers(1,1);
D_running = {A_D_load{1},R_D_load{1},L_D_load{1}};
corresponding_A_W_running = {true_W_center_load{1},true_W_center_load{1},true_W_center_load{1}};
anchor_str = {anchor_str_load{1},anchor_str_load{1},anchor_str_load{1}};
append_str = {"A","R","L"};
direction_index = {0,1,2};

f = S(1).f;

count = size(D_running);
count_2 = count(2);
new_centers = {};
r = centers_load(1).r;

for c = 1:count_2;
	new_W_ = eval(f,double(D_running{c}));
	f_ = center(f,new_W_);
	A_D_from_ = evalinv(f_,double(corresponding_A_W_running{c}));
	A_D_from_ = A_D_from_*r/abs(A_D_from_);
	p = polygon(f_);
	w = vertex(p);
	beta = angle(p) - 1;
	z = get(f_,'pre');
	new_centers{c,1} = strcat(anchor_str{c},append_str{c});
	new_centers{c,2} = new_W_;
	new_centers{c,3} = corresponding_A_W_running{c};
	new_centers{c,4} = direction_index{c};
	new_centers{c,5} = A_D_from_;
	new_centers{c,6} = A_D_from_*exp(pi*2/3*1i);
	new_centers{c,7} = A_D_from_*exp(-pi*2/3*1i);
	new_centers{c,8} = w;
	new_centers{c,9} = z;
	new_centers{c,10} = beta;
	new_centers{c,11} = f_;
	new_centers{c,12} = eval(f_,A_D_from_*exp(pi*2/3*1i));
	new_centers{c,13} = eval(f_,A_D_from_*exp(-pi*2/3*1i));
end;
clearvars S centers_load;
save("-v7", outfile_level_computed);
result = 0;
end