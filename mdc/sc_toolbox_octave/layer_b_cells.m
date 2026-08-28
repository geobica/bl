function new_centers = layer_b_cells(lo,hi,D_running,W_running,corresponding_A_W_running,anchor_str,append_str,direction_index,f_set,r)
new_centers = {};
for c = lo:hi;
	f = f_set{c};
	new_W_ = double(W_running{c});
	[f_,A_D_from_] = center_fast(f,new_W_,double(D_running{c}));
	A_D_from_ = A_D_from_*r/abs(A_D_from_);
	p = polygon(f_);
	w = vertex(p);
	beta = angle(p) - 1;
	z = get(f_,'pre');
	k = c-lo+1;
	new_centers{k,1} = strcat(anchor_str{c},append_str{c});
	new_centers{k,2} = new_W_;
	new_centers{k,3} = corresponding_A_W_running{c};
	new_centers{k,4} = direction_index{c};
	new_centers{k,5} = A_D_from_;
	new_centers{k,6} = A_D_from_*exp(pi*2/3*1i);
	new_centers{k,7} = A_D_from_*exp(-pi*2/3*1i);
	new_centers{k,8} = w;
	new_centers{k,9} = z;
	new_centers{k,10} = beta;
	new_centers{k,11} = f_;
	new_centers{k,12} = eval(f_,A_D_from_*exp(pi*2/3*1i));
	new_centers{k,13} = eval(f_,A_D_from_*exp(-pi*2/3*1i));
end;
end
