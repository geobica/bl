function [D_running,W_running,corresponding_A_W_running,anchor_str,append_str,direction_index,f_set,r] = layer_b_setup(infile_level,infile_level_computed)
centers_load = load(infile_level);
centers = centers_load(1).centers;
current_level = centers_load(1).current_level;
index_within_last = centers_load(1).index_within_last;
count_current_level_2 = numel(current_level);
prev_level_f_load = load(infile_level_computed);
prev_level_centers = prev_level_f_load(1).new_centers;
r = centers_load(1).r;

D_running = {};
W_running = {};
corresponding_A_W_running = {};
anchor_str = {};
append_str = {};
direction_index = {};
f_set = {};
for c_i = 1:count_current_level_2;
	c = current_level(c_i)+1;

	R_D_load = centers(c,6);
	D_running{1,2*c_i-1} = R_D_load{1};
	L_D_load = centers(c,7);
	D_running{1,2*c_i} = L_D_load{1};
	R_W_load = centers(c,12);
	W_running{1,2*c_i-1} = R_W_load{1};
	L_W_load = centers(c,13);
	W_running{1,2*c_i} = L_W_load{1};
	true_W_center_load = centers(c,2);
	corresponding_A_W_running{1,2*c_i-1} = true_W_center_load{1};
	corresponding_A_W_running{1,2*c_i} = true_W_center_load{1};
	anchor_str{1,2*c_i-1} = "";
	anchor_str{1,2*c_i} = "";
	append_str{1,2*c_i-1} = "R";
	append_str{1,2*c_i} = "L";
	direction_index{1,2*c_i-1} = 1;
	direction_index{1,2*c_i} = 2;
	prev_f = prev_level_centers(index_within_last(c_i)+1,11);
	f_set{1,2*c_i-1} = prev_f{1};
	f_set{1,2*c_i} = prev_f{1};

end;
end
