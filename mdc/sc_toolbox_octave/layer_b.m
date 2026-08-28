function result = layer_b(infile_sc,infile_level,infile_level_computed,outfile_level_computed)
warning('off', 'Octave:possible-matlab-short-circuit-operator');
warning('off', 'Octave:legacy-function');
[D_running,W_running,corresponding_A_W_running,anchor_str,append_str,direction_index,f_set,r] = ...
	layer_b_setup(infile_level,infile_level_computed);
new_centers = layer_b_cells(1,numel(D_running),D_running,W_running,corresponding_A_W_running, ...
	anchor_str,append_str,direction_index,f_set,r);
save("-v7", outfile_level_computed, "new_centers");
result = 0;
end