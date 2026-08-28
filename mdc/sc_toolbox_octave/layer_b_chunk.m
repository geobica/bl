function result = layer_b_chunk(infile_sc,infile_level,infile_level_computed,worker,nworkers,prefix)
warning('off', 'Octave:possible-matlab-short-circuit-operator');
warning('off', 'Octave:legacy-function');
[D_running,W_running,corresponding_A_W_running,anchor_str,append_str,direction_index,f_set,r] = ...
	layer_b_setup(infile_level,infile_level_computed);
count_2 = numel(D_running);
lo = floor((worker-1)*count_2/nworkers)+1;
hi = floor(worker*count_2/nworkers);
new_centers = layer_b_cells(lo,hi,D_running,W_running,corresponding_A_W_running, ...
	anchor_str,append_str,direction_index,f_set,r);
save("-v7", sprintf('%schunk_%d.mat',prefix,worker), "new_centers");
result = 0;
end
