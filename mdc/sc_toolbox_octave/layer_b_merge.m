function result = layer_b_merge(prefix,nworkers,outfile_level_computed)
warning('off', 'Octave:possible-matlab-short-circuit-operator');
new_centers = {};
for worker = 1:nworkers;
	name = sprintf('%schunk_%d.mat',prefix,worker);
	part = load(name);
	new_centers = [new_centers; part(1).new_centers];
	delete(name);
end;
save("-v7", outfile_level_computed, "new_centers");
result = 0;
end
