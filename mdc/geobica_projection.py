import numpy as np
from math import e,pi,sqrt,sin,cos,tan,atan
import os
import matplotlib.pyplot as plt
from scipy import io
import scipy.special
import copy
from shapely.geometry import Polygon, MultiPolygon, Point
import scipy.integrate as integrate
import pickle
from scipy.optimize import minimize
import time
import sys
import oct2py
from oct2py import octave
oc = oct2py.Oct2Py()
D = 12756274

def f_W_ring(w,z,beta,power=8):
	theta = np.linspace(0,2*pi,2**power)
	eth_k = np.remainder(np.imag(np.log(z)),2*pi)
	eth_k[eth_k>2*pi-0.0000000000001]=0
	sums_d2 = np.prod((2-2*np.cos(theta[:,None]-eth_k[None,:]))**(beta_raw[None,:]/2),axis=1)
	prev_eth = np.zeros_like(theta)
	w_pos = np.zeros_like(theta+0j)
	for i in range(z.shape[0]):
		if eth_k[i]<eth_k[(i+1)%z.shape[0]]:
			prev_eth[np.logical_and(theta>=eth_k[i],theta<eth_k[(i+1)%z.shape[0]])] = i
	for i in range(z.shape[0]):
		wheres = np.where(prev_eth==i)[0]
		if wheres.shape[0]==1:
			total_len = eth_k[(i+1)%z.shape[0]]-eth_k[i]
			if total_len<=theta[1]:
				w_pos[wheres[0]] = w[i]
			else:
				lambd = np.nan_to_num((theta[wheres[0]]-eth_k[i])/total_len)
				w_pos[wheres[0]] = w[i]*(1-lambd)+w[(i+1)%z.shape[0]]*lambd
		elif wheres.shape[0]>1:	
			test = np.prod((2-2*np.cos(eth_k[i]-eth_k[:i]))**(beta_raw[:i]/2))*np.prod((2-2*np.cos(eth_k[i]-eth_k[i+1:]))**(beta_raw[i+1:]/2))
			test = np.nan_to_num(test)
			test3 = np.prod((2-2*np.cos(eth_k[(i+1)%z.shape[0]]-eth_k[:(i+1)%z.shape[0]]))**(beta_raw[:(i+1)%z.shape[0]]/2))*np.prod((2-2*np.cos(eth_k[(i+1)%z.shape[0]]-eth_k[(i+1)%z.shape[0]+1:]))**(beta_raw[(i+1)%z.shape[0]+1:]/2))
			test3 = np.nan_to_num(test3)

			#integral from eth_k[i] to theta[wheres[0]]:
			b = beta_raw[i]/2
			v = eth_k[i]
			br = beta_raw[(i+1)%z.shape[0]]/2
			vr = eth_k[(i+1)%z.shape[0]]

			integral_left = lambda x:-(np.sqrt(2)*np.sqrt(np.cos(v - x) + 1)*np.tan((v - x)/2)*(2 - 2*np.cos(v - x))**b*scipy.special.hyp2f1(1/2, b + 1/2, b + 3/2, np.sin((v - x)/2)**2))/(2*b + 1)
			integral_right = lambda x:-(np.sqrt(2)*np.sqrt(np.cos(vr - x) + 1)*np.tan((vr - x)/2)*(2 - 2*np.cos(vr - x))**br*scipy.special.hyp2f1(1/2, br + 1/2, br + 3/2, np.sin((vr - x)/2)**2))/(2*br + 1)
			
			integral_dtheta = np.zeros_like(theta[wheres]+0j)
			int_left_eval = integral_left(theta[wheres[0]])
			int_right_eval = integral_right(theta[wheres[-1]])
			if not np.any(np.isinf(int_left_eval)):
				integral_dtheta += test*int_left_eval
			sums_here = sums_d2[wheres]
			if np.any(np.isinf(sums_here[0])):
				sums_here[0] = sums_here[1]
			
			integral_dtheta[1:]+=np.cumsum((sums_here[:-1]+sums_here[1:])/2)*theta[1]

			total_len = integral_dtheta[-1]
			if not np.any(np.isinf(int_right_eval)):
				total_len-=test3*integral_right(theta[wheres[-1]]) #this is correct
			if total_len<=0:
				w_pos[wheres] = w[i]
			else:
				lambd = np.nan_to_num(integral_dtheta/total_len)
				w_pos[wheres] = w[i]*(1-lambd)+w[(i+1)%z.shape[0]]*(lambd)
				if np.any(np.isnan(w_pos[wheres])):
					w_pos[wheres] = w[i]

	return w_pos

def prepare_coefficients(w,z,beta,power=8):
	sample_density = 2**power
	f_M_z_in_pre = np.linspace(0,2*pi,sample_density,endpoint=False)

	f_M_z_loopup_table = f_W_ring(w,z,beta,power=power)

	t1_ln_tlu = np.log(np.absolute(f_M_z_loopup_table/D)**2+1)

	fft_get = -np.fft.fft(t1_ln_tlu)/2**power

	return fft_get

def sc_map_continuous(zp,z,beta,p):
	k = np.arange(p.shape[0])
	f_prime_W = lambda zeta : np.prod((1-zeta/z)**(beta),axis=0)*np.exp(np.sum(p*zeta**k))
	f_prime_W_inv = lambda zeta : np.prod((1-(1/zeta)/z)**(beta),axis=0)*np.exp(np.sum(p*np.conj(zeta)**(-k)))
	f_prime_W_arr = lambda zeta : np.prod((1-zeta[None,:]/z[:,None])**(beta[:,None]),axis=0)*np.exp(np.sum(p[:,None]*zeta[None,:]**k[:,None],axis=0))
	f_prime_W_arr_inv = lambda zeta : np.nan_to_num(np.prod((1-1/zeta[None,:]/z[:,None])**(beta[:,None]),axis=0)*np.exp(np.sum(p[:,None]*zeta[None,:]**k[:,None],axis=0)))/zeta


	wp = np.zeros_like(zp).astype('complex')

	ii = 0
	for i in range(wp.shape[0]):
		ii+=1
		if ii%1000==0:
			print(ii)
		wreal = integrate.quad(lambda zeta : np.real(f_prime_W(zp[i]/np.absolute(zp[i])*zeta)*zp[i]/np.absolute(zp[i])), 0, np.absolute(zp[i]))[0]
		wimag = integrate.quad(lambda zeta : np.imag(f_prime_W(zp[i]/np.absolute(zp[i])*zeta)*zp[i]/np.absolute(zp[i])), 0, np.absolute(zp[i]))[0]
		wp[i] = (wreal+1j*wimag)
	wp[zp==0] = 0
	return wp

def tlu_M_over_tlu_W(zp,z,beta,p):
	k = np.arange(p.shape[0])
	return np.absolute(np.exp(np.sum(p[:,None]*zp[None,:]**k[:,None],axis=0)))

def explore_polygon(poly_coords,preceeder,centers_pickled=0,inverted_buffer_metric=False):
	poly_string = 'A = ['
	x = np.real(poly_coords)
	y = np.imag(poly_coords)
	for i in range(x.shape[0]):
		poly_string = poly_string+str(x[i])
		if y[i]>=0:
			poly_string = poly_string+'+'+str(y[i])+'i'
		else:
			poly_string = poly_string+str(y[i])+'i'
		if i<x.shape[0]-1:
			poly_string = poly_string+','
		else:
			poly_string = poly_string+'];'
	#print(poly_string)
	sides = 7
	r = sqrt(-5 + sqrt(3*sides)*cos(1/3*atan(1/(3*sqrt(3)))) + 3*sqrt(sides)*sin(1/3*atan(1/(3*sqrt(3)))))

	if not os.path.isfile(f'pickle/matlab_saves/{preceeder}sc.mat'):
		oc.eval(poly_string)
		oc.eval(f'A_D = {r};')
		oc.eval(f'outfile = strcat(\'pickle/matlab_saves/{preceeder}sc.mat\');')
		oc.feval('sc_toolbox_octave/layer_0.m',oc.pull('A'),oc.pull('A_D'),oc.pull('outfile'))
	#data = np.loadtxt(f'pickle/matlab_saves/{preceeder}sc.mat')
	loaded = io.loadmat(f'pickle/matlab_saves/{preceeder}sc', verify_compressed_data_integrity=False)
	w_raw = loaded['w'][:,0]
	beta_raw = loaded['beta'][:,0]
	z_raw = loaded['z'][:,0]

	centers = []
	centers.append(['',0,loaded['A_W'][:,0],0,loaded['A_D'][:,0],loaded['R_D'][:,0],loaded['L_D'][:,0],w_raw,z_raw,beta_raw,'',loaded['R_W'][:,0],loaded['L_W'][:,0]])
	scipy.io.savemat(f'pickle/matlab_saves/{preceeder}level_0.mat', mdict={'centers': centers,'r': r})
	if not os.path.isfile(f'pickle/matlab_saves/{preceeder}level_0_computed.mat'):
		oc.eval(f'infile_sc = \'pickle/matlab_saves/{preceeder}sc.mat\';')
		oc.eval(f'infile_level = \'pickle/matlab_saves/{preceeder}level_0.mat\';')
		oc.eval(f'outfile_level_computed = strcat(\'pickle/matlab_saves/{preceeder}level_0_computed.mat\');')
		oc.feval('sc_toolbox_octave/layer_a.m',oc.pull('infile_sc'),oc.pull('infile_level'),oc.pull('outfile_level_computed'))

	loaded = io.loadmat(f'pickle/matlab_saves/{preceeder}level_0_computed', verify_compressed_data_integrity=False)
	new_centers = loaded['new_centers']
	current_level = []
	index_within_last = []
	for j in range(len(new_centers)):
		new_center_to_add = ['']
		for i in range(1,10):
			new_center_to_add.append(copy.copy(new_centers[j,i][:,0]))
		new_center_to_add.append('')
		for i in range(11,13):
			new_center_to_add.append(copy.copy(new_centers[j,i][:,0]))
		w_ = new_center_to_add[7]
		w_cen = new_center_to_add[1]
		w_point = Point(np.real(w_cen),np.imag(w_cen))
		current_level.append(len(centers))
		centers.append(new_center_to_add)
		index_within_last.append(j)

	scipy.io.savemat(f'pickle/matlab_saves/{preceeder}level_1.mat', mdict={'centers': centers,'r': r,'current_level': current_level,'index_within_last': index_within_last})

	done_exploring = False
	level_i = 0
	while not done_exploring:
		level_i += 1
		print(level_i)
		if level_i>=centers_pickled:
			oc.eval(f'infile_sc = \'pickle/matlab_saves/{preceeder}sc.mat\';')
			oc.eval(f'infile_level = \'pickle/matlab_saves/{preceeder}level_{level_i}.mat\';')
			oc.eval(f'infile_level_computed = \'pickle/matlab_saves/{preceeder}level_{level_i-1}_computed.mat\';')
			oc.eval(f'outfile_level_computed = strcat(\'pickle/matlab_saves/{preceeder}level_{level_i}_computed.mat\');')
			oc.feval('sc_toolbox_octave/layer_b.m',oc.pull('infile_sc'),oc.pull('infile_level'),oc.pull('infile_level_computed'),oc.pull('outfile_level_computed'))

		if level_i>=centers_pickled:
			loaded = io.loadmat(f'pickle/matlab_saves/{preceeder}level_{level_i}_computed', verify_compressed_data_integrity=False)
			new_centers = loaded['new_centers']
			current_level = []
			index_within_last = []
			#print(len(centers))
			for j in range(len(new_centers)):
				new_center_to_add = ['']
				for i in range(1,10):
					new_center_to_add.append(copy.copy(new_centers[j,i][:,0]))
				new_center_to_add.append('')
				for i in range(11,13):
					new_center_to_add.append(copy.copy(new_centers[j,i][:,0]))
				w_of_already_placed = np.array([center[1] for center in centers])
				w_cen = new_center_to_add[1]
				min_val = np.amin(np.hstack(np.absolute(w_of_already_placed-w_cen)))
				in_range = False
				if inverted_buffer_metric:
					if w_cen == 0:
						in_range = True
					else:
						buffering = np.absolute(1/w_cen)/80
						w_point = Point(np.real(1/w_cen),np.imag(1/w_cen))
						w_polygon = Polygon(np.concatenate((np.real(1/w_raw).reshape([w_raw.shape[0],1]),np.imag(1/w_raw).reshape([w_raw.shape[0],1])),axis=1)).buffer(buffering)
						in_range = w_polygon.contains(w_point)
				else:
					buffering = (np.absolute(w_cen/D)**2+1)*D/240
					w_point = Point(np.real(w_cen),np.imag(w_cen))
					w_polygon = Polygon(np.concatenate((np.real(w_raw).reshape([w_raw.shape[0],1]),np.imag(w_raw).reshape([w_raw.shape[0],1])),axis=1)).buffer(-buffering)
					in_range = w_polygon.contains(w_point)
				if min_val>0.01*np.absolute(new_center_to_add[11]-new_center_to_add[1]) and in_range and np.absolute(new_center_to_add[2]-new_center_to_add[1])<3*np.absolute(new_center_to_add[11]-new_center_to_add[1]):
					current_level.append(len(centers))
					centers.append(new_center_to_add)
					index_within_last.append(j)

			print(len(centers))
			with open(f'pickle/matlab_saves/{preceeder}level_{level_i}_centers.txt', 'wb') as fh:
				pickle.dump(centers, fh)
			with open(f'pickle/matlab_saves/{preceeder}level_{level_i}_current_level.txt', 'wb') as fh:
				pickle.dump(current_level, fh)
			with open(f'pickle/matlab_saves/{preceeder}level_{level_i}_index_within_last.txt', 'wb') as fh:
				pickle.dump(index_within_last, fh)
		elif level_i==centers_pickled-1:
			pickle_off = open(f"pickle/matlab_saves/{preceeder}level_{level_i}_centers.txt", "rb")
			centers = pickle.load(pickle_off)
			pickle_off = open(f"pickle/matlab_saves/{preceeder}level_{level_i}_current_level.txt", "rb")
			current_level = pickle.load(pickle_off)
			pickle_off = open(f"pickle/matlab_saves/{preceeder}level_{level_i}_index_within_last.txt", "rb")
			index_within_last = pickle.load(pickle_off)

		if level_i>=centers_pickled-1:
			if len(current_level)==0:
				done_exploring = True
			else:
				scipy.io.savemat(f'pickle/matlab_saves/{preceeder}level_{level_i+1}.mat', mdict={'centers': centers,'r': r,'current_level': current_level,'index_within_last': index_within_last})

	with open(f'pickle/matlab_saves/{preceeder}centers_file.txt', 'wb') as fh:
		 pickle.dump(centers, fh)

	return w_raw,beta_raw,z_raw,centers


def M_centers_from_W(centers_W,w_raw,beta_raw,z_raw):
	line_placements = []
	line_posers = []
	line_placements_w = []
	line_posers_w = []

	line_placements_raw_pick = np.array([])
	line_placements_raw = []
	line_placements_raw_w = []

	ii=0
	w_center_array = np.hstack(np.array([cent[1] for cent in centers_W]))
	w_anchor_array = np.hstack(np.array([cent[2] for cent in centers_W]))
	w_right_array = np.hstack(np.array([cent[11] for cent in centers_W]))
	w_left_array = np.hstack(np.array([cent[12] for cent in centers_W]))
	power = 8
	w_samples = []
	isolated_points = []
	isolated_points_w = []


	important = []
	target_resize = 0
	resize_factor = []
	previous_index = []
	exp_p_0 = []
	associates = []
	tlu_M_list = []
	extra_tlu_M_list = []

	extra_tlu_M_over_tlu_W_list = []

	for center in centers_W:
		if ii%200==0:
			print(f'{ii}/{len(centers_W)}')
		important.append(0)
		previous_index.append(0)

		associates.append(np.absolute(center[8]-np.roll(center[8],1))/np.sum(np.absolute(center[8]-np.roll(center[8],1))))
		p = prepare_coefficients(center[7],center[8],center[9],power)
		detail_w = np.ndarray.flatten(np.array(center[7])[:,None]*np.linspace(0,1,10,endpoint=False)[None,:]+np.roll(np.array(center[7]),1)[:,None]*(1-np.linspace(0,1,10,endpoint=False))[None,:])
		detail_w/=D
		a = center[1]/D
		detail_w = (detail_w-a)/(1+detail_w*np.conj(a))
		stereo_dist = np.amin(np.absolute(detail_w))

		exp_p_0 = sc_map_continuous(np.array([center[4]+0j]),center[8],center[9],np.ones((1))*p[0])[0]/(center[2]-center[1])
		
		extra_D = np.array([center[4],center[5],center[6],center[4],center[5],center[6],center[4],center[5],center[6],center[4],center[5],center[6],center[4],center[5],center[6]])[:,0]
		extra_D[3:6] *= 1/3
		extra_D[6:9] *= e**(pi*1j/3)/2
		extra_D[9:12] *= e**(pi*1j/3)/4
		extra_D[12:] *= e**(pi*1j/3)*3/4
		M_vals = sc_map_continuous(extra_D,center[8],center[9],p*2)
		W_vals = sc_map_continuous(extra_D,center[8],center[9],np.ones((1)))

		tlu_W = np.absolute(center[1]/D)**2+1
		tlu_M = np.absolute(tlu_W/np.exp(np.absolute(p[0])))
		exp_p_0 /= (1/(np.absolute(stereo_dist)**2+1)/tlu_M)
		tlu_M = 1/(np.absolute(stereo_dist)**2+1)

		extra_tlu_M_over_tlu_W = np.hstack([p[0]*2,*list(tlu_M_over_tlu_W(extra_D,center[8],center[9],p*2))])
		extra_tlu_M_over_tlu_W *= tlu_W/tlu_M
		
		M_vals = np.nan_to_num(np.hstack(M_vals))

		ii+= 1
		if len(line_placements)==0:
			line_placements.append([0,*list(M_vals)])
			line_placements_raw.append([0,*list(M_vals)])

			line_posers.append(1+0j)
			resize_factor.append(exp_p_0)
			w_samples.append(0+0j)
		else:
			# find center[2] in centers[:,1]
			anchor_i = np.argmin(np.absolute(w_center_array-center[2]))
			previous_index[len(previous_index)-1] = anchor_i
			# so this one's 0 corresponds to 
			this_cent = line_placements[anchor_i][int(1+center[3][0])]


			this_anchor = line_placements[anchor_i][0]
			w_samples.append(center[1])
			resize_factor.append(exp_p_0)

			line_placements_raw.append([0,*list(M_vals)])
			line_posers.append((this_anchor-this_cent)/M_vals[0])
			M_vals = this_cent+M_vals*(this_anchor-this_cent)/M_vals[0]

			line_placements.append([this_cent,*list(M_vals)])
		tlu_M_list.append(tlu_M)
		extra_tlu_M_over_tlu_W_list.append(extra_tlu_M_over_tlu_W)

		line_placements_w.append([0,*list(W_vals)])
		line_placements_raw_w.append([0,*list(W_vals)])

	with open(f'pickle/matlab_saves/{preceeder}line_placements.txt', 'wb') as fh:
		 pickle.dump(line_placements, fh)
	with open(f'pickle/matlab_saves/{preceeder}line_posers.txt', 'wb') as fh:
		 pickle.dump(line_posers, fh)
	with open(f'pickle/matlab_saves/{preceeder}line_placements_raw.txt', 'wb') as fh:
		 pickle.dump(line_placements_raw, fh)
	with open(f'pickle/matlab_saves/{preceeder}resize_factor.txt', 'wb') as fh:
		 pickle.dump(resize_factor, fh)
	rescaled_raw = np.absolute(np.hstack(resize_factor))

	line_posers = line_posers/np.absolute(line_posers)/rescaled_raw

	line_placements = []
	line_placements_w = []
	line_placements.append(list(np.hstack(line_placements_raw[0])*line_posers[0]))
	centers_M = []

	center = centers_W[0]
	detail_w = np.ndarray.flatten(np.array(center[7])[:,None]*np.linspace(0,1,10,endpoint=False)[None,:]+np.roll(np.array(center[7]),1)[:,None]*(1-np.linspace(0,1,10,endpoint=False))[None,:])
	detail_w/=D

	for i in range(0,len(line_placements_raw)):
		line_placements_w.append(w_center_array[i]+(w_anchor_array[i]-w_center_array[i])*np.hstack(line_placements_raw_w[i])/np.hstack(line_placements_raw_w[i])[1])
		extra_tlu_M_list.append(extra_tlu_M_over_tlu_W_list[i]*(np.absolute(line_placements_w[i]/D)**2+1))
		for ii in range(len(line_placements_w[i])):
			a = line_placements_w[i][ii]/D
			detail_w_mod = (detail_w-a)/(1+detail_w*np.conj(a))
			#detail_w*=D
			stereo_dist = np.amin(np.absolute(detail_w_mod))
			extra_tlu_M_list[i][ii] = 1/(np.absolute(stereo_dist)**2+1)


		if i > 0:
			if np.amin(np.absolute(w_center_array[:i]-w_right_array[i]))<0.01*np.absolute(w_center_array[i]-w_right_array[i]):
				print(np.amin(np.absolute(w_center_array[:i]-w_right_array[i])))
				line_placements.append((np.hstack(line_placements[np.argmin(np.absolute(w_center_array[:i]-w_right_array[i]))])[0]+line_posers[i]*(np.hstack(line_placements_raw[i])-np.hstack(line_placements_raw[i])[2])+\
					np.hstack(line_placements[previous_index[i]])[0]+line_posers[i]*(np.hstack(line_placements_raw[i])-np.hstack(line_placements_raw[i])[1]))/2)
			else:
				line_placements.append(np.hstack(line_placements[previous_index[i]])[0]+line_posers[i]*(np.hstack(line_placements_raw[i])-np.hstack(line_placements_raw[i])[1]))
		centers_M.append(['',line_placements[len(line_placements)-1][0],line_placements[len(line_placements)-1][1],line_posers[i],0,0,0,0,0,0,0,line_placements[len(line_placements)-1][2],line_placements[len(line_placements)-1][3]])

	return centers_M,associates,tlu_M_list,line_placements,line_placements_w,extra_tlu_M_list#,line_placements_W_inv,line_placements_M_inv

'''
Create correspondence between:
string representation of heptagon vertex
position in W
position of its anchor in W (backward in tree), already known, store as A_W
whether this corresponds to the A (0), R (1), or L (2) of its anchor
position of its anchor in D calculated with A_D = evalinv(f,A_W)
position of its right branch in D
position of its left branch in D
w, z, and beta
'''

from_start = True


preceeder = f'{sys.argv[1]}_'

pickle_off = open(f"pickle/bubble_wrap/{preceeder}full_wrapping.txt", "rb")
loaded_pickle = pickle.load(pickle_off)
full_wrapping = np.array(loaded_pickle[0])[::-10]

nemo_point = loaded_pickle[1]-pi

opposing_point = np.conj(nemo_point)+pi
opposing_point = np.conj(nemo_point)
k = D/(1+sin(opposing_point.imag)*np.sin(np.imag(full_wrapping))+cos(opposing_point.imag)*np.cos(np.imag(full_wrapping))*np.cos(np.real(full_wrapping)-opposing_point.real))
x = k*np.cos(np.imag(full_wrapping))*np.sin(np.real(full_wrapping)-opposing_point.real)
y = k*(cos(opposing_point.imag)*np.sin(np.imag(full_wrapping))-sin(opposing_point.imag)*np.cos(np.imag(full_wrapping))*np.cos(np.real(full_wrapping)-opposing_point.real))
centers_pickled = 0
if len(sys.argv)>2:
	centers_pickled = int(sys.argv[2])
if from_start:
	poly_use = (x+y*1j)

	w_raw,beta_raw,z_raw,centers = explore_polygon(poly_use,preceeder,centers_pickled=centers_pickled,inverted_buffer_metric=False)

	centers_M,associates,tlu_M_list,line_placements_m,line_placements_w,extra_tlu_M_list = M_centers_from_W(centers,w_raw,beta_raw,z_raw)
	with open(f'pickle/matlab_saves/{preceeder}centers_M.txt', 'wb') as fh:
		 pickle.dump(centers_M, fh)
	with open(f'pickle/matlab_saves/{preceeder}associates.txt', 'wb') as fh:
		 pickle.dump(associates, fh)
	with open(f'pickle/matlab_saves/{preceeder}tlu_M_list.txt', 'wb') as fh:
		 pickle.dump(tlu_M_list, fh)
	with open(f'pickle/matlab_saves/{preceeder}line_placements_m.txt', 'wb') as fh:
		 pickle.dump(line_placements_m, fh)
	with open(f'pickle/matlab_saves/{preceeder}line_placements_w.txt', 'wb') as fh:
		 pickle.dump(line_placements_w, fh)
	with open(f'pickle/matlab_saves/{preceeder}extra_tlu_M_list.txt', 'wb') as fh:
		 pickle.dump(extra_tlu_M_list, fh)
else:
	poly_use = (x+y*1j)

	w_raw,beta_raw,z_raw,centers = explore_polygon(poly_use,preceeder,centers_pickled=centers_pickled,inverted_buffer_metric=False)

	pickle_off = open(f"pickle/matlab_saves/{preceeder}centers_M.txt", "rb")
	centers_M = pickle.load(pickle_off)
	pickle_off = open(f"pickle/matlab_saves/{preceeder}associates.txt", "rb")
	associates = pickle.load(pickle_off)
	pickle_off = open(f"pickle/matlab_saves/{preceeder}tlu_M_list.txt", "rb")
	tlu_M_list = pickle.load(pickle_off)
	pickle_off = open(f"pickle/matlab_saves/{preceeder}line_placements_m.txt", "rb")
	line_placements_m = pickle.load(pickle_off)
	pickle_off = open(f"pickle/matlab_saves/{preceeder}line_placements_w.txt", "rb")
	line_placements_w = pickle.load(pickle_off)
	pickle_off = open(f"pickle/matlab_saves/{preceeder}extra_tlu_M_list.txt", "rb")
	extra_tlu_M_list = pickle.load(pickle_off)

fig = plt.figure()
ax = fig.add_subplot(111)


power = 8
w_center_array = np.hstack(np.array([cent[1] for cent in centers]))
w_anchor_array = np.hstack(np.array([cent[2] for cent in centers]))
m_center_array = np.hstack(np.array([cent[1] for cent in centers_M]))
m_anchor_array = np.hstack(np.array([cent[2] for cent in centers_M]))
m_posers_array = np.hstack(np.array([cent[3] for cent in centers_M]))


nearest_to_W_vert = np.argmin(np.absolute(w_center_array[:,None]-w_raw[None,:]),axis=0)

M_vertices_est = []
ii = 0
for i in range(nearest_to_W_vert.shape[0]):
	if np.absolute(w_raw[i]-poly_use[ii])<0.00000001:
		center = centers[nearest_to_W_vert[i]]
		p = prepare_coefficients(center[7],center[8],center[9],power)
		M_vals = sc_map_continuous(np.array([center[8][i]]),center[8],center[9],p)
		M_vertices_est.append(m_center_array[nearest_to_W_vert[i]]+m_posers_array[nearest_to_W_vert[i]]*M_vals[0])
		ii+=1

np.save(f'pickle/interpolation_points/{preceeder}W.npy',np.hstack([w_center_array,poly_use,np.ndarray.flatten(np.array(line_placements_w)[:,4:])]))
np.save(f'pickle/interpolation_points/{preceeder}M.npy',np.hstack([m_center_array,M_vertices_est,np.ndarray.flatten(np.array(line_placements_m)[:,4:])]))
np.save(f'pickle/interpolation_points/{preceeder}tlu_M.npy',np.hstack([np.array(extra_tlu_M_list)[:,0],np.ones_like(M_vertices_est),np.ndarray.flatten(np.array(extra_tlu_M_list)[:,4:])]))