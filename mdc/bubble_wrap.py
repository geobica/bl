import numpy as np
import csv
import sys
import shapely.wkt
import shapely.ops
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import minimum_spanning_tree
import matplotlib.pyplot as plt
import math
from math import pi,e,sin,cos,asin,acos,atan2
import scipy.linalg
from scipy.spatial.transform import Rotation as R
from shapely.geometry import Polygon, Point
from shapely.geometry import shape
import pickle
import os
import sys
import fiona
import geopandas as gpd
from tqdm import tqdm

import auto_center

WRAP_BUFFER = 0.005
BRIDGE_BUFFER = 0.01
SIMPLIFY = 2e-3
PRE_SIMPLIFY = 1e-3
UNION_CHUNK = 4000

def polygons_of(geom):
	if geom.geom_type=='Polygon':
		return [geom]
	if geom.geom_type=='MultiPolygon':
		return list(geom.geoms)
	return [g for g in geom.geoms if g.geom_type=='Polygon']

def rotation_matrix(axis,angle):
	axis /= np.linalg.norm(axis)
	return np.array([[np.cos(angle)+axis[0]**2*(1-np.cos(angle)), axis[0]*axis[1]*(1-np.cos(angle))-axis[2]*np.sin(angle), axis[0]*axis[2]*(1-np.cos(angle))+axis[1]*np.sin(angle)],
						[axis[0]*axis[1]*(1-np.cos(angle))+axis[2]*np.sin(angle), np.cos(angle)+axis[1]**2*(1-np.cos(angle)), axis[1]*axis[2]*(1-np.cos(angle))-axis[0]*np.sin(angle)],
						[axis[0]*axis[2]*(1-np.cos(angle))-axis[1]*np.sin(angle), axis[1]*axis[2]*(1-np.cos(angle))+axis[0]*np.sin(angle), np.cos(angle)+axis[2]**2*(1-np.cos(angle))]])

def np_equi_to_coord(equi_com):
	return np.array([np.cos(np.imag(equi_com))*np.cos(np.real(equi_com)),np.cos(np.imag(equi_com))*np.sin(np.real(equi_com)),np.sin(np.imag(equi_com))]).T
def equi_to_coord(equi_com):
	return np.array([cos(equi_com.imag)*cos(equi_com.real),cos(equi_com.imag)*sin(equi_com.real),sin(equi_com.imag)])
def coord_to_equi(coord_vec):
	return atan2(coord_vec[1],coord_vec[0])+asin(coord_vec[2])*1j
def np_coord_to_equi(coord_vec):
	return np.arctan2(coord_vec[:,1],coord_vec[:,0])+np.arcsin(coord_vec[:,2])*1j

def equi_to_stereo(equi,center):
	#equi *= pi/180
	k = 2/(1+sin(center.imag)*np.sin(np.imag(equi))+cos(center.imag)*np.cos(np.imag(equi))*np.cos(np.real(equi)-center.real))
	x = k*np.cos(np.imag(equi))*np.sin(np.real(equi)-center.real)
	y = k*(cos(center.imag)*np.sin(np.imag(equi))-sin(center.imag)*np.cos(np.imag(equi))*np.cos(np.real(equi)-center.real))

	x[np.isnan(x)] = 0
	y[np.isnan(y)] = 0
	x = np.minimum(x,10000000000000)
	x = np.maximum(x,-10000000000000)
	y = np.minimum(y,10000000000000)
	y = np.maximum(y,-10000000000000)

	return (x+y*1j)

def stereo_to_equi(stereo,center):
	x = np.real(stereo)
	y = np.imag(stereo)
	lam_0 = center.real
	phi_0 = center.imag
	rho = np.absolute(stereo)
	c = 2*np.arctan(rho/2)
	phi = np.arcsin(np.cos(c)*math.sin(phi_0)+y*np.sin(c)*math.cos(phi_0)/rho)
	lam = lam_0+np.arctan2(x*np.sin(c),(rho*math.cos(phi_0)*np.cos(c)-y*math.sin(phi_0)*np.sin(c)))
	equi_com = lam+1j*phi
	return equi_com

def haversine(cent,arr):
	dlon = np.real(arr - cent)
	dlat = np.imag(arr - cent)
	a = np.sin(dlat/2)**2 + np.cos(np.imag(arr)) * np.cos(np.imag(cent)) * np.sin(dlon/2)**2
	return np.arcsin(np.sqrt(a))

def haversine_arr(arr_0,arr_1):
	if len(arr_0)>1000:
		arr_0 = arr_0[::int(len(arr_0)/1000)]
	if len(arr_1)>1000:
		arr_1 = arr_1[::int(len(arr_1)/1000)]
	dlon = np.real(arr_1[:,None] - arr_0[None,:])
	dlat = np.imag(arr_1[:,None] - arr_0[None,:])
	a = np.sin(dlat/2)**2 + np.cos(np.imag(arr_1[:,None])) * np.cos(np.imag(arr_0[None,:])) * np.sin(dlon/2)**2
	return np.arcsin(np.sqrt(a))

def haversine_2p(p_0,p_1):
	return 2*np.arcsin(np.sqrt(np.sin((np.imag(p_1)-np.imag(p_0))/2)**2+np.cos(np.imag(p_0))*np.cos(np.imag(p_1))*np.sin((np.real(p_1)-np.real(p_0))/2)**2))

def nearest_pair(arr_0,arr_1):
	tree = cKDTree(np_equi_to_coord(arr_0))
	chord,idx = tree.query(np_equi_to_coord(arr_1))
	j = int(np.argmin(chord))
	return arr_0[idx[j]],arr_1[j],2*np.arcsin(min(1,chord[j]/2))

def dash_between(arr_0,arr_1,step,radius):
	point_0,point_1,gap = nearest_pair(arr_0,arr_1)
	direction = direction_between(point_0,point_1)
	dash_instances = 2*np.linspace(-radius,haversine(point_0,point_1)+radius,num=int((haversine(point_0,point_1)+2*radius)/step),endpoint=True)
	dashes = []
	for i in range(dash_instances.shape[0]):
		dashes.append(move(point_0,direction,dash_instances[i]))
	return np.array(dashes)

def direction_between(start,end):
	a = equi_to_coord(start)
	b = equi_to_coord(end)
	v = np.cross(a, b)
	c = np.dot(a, b)
	s = np.linalg.norm(v)
	kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
	rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
	angle = np.arccos(np.trace(rotation_matrix)/2-.5)
	return scipy.linalg.fractional_matrix_power(rotation_matrix,1/angle)

def move(start,direction,angle):
	return coord_to_equi(np.matmul(scipy.linalg.fractional_matrix_power(direction,angle),equi_to_coord(start)))

def buffer_3d_full_poly(v_0,v_1,dist,segments,nemo_point):
	u_0 = np_equi_to_coord(v_0)
	u_1 = np_equi_to_coord(v_1)
	full_lengths = haversine_2p(v_0,v_1)
	steps_along = (full_lengths/(dist*2*pi/segments)+1).astype('int')

	segment_polys = []
	for i in range(steps_along.shape[0]):
		portions = [[],[],[],[]]
		all_portions = []
		axis = np.cross(u_0[i],u_1[i])
		for j in range(steps_along[i]):
			along_vecs = np.matmul(rotation_matrix(axis,full_lengths[i]*j/steps_along[i]),u_0[i])
			away = np.cross(along_vecs,axis)
			fixer = rotation_matrix(away,dist)
			portions[1].append(np.matmul(fixer,along_vecs))
			fixer_inv = rotation_matrix(away,-dist)
			portions[3].append(np.matmul(fixer_inv,along_vecs))
		away = np.cross(u_0[i],axis)
		fixer = rotation_matrix(away,dist)
		beside = np.matmul(fixer,u_0[i])
		for j in range(segments):
			stirrer = rotation_matrix(u_0[i],pi/segments*j)
			aroundy = np.matmul(stirrer,beside)
			portions[0].append(aroundy)
		away = np.cross(u_1[i],axis)
		fixer = rotation_matrix(away,dist)
		beside = np.matmul(fixer,u_1[i])
		for j in range(segments):
			stirrer = rotation_matrix(u_1[i],pi+pi/segments*j)
			aroundy = np.matmul(stirrer,beside)
			portions[2].append(aroundy)

		all_portions = portions[0][::-1]+portions[1]+portions[2][::-1]+portions[3][::-1]
		stereo = equi_to_stereo(np_coord_to_equi(np.array(all_portions)),np.conj(nemo_point)+pi)

		segment_polys.append(Polygon(np.array([np.real(stereo),np.imag(stereo)]).T))
	return shapely.ops.unary_union(segment_polys)

def buffer_3d(v_0,v_1,dist,segments,nemo_point):
	full_poly = buffer_3d_full_poly(v_0,v_1,dist,segments,nemo_point)
	x,y = max(polygons_of(full_poly),key=lambda part:part.area).exterior.coords.xy
	equi_poly = stereo_to_equi(np.array(x)+1j*np.array(y),np.conj(nemo_point)+pi)

	return equi_poly

def comp_to_poly(comp):
	return Polygon(np.array([np.real(comp),np.imag(comp)]).T)

def poly_to_comp(poly):
	x,y = poly.exterior.coords.xy
	return np.array(x)+1j*np.array(y)

def reduce_line(line_st,dist_needed=0.005):
	included_points = [line_st[0]]
	for i in range(line_st.shape[0]-1):
		if haversine_2p(included_points[-1],line_st[i])>=dist_needed:
			included_points.append(line_st[i])
	included_points.append(line_st[-1])
	included_points = np.array(included_points)
	return included_points

def lonlat_to_nemo(lon_in,lat_in):
	# the point opposite the given input point, which should be outside the region
	return (float((lon_in+360)%360-180)-lat_in*1j)*pi/180

def ring_spanning_tree(rings):
	n = len(rings)
	distances = np.zeros((n,n))
	for i in tqdm(range(n-1)):
		for j in range(i+1,n):
			_,_,gap = nearest_pair(rings[i],rings[j])
			distances[i,j] = gap
			distances[j,i] = gap
	return minimum_spanning_tree(distances).tocoo()

def wrap_lon(comp):
	return np.remainder(np.real(comp)*180/pi+180,360)-180

def run_bubble_wrap(name,lon_in=None,lat_in=None,input_path='input_sample',pickle_path='pickle/bubble_wrap'):
	if lon_in is None or lat_in is None:
		gpkg = os.path.join(f'{input_path}',f'{name}.gpkg')
		with fiona.open(gpkg, layer=fiona.listlayers(gpkg)[0]) as layer:
			geoms = [shape(feature['geometry']).buffer(0) for feature in layer]
		region = shapely.ops.unary_union(geoms)
		lon_in,lat_in,score = auto_center.find_center(region)
		print(f'auto-selected center: lon={lon_in:.3f} lat={lat_in:.3f} (min boundary distance {score:.3f})')

	print("Starting bubble_wrap routine...")

	nemo_point = lonlat_to_nemo(lon_in,lat_in)

	stereo_poly_combine = None
	stripped_equis = []

	gpkg = os.path.join(f'{input_path}',f'{name}.gpkg')
	with fiona.open(gpkg, layer=fiona.listlayers(gpkg)[0]) as layer:
		pieces = []
		pre_verts_before = pre_verts_after = 0
		for feature_i,feature in enumerate(layer):
			print(f"Reprojecting every polygon in feature {feature_i}/{len(layer)}:")
			for poly in tqdm(polygons_of(shape(feature['geometry']))):
				pre_verts_before += len(poly.exterior.coords)
				if PRE_SIMPLIFY>0:
					poly = poly.simplify(PRE_SIMPLIFY,preserve_topology=True)
				for simple_poly in polygons_of(poly):
					pre_verts_after += len(simple_poly.exterior.coords)
					pieces.append(comp_to_poly(equi_to_stereo(poly_to_comp(simple_poly)*pi/180,np.conj(nemo_point)+pi)).buffer(WRAP_BUFFER))
		if PRE_SIMPLIFY>0:
			print(f'pre-simplify {PRE_SIMPLIFY}deg: {pre_verts_before} -> {pre_verts_after} vertices ({len(pieces)} parts)')
		stereo_poly_combine = shapely.ops.unary_union([shapely.ops.unary_union(pieces[i:i+UNION_CHUNK]) for i in range(0,len(pieces),UNION_CHUNK)])
	if SIMPLIFY>0:
		before = sum(len(p.exterior.coords) for p in polygons_of(stereo_poly_combine))
		stereo_poly_combine = stereo_poly_combine.simplify(SIMPLIFY,preserve_topology=True)
		after = sum(len(p.exterior.coords) for p in polygons_of(stereo_poly_combine))
		print(f'simplify {SIMPLIFY}: {before} -> {after} vertices ({len(polygons_of(stereo_poly_combine))} parts)')
	for poly in polygons_of(stereo_poly_combine):
		x,y = poly.exterior.coords.xy
		stripped_equis.append(stereo_to_equi(np.array(x)+1j*np.array(y),np.conj(nemo_point)+pi))

	step = 0.001
	radius = 0.005

	print(f'{len(stripped_equis)} rings to join')
	included_lines = list(stripped_equis)
	bridge_log = []
	if len(stripped_equis)>1:
		mst = ring_spanning_tree(stripped_equis)
		for i,j in tqdm(list(zip(mst.row,mst.col))):
			point_0,point_1,gap = nearest_pair(stripped_equis[i],stripped_equis[j])
			bridge_log.append(f'{i}\t{j}\t{np.real(point_0)*180/pi:.5f},{np.imag(point_0)*180/pi:.5f}'
				f'\t{np.real(point_1)*180/pi:.5f},{np.imag(point_1)*180/pi:.5f}\t{gap*180/pi:.5f}')
			included_lines.append(stereo_to_equi(poly_to_comp(comp_to_poly(equi_to_stereo(dash_between(stripped_equis[i],stripped_equis[j],step,radius),np.conj(nemo_point)+pi)).buffer(BRIDGE_BUFFER)),np.conj(nemo_point)+pi))
	with open(os.path.join(f'{pickle_path}',f'{name}_bridges.txt'),'w') as fh:
		fh.write('ring_a\tring_b\tfrom_lon,lat\tto_lon,lat\tdegrees\n')
		fh.write('\n'.join(bridge_log))

	v_0 = []
	v_1 = []
	for st in tqdm(included_lines):
		st_red = reduce_line(st)
		v_0 = v_0+list(st_red[:-1])
		v_1 = v_1+list(st_red[1:])
	v_0 = np.array(v_0)
	v_1 = np.array(v_1)
	double_buffed = buffer_3d(v_0,v_1,WRAP_BUFFER*2,10,nemo_point)
	double_buffed = reduce_line(double_buffed)
	third_poly = buffer_3d_full_poly(double_buffed[:-1],double_buffed[1:],WRAP_BUFFER,10,nemo_point)

	final_stereo_poly = comp_to_poly(equi_to_stereo(double_buffed,np.conj(nemo_point)+pi)).buffer(0).difference(third_poly.buffer(0))
	parts = polygons_of(final_stereo_poly)
	biggest_land = comp_to_poly(equi_to_stereo(max(stripped_equis,key=len),np.conj(nemo_point)+pi))
	anchor = biggest_land.buffer(0).representative_point()
	holding = [part for part in parts if part.contains(anchor)]
	if len(parts)>1:
		print(f'difference left {len(parts)} pieces; keeping the largest')
	poly = holding[0] if holding else max(parts,key=lambda part:part.area)
	for hole in poly.interiors:
		print(f'enclosed ocean: ring of {len(hole.coords)} points, area {Polygon(hole).area:.5g}')
	x,y = poly.exterior.coords.xy
	comp_xy = stereo_to_equi((np.array(x)+1j*np.array(y)),np.conj(nemo_point)+pi)

	included_points = [comp_xy[0]]
	dist_needed = 0.002
	for i in range(comp_xy.shape[0]):
		if haversine_2p(included_points[-1],comp_xy[i])>=dist_needed:
			included_points.append(comp_xy[i])
	included_points = np.array(included_points)

	with open(os.path.join(f'{pickle_path}',f'{name}_full_wrapping.txt'), 'wb') as fh:
		pickle.dump([included_points,nemo_point], fh)

	fig,ax = plt.subplots(figsize=(16,8))
	for ring in stripped_equis:
		ax.plot(wrap_lon(ring),np.imag(ring)*180/pi,'.',color='#999999',markersize=0.5)
	ax.plot(wrap_lon(included_points),np.imag(included_points)*180/pi,'.',color='#cc5511',markersize=1.2)
	ax.set_xlim(-180,180)
	ax.set_ylim(-90,90)
	ax.set_aspect('equal')
	fig.savefig(os.path.join(f'{pickle_path}',f'{name}_boundary.png'),dpi=140,bbox_inches='tight')

	return included_points,nemo_point

if __name__=='__main__':
	run_bubble_wrap(sys.argv[1])