import numpy as np
import math
import scipy as sp
from scipy.special import gamma, factorial
from scipy import io
from scipy.integrate import solve_ivp
import fiona
from shapely.geometry import Polygon, MultiPolygon, Point
import sys
import shapely.wkt
import random
import scipy.interpolate
from shapely.geometry import mapping, Polygon, LineString
import fiona
import matplotlib.pyplot as plt
import cplot
import pickle
from math import pi,e,sin,cos,tan
import os
from shapely.geometry import shape


grid_sample_w = np.load(f'pickle/interpolation_points/{sys.argv[1]}_W.npy')/12756274#*12756274
dx = 0.00001
slightly_to_the_right = grid_sample_w+dx
mapped_sample_points = np.load(f'pickle/interpolation_points/{sys.argv[1]}_M.npy')/12756274


pickle_off = open(f"pickle/bubble_wrap/{sys.argv[1]}_full_wrapping.txt", "rb")
loaded_pickle = pickle.load(pickle_off)
nemo_point = loaded_pickle[1]
opposing_point = pi+np.conj(nemo_point)

outline_point_location = [int((grid_sample_w.shape[0]-math.ceil(loaded_pickle[0].shape[0]/10))/13),int((grid_sample_w.shape[0]-math.ceil(loaded_pickle[0].shape[0]/10))/13+math.ceil(loaded_pickle[0].shape[0]/10))]

grid_sample_w = np.concatenate((grid_sample_w[:outline_point_location[0]],grid_sample_w[outline_point_location[1]:]))
mapped_sample_points = np.concatenate((mapped_sample_points[:outline_point_location[0]],mapped_sample_points[outline_point_location[1]:]))

# real_rbf = sp.interpolate.LinearNDInterpolator(list(zip(np.real(grid_sample_w),np.imag(grid_sample_w))),np.real(mapped_sample_points))
# imag_rbf = sp.interpolate.LinearNDInterpolator(list(zip(np.real(grid_sample_w),np.imag(grid_sample_w))),np.imag(mapped_sample_points))
#needed for normal vv
real_rbf = sp.interpolate.CloughTocher2DInterpolator(list(zip(np.real(grid_sample_w),np.imag(grid_sample_w))),np.real(mapped_sample_points),tol=0.1)
imag_rbf = sp.interpolate.CloughTocher2DInterpolator(list(zip(np.real(grid_sample_w),np.imag(grid_sample_w))),np.imag(mapped_sample_points),tol=0.1)
#needed for normal ^^
# real_clough = sp.interpolate.CloughTocher2DInterpolator(list(zip(np.real(grid_sample_w),np.imag(grid_sample_w))),np.real(mapped_sample_points))
# imag_clough = sp.interpolate.CloughTocher2DInterpolator(list(zip(np.real(grid_sample_w),np.imag(grid_sample_w))),np.imag(mapped_sample_points))
#needed for normal vv
real_nearest = sp.interpolate.NearestNDInterpolator(list(zip(np.real(grid_sample_w),np.imag(grid_sample_w))),np.real(mapped_sample_points))
imag_nearest = sp.interpolate.NearestNDInterpolator(list(zip(np.real(grid_sample_w),np.imag(grid_sample_w))),np.imag(mapped_sample_points))
#needed for normal ^^

def equi_to_stereo(equi):
	equi *= pi/180
	k = 12756274/(1+sin(opposing_point.imag)*np.sin(np.imag(equi))+cos(opposing_point.imag)*np.cos(np.imag(equi))*np.cos(np.real(equi)-opposing_point.real))
	x = k*np.cos(np.imag(equi))*np.sin(np.real(equi)-opposing_point.real)
	y = k*(cos(opposing_point.imag)*np.sin(np.imag(equi))-sin(opposing_point.imag)*np.cos(np.imag(equi))*np.cos(np.real(equi)-opposing_point.real))

	x[np.isnan(x)] = 0
	y[np.isnan(y)] = 0
	x = np.minimum(x,10000000000000)
	x = np.maximum(x,-10000000000000)
	y = np.minimum(y,10000000000000)
	y = np.maximum(y,-10000000000000)

	return (x+y*1j)

def reproj(equi,stereoproj=True):
	if stereoproj:
		stereo = equi_to_stereo(equi)
	else:
		stereo = equi

	# #return stereo
	out_real = real_rbf(np.real(stereo/12756274),np.imag(stereo/12756274))
	out_imag = imag_rbf(np.real(stereo/12756274),np.imag(stereo/12756274))
	out_real[np.isnan(out_real)] = real_nearest(np.real(stereo[np.isnan(out_real)]/12756274),np.imag(stereo[np.isnan(out_real)]/12756274))
	out_imag[np.isnan(out_imag)] = imag_nearest(np.real(stereo[np.isnan(out_imag)]/12756274),np.imag(stereo[np.isnan(out_imag)]/12756274))
	return out_real+1j*out_imag

def multi_to_file(filename,multi):
	schema = {
			'geometry': 'Polygon',
	}
	with fiona.open(filename, 'w', 'ESRI Shapefile', schema) as c:
			for poly in list(multi):
				c.write({
						'geometry': mapping(MultiPolygon(poly)),
				})

def lines_to_file(filename,lines):
	schema = {
			'geometry': 'LineString',
	}
	with fiona.open(filename, 'w', 'ESRI Shapefile', schema) as c:
			for line in list(lines):
				c.write({
						'geometry': mapping(LineString(line)),
				})


##outline files
polyarr = np.zeros((loaded_pickle[0].shape[0],2))
polyarr[:,0] = np.real(loaded_pickle[0])
polyarr[:,1] = np.imag(loaded_pickle[0])
polyarr*=180/pi
poly = Polygon(polyarr)
poly_ext_out = np.zeros_like(np.array(poly.exterior.coords))
poly_ext = equi_to_stereo(np.array(poly.exterior.coords)[:,0]+1j*np.array(poly.exterior.coords)[:,1])
poly_ext_out[:,0] = np.real(poly_ext)
poly_ext_out[:,1] = np.imag(poly_ext)
multipolys = [[[poly_ext_out,[]]]]
original_outline = Polygon(poly_ext_out)
multi_to_file(f'maps_projected/original_outline_{sys.argv[1]}.shp',multipolys)
poly_ext = reproj(np.array(poly.exterior.coords)[:,0]+1j*np.array(poly.exterior.coords)[:,1])
poly_ext_out[:,0] = np.real(poly_ext)
poly_ext_out[:,1] = np.imag(poly_ext)
multipolys = [[[poly_ext_out,[]]]]
multi_to_file(f'maps_projected/outline_{sys.argv[1]}.shp',multipolys)


##gridlines files
gridlines_equi_long = np.arange(-180,180,5)[:,None]+1j*np.arange(-90,90.1,0.1)[None,:]
gridlines_equi_lat = (np.arange(-180,180.1,0.1)[:,None]+1j*np.arange(-85,90,5)[None,:]).T
line_list = []
for gridline in gridlines_equi_long:
	line_comp = equi_to_stereo(gridline)
	line_split = np.zeros((line_comp.shape[0],2))
	line_split[:,0] = np.real(line_comp)
	line_split[:,1] = np.imag(line_comp)
	line_list.append(line_split)
for gridline in gridlines_equi_lat:
	line_comp = equi_to_stereo(gridline)
	line_split = np.zeros((line_comp.shape[0],2))
	line_split[:,0] = np.real(line_comp)
	line_split[:,1] = np.imag(line_comp)
	line_list.append(line_split)
lines_to_file(f'maps_projected/original_gridline_{sys.argv[1]}.shp',line_list)
gridlines_equi_long = np.arange(-180,180,5)[:,None]+1j*np.arange(-90,90.1,0.1)[None,:]
gridlines_equi_lat = (np.arange(-180,180.1,0.1)[:,None]+1j*np.arange(-85,90,5)[None,:]).T
line_list = []
from shapely.validation import make_valid
for gridline in gridlines_equi_long:
	line_comp = equi_to_stereo(gridline)
	line_split = np.zeros((line_comp.shape[0],2))
	line_split[:,0] = np.real(line_comp)
	line_split[:,1] = np.imag(line_comp)
	intersected = LineString(line_split).intersection(make_valid(original_outline))
	if intersected.type=='LineString':
		if np.array(intersected.coords).shape[0]>0:
			reprojed = reproj(np.array(intersected.coords)[:,0]+1j*np.array(intersected.coords)[:,1],False)
			reprojed_split = np.zeros((reprojed.shape[0],2))
			reprojed_split[:,0] = np.real(reprojed)
			reprojed_split[:,1] = np.imag(reprojed)
			line_list.append(reprojed_split)
	elif intersected.type=='MultiLineString':
		for linest in LineString(line_split).intersection(make_valid(original_outline)):
			reprojed = reproj(np.array(linest.coords)[:,0]+1j*np.array(linest.coords)[:,1],False)
			reprojed_split = np.zeros((reprojed.shape[0],2))
			reprojed_split[:,0] = np.real(reprojed)
			reprojed_split[:,1] = np.imag(reprojed)
			line_list.append(reprojed_split)

for gridline in gridlines_equi_lat:
	line_comp = equi_to_stereo(gridline)
	line_split = np.zeros((line_comp.shape[0],2))
	line_split[:,0] = np.real(line_comp)
	line_split[:,1] = np.imag(line_comp)
	intersected = LineString(line_split).intersection(make_valid(original_outline))
	if intersected.type=='LineString':
		if np.array(intersected.coords).shape[0]>0:
			reprojed = reproj(np.array(intersected.coords)[:,0]+1j*np.array(intersected.coords)[:,1],False)
			reprojed_split = np.zeros((reprojed.shape[0],2))
			reprojed_split[:,0] = np.real(reprojed)
			reprojed_split[:,1] = np.imag(reprojed)
			line_list.append(reprojed_split)
	elif intersected.type=='MultiLineString':
		for linest in LineString(line_split).intersection(make_valid(original_outline)):
			reprojed = reproj(np.array(linest.coords)[:,0]+1j*np.array(linest.coords)[:,1],False)
			reprojed_split = np.zeros((reprojed.shape[0],2))
			reprojed_split[:,0] = np.real(reprojed)
			reprojed_split[:,1] = np.imag(reprojed)
			line_list.append(reprojed_split)
lines_to_file(f'maps_projected/gridline_{sys.argv[1]}.shp',line_list)


def read_in_file(filename,transform):
	multipolys = []

	all_points = []
	with fiona.open(filename, layer=fiona.listlayers(filename)[0]) as layer:
		for feature in layer:
			multipolys.append([])
			if shape(feature['geometry']).type=='MultiPolygon':
				for poly in shape(feature['geometry']):
					ext_coords = poly.exterior.coords
					poly_ext_out = np.zeros_like(np.array(poly.exterior.coords))
					poly_int_out = [np.zeros_like(np.array(b.coords)) for b in list(poly.interiors)]
					if len(ext_coords)>0:
						poly_ext = transform(np.array(ext_coords)[:,0]+1j*np.array(ext_coords)[:,1])

						poly_ext_out[:,0] = np.real(poly_ext)
						poly_ext_out[:,1] = np.imag(poly_ext)
						poly_int = [transform(np.array(b.coords)[:,0]+1j*np.array(b.coords)[:,1]) for b in list(poly.interiors)]
						for i in range(len(poly_int)):
							poly_int_out[i][:,0] = np.real(poly_int[i])
							poly_int_out[i][:,1] = np.imag(poly_int[i])
						multipolys[-1].append([poly_ext_out,poly_int_out])
			else:
				poly = shape(feature['geometry'])
				ext_coords = poly.exterior.coords
				poly_ext_out = np.zeros_like(np.array(poly.exterior.coords))
				poly_int_out = [np.zeros_like(np.array(b.coords)) for b in list(poly.interiors)]
				if len(ext_coords)>0:
					poly_ext = transform(np.array(ext_coords)[:,0]+1j*np.array(ext_coords)[:,1])

					poly_ext_out[:,0] = np.real(poly_ext)
					poly_ext_out[:,1] = np.imag(poly_ext)
					poly_int = [transform(np.array(b.coords)[:,0]+1j*np.array(b.coords)[:,1]) for b in list(poly.interiors)]
					for i in range(len(poly_int)):
						poly_int_out[i][:,0] = np.real(poly_int[i])
						poly_int_out[i][:,1] = np.imag(poly_int[i])
					multipolys[-1].append([poly_ext_out,poly_int_out])
	return multipolys


multipolys = read_in_file(sys.argv[2],equi_to_stereo)
multi_to_file(f'maps_projected/original_{sys.argv[1]}.shp',multipolys)

multipolys = read_in_file(sys.argv[2],reproj)
multi_to_file(f'maps_projected/{sys.argv[1]}.shp',multipolys)