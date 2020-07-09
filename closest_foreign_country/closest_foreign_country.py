from PIL import Image, ImageDraw, ImageFont, ImageColor
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import gc
from matplotlib import pyplot as plt
import cv2
import os
import PIL
from PIL import Image
import random
import time
from multiprocessing import Process

def boundries_of(array_in,search_for): # returns [top,bottom,left,right]
	return [np.min(array_in.shape[0]*(1-np.max(array_in==search_for,axis=0))+np.argmax(array_in==search_for,axis=0))
	,array_in.shape[0]-np.min(array_in.shape[0]*(1-np.max(np.flip(array_in,0)==search_for,axis=0))+np.argmax(np.flip(array_in,0)==search_for,axis=0))
	,np.min(array_in.shape[1]*(1-np.max(array_in==search_for,axis=1))+np.argmax(array_in==search_for,axis=1))
	,array_in.shape[1]-np.min(array_in.shape[1]*(1-np.max(np.flip(array_in,1)==search_for,axis=1))+np.argmax(np.flip(array_in,1)==search_for,axis=1))]

code_map = cv2.imread('code_map.png')
code_map[:,:,0] = code_map[:,:,2]
size = code_map.shape[0]/2

output_map = np.zeros(code_map.shape, dtype=np.uint32)

grid_size = 1
flags_to_place = []
flags_to_place_codes = []
search_scale = 4

simple_code_map = code_map[:,:,0].astype(np.uint32)*50+code_map[:,:,1].astype(np.uint32)
code_map_small = np.zeros((int(3600/search_scale),int(7200/search_scale),3))
code_map = cv2.resize(code_map,(7200,3600),interpolation = cv2.INTER_NEAREST)
closest_codes = np.zeros((int(code_map.shape[0]/grid_size),int(code_map.shape[1]/grid_size)))

all_lat = np.tile(np.reshape(np.arange(code_map_small.shape[0]),(code_map_small.shape[0],1)),(1,code_map_small.shape[1]))*math.pi/code_map_small.shape[0]
all_long = np.tile(np.reshape(np.arange(code_map_small.shape[1]),(1,code_map_small.shape[1])),(code_map_small.shape[0],1))*math.pi*2/code_map_small.shape[1]
all_x = np.sin(all_lat)*np.cos(all_long)
all_y = np.sin(all_lat)*np.sin(all_long)
all_z = np.cos(all_lat)

country_list = []
bool_map_list = []
code_map_small_list = []
bool_map = (code_map[:,:,0] != 0)

def run_latitude(latitude_pixel_base):
	for latitude_pixel in range(latitude_pixel_base,latitude_pixel_base+100,grid_size):
		if latitude_pixel%100 == 0:
			country_list = []
			bool_map_list = []
			code_map_small_list = []
		print(latitude_pixel)
		start_time = time.time()
		here_z = math.cos(latitude_pixel*math.pi/code_map.shape[0])
		long_z = (all_z-here_z)**2
		for longitude_pixel in range(0,code_map.shape[1],grid_size):
			original_countries = []
			if code_map[latitude_pixel,longitude_pixel,0] != 0:
				starting_country = code_map[latitude_pixel,longitude_pixel]
				if (50*starting_country[0]+starting_country[1]) in country_list:
					bool_map = bool_map_list[country_list.index((50*starting_country[0]+starting_country[1]))]
					code_map_small = code_map_small_list[country_list.index((50*starting_country[0]+starting_country[1]))]
				else:
					country_list.append((50*starting_country[0]+starting_country[1]))
					bool_map = (code_map[:,:,0] != 0) * (1 - (code_map[:,:,0] == code_map[latitude_pixel,longitude_pixel,0])*(code_map[:,:,1] == code_map[latitude_pixel,longitude_pixel,1])*(code_map[:,:,2] == code_map[latitude_pixel,longitude_pixel,2]))
					simple_code_map_here = simple_code_map*bool_map.astype(np.uint32)
					code_map_small = np.zeros((int(3600/search_scale),int(7200/search_scale),3))
					bool_map = cv2.resize(bool_map,(int(7200/search_scale),int(3600/search_scale)),interpolation = cv2.INTER_NEAREST).astype(np.uint32)
					bool_map_list.append(bool_map.astype(np.uint32))
					for i in range(0,int(3600/search_scale)):
						for j in range(0,int(7200/search_scale)):
							if code_map_small[i,j,0] == 0:
								code_map_small[i,j,0] = int(np.max(simple_code_map_here[search_scale*i:search_scale*i+search_scale,search_scale*j:search_scale*j+search_scale])/50)
								code_map_small[i,j,1] = int(np.max(simple_code_map_here[search_scale*i:search_scale*i+search_scale,search_scale*j:search_scale*j+search_scale])%50)
					code_map_small_list.append(code_map_small.astype(np.uint32))
				here_x = math.sin(latitude_pixel*math.pi/code_map.shape[0])*math.cos(longitude_pixel*math.pi*2/code_map.shape[1])
				here_y = math.sin(latitude_pixel*math.pi/code_map.shape[0])*math.sin(longitude_pixel*math.pi*2/code_map.shape[1])
				all_distances = (all_x-here_x)**2+(all_y-here_y)**2+long_z-4
				all_distances *= bool_map
				nearest_point = np.where(all_distances == np.amin(all_distances))
				first_landmass = code_map_small[nearest_point[0][0],nearest_point[1][0]]
				closest_codes[int(latitude_pixel/grid_size),int(longitude_pixel/grid_size)] = 50*50*(50*starting_country[0]+starting_country[1])+(50*first_landmass[0]+first_landmass[1])
				if 50*50*(50*starting_country[0]+starting_country[1])+(50*first_landmass[0]+first_landmass[1]) not in flags_to_place_codes:
					flags_to_place_codes.append(50*50*(50*starting_country[0]+starting_country[1])+(50*first_landmass[0]+first_landmass[1]))
					flags_to_place.append(chr(int(64+first_landmass[0]))+chr(int(64+first_landmass[1])))
		print(time.time()-start_time)
	np.save('npy/'+str(latitude_pixel_base)+'_closest_codes.npy', closest_codes)
	np.save('npy/'+str(latitude_pixel_base)+'_flags_to_place.npy', flags_to_place)
	np.save('npy/'+str(latitude_pixel_base)+'_flags_to_place_codes.npy', flags_to_place_codes)

proc = []
parallel = 3
index = 0
for latitude_pixel_base_2 in range(0,3600,100*parallel):
	if latitude_pixel_base_2 == 1200:
		latitude_pixel_base = 900
		print('Starting: '+str(latitude_pixel_base))
		p = Process(target=run_latitude,args=(latitude_pixel_base,))
		index += 1
		p.start()
		proc.append(p)
	for latitude_pixel_base in range(latitude_pixel_base_2,latitude_pixel_base_2+100*parallel,100):
		print('Starting: '+str(latitude_pixel_base))
		p = Process(target=run_latitude,args=(latitude_pixel_base,))
		index += 1
		p.start()
		proc.append(p)
	for p in proc:
		p.join()

closest_codes = np.zeros((int(code_map.shape[0]/grid_size),int(code_map.shape[1]/grid_size)))
output_map = np.zeros(code_map.shape, dtype=np.uint32)
flags_to_place = []
flags_to_place_codes = []
flags_to_place = np.array(flags_to_place)
for row in range(0,36):
    closest_codes += np.load('npy/'+str(100*row)+'_closest_codes.npy',allow_pickle=True)
    flags_to_place = np.concatenate((flags_to_place,np.load('npy/'+str(100*row)+'_flags_to_place.npy')))
    flags_to_place_codes = np.concatenate((flags_to_place_codes,np.load('npy/'+str(100*row)+'_flags_to_place_codes.npy')))

done_codes = []

for i in range(0,len(flags_to_place_codes)):
    if flags_to_place_codes[i] not in done_codes:
        done_codes.append(flags_to_place_codes[i])
        print(flags_to_place[i])
        current_flag = cv2.imread('flags/'+flags_to_place[i]+'.png')
        boundaries = boundries_of(closest_codes,flags_to_place_codes[i])
        scale_flag = max((boundaries[1]-boundaries[0])/current_flag.shape[0],(boundaries[3]-boundaries[2])/current_flag.shape[1])
        if int(current_flag.shape[0]*scale_flag) > 0 and int(current_flag.shape[1]*scale_flag) > 0:
            output_map[boundaries[0]:boundaries[1],boundaries[2]:boundaries[3],:] += (np.tile(np.reshape((closest_codes==flags_to_place_codes[i]).astype(float)[boundaries[0]:boundaries[1],boundaries[2]:boundaries[3]],(boundaries[1]-boundaries[0],boundaries[3]-boundaries[2],1)),(1,1,3))*\
            cv2.resize(current_flag,(int(1+current_flag.shape[1]*scale_flag),int(1+current_flag.shape[0]*scale_flag)))[\
            int((int(1+current_flag.shape[0]*scale_flag))/2-(boundaries[1]-boundaries[0])/2):\
            int((int(1+current_flag.shape[0]*scale_flag))/2-(boundaries[1]-boundaries[0])/2+boundaries[1]-boundaries[0]),\
            int((int(1+current_flag.shape[1]*scale_flag))/2-(boundaries[3]-boundaries[2])/2):\
            int((int(1+current_flag.shape[1]*scale_flag))/2-(boundaries[3]-boundaries[2])/2+boundaries[3]-boundaries[2]),:\
            ]).astype(np.uint32)
output_map = output_map.astype(np.uint32)
output_map[:,:,0] += 127*(code_map[:,:,0] == 0).astype(np.uint32)
output_map[:,:,1] += 127*(code_map[:,:,0] == 0).astype(np.uint32)
output_map[:,:,2] += 127*(code_map[:,:,0] == 0).astype(np.uint32)
output_map = output_map.astype(np.uint8)
output_map[:closest_codes.shape[0],:closest_codes.shape[1],0] *= (closest_codes == np.roll(closest_codes,1,axis=0)).astype(np.uint32)
output_map[:closest_codes.shape[0],:closest_codes.shape[1],1] *= (closest_codes == np.roll(closest_codes,1,axis=0)).astype(np.uint32)
output_map[:closest_codes.shape[0],:closest_codes.shape[1],2] *= (closest_codes == np.roll(closest_codes,1,axis=0)).astype(np.uint32)
output_map[:closest_codes.shape[0],:closest_codes.shape[1],0] *= (closest_codes == np.roll(closest_codes,1,axis=1)).astype(np.uint32)
output_map[:closest_codes.shape[0],:closest_codes.shape[1],1] *= (closest_codes == np.roll(closest_codes,1,axis=1)).astype(np.uint32)
output_map[:closest_codes.shape[0],:closest_codes.shape[1],2] *= (closest_codes == np.roll(closest_codes,1,axis=1)).astype(np.uint32)
im = Image.fromarray(np.flip(output_map,2))
im.save('output.png')
