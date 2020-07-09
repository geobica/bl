import numpy as np
import math
import cv2
from PIL import Image
import sys
import time

def hsv_to_rgb(hsv):
	# Does not use the actual HSV color space.
	# Hues are spaced evenly between the red, yellow, green, and blue defined below.
	primary_colors = np.array([[223,0,52],[255,219,0],[0,178,103],[0,133,191]])
	rgb = np.zeros(hsv.shape)
	rgb[..., 3:] = hsv[..., 3:]
	h, s, v = hsv[..., 0], hsv[..., 1], hsv[..., 2]
	for pr in range(0,primary_colors.shape[0]):
		rgb += ((s*v*max(0,min(1,1-min(abs(h-pr/primary_colors.shape[0]),1-abs(h-pr/primary_colors.shape[0]))*(primary_colors.shape[0]))))*primary_colors[pr,:]*primary_colors[pr,:])
	rgb += (v-s*v)*np.array([255,255,255])
	return np.sqrt(rgb).astype('uint8')

folder = 'pop_arr/'

land_map = []
side_map = []
composite = np.zeros((1800,3600,3))
size = composite.shape[0]
amount_of_divisions = int(sys.argv[2])

for i in range(0,amount_of_divisions):
	color_id = ((i)%amount_of_divisions)
	land_map.append(np.load(folder+sys.argv[1]+str(i)+'_land.npy'))
	side_map.append(np.load(folder+sys.argv[1]+str(i)+'_side.npy'))
	land_map[i] = cv2.resize(np.float32(land_map[i]), dsize=(int(size*2),int(size)),interpolation = cv2.INTER_NEAREST)
	side_map[i] = cv2.resize(np.float32(side_map[i]), dsize=(int(size*2),int(size)),interpolation = cv2.INTER_NEAREST)
	composite[:,:,:] += 0.75*np.tile(land_map[i].astype('uint8').reshape((composite.shape[0],composite.shape[1],1)),(3))*np.tile(hsv_to_rgb(np.array([color_id/amount_of_divisions,1,0.5+0.5*((color_id)%2)])),((composite.shape[0],composite.shape[1],1)))
	composite[:,:,:] += 0.25*np.tile(side_map[i].astype('uint8').reshape((composite.shape[0],composite.shape[1],1)),(3))*np.tile(hsv_to_rgb(np.array([color_id/amount_of_divisions,1,0.5+0.5*((color_id)%2)])),((composite.shape[0],composite.shape[1],1)))

del land_map
del side_map
equi_log_pop = np.load('equi_log_pop.npy')/255

print(np.max(0.25+0.75*np.maximum(np.maximum(np.maximum(np.abs(np.sign(np.roll(composite,1,1)-composite)),np.abs(np.sign(np.roll(composite,1,0)-composite))),np.maximum(np.abs(np.sign(np.roll(composite,-1,1)-composite)),np.abs(np.sign(np.roll(composite,-1,0)-composite)))),np.tile(np.reshape(equi_log_pop,(equi_log_pop.shape[0],equi_log_pop.shape[1],1)),(1,1,3)))))
composite *= 0.25+0.75*np.maximum(np.maximum(np.maximum(np.abs(np.sign(np.roll(composite,1,1)-composite)),np.abs(np.sign(np.roll(composite,1,0)-composite))),np.maximum(np.abs(np.sign(np.roll(composite,-1,1)-composite)),np.abs(np.sign(np.roll(composite,-1,0)-composite)))),np.tile(np.reshape(equi_log_pop,(equi_log_pop.shape[0],equi_log_pop.shape[1],1)),(1,1,3)))

img = Image.fromarray(composite.astype('uint8'), 'RGB')
img.save('output_map.png')
