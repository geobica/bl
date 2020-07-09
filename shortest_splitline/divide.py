import numpy as np
import math
from PIL import Image
import sys
import time

folder = 'pop_arr/'

n_pop_ori = np.load(folder+'n_array.npy')
s_pop_ori = np.load(folder+'s_array.npy')
a_pop_ori = np.load(folder+'a_array.npy')
w_pop_ori = np.load(folder+'w_array.npy')
m_pop_ori = np.load(folder+'m_array.npy')
e_pop_ori = np.load(folder+'e_array.npy')

n_pop = np.sign(np.load(folder+sys.argv[1]+'n_array.npy'))*n_pop_ori
s_pop = np.sign(np.load(folder+sys.argv[1]+'s_array.npy'))*s_pop_ori
a_pop = np.sign(np.load(folder+sys.argv[1]+'a_array.npy'))*a_pop_ori
w_pop = np.sign(np.load(folder+sys.argv[1]+'w_array.npy'))*w_pop_ori
m_pop = np.sign(np.load(folder+sys.argv[1]+'m_array.npy'))*m_pop_ori
e_pop = np.sign(np.load(folder+sys.argv[1]+'e_array.npy'))*e_pop_ori

land_map = np.load(folder+sys.argv[1]+'land.npy')
side_map = np.load(folder+sys.argv[1]+'side.npy')
np.save(folder+'side.npy',side_map.astype('uint8'))

divide = 0.5
if len(sys.argv) > 4:
	divide = float(sys.argv[4])

total_pop = np.sum(n_pop)+np.sum(s_pop)+np.sum(a_pop)+np.sum(w_pop)+np.sum(m_pop)+np.sum(e_pop)
print(total_pop)

output_radius = n_pop.shape[0]/2

scale = 0.2/0.00015707963397141786
big_number = 1000000000000

def great_circle_ratio(latitude,longitude,save):
	x = math.cos(longitude%(math.pi/2))*math.sin(latitude)
	y = math.sin(longitude%(math.pi/2))*math.sin(latitude)
	z = math.cos(latitude)
	if longitude < math.pi*2:
		a_dist = -x/math.sqrt(1-x*x)*scale
		w_dist = y/math.sqrt(1-y*y)*scale
		n_dist = -z/math.sqrt(1-z*z)*scale
		m_dist = a_dist
		e_dist = w_dist
		s_dist = n_dist
		a_slope = z/y
		w_slope = -z/x
		n_slope = y/x
		m_slope = -a_slope
		e_slope = -w_slope
		s_slope = n_slope

		needed_n = (np.linspace(np.linspace(1+(-output_radius*n_slope-output_radius-math.sqrt(1+n_slope*n_slope)*n_dist)/big_number,1+(output_radius*n_slope-output_radius-math.sqrt(1+n_slope*n_slope)*n_dist)/big_number,int(2*output_radius),endpoint=False),np.linspace(1+(-output_radius*n_slope+output_radius-math.sqrt(1+n_slope*n_slope)*n_dist)/big_number,1+(output_radius*n_slope+output_radius-math.sqrt(1+n_slope*n_slope)*n_dist)/big_number,int(2*output_radius),endpoint=False),int(2*output_radius),endpoint=False).astype(int))

		needed_s = (1-needed_n)

		needed_a = (np.linspace(np.linspace(1+(-output_radius*a_slope-output_radius-math.sqrt(1+a_slope*a_slope)*a_dist)/big_number,1+(output_radius*a_slope-output_radius-math.sqrt(1+a_slope*a_slope)*a_dist)/big_number,int(2*output_radius),endpoint=False),np.linspace(1+(-output_radius*a_slope+output_radius-math.sqrt(1+a_slope*a_slope)*a_dist)/big_number,1+(output_radius*a_slope+output_radius-math.sqrt(1+a_slope*a_slope)*a_dist)/big_number,int(2*output_radius),endpoint=False),int(2*output_radius),endpoint=False).astype(int))

		needed_w = (np.linspace(np.linspace(1-(-output_radius*w_slope-output_radius-math.sqrt(1+w_slope*w_slope)*w_dist)/big_number,1-(output_radius*w_slope-output_radius-math.sqrt(1+w_slope*w_slope)*w_dist)/big_number,int(2*output_radius),endpoint=False),np.linspace(1-(-output_radius*w_slope+output_radius-math.sqrt(1+w_slope*w_slope)*w_dist)/big_number,1-(output_radius*w_slope+output_radius-math.sqrt(1+w_slope*w_slope)*w_dist)/big_number,int(2*output_radius),endpoint=False),int(2*output_radius),endpoint=False).astype(int))

		needed_m = (np.linspace(np.linspace(1-(-output_radius*m_slope-output_radius-math.sqrt(1+m_slope*m_slope)*m_dist)/big_number,1-(output_radius*m_slope-output_radius-math.sqrt(1+m_slope*m_slope)*m_dist)/big_number,int(2*output_radius),endpoint=False),np.linspace(1-(-output_radius*m_slope+output_radius-math.sqrt(1+m_slope*m_slope)*m_dist)/big_number,1-(output_radius*m_slope+output_radius-math.sqrt(1+m_slope*m_slope)*m_dist)/big_number,int(2*output_radius),endpoint=False),int(2*output_radius),endpoint=False).astype(int))

		needed_e = (np.linspace(np.linspace(1+(-output_radius*e_slope-output_radius-math.sqrt(1+e_slope*e_slope)*e_dist)/big_number,1+(output_radius*e_slope-output_radius-math.sqrt(1+e_slope*e_slope)*e_dist)/big_number,int(2*output_radius),endpoint=False),np.linspace(1+(-output_radius*e_slope+output_radius-math.sqrt(1+e_slope*e_slope)*e_dist)/big_number,1+(output_radius*e_slope+output_radius-math.sqrt(1+e_slope*e_slope)*e_dist)/big_number,int(2*output_radius),endpoint=False),int(2*output_radius),endpoint=False).astype(int))

		if longitude < math.pi/2:
			n_pop_side = n_pop*needed_n
			s_pop_side = s_pop*needed_s
			a_pop_side = a_pop*needed_a
			w_pop_side = w_pop*needed_w
			m_pop_side = m_pop*needed_m
			e_pop_side = e_pop*needed_e
		elif longitude < math.pi:
			n_pop_side = n_pop*np.rot90(needed_n)
			s_pop_side = s_pop*np.rot90(needed_s)
			w_pop_side = w_pop*needed_a
			m_pop_side = m_pop*needed_w
			e_pop_side = e_pop*needed_m
			a_pop_side = a_pop*needed_e
		elif longitude < 3*math.pi/2:
			n_pop_side = n_pop*np.rot90(needed_n,2)
			s_pop_side = s_pop*np.rot90(needed_s,2)
			m_pop_side = m_pop*needed_a
			e_pop_side = e_pop*needed_w
			a_pop_side = a_pop*needed_m
			w_pop_side = w_pop*needed_e
		else:
			n_pop_side = n_pop*np.rot90(needed_n,3)
			s_pop_side = s_pop*np.rot90(needed_s,3)
			e_pop_side = e_pop*needed_a
			a_pop_side = a_pop*needed_w
			w_pop_side = w_pop*needed_m
			m_pop_side = m_pop*needed_e

		if save == 1:
			x = math.cos(longitude)*math.sin(latitude)
			y = math.sin(longitude)*math.sin(latitude)
			z = math.cos(latitude)
			size = land_map.shape[0]
			latitudes = np.tile(np.reshape(np.arange(0,math.pi,math.pi/2/(size/2)),(int(size),1)),(1,int(size*2)))
			longitudes = np.tile(np.arange(0,math.pi*2,math.pi*2/(size*2)),(int(size),1))
			coords = np.zeros((int(size),int(size*2),3))
			coords[:,:,0] = np.cos(longitudes)*np.sin(latitudes)
			coords[:,:,1] = np.sin(longitudes)*np.sin(latitudes)
			coords[:,:,2] = np.cos(latitudes)
			del latitudes
			del longitudes
			coords[:,:,0] *= x
			coords[:,:,1] *= y
			coords[:,:,2] *= z
			side_array = coords[:,:,0]+coords[:,:,1]+coords[:,:,2]

			np.save(folder + sys.argv[2] + 'land.npy',(land_map*(1+np.sign(side_array))/2).astype('uint8'))
			np.save(folder + sys.argv[3] + 'land.npy',(land_map*(1-np.sign(side_array))/2).astype('uint8'))

			np.save(folder + sys.argv[2] + 'side.npy',(side_map*(1+np.sign(side_array))/2).astype('uint8'))
			np.save(folder + sys.argv[3] + 'side.npy',(side_map*(1-np.sign(side_array))/2).astype('uint8'))

			np.save(folder + sys.argv[2] + 'n_array.npy',np.sign(n_pop_side).astype('uint8'))
			np.save(folder + sys.argv[2] + 's_array.npy',np.sign(s_pop_side).astype('uint8'))
			np.save(folder + sys.argv[2] + 'a_array.npy',np.sign(a_pop_side).astype('uint8'))
			np.save(folder + sys.argv[2] + 'w_array.npy',np.sign(w_pop_side).astype('uint8'))
			np.save(folder + sys.argv[2] + 'm_array.npy',np.sign(m_pop_side).astype('uint8'))
			np.save(folder + sys.argv[2] + 'e_array.npy',np.sign(e_pop_side).astype('uint8'))
			np.save(folder + sys.argv[3] + 'n_array.npy',np.sign(n_pop-n_pop_side).astype('uint8'))
			np.save(folder + sys.argv[3] + 's_array.npy',np.sign(s_pop-s_pop_side).astype('uint8'))
			np.save(folder + sys.argv[3] + 'a_array.npy',np.sign(a_pop-a_pop_side).astype('uint8'))
			np.save(folder + sys.argv[3] + 'w_array.npy',np.sign(w_pop-w_pop_side).astype('uint8'))
			np.save(folder + sys.argv[3] + 'm_array.npy',np.sign(m_pop-m_pop_side).astype('uint8'))
			np.save(folder + sys.argv[3] + 'e_array.npy',np.sign(e_pop-e_pop_side).astype('uint8'))

		return ((np.sum(n_pop_side)+np.sum(s_pop_side)+np.sum(a_pop_side)+np.sum(w_pop_side)+np.sum(m_pop_side)+np.sum(e_pop_side))/total_pop)

def rect_from_sph(latitude_in,longitude_in,altitude_in):
	return altitude_in*np.array([math.cos(latitude_in)*math.cos(longitude_in),math.cos(latitude_in)*math.sin(longitude_in),math.sin(latitude_in)])

def sph_from_rect(x_in,y_in,z_in):
	return [math.atan2(y_in,x_in),math.pi/2-math.atan2(math.sqrt(x_in*x_in+y_in*y_in),z_in),math.sqrt(x_in*x_in+y_in*y_in+z_in*z_in)]

ring_size = 1000

def great_circle_length(latitude,longitude):
	x = math.cos(longitude)*math.sin(latitude)
	y = math.sin(longitude)*math.sin(latitude)
	z = math.cos(latitude)
	coord_0 = np.array([x,y,z])
	coord_1 = np.cross([0,0,1],coord_0)
	coord_1 /= np.dot(coord_1,coord_1)
	coord_2 = np.cross(coord_0,coord_1)
	theta_ring = np.linspace(0,math.pi*2,ring_size,endpoint=False)
	coord_ring = np.matmul(np.array([coord_1]).T,np.array([np.cos(theta_ring)]))+np.matmul(np.array([coord_2]).T,np.array([np.sin(theta_ring)]))
	sph_ring_lat = math.pi/2-np.arctan2(np.sqrt(coord_ring[0,:]*coord_ring[0,:]+coord_ring[1,:]*coord_ring[1,:]),coord_ring[2,:])
	sph_ring_long = np.arctan2(coord_ring[1,:],coord_ring[0,:])
	ring_land = 0
	for i in range(0,len(sph_ring_lat)-1):
		if land_map[int(land_map.shape[0]/2-land_map.shape[0]/2*sph_ring_lat[i]/(math.pi/2)),int(land_map.shape[1]/2*sph_ring_long[i]/math.pi)] == 1:
			for j in range(i+1,len(sph_ring_lat)):
				if land_map[int(land_map.shape[0]/2-land_map.shape[0]/2*sph_ring_lat[j]/(math.pi/2)),int(land_map.shape[1]/2*sph_ring_long[j]/math.pi)] == 1:
					if math.pi-abs(abs(theta_ring[i]-theta_ring[j])-math.pi) > ring_land:
						ring_land = math.pi-abs(abs(theta_ring[i]-theta_ring[j])-math.pi)
	return ring_land

minimum_step_i = math.pi*2/400000
minimum_step_j = math.pi*2/4000
j_step = math.pi*2/10
checked_i = [0]
checked_j = [math.pi]
checked_lengths = [ring_size]
check_radius = math.pi
start_time = time.time()
shrink_factor = 4
while abs(j_step) > minimum_step_j:
	check_upto = len(checked_lengths)
	for check in range(0,check_upto):
		if checked_lengths[check] == np.min(np.array(checked_lengths)[:check_upto]):
			for j in np.arange(checked_j[check]-check_radius+minimum_step_j,checked_j[check]+check_radius-minimum_step_j,j_step):
				print('Searching: '+str(j))
				i = minimum_step_i
				step = math.pi/shrink_factor
				last = -1
				here = -1
				while i < math.pi/2 + math.pi/2*(divide != 0.5):
					here = great_circle_ratio(i,(2*math.pi+j)%(2*math.pi),0)
					if last != -1 and (last-divide)*(here-divide) < 0:
						step *= -0.5
					i += step
					if i > (math.pi/2 + math.pi/2*(divide != 0.5)):
						step *= 0.5
						i -= step
					last = here
					if abs(step) < minimum_step_i:
						step = math.pi/shrink_factor
						if i < (math.pi/2 + math.pi/2*(divide != 0.5))-minimum_step_i:
							checked_i.append(i)
							checked_j.append((2*math.pi+j)%(2*math.pi))
							checked_lengths.append(great_circle_length(i,(2*math.pi+j)%(2*math.pi)))
						last = -1
			break
	j_step /= shrink_factor
	check_radius /= shrink_factor
for check in range(0,len(checked_lengths)):
	if checked_lengths[len(checked_lengths)-1-check] == np.min(np.array(checked_lengths)):
		print(str(checked_i[len(checked_lengths)-1-check])+','+str(checked_j[len(checked_lengths)-1-check])+' ~ ratio='+str(int(1000000*great_circle_ratio(checked_i[len(checked_lengths)-1-check],checked_j[len(checked_lengths)-1-check],1))/1000000)+' length='+str((great_circle_length(checked_i[len(checked_lengths)-1-check],checked_j[len(checked_lengths)-1-check]))))
		print(sys.argv[2] + ', ' + sys.argv[3] + ' done.')
		print(time.time()-start_time)
		break
