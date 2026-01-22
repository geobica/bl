import hsluv
import math
import numpy as np
import os
import math
from matplotlib import pyplot as plt
from PIL import Image



# these files are large and I will not be putting them in the repo
pr_by_year = np.load(f"../../../../gravel/unorganized_projects/april_showers/by_year_pr_(40, 12, 522, 1080).npy")[1:]
tas_by_year = np.load(f"../../../../gravel/unorganized_projects/april_showers/by_year_tas_(24, 12, 522, 1080).npy") / 10

veg_total_by_month = np.load("../../../../gravel/unorganized_projects/vegetation/1k_NDVI_total_by_month.npy")
veg_count_by_month = np.load("../../../../gravel/unorganized_projects/vegetation/1k_NDVI_count_by_month.npy")
# pr_by_year = np.load(f"downsampled_chelsa_pr.npy")[1:]
# tas_by_year = np.load(f"downsampled_chelsa_tas.npy") / 10

# veg_total_by_month = np.load("1k_NDVI_total_by_month.npy")
# veg_count_by_month = np.load("1k_NDVI_count_by_month.npy")
veg_mean_by_month = (veg_total_by_month/veg_count_by_month).astype(float)

def equi_to_merc(equi_arr:np.ndarray,file_code:str="pr"):
    # equi_arr should be shape (n, h, w)
    # n is the number of channels, based on what data is being used
    # for CHELSA pr for example, n=12 because it's the averaged data for each of the 12 months
    # h = the height of the map
    # w = the width of the map

    if os.path.exists(f"{file_code}_mercator.npy"):
        return np.load(f"{file_code}_mercator.npy")

    # because the vegetation index had no data for some months because I ignored places with too little to avoid cloud cover issues
    # for the sake of /bl/msc I am plotting it as log(max(veg,1000)/1000)
    equi_arr = np.nan_to_num(equi_arr,nan=0)

    if file_code in ["pr","tas"]:
        # the CHELSA data excludes some land around the north pole
        lon_min, lat_max, lon_max, lat_min = -180, 84, 180, -90
    else:
        lon_min, lat_max, lon_max, lat_min = -180, 90, 180, -90
    mercator_x = np.linspace(0,1080,num=1080,endpoint=False).astype(int)
    mercator_y = np.linspace(0,1080,num=1080,endpoint=False)
    xx,yy = np.meshgrid(mercator_x,mercator_y)
    # convert these coordinates to equirectangular to sample from equi_arr with
    lon = np.floor(equi_arr.shape[2]/2+0.5+equi_arr.shape[2]/2/(1080/2)*(xx-1080/2)).astype(int)
    lat = np.arctan(np.sinh(math.pi/(1080/2)*(yy-1080/2)))*180/math.pi
    lat_j = equi_arr.shape[1]*(lat+lat_max)/(lat_max-lat_min)
    lat_j_floor = np.floor(lat_j).astype(int)
    # progress from row lat_j_floor to row lat_j_floor+1
    lat_j_p = lat_j-lat_j_floor

    as_mercator = np.zeros((equi_arr.shape[0],1080,1080))
    for channel_i in range(equi_arr.shape[0]):
        before = equi_arr[channel_i,np.minimum(equi_arr.shape[1]-1,np.maximum(0,lat_j_floor)),np.minimum(equi_arr.shape[2]-1,np.maximum(0,lon))]
        after = equi_arr[channel_i,np.minimum(equi_arr.shape[1]-1,np.maximum(0,lat_j_floor+1)),np.minimum(equi_arr.shape[2]-1,np.maximum(0,lon))]
        as_mercator[channel_i] = after*lat_j_p+before*(1-lat_j_p)
    np.save(f"{file_code}_mercator.npy",as_mercator)
    return as_mercator

pr_mercator = equi_to_merc(np.mean(pr_by_year,axis=0),file_code="pr")
tas_mercator = equi_to_merc(np.mean(tas_by_year,axis=0),file_code="tas")
veg_mercator = equi_to_merc(veg_mean_by_month,file_code="veg")

def split_raster_into_rgb_pngs(arr:np.ndarray,js_name:str="month_thresholds",file_code:str="pr"):
    # arr should be shape (n, h, w)
    # where n is the number of rasters
    # in standard usage for /bl/msc, n=12, as there are 12 rasters, one for each month of the year
    # h = the height of the map
    # w = the width of the map
    # for /bl/msc, these should be mercator maps stretching from 85.05113°N, 180°W to 85.05113°S, 180°E
    percentile_rgb = np.zeros((arr.shape[1],arr.shape[2],3))
    residual_rgb = np.zeros((arr.shape[1],arr.shape[2],3))
    js_strings = []
    for rgb_set_i in range(int(math.ceil(arr.shape[0]/3))):
        for channel_i in range(3):
            n_i = 3*rgb_set_i+channel_i

            by_percentile = np.zeros((arr.shape[1],arr.shape[2]))
            residual = np.zeros((arr.shape[1],arr.shape[2]))
            all_thresholds = []
            for i in range(255):
                threshold = np.percentile(arr[n_i],i/254*100)
                all_thresholds.append(threshold)
            all_thresholds.append(np.max(arr[n_i]))
            js_strings.append(f"{js_name}[{n_i}] = ["+",".join([str(el) for el in all_thresholds])+"];")
            for i in range(256):
                threshold = all_thresholds[i]
                by_percentile += (arr[n_i]>=threshold).astype(int)
            for i in range(255):
                threshold = all_thresholds[i]
                if all_thresholds[i+1]>threshold:
                    residual += (by_percentile==i+1).astype(float)*(arr[n_i]-threshold)/(all_thresholds[i+1]-threshold)
            percentile_rgb[:,:,channel_i] = by_percentile-1
            residual_rgb[:,:,channel_i] = residual
        percentile_rgb_img = Image.fromarray(np.maximum(0,np.minimum(255,percentile_rgb)).astype('uint8'))
        percentile_rgb_img.save(f"percentile_{file_code}_{rgb_set_i}.png")
        residual_rgb_img = Image.fromarray(np.maximum(0,np.minimum(255,residual_rgb*256)).astype('uint8'))
        residual_rgb_img.save(f"residual_{file_code}_{rgb_set_i}.png")
    return f"""
var {js_name} = [];
for(var i=0;i<12;i+=1){{
    {js_name}.push([]);
}}
"""+"\n".join(js_strings)
with open("month_thresholds.js","w") as f:
    f.write(split_raster_into_rgb_pngs(pr_mercator,js_name="month_thresholds",file_code="pr"))
    f.write(split_raster_into_rgb_pngs(tas_mercator,js_name="tas_month_thresholds",file_code="tas"))
    f.write(split_raster_into_rgb_pngs(veg_mercator,js_name="veg_month_thresholds",file_code="veg"))
