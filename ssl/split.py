from PIL import Image, ImageDraw, ImageFont, ImageColor
import random
import math
import numpy as np
from osgeo import gdal
import matplotlib.pyplot as plt
import gc
import os

#densities=[60000,60000,60000,60000,60000,60000,37500,37500,37500,37500,37500,37500,37500,37500,37500,37500,37500,37500,37500,37500,17500,17500,17500,17500,17500,17500,17500,17500,17500,17500,17500,17500,9000,9000,9000,9000,9000,9000,9000,9000,9000,9000,9000,9000,9000,7000,7000,7000,7000,7000,7000,7000,5000,5000,5000,5000,5000,5000,3500,3500,3500,3500,3500,3500,2750,2750,2750,2750,2750,2750,2250,2250,2250,2250,2250,2250,1875,1875,1875,1875,1875,1875,1875,1625,1625,1625,1625,1625,1375,1375,1375,1375,1375,1375,1375,1125,1125,1125,1125,1125,900,900,900,900,900,900,900,700,700,700,700,700,700,550,550,550,550,550,550,550,450,450,450,450,450,450,450,450,350,350,350,350,350,350,250,250,250,250,250,250,250,250,175,175,175,175,175,175,175,125,125,125,125,125,125,125,125,75,75,75,75,75,75,75,75,37.5,37.5,37.5,37.5,37.5,37.5,37.5,37.5,17.5,17.5,17.5,17.5,17.5,17.5,17.5,17.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,4,4,4,4,4,4,4,4,2,2,2,2,2,2,2,2,2,2,2,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
densities=[]
densities.append(0.0001)
for i in range(1,254):
    densities.append(i)
densities.append(0)
densities.append(0)

size=2160.0
scal=1.0/(21600.0/size*21600.0/size)

scale=int(21600.0/size)
isize=int(size)

full_width=2*isize
full_height=isize
origin_x=size
origin_y=size/2.0
scale_x=size
scale_y=size/2.0

width=2*isize
height=isize
start_part=200
end_part=1200
graft=50
greatest_density=10000000
histy=[]
index=[]
binns=100.0
for i in range(0,int(binns)):
    histy.append(0)
    index.append(7.0/binns*i)
    
if os.path.exists('equi.tif'):
    ds = gdal.Open('equi.tif', gdal.GA_ReadOnly)
    rb = ds.GetRasterBand(1)
    img_array2 = rb.ReadAsArray()
else:
    print('File `equi.tif` does not exist, which is needed to run the program. I decided not to store this file on github as it\'s very large, but it consists of a 21600x43200 raster file of global population density on an equirectangular map.')
del ds
del rb
while(end_part<20000):

    img_array=img_array2[start_part-200:end_part+200,:]
    print(img_array.shape)
    gc.collect()
    no_ocean=np.maximum(img_array,0)
    rows=np.zeros((1,2*21600))
    rows_stand=np.zeros((1,2*21600))

    pix_at_eq=10.0/0.92592592593*10.0/0.92592592593
    #map_shifts=np.zeros((40,40,21600,21600*2))
    for i in range(start_part/graft,end_part/graft):
        print(i*graft)
        band=no_ocean[i*graft-start_part+200-20:i*graft-start_part+200+20+graft,:]
        print(band)
        row_add=no_ocean[i*graft-start_part+200:i*graft-start_part+200+graft,:]
        rows_stand=np.concatenate((rows_stand,row_add))
        for a in range(10,30):
            for b in range(0,40):
                if (a-20)*(a-20)+(b-20)*(b-20)*math.sin(i*graft*math.pi/21600.0)*math.sin(i*graft*math.pi/21600.0)<=pix_at_eq:
                    band_alt=np.roll(np.roll(band,a-20,axis=0),b-20,axis=1)
                    row_add=row_add+band_alt[20:20+graft,:]
        row_add=row_add*math.sin(i*graft*math.pi/21600.0)
        #print(rows.shape)
        #print(row_add.shape)
        rows=np.concatenate((rows,row_add))
    print(np.histogram(rows_stand,bins=[0,1000,2000,3000,4000,5000,6000,1000000]))
    print(np.histogram(rows,bins=[0,1000,2000,3000,4000,5000,6000,1000000]))
    for i in range(0,int(binns)):
        what=np.sign(np.maximum(np.abs(rows-(math.pow(10.0,7.0*i/binns)+math.pow(10.0,7.0*(i+1.0)/binns))*0.5)*(-1.0)+(math.pow(10.0,7.0*(i+1.0)/binns)-math.pow(10.0,7.0*i/binns))*0.5,0))
        print(math.pow(10.0,7.0*(i+1.0)/binns))
        print(math.pow(10.0,7.0*i/binns))
        print(np.average(what))
        histy[i]+=np.sum(what*rows_stand)
    start_part+=1000
    end_part+=1000
print(greatest_density)
print(histy)
_ = plt.bar(index,histy, width = greatest_density/binns/2.0, align='center', alpha=0.5)  # arguments are passed to np.histogram
plt.title("histogram")
plt.savefig('histogram.png')
plt.show()

rows=rows*rows

popmap_array=np.zeros((height,width))
    
shrunk=np.zeros((isize,isize*2))
land=np.zeros((isize,isize*2))
    
for i in range(0,isize):
    if i%100==0:
        print(i)
    for j in range(0,isize*2):
        popmap_array[i,j]=np.sum(np.maximum(img_array[scale*i:scale*i+scale,scale*j:scale*j+scale],0))*scal #so ocean isn't -200
        shrunk[i,j]=np.sum(img_array[scale*i:scale*i+scale,scale*j:scale*j+scale])*scal
#for i in range(200,21600-200):
#    for i in 
#full_width=1383#6911/5
#full_height=549#2745/5
#origin_x=3464.0/5.0
#origin_y=1632.0/5.0
#scale_x=3464.0/5.0
#scale_y=1728.9305949/5.0
# Open an Image
def open_image(path):
  newImage = Image.open(path)
  return newImage

# Save Image
def save_image(image, path):
  image.save(path, 'png')

# Get the pixel from the given image
def get_pixel(image, i, j):
  # Inside image bounds?
  width, height = image.size
  if i > width or j > height:
    return None

  # Get Pixel
  pixel = image.getpixel((i, j))
  return pixel

# Create a new image with the given size
def create_image(i, j):
  image = Image.new("RGBA", (i, j), "white")
  return image

def fractov(n):
    if(n==0):
        return 0.0
    if(n==1):
        return 0.5
    return fractov(n-2.0**(math.floor(math.log(n)/math.log(2))))+1.0/2.0**(1.0+math.floor(math.log(n)/math.log(2.0)))

def draw_map(image,divisions,region_data):
  new = create_image(2*isize,isize)
  pixels = new.load()
  for i in range(width):
    for j in range(height):
      # Get Land Value
      blue =  image[j,i]
      if(blue<0.5):
          pixels[i, j] =(255,255,255,255)
      else:
          hue=fractov(region_data[i][j]-(2.0**divisions-1.0))
          red1=0
          if (hue < 1 / 6.0):
              red1 = 1.0
          elif (hue < 1 / 3.0):
              red1 = 1 - 6.0 * (hue - 1.0 / 6.0)
          elif (hue < 2 / 3.0):
              red1 = 0
          elif (hue < 5 / 6.0):
              red1 = 6.0 * (hue - 2.0 / 3.0)
          else:
              red1 = 1
          green1 = 0.0
          if (hue < 1 / 6.0):
              green1 = 6.0 * hue
          elif (hue < 1.0 / 2.0):
              green1 = 1
          elif (hue < 2.0 / 3.0):
              green1 = 1 - 6.0 * (hue - 1.0 / 2.0)
          else:
              green1 = 0
          blue1 = 0.0
          if (hue < 1 / 3.0):
              blue1 = 0.0
          elif (hue < 1 / 2.0):
              blue1 = 6.0 * (hue - 1.0 / 3.0)
          elif (hue < 5 / 6.0):
              blue1 = 1.0
          else:
              blue1 = 1 - 6.0 * (hue - 5.0 / 6.0)
          pixels[i, j] = (int(255*red1),int(255*green1),int(255*blue1),255)
      if(i>0):
        if(region_data[i][j]!=region_data[i-1][j] or (blue>0.5 and image[j,i-1]<0.5) or (blue<=0.5 and image[j,i-1]>0.5)):
            pixels[i, j] = (0,0,0, 255)
      if(j>0):
        if(region_data[i][j]!=region_data[i][j-1] or (blue>0.5 and image[j-1,i]<=0.5) or (blue<=0.5 and image[j-1,i]>0.5)):
            pixels[i, j] = (0,0,0, 255)


  # Return new image
  return new

def find_average_point(image,region_data,id):
  # Get size
  width, height = image.size

  # Create new Image and a Pixel Map
  new = create_image(width, height)
  pixels = new.load()
  average_x=0
  average_y=0
  average_z=0
  total_density=0
  for i in range(width):
    for j in range(height):
      # Get Pixel
      if(region_data[i][j]==id):
          pixel = get_pixel(image, i, j)

          # Get R, G, B values (This are int from 0 to 255)
          red =   pixel[0]
          green = pixel[1]
          blue =  pixel[2]
          value = (red+green+blue)/3.0
          t=math.pi*float((i-origin_x)/scale_x)
          p=-math.pi/2.0*float((j-origin_y)/scale_y)
          average_x+=math.cos(p)*densities[int(value)]*math.cos(t)*math.cos(p)
          average_y+=math.cos(p)*densities[int(value)]*math.sin(t)*math.cos(p)
          average_z+=math.cos(p)*densities[int(value)]*math.sin(p)
          total_density+=math.cos(p)*densities[int(value)]
  average_x /= total_density
  average_y /= total_density
  average_z /= total_density
  radius=math.sqrt(average_x*average_x+average_y*average_y+average_z*average_z)
  coord=[math.atan2(average_y,average_x),math.asin(average_z/radius)]
  return coord

def measure(average_point,angle,id,region_data,original):
    v_0=[math.cos(average_point[0])*math.cos(average_point[1]),math.sin(average_point[0])*math.cos(average_point[1]),math.sin(average_point[1])]
    v_1=[v_0[1]/math.sqrt(v_0[0]*v_0[0]+v_0[1]*v_0[1]),-v_0[0]/math.sqrt(v_0[0]*v_0[0]+v_0[1]*v_0[1]),0]
    v_2=[v_0[1]*v_1[2]-v_0[2]*v_1[1],v_0[2]*v_1[0]-v_0[0]*v_1[2],v_0[0]*v_1[1]-v_0[1]*v_1[0]]
    v_t=[v_1[0]*math.cos(angle*2.0*math.pi/360.0)+v_2[0]*math.sin(angle*2.0*math.pi/360.0),v_1[1]*math.cos(angle*2.0*math.pi/360.0)+v_2[1]*math.sin(angle*2.0*math.pi/360.0),v_1[2]*math.cos(angle*2.0*math.pi/360.0)+v_2[2]*math.sin(angle*2.0*math.pi/360.0)]
    leng=0
    for i in range(0,1000):
        v_c=[v_1[0]*math.cos(i*2.0*math.pi/1000.0)+v_2[0]*math.sin(i*2.0*math.pi/1000.0),v_1[1]*math.cos(i*2.0*math.pi/1000.0)+v_2[1]*math.sin(i*2.0*math.pi/1000.0),v_1[2]*math.cos(i*2.0*math.pi/1000.0)+v_2[2]*math.sin(i*2.0*math.pi/1000.0)]
        if (origin_y - scale_y * math.asin(v_c[2])/(math.pi/2) >= full_height or origin_y - scale_y * math.asin(v_c[2])/(math.pi/2) <0):
            red = 255
            green = 255
            blue = 255
        else:
            value = original[int(origin_y-scale_y*math.asin(v_c[2])/(math.pi/2)),int((10*full_width+origin_x+scale_x*math.atan2(v_c[1],v_c[0])/(math.pi))%full_width)]
            if (region_data[int(origin_x+scale_x*math.atan2(v_c[1],v_c[0])/(math.pi))%full_width][int(origin_y-scale_y*math.asin(v_c[2])/(math.pi/2))] == id and value>0.5):
                leng+=1
    return leng

def split_region(region_data,average_point,initial,primary,secondary):
    print(average_point)
    v_0 = [math.cos((average_point[0] - origin_x) * (math.pi) / scale_x) * math.cos(
                    -(average_point[1] - origin_y) * (math.pi / 2) / scale_y),
                       math.sin((average_point[0] - origin_x) * (math.pi) / scale_x) * math.cos(
                           -(average_point[1] - origin_y) * (math.pi / 2) / scale_y),
                       math.sin(-(average_point[1] - origin_y) * (math.pi / 2) / scale_y)]
    for i in range(0,full_width):
        for j in range(0,full_height):
            if(region_data[i][j]==initial):
                v_1 = [math.cos((i - origin_x) * (math.pi) / scale_x) * math.cos(
                    -(j - origin_y) * (math.pi / 2) / scale_y),
                       math.sin((i - origin_x) * (math.pi) / scale_x) * math.cos(
                           -(j - origin_y) * (math.pi / 2) / scale_y),
                       math.sin(-(j - origin_y) * (math.pi / 2) / scale_y)]
                if(v_0[0] * v_1[0] + v_0[1] * v_1[1] + v_0[2] * v_1[2]<0):
                    region_data[i][j]=primary
                else:
                    region_data[i][j]=secondary
    return region_data
def y(x,t,p):
    if p<=-math.pi/2 or p>=math.pi/2:
        return 0
    bar=math.cos(t)*math.cos(p)*math.cos(x)+math.sin(t)*math.cos(p)*math.sin(x)
    bumpy=math.acos(math.sin(p)/math.sqrt(bar*bar+math.sin(p)*math.sin(p)))
    if p>0:
        if abs(t-x)<math.pi/2 or abs(t-x)>3*math.pi/2:
            return -bumpy
        else:
            return bumpy
    else:
        if abs(t-x)<math.pi/2 or abs(t-x)>3*math.pi/2:
            return math.pi/2-bumpy
        else:
            return bumpy-math.pi
    
def sideRatio(density2,zoning,zone,t,p):
    restrictedDensity=density2*(1-np.maximum(np.minimum(np.abs(zoning-zone),1),0))
    if p<=0:
        return np.sum(restrictedDensity[:isize/2,:])/np.sum(restrictedDensity[isize/2:,:])
    if p>=size-1:
        return np.sum(restrictedDensity[isize/2:,:])/np.sum(restrictedDensity[:isize/2,:])
    if p==size/2:
        if t<size/2:
            return (np.sum(restrictedDensity[:,:t+isize/2])+np.sum(restrictedDensity[:,t+3*isize/2:]))/np.sum(restrictedDensity[:,t+isize/2:t+3*isize/2])
        if t>3*size/2:
            return (np.sum(restrictedDensity[:,:t-3*isize/2])+np.sum(restrictedDensity[:,t-isize/2:]))/np.sum(restrictedDensity[:,t-3*isize/2:t-3*isize/2])
        return np.sum(restrictedDensity[:,t-isize/2:t+3*isize/2])/(np.sum(restrictedDensity[:,:t-isize/2])+np.sum(restrictedDensity[:,t+isize/2:]))
    top=0
    bottom=0
    for j in range(0,isize*2):
        top+=np.sum(restrictedDensity[:int(size/math.pi*y(math.pi*j/size-math.pi,math.pi*t/size-math.pi,math.pi*p/size-math.pi/2.0)+size/2.0),j])
        bottom+=np.sum(restrictedDensity[int(size/math.pi*y(math.pi*j/size-math.pi,math.pi*t/size-math.pi,math.pi*p/size-math.pi/2.0)+size/2.0):,j])
    if p<size/2:
        return top/bottom
    return bottom/top
    
    
def split_difference(region_data,average_point,minimum_angle,initial):
    #print(str(int(origin_x+scale_x*average_point[0]/math.pi))+", "+str(int(origin_y-scale_y*average_point[1]/math.pi*2.0)))
    return sideRatio(popmap_array,np.swapaxes(np.array(region_data),0,1),initial,int(origin_x+scale_x*average_point[0]/math.pi),int(origin_y-scale_y*average_point[1]/math.pi*2.0))
    #print(average_point)
    #print(minimum_angle)
    primary_total=0
    secondary_total=0
    v_0 = [math.cos(average_point[0]) * math.cos(average_point[1]),
           math.sin(average_point[0]) * math.cos(average_point[1]), math.sin(average_point[1])]

    dot_products_X2=[]
    for i in range(0,width/100):
        row = []
        for j in range(0,height/100):
            v_1=[math.cos((100*i-origin_x)*(math.pi)/scale_x) * math.cos(-(100*j-origin_y)*(math.pi/2)/scale_y),
           math.sin((100*i-origin_x)*(math.pi)/scale_x) * math.cos(-(100*j-origin_y)*(math.pi/2)/scale_y), math.sin(-(100*j-origin_y)*(math.pi/2)/scale_y)]
            row.append(v_0[0]*v_1[0]+v_0[1]*v_1[1]+v_0[2]*v_1[2])
        dot_products_X2.append(row)
    for i in range(0,width/100-1):
        for j in range(0,height/100-1):
            any_in_region=(region_data[100*i][100*j]==initial or region_data[100*(i+1)][100*j]==initial or region_data[100*(i)][100*(j+1)]==initial or region_data[100*(i+1)][100*(j+1)]==initial)
            all_in_region=(region_data[100*i][100*j]==initial and region_data[100*(i+1)][100*j]==initial and region_data[100*(i)][100*(j+1)]==initial and region_data[100*(i+1)][100*(j+1)]==initial)
            if(dot_products_X2[i][j]*dot_products_X2[i+1][j]<0 or dot_products_X2[i][j]*dot_products_X2[i][j+1]<0 or (any_in_region and not all_in_region)):
                dot_products_X1 = []
                for k in range(0, 10):
                    row = []
                    for l in range(0, 10):
                        v_1 = [math.cos((100 * i+10*k - origin_x) * (math.pi) / scale_x) * math.cos(
                            -(100 * j +10*l- origin_y) * (math.pi / 2) / scale_y),
                               math.sin((100 * i+10*k - origin_x) * (math.pi) / scale_x) * math.cos(
                                   -(100 * j+10*l - origin_y) * (math.pi / 2) / scale_y),
                               math.sin(-(100 * j+10*l - origin_y) * (math.pi / 2) / scale_y)]
                        row.append(v_0[0] * v_1[0] + v_0[1] * v_1[1] + v_0[2] * v_1[2])
                    dot_products_X1.append(row)
                for k in range(0, 9):
                    for l in range(0, 9):
                        any_in_region = (
                        region_data[100 * i+10*k][100 * j+10*l] == initial or region_data[100 * (i)+10*(k+1)][100 * j+10*l] == initial or
                        region_data[100 * (i)][100 * j + 10*(l+1)] == initial or region_data[100 * (i)+10*(k+1)][
                            100 * j + 10*(l+1)] == initial)
                        all_in_region = (
                        region_data[100 * i+10*k][100 * j+10*l] == initial and region_data[100 * (i)+10*(k+1)][100 * j+10*l] == initial and
                        region_data[100 * (i)][100 * j + 10*(l+1)] == initial and region_data[100 * (i)+10*(k+1)][
                            100 * j + 10*(l+1)] == initial)

                        if (dot_products_X1[k][l] * dot_products_X1[k + 1][l] < 0 or dot_products_X1[k][l] *
                            dot_products_X1[k][l + 1] < 0 or (any_in_region and not all_in_region)):
                            dot_products_X0 = []
                            for m in range(0, 10):
                                row = []
                                for n in range(0, 10):
                                    v_1 = [math.cos((100 * i + 10 * k+m - origin_x) * (math.pi) / scale_x) * math.cos(
                                        -(100 * j + 10 * l+n - origin_y) * (math.pi / 2) / scale_y),
                                           math.sin((100 * i + 10 * k+m - origin_x) * (math.pi) / scale_x) * math.cos(
                                               -(100 * j + 10 * l+n - origin_y) * (math.pi / 2) / scale_y),
                                           math.sin(-(100 * j + 10 * l+n - origin_y) * (math.pi / 2) / scale_y)]
                                    row.append(v_0[0] * v_1[0] + v_0[1] * v_1[1] + v_0[2] * v_1[2])
                                dot_products_X0.append(row)
                            for m in range(0, 10):
                                for n in range(0, 10):
                                    if(region_data[100 * i + 10*k+m][100 * j + 10*l+n]==initial):
                                        if (dot_products_X0[m][n] > 0):
                                            primary_total += popmap_array[100 * j + 10*l+n,100 * i + 10*k+m]
                                        else:
                                            secondary_total += popmap_array[100 * j + 10*l+n,100 * i + 10*k+m]
                        else:
                            if(all_in_region):
                                if (dot_products_X1[k][l] > 0):
                                    primary_total += popmap_array_X1[10*i+k][10*j+l]
                                else:
                                    secondary_total += popmap_array_X1[10*i+k][10*j+l]
            else:
                if(all_in_region):
                    if(dot_products_X2[i][j]>0):
                        primary_total+=popmap_array_X2[i][j]
                    else:
                        secondary_total+=popmap_array_X2[i][j]
    if(primary_total==0):
        return 9999
    else:
        return secondary_total/primary_total

if 0==0:
    #original = open_image('truepop_equi_2.png')
    #original10 = open_image('popmap_10.png')
    
    land=0.5*np.sign(np.abs(shrunk+200)-0.00001)+0.5 #so land is 1 and ocean is 0
    
    popmap_array=popmap_array-1+land #so ocean is -1
    
    popmap_array_X1=[]
    for i in range(0,width/10):
        row = []
        for j in range(0,height/10):
            total=0
            for a in range(0, 10):
                for b in range(0, 10):
                    total+=popmap_array[10*j+b,10*i+a]
            row.append(total)
        popmap_array_X1.append(row)
    popmap_array_X2=[]
    for i in range(0,width/100):
        row = []
        for j in range(0,height/100):
            total=0
            for a in range(0, 10):
                for b in range(0, 10):
                    total+=popmap_array_X1[10*i+a][10*j+b]
            row.append(total)
        popmap_array_X2.append(row)
  # Example Pixel Color
    rrr=[[random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5],[random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5],[random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5],[random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5],[random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5],[random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5],[random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5],[random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5],[random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5],[random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5,random.random()-0.5]]
  # Convert to Grayscale and save
    divisions=2

    region_data = []
    for i in range(0,full_width):
        row=[]
        for j in range(0,full_height):
            row.append(0)
        region_data.append(row)

    region_centers = []
    for i in range(0,2**(divisions+1)-1):
        coord=[-1,-1]
        region_centers.append(coord)
    for i in range(0,divisions):
        for j in range(0, 2**i):
            pole_array_X2=[]
            for pole_long in range(0, width / 100):
                row = []
                for pole_lat in range(0, int(origin_y / 100)):
                    row.append(split_difference(region_data,[math.pi*float((100*pole_long-origin_x)/scale_x),-math.pi/2.0*float((100*pole_lat-origin_y)/scale_y)],99,2 ** i - 1 + j))
                row.append(split_difference(region_data,[math.pi*float((100*pole_long-origin_x)/scale_x),0],99,2 ** i - 1 + j))
                pole_array_X2.append(row)
            print(pole_array_X2)
            boxes_to_search=[]
            for pole_long_X2 in range(0, width / 100):
                for pole_lat_X2 in range(0, int(origin_y / 100)):
                    some_within=pole_array_X2[pole_long_X2][pole_lat_X2]<1 or pole_array_X2[(pole_long_X2+1)%(width/100)][pole_lat_X2]<1 or pole_array_X2[pole_long_X2][pole_lat_X2+1]<1 or pole_array_X2[(pole_long_X2+1)%(width/100)][pole_lat_X2+1]<1
                    within=pole_array_X2[pole_long_X2][pole_lat_X2]<1 and pole_array_X2[(pole_long_X2+1)%(width/100)][pole_lat_X2]<1 and pole_array_X2[pole_long_X2][pole_lat_X2+1]<1 and pole_array_X2[(pole_long_X2+1)%(width/100)][pole_lat_X2+1]<1
                    if(some_within and not within):
                        boxes_to_search.append([pole_long_X2,pole_lat_X2])
            boxes_to_search_X1=[]
            for box in boxes_to_search:
                pole_array_X1=[]
                for pole_long in range(10*box[0],10*box[0]+11):
                    row = []
                    for pole_lat in range(10*box[1],10*box[1]+11):
                        row.append(split_difference(region_data,[math.pi*float((10*pole_long-origin_x)/scale_x),-math.pi/2.0*float((10*pole_lat-origin_y)/scale_y)],99,2 ** i - 1 + j))
                    pole_array_X1.append(row)
                for pole_long_X1 in range(0,10):
                    for pole_lat_X1 in range(0,10):
                        some_within = pole_array_X1[pole_long_X1][pole_lat_X1] < 1 or \
                                      pole_array_X1[(pole_long_X1 + 1)][pole_lat_X1] < 1 or \
                                      pole_array_X1[pole_long_X1][pole_lat_X1 + 1] < 1 or \
                                      pole_array_X1[(pole_long_X1 + 1)][pole_lat_X1 + 1] < 1
                        within = pole_array_X1[pole_long_X1][pole_lat_X1] < 1 and \
                                      pole_array_X1[(pole_long_X1 + 1)][pole_lat_X1] < 1 and \
                                      pole_array_X1[pole_long_X1][pole_lat_X1 + 1] < 1 and \
                                      pole_array_X1[(pole_long_X1 + 1)][pole_lat_X1 + 1] < 1
                        if (some_within and not within):
                            boxes_to_search_X1.append([10*box[0]+pole_long_X1, 10*box[1]+pole_lat_X1])
            print(str(2**i-1+j)+" X1 "+str(len(boxes_to_search_X1)))
            boxes_to_search_X2 = []
            for box in boxes_to_search_X1:
                pole_array_X2 = []
                for pole_long in range(10 * box[0], 10 * box[0] + 11):
                    row = []
                    for pole_lat in range(10 * box[1], 10 * box[1] + 11):
                        row.append(split_difference(region_data,
                                                    [math.pi * float((pole_long - origin_x) / scale_x),
                                                     -math.pi / 2.0 * float((pole_lat - origin_y) / scale_y)], 99,
                                                    2 ** i - 1 + j))
                    pole_array_X2.append(row)
                for pole_long_X2 in range(0, 10):
                    for pole_lat_X2 in range(0, 10):
                        some_within = pole_array_X2[pole_long_X2][pole_lat_X2] < 1 or \
                                      pole_array_X2[(pole_long_X2 + 1)][pole_lat_X2] < 1 or \
                                      pole_array_X2[pole_long_X2][pole_lat_X2 + 1] < 1 or \
                                      pole_array_X2[(pole_long_X2 + 1)][pole_lat_X2 + 1] < 1
                        within = pole_array_X2[pole_long_X2][pole_lat_X2] < 1 and \
                                 pole_array_X2[(pole_long_X2 + 1)][pole_lat_X2] < 1 and \
                                 pole_array_X2[pole_long_X2][pole_lat_X2 + 1] < 1 and \
                                 pole_array_X2[(pole_long_X2 + 1)][pole_lat_X2 + 1] < 1
                        if (some_within and not within):
                            boxes_to_search_X2.append([10 * box[0] + pole_long_X2, 10 * box[1] + pole_lat_X2])
            minimum_point=[0,origin_y]
            minimum_measurement=full_width
            print(str(2**i-1+j)+" X2 "+str(len(boxes_to_search_X2)))
            for box in boxes_to_search_X2:
                measurement = measure(
                    [math.pi * float((box[0] - origin_x) / scale_x), -math.pi / 2.0 * float((box[1] - origin_y) / scale_y)],
                    a, 2 ** i - 1 + j, region_data, land)
                if(measurement<minimum_measurement):
                    minimum_point=box
                    minimum_measurement=measurement
            region_data=split_region(region_data,minimum_point,2**i-1+j,2**(i+1)-1+2*j,2**(i+1)+2*j)
            print(""+str(2**(i+1)-1+2*j)+" "+str(region_data[0][0])+" "+str(region_data[100][40]))

    new = draw_map(land,divisions,region_data)
    save_image(new, 'shortest.png')
