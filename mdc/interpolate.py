import sys

import numpy as np
import time
import scipy as sp
import fiona
import shapely
from shapely.geometry import MultiPolygon, Point
import scipy.interpolate
from scipy.spatial import Delaunay, cKDTree
from scipy import ndimage
from shapely.geometry import mapping, Polygon, LineString
import pickle
from math import pi,e,sin,cos
import os
from shapely.geometry import shape
from shapely.validation import make_valid
from utils import *

import matplotlib.tri as mtri


class ConstrainedTriangulation:
    def __init__(self,points,triangles):
        self.points = np.asarray(points,dtype=float)
        self.simplices = np.asarray(triangles,dtype=np.int64)

        if self.points.ndim!=2 or self.points.shape[1]!=2:
            raise ValueError('points must have shape (N, 2)')
        if self.simplices.ndim!=2 or self.simplices.shape[1]!=3:
            raise ValueError('triangles must have shape (M, 3)')

        self.tri_points = self.points[self.simplices]
        p0,p1,p2 = self.tri_points[:,0],self.tri_points[:,1],self.tri_points[:,2]
        self.v0 = p1-p0
        self.v1 = p2-p0
        self.den = self.v0[:,0]*self.v1[:,1]-self.v1[:,0]*self.v0[:,1]
        self.valid_triangles = np.abs(self.den)>1e-14

        self._triangulation = mtri.Triangulation(self.points[:,0],self.points[:,1],self.simplices)
        self._finder = self._triangulation.get_trifinder()

    def find_simplex(self,xy):
        xy = np.asarray(xy,dtype=float)
        if xy.shape[-1]!=2:
            raise ValueError('points must have shape (..., 2)')
        flat = xy.reshape(-1,2)
        found = self._finder(flat[:,0],flat[:,1]).astype(np.int64)
        return found.reshape(xy.shape[:-1])

    def barycentric(self,xy,simplex):
        tri = np.maximum(simplex,0)
        offset = xy-self.tri_points[tri,0]
        v0,v1,den = self.v0[tri],self.v1[tri],self.den[tri]
        u = (offset[:,0]*v1[:,1]-v1[:,0]*offset[:,1])/den
        v = (v0[:,0]*offset[:,1]-offset[:,0]*v0[:,1])/den
        return 1.0-u-v,u,v

    def interpolate(self,xy,values):
        xy = np.asarray(xy,dtype=float)
        values = np.asarray(values)
        flat = xy.reshape(-1,2)
        out = np.full(len(flat),np.nan,dtype=np.result_type(values,float))
        simplex = self.find_simplex(flat)
        hit = simplex>=0
        if np.any(hit):
            w,u,v = self.barycentric(flat[hit],simplex[hit])
            corners = self.simplices[simplex[hit]]
            out[hit] = w*values[corners[:,0]]+u*values[corners[:,1]]+v*values[corners[:,2]]
        return out.reshape(xy.shape[:-1])

    def linear_interpolator(self,values):
        interp = mtri.LinearTriInterpolator(self._triangulation,np.asarray(values,dtype=float))
        def evaluate(x,y):
            return np.ma.filled(interp(x,y),np.nan)
        return evaluate


D = 12756274
LOCALITY = 4.0
max_neighbor_dist = None

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

def points_to_file(filename,points):
    schema = {
            'geometry': 'Point',
    }
    with fiona.open(filename, 'w', 'ESRI Shapefile', schema) as c:
            for p in points:
                c.write({
                        'geometry': mapping(Point(np.real(p),np.imag(p))),
                })

def points_with_attrs_to_file(filename,points,attrs):
    schema = {
            'geometry': 'Point',
            'properties': {name: (type(vals[0]).__name__ if vals else 'str') for name,vals in attrs.items()},
    }
    with fiona.open(filename, 'w', 'ESRI Shapefile', schema) as c:
            for i,p in enumerate(points):
                c.write({
                        'geometry': mapping(Point(np.real(p),np.imag(p))),
                        'properties': {name: vals[i] for name,vals in attrs.items()},
                })

def polygons_with_attrs_to_file(filename,polys_xy,attrs):
    schema = {
            'geometry': 'Polygon',
            'properties': {name: (type(vals[0]).__name__ if vals else 'str') for name,vals in attrs.items()},
    }
    with fiona.open(filename, 'w', 'ESRI Shapefile', schema) as c:
            for i,poly_xy in enumerate(polys_xy):
                c.write({
                        'geometry': mapping(Polygon(poly_xy)),
                        'properties': {name: vals[i] for name,vals in attrs.items()},
                })

DBF_TEXT_LIMIT = 254

def named_points_to_file(filename,points,names,depths,arms=None):
    props = {'name': f'str:{DBF_TEXT_LIMIT}','name_ext': f'str:{DBF_TEXT_LIMIT}','depth': 'int'}
    if arms is not None:
        props['arm'] = 'str:3'
    schema = {'geometry': 'Point','properties': props}
    with fiona.open(filename, 'w', 'ESRI Shapefile', schema) as c:
            for i,(p,name,depth) in enumerate(zip(points,names,depths)):
                row = {'name': name[:DBF_TEXT_LIMIT],'name_ext': name[DBF_TEXT_LIMIT:],
                    'depth': int(depth)}
                if arms is not None:
                    row['arm'] = arms[i]
                c.write({
                        'geometry': mapping(Point(np.real(p),np.imag(p))),
                        'properties': row,
                })

EXTRA_POINT_ARMS = ('h','r','l','-l','-h','-r','-2l','-2h','-2r')

def mesh_point_names(tree_names,tree_depths,cell_count,point_count):
    extra = point_count-cell_count
    if extra<0 or cell_count==0 or extra%cell_count!=0:
        return None,None,None
    if extra//cell_count!=len(EXTRA_POINT_ARMS):
        return None,None,None

    names = list(tree_names)
    depths = list(tree_depths)
    arms = ['']*cell_count
    for i in range(cell_count):
        for arm in EXTRA_POINT_ARMS:
            names.append(tree_names[i]+arm)
            depths.append(tree_depths[i])
            arms.append(arm)
    return names,depths,arms

def lattice_structure(projection_name,cell_count):
    path = f'pickle/matlab_saves/{projection_name}_line_placements_w.txt'
    if not os.path.isfile(path):
        return None
    with open(path,'rb') as fh:
        line_placements_w = pickle.load(fh)
    if len(line_placements_w)!=cell_count:
        return None

    arms = np.array([np.hstack(row)[:4] for row in line_placements_w])
    self_w,anchor_w,right_w,left_w = arms[:,0],arms[:,1],arms[:,2],arms[:,3]

    w_tree = cKDTree(np.stack([self_w.real,self_w.imag],axis=1))
    _,parent = w_tree.query(np.stack([anchor_w.real,anchor_w.imag],axis=1))
    parent = parent.astype(int)
    parent[0] = 0

    letter = np.empty(cell_count,dtype='<U1')
    for i in range(1,cell_count):
        p = parent[i]
        reach_to = (anchor_w[p],right_w[p],left_w[p])
        letter[i] = 'HRL'[int(np.argmin([abs(self_w[i]-arm) for arm in reach_to]))]

    reach = np.abs(right_w-self_w)
    own = np.arange(cell_count)
    dist_r,near_r = w_tree.query(np.stack([right_w.real,right_w.imag],axis=1))
    dist_l,near_l = w_tree.query(np.stack([left_w.real,left_w.imag],axis=1))
    right_of = np.where((near_r!=own)&(dist_r<=0.01*reach),near_r,-1)
    left_of = np.where((near_l!=own)&(dist_l<=0.01*reach),near_l,-1)

    children = [[] for _ in range(cell_count)]
    for i in range(1,cell_count):
        children[parent[i]].append(i)

    depth = np.zeros(cell_count,dtype=int)
    stack = [(child,1) for child in children[0]]
    while stack:
        i,d = stack.pop()
        depth[i] = d
        for child in children[i]:
            stack.append((child,d+1))

    return {'parent':parent,'letter':letter,'children':children,'depth':depth,
        'right_of':right_of,'left_of':left_of}

def lattice_tree_names(projection_name,cell_count):
    structure = lattice_structure(projection_name,cell_count)
    if structure is None:
        return None,None
    letter,children = structure['letter'],structure['children']

    names = ['']*cell_count
    depths = np.zeros(cell_count,dtype=int)
    names[0] = 'O'
    stack = [(child,'',1) for child in children[0]]
    while stack:
        i,prefix,depth = stack.pop()
        names[i] = prefix+letter[i]
        depths[i] = depth
        for child in children[i]:
            stack.append((child,names[i],depth+1))
    return names,depths

def tree_links_to_file(filename,points,links):
    schema = {
            'geometry': 'LineString',
            'properties': {'kind': 'str:8','hops': 'int','a': 'int','b': 'int',
                'depth_a': 'int','depth_b': 'int'},
    }
    with fiona.open(filename, 'w', 'ESRI Shapefile', schema) as c:
            for a,b,kind,hops,depth_a,depth_b in links:
                c.write({
                        'geometry': mapping(LineString([
                            (float(np.real(points[a])),float(np.imag(points[a]))),
                            (float(np.real(points[b])),float(np.imag(points[b])))])),
                        'properties': {'kind': kind,'hops': int(hops),'a': int(a),'b': int(b),
                            'depth_a': int(depth_a),'depth_b': int(depth_b)},
                })

def _tree_hop_count(a,b,parent,depth):
    hops = 0
    while depth[a]>depth[b]:
        a = parent[a]; hops += 1
    while depth[b]>depth[a]:
        b = parent[b]; hops += 1
    while a!=b:
        a = parent[a]; b = parent[b]; hops += 2
    return hops

def lattice_tree_links(projection_name,cell_count):
    structure = lattice_structure(projection_name,cell_count)
    if structure is None:
        return None
    parent,depth = structure['parent'],structure['depth']
    right_of,left_of = structure['right_of'],structure['left_of']

    tree_pairs = {(min(parent[i],i),max(parent[i],i)) for i in range(1,cell_count)}
    neighbour_pairs = set()
    for i in range(cell_count):
        for j in (right_of[i],left_of[i]):
            if j>=0:
                neighbour_pairs.add((min(i,j),max(i,j)))

    links = []
    for a,b in sorted(tree_pairs|neighbour_pairs):
        in_tree = (a,b) in tree_pairs
        in_neighbour = (a,b) in neighbour_pairs
        kind = 'both' if (in_tree and in_neighbour) else ('tree' if in_tree else 'neighbor')
        links.append((a,b,kind,_tree_hop_count(a,b,parent,depth),depth[a],depth[b]))
    return links

ARM_INDEX = {arm:i for i,arm in enumerate(EXTRA_POINT_ARMS)}

DIRECTION_CYCLE = ('h','-l','r','-h','l','-r')
ARM_POSITION = {'h':0,'r':2,'l':4}

def seam_pairs(arm_here,arm_there):
    here,there = ARM_POSITION[arm_here],ARM_POSITION[arm_there]
    return tuple(('-2'+DIRECTION_CYCLE[(here+turn)%6][1],'-2'+DIRECTION_CYCLE[(there-turn)%6][1])
        for turn in (1,-1))

def lattice_mesh_links(projection_name,cell_count,point_count):
    structure = lattice_structure(projection_name,cell_count)
    if structure is None:
        return None
    extra = point_count-cell_count
    if extra<0 or cell_count==0 or extra%cell_count!=0 or extra//cell_count!=len(EXTRA_POINT_ARMS):
        return None
    parent,letter = structure['parent'],structure['letter']
    right_of,left_of = structure['right_of'],structure['left_of']

    anchor_of = np.array(parent,dtype=int)
    anchor_of[0] = -1
    for child in structure['children'][0]:
        if letter[child]=='H':
            anchor_of[0] = child

    def sample(cell,arm):
        return cell_count+cell*len(EXTRA_POINT_ARMS)+ARM_INDEX[arm]

    links = []
    seen = set()
    def add(a,b,kind):
        key = (min(a,b),max(a,b))
        if key not in seen:
            seen.add(key)
            links.append((key[0],key[1],kind))

    for i in range(cell_count):
        for arm in ('h','r','l'):
            add(i,sample(i,arm),'spoke')
        for arm in ('-h','-r','-l'):
            add(i,sample(i,arm),'spoke')
            add(sample(i,arm),sample(i,'-2'+arm[1]),'radial')
            for other in ('h','r','l'):
                if other!=arm[1]:
                    add(sample(i,arm),sample(i,other),'ring')

        for arm in [['r','-2l'],['r','-2h'],['l','-2r'],['l','-2h'],['h','-2r'],['h','-2l']]:
            add(sample(i,arm[0]),sample(i,arm[1]),'rectlong')

    heptagonal_name_dict = {}
    biggest_heptagonal_name = 0
    heptagonal_link_list = []

    for i in range(cell_count):
        for arm_here,neighbour in (('h',anchor_of[i]),('r',right_of[i]),('l',left_of[i])):
            if neighbour<0:
                continue
            if anchor_of[neighbour]==i:
                arm_there = 'h'
            elif right_of[neighbour]==i:
                arm_there = 'r'
            elif left_of[neighbour]==i:
                arm_there = 'l'
            else:
                continue
            add(i,sample(i,arm_here),'spoke')
            add(sample(i,arm_here),sample(neighbour,arm_there),'span')
            add(sample(neighbour,arm_there),neighbour,'spoke')

            if sample(i,arm_here)<sample(neighbour,arm_there):
                for other_arm_there in ['h','l','r']:
                    if other_arm_there!=arm_there:
                        add(sample(i,arm_here),sample(neighbour,'-2'+other_arm_there),'rectdiag')
            else:
                for other_arm_here in ['h','l','r']:
                    if other_arm_here!=arm_here:
                        add(sample(i,'-2'+other_arm_here),sample(neighbour,arm_there),'rectdiag')


            for mine,theirs in seam_pairs(arm_here,arm_there):
                add(sample(i,mine),sample(neighbour,theirs),'seam')

                heptagonal_link_list.append([int(sample(i,mine)),int(sample(neighbour,theirs))])
                if sample(i,mine) in heptagonal_name_dict:
                    if sample(neighbour,theirs) in heptagonal_name_dict:
                        to_convert = heptagonal_name_dict[sample(neighbour,theirs)]
                        if to_convert != heptagonal_name_dict[sample(i,mine)]:
                            for key in heptagonal_name_dict:
                                if heptagonal_name_dict[key] == to_convert:
                                    heptagonal_name_dict[key] = heptagonal_name_dict[sample(i,mine)]
                    else:
                        heptagonal_name_dict[sample(neighbour,theirs)] = heptagonal_name_dict[sample(i,mine)]
                elif sample(neighbour,theirs) in heptagonal_name_dict:
                    heptagonal_name_dict[sample(i,mine)] = heptagonal_name_dict[sample(neighbour,theirs)]
                else:
                    heptagonal_name_dict[sample(i,mine)] = biggest_heptagonal_name
                    heptagonal_name_dict[sample(neighbour,theirs)] = biggest_heptagonal_name
                    biggest_heptagonal_name += 1

    heptagon_names_used = []
    for key in heptagonal_name_dict:
         if heptagonal_name_dict[key] not in heptagon_names_used:
            heptagon_names_used.append(heptagonal_name_dict[key])
    points_by_heptagon = {name: [] for name in heptagon_names_used}

    for key, name in heptagonal_name_dict.items():
        if name in points_by_heptagon:
            points_by_heptagon[name].append(key)

    for name in tqdm(heptagon_names_used):
        if len(points_by_heptagon[name]) == 7:
            A = np.min(points_by_heptagon[name])

            A_links_to = []
            for link in heptagonal_link_list:
                if A in link:
                    if link[link.index(A)^1] not in A_links_to:
                        A_links_to.append(link[link.index(A)^1])
            B = np.min(A_links_to)
            G = np.max(A_links_to)

            B_links_to = []
            for link in heptagonal_link_list:
                if B in link:
                    if link[link.index(B)^1] not in B_links_to:
                        B_links_to.append(link[link.index(B)^1])
            C = B_links_to[B_links_to.index(A)^1]

            C_links_to = []
            for link in heptagonal_link_list:
                if C in link:
                    if link[link.index(C)^1] not in C_links_to:
                        C_links_to.append(link[link.index(C)^1])
            D = C_links_to[C_links_to.index(B)^1]

            D_links_to = []
            for link in heptagonal_link_list:
                if D in link:
                    if link[link.index(D)^1] not in D_links_to:
                        D_links_to.append(link[link.index(D)^1])
            E = D_links_to[D_links_to.index(C)^1]

            E_links_to = []
            for link in heptagonal_link_list:
                if E in link:
                    if link[link.index(E)^1] not in E_links_to:
                        E_links_to.append(link[link.index(E)^1])
            F = E_links_to[E_links_to.index(D)^1]

            add(A,C,'heptAC')
            add(A,D,'heptAD')
            add(G,D,'heptGD')
            add(G,E,'heptGE')

    return links

_lattice_triangle_cache = {}

def lattice_triangles(projection_name,cell_count,point_count):
    key = (projection_name,cell_count,point_count)
    if key in _lattice_triangle_cache:
        return _lattice_triangle_cache[key]

    links = lattice_mesh_links(projection_name,cell_count,point_count)
    if links is None:
        return None

    neighbours = [set() for _ in range(point_count)]
    for a,b,_kind in links:
        neighbours[a].add(b)
        neighbours[b].add(a)

    triangles = []
    for a in range(point_count):
        near_a = neighbours[a]
        for b in near_a:
            if b<=a:
                continue
            for c in near_a&neighbours[b]:
                if c>b:
                    triangles.append((a,b,c))

    triangles = np.asarray(triangles,dtype=np.int64)
    _lattice_triangle_cache[key] = triangles
    return triangles

def lattice_triangulation(projection_name,points_complex,cell_count):
    triangles = lattice_triangles(projection_name,cell_count,points_complex.shape[0])
    if triangles is None or triangles.shape[0]==0:
        return None
    node_xy = np.stack([np.real(points_complex),np.imag(points_complex)],axis=1)
    return ConstrainedTriangulation(node_xy,triangles)

def mesh_links_to_file(filename,points,links):
    schema = {
            'geometry': 'LineString',
            'properties': {'kind': 'str:8','a': 'int','b': 'int'},
    }
    with fiona.open(filename, 'w', 'ESRI Shapefile', schema) as c:
            for a,b,kind in links:
                c.write({
                        'geometry': mapping(LineString([
                            (float(np.real(points[a])),float(np.imag(points[a]))),
                            (float(np.real(points[b])),float(np.imag(points[b])))])),
                        'properties': {'kind': kind,'a': int(a),'b': int(b)},
                })

def read_in_file(filename,transform,clip_poly=None,clip_forward=None,clip_inverse=None,progress_label=None):
    multipolys = []
    with fiona.open(filename, layer=fiona.listlayers(filename)[0]) as layer:
        total_features = len(layer) if progress_label else 0
        start_time = time.time()
        for feature_i,feature in enumerate(layer):
            print(feature_i)
            if progress_label:
                progress_bar(progress_label,feature_i/max(total_features,1),start_time)
            geom = shape(feature['geometry'])
            if clip_poly is not None:
                if clip_forward is not None:
                    geom = clip_inverse(clip_forward(geom).intersection(clip_poly))
                else:
                    geom = geom.intersection(clip_poly)
                if geom.is_empty:
                    continue
            multipolys.append([])
            polys = geom.geoms if geom.geom_type in ('MultiPolygon','GeometryCollection') else [geom]
            for poly in polys:
                if poly.geom_type!='Polygon' or poly.is_empty:
                    continue
                ext_coords = poly.exterior.coords
                if len(ext_coords)==0:
                    continue
                poly_ext = transform(np.array(ext_coords)[:,0]+1j*np.array(ext_coords)[:,1])
                poly_ext_out = np.stack([np.real(poly_ext),np.imag(poly_ext)],axis=1)
                poly_ext_out = poly_ext_out[np.isfinite(poly_ext_out).all(axis=1)]
                if poly_ext_out.shape[0]<3:
                    continue
                poly_int_out = []
                for b in list(poly.interiors):
                    poly_int = transform(np.array(b.coords)[:,0]+1j*np.array(b.coords)[:,1])
                    b_out = np.stack([np.real(poly_int),np.imag(poly_int)],axis=1)
                    b_out = b_out[np.isfinite(b_out).all(axis=1)]
                    if b_out.shape[0]>=3:
                        poly_int_out.append(b_out)
                multipolys[-1].append([poly_ext_out,poly_int_out])
            if not multipolys[-1]:
                multipolys.pop()
        if progress_label:
            progress_bar(progress_label,1.0,start_time,done=True)
    return multipolys

def _compute_rotation_phase(name):
    grid_sample_w = np.load(f'pickle/interpolation_points/{name}_W.npy')/D
    mapped_sample_points = np.load(f'pickle/interpolation_points/{name}_M.npy')/D

    with open(f"pickle/bubble_wrap/{name}_full_wrapping.txt","rb") as fh:
        loaded_pickle = pickle.load(fh)
    nemo_point = loaded_pickle[1]
    opposing_point = pi+np.conj(nemo_point)

    center_deg = (opposing_point.real+1j*opposing_point.imag)*180/pi
    north_eps_rad = 1e-4
    north_deg = (opposing_point.real+1j*(opposing_point.imag+north_eps_rad))*180/pi
    w_center = stereo_w(np.array([center_deg]),opposing_point)[0]
    w_north = stereo_w(np.array([north_deg]),opposing_point)[0]

    node_xy_raw = np.stack([np.real(grid_sample_w),np.imag(grid_sample_w)],axis=1)
    _,knn = cKDTree(node_xy_raw).query([w_center.real,w_center.imag],k=10)
    Wc = grid_sample_w[knn]-w_center
    Mk = mapped_sample_points[knn]
    Asys = np.stack([np.ones_like(Wc),Wc,np.conj(Wc)],axis=1)
    sol,_,_,_ = np.linalg.lstsq(Asys,Mk,rcond=None)
    a,b = sol[1],sol[2]
    north_dir_m = a*(w_north-w_center)+b*np.conj(w_north-w_center)
    rotation_phase = pi/2-np.angle(north_dir_m)
    return rotation_phase,opposing_point,loaded_pickle

def _build_optimized_mesh(name):
    out_path = f'pickle/interpolation_points/{name}_optimized_M.npy'
    if os.path.isfile(out_path):
        return out_path

    from scipy.optimize import minimize as _minimize

    with open(f'pickle/matlab_saves/{name}_centers_file.txt','rb') as fh:
        centers = pickle.load(fh)
    with open(f'pickle/matlab_saves/{name}_centers_M.txt','rb') as fh:
        centers_M = pickle.load(fh)

    N = len(centers)
    w_pos = np.array([complex(np.ravel(c[1])[0]) for c in centers])
    w_anchor = np.array([complex(np.ravel(c[2])[0]) for c in centers])
    w_right = np.array([complex(np.ravel(c[11])[0]) for c in centers])
    w_left = np.array([complex(np.ravel(c[12])[0]) for c in centers])
    reach = np.abs(w_right-w_pos)

    w_tree = cKDTree(np.stack([w_pos.real,w_pos.imag],axis=1))
    _,parent = w_tree.query(np.stack([w_anchor.real,w_anchor.imag],axis=1))
    parent = parent.astype(int)

    m_pos = np.array([complex(np.ravel(c[1])[0]) for c in centers_M])
    rotation_phase,_,_ = _compute_rotation_phase(name)
    pos_orig = m_pos/D*np.exp(1j*rotation_phase)

    def find_RL(i):
        dR,kR = w_tree.query([w_right[i].real,w_right[i].imag])
        dL,kL = w_tree.query([w_left[i].real,w_left[i].imag])
        R = int(kR) if (kR!=i and dR<=0.01*reach[i]) else -1
        L = int(kL) if (kL!=i and dL<=0.01*reach[i]) else -1
        return R,L
    R_of = np.full(N,-1,dtype=int); L_of = np.full(N,-1,dtype=int)
    for i in range(N):
        R_of[i],L_of[i] = find_RL(i)

    def neighbors3(x):
        H = -1 if x==0 else int(parent[x])
        return H,int(R_of[x]),int(L_of[x])
    def others_excluding(x,p):
        H,R,L = neighbors3(x)
        return [n for n in (H,R,L) if n>=0 and n!=p]

    full_P = []
    neighborhood_of = {}
    for P in range(1,N):
        H,R,L = neighbors3(P)
        if R<0 or L<0:
            continue
        if H==0 or R==0 or L==0:
            continue
        ok = True
        nb = []
        for X in (H,R,L):
            oth = others_excluding(X,P)
            if len(oth)!=2:
                ok = False
                break
            nb.extend([X]+oth)
        if not ok:
            continue
        full_P.append(P)
        neighborhood_of[P] = nb

    full_P_arr = np.array(full_P)
    nb_arr = np.array([neighborhood_of[P] for P in full_P])
    V0 = pos_orig[nb_arr]
    Mcoef,_,_,_ = np.linalg.lstsq(V0,pos_orig[full_P_arr],rcond=None)

    tree_child = np.arange(1,N)
    tree_par = parent[1:]
    orig_edge = pos_orig[tree_child]-pos_orig[tree_par]

    def unpack(x):
        pos = pos_orig.copy()
        pos[1:] = x[:N-1]+1j*x[N-1:]
        return pos
    def pack(pos):
        return np.concatenate([pos[1:].real,pos[1:].imag])

    W2 = 0.001

    def loss_and_grad(x):
        pos = unpack(x)
        G = np.zeros(N,dtype=complex)

        v = pos[nb_arr]
        resid = v@Mcoef-pos[full_P_arr]
        loss = float(np.sum(np.abs(resid)**2))
        coeff = 2*resid
        for k in range(9):
            np.add.at(G,nb_arr[:,k],np.conj(Mcoef[k])*coeff)
        np.add.at(G,full_P_arr,-coeff)

        new_edge = pos[tree_child]-pos[tree_par]
        ratio = new_edge/orig_edge
        r = np.abs(ratio); theta = np.angle(ratio)
        loss += W2*float(np.sum((r-1)**2+theta**2))
        G_ratio = (ratio/np.maximum(r,1e-300))*(2*(r-1)+1j*2*theta/np.maximum(r,1e-300))
        G_ratio *= W2
        c = 1.0/orig_edge
        G_newedge = np.conj(c)*G_ratio
        np.add.at(G,tree_child,G_newedge)
        np.add.at(G,tree_par,-G_newedge)

        grad = np.concatenate([G[1:].real,G[1:].imag])
        return loss,grad

    MAXITER = 15000
    _opt_start = time.time()
    _opt_iter = [0]
    def _opt_callback(xk):
        _opt_iter[0] += 1
        progress_bar('Optimizing lattice positions',_opt_iter[0]/MAXITER,_opt_start)

    x0 = pack(pos_orig)
    res = _minimize(loss_and_grad,x0,jac=True,method='L-BFGS-B',
        options={'maxiter':MAXITER,'maxfun':10*MAXITER,'ftol':0.0,'gtol':1e-16},
        callback=_opt_callback)
    progress_bar('Optimizing lattice positions',1.0,_opt_start,done=True)

    pos_opt = unpack(res.x)
    delta_out = pos_opt-pos_orig
    write_shifted_mesh(name,N,w_pos,delta_out,rotation_phase,'_optimized')
    return out_path

ANGLE_ERROR_BUCKETS_DEG = [0,1,5,15,45,180]
LENGTH_ERROR_BUCKETS = [0,0.05,0.2,0.5,1.0,float('inf')]

def _bucket_percentages(values,edges):
    values = np.asarray(values)
    n = values.shape[0]
    out = []
    for i in range(len(edges)-1):
        lo,hi = edges[i],edges[i+1]
        if hi==float('inf'):
            mask = values>=lo
            label = f'>= {lo:g}'
        else:
            mask = (values>=lo)&(values<hi)
            label = f'{lo:g} - {hi:g}'
        count = int(np.sum(mask))
        pct = 100*count/n if n else 0.0
        out.append((label,count,pct))
    return out

def stereo_w(equi_deg,opposing_point):
    equi_rad = equi_deg*pi/180
    k = D/(1+sin(opposing_point.imag)*np.sin(np.imag(equi_rad))+cos(opposing_point.imag)*np.cos(np.imag(equi_rad))*np.cos(np.real(equi_rad)-opposing_point.real))
    x = k*np.cos(np.imag(equi_rad))*np.sin(np.real(equi_rad)-opposing_point.real)
    y = k*(cos(opposing_point.imag)*np.sin(np.imag(equi_rad))-sin(opposing_point.imag)*np.cos(np.imag(equi_rad))*np.cos(np.real(equi_rad)-opposing_point.real))
    return (x+y*1j)/D

def stereo_to_equi_deg(w,opposing_point):
    x = np.real(w)*D
    y = np.imag(w)*D
    rho = np.sqrt(x**2+y**2)
    c = 2*np.arctan2(rho,D)
    op_lat,op_lon = opposing_point.imag,opposing_point.real
    rho_safe = np.where(rho>0,rho,1.0)
    lat = np.arcsin(np.clip(np.cos(c)*np.sin(op_lat)+y*np.sin(c)*np.cos(op_lat)/rho_safe,-1,1))
    lon = op_lon+np.arctan2(x*np.sin(c),rho_safe*np.cos(op_lat)*np.cos(c)-y*np.sin(op_lat)*np.sin(c))
    at_center = rho==0
    lat = np.where(at_center,op_lat,lat)
    lon = np.where(at_center,op_lon,lon)
    lon = (lon+pi)%(2*pi)-pi
    return np.degrees(lon)+1j*np.degrees(lat)

EARTH_RADIUS = D/2

def great_circle_distance(lon_a,lat_a,lon_b,lat_b):
    lat_a,lat_b = np.radians(lat_a),np.radians(lat_b)
    half_chord = (np.sin((lat_b-lat_a)/2)**2
        +np.cos(lat_a)*np.cos(lat_b)*np.sin(np.radians(lon_b-lon_a)/2)**2)
    return 2*EARTH_RADIUS*np.arcsin(np.sqrt(np.clip(half_chord,0.0,1.0)))

def triangle_affine_coefficients(triangulation,source,target):
    verts = triangulation.simplices
    corners = source[verts]
    system = np.empty((verts.shape[0],3,3),dtype=complex)
    system[:,:,0] = corners; system[:,:,1] = np.conj(corners); system[:,:,2] = 1.0
    system_inv = np.zeros_like(system)
    nondegenerate = np.abs(np.linalg.det(system))>1e-300
    system_inv[nondegenerate] = np.linalg.inv(system[nondegenerate])
    mapped = target[verts]
    return (np.sum(system_inv[:,0,:]*mapped,axis=1),
        np.sum(system_inv[:,1,:]*mapped,axis=1),
        np.sum(system_inv[:,2,:]*mapped,axis=1),
        nondegenerate)

def map_through_nearest_triangle(triangulation,source,target,query,skip=None):
    coef_a,coef_b,coef_c,nondegenerate = triangle_affine_coefficients(triangulation,source,target)
    usable = nondegenerate if skip is None else (nondegenerate&~skip)
    usable_idx = np.where(usable)[0]
    corners = source[triangulation.simplices]
    centroids = np.stack([np.real(corners).mean(axis=1),np.imag(corners).mean(axis=1)],axis=1)
    _,nearest = cKDTree(centroids[usable_idx]).query(
        np.stack([np.real(query),np.imag(query)],axis=1))
    chosen = usable_idx[nearest]
    return coef_a[chosen]*query+coef_b[chosen]*np.conj(query)+coef_c[chosen]

def reproj_boundary_ring_m(ring_equi_deg,opposing_point,boundary_W,boundary_M,tol=1e-6):
    if boundary_W.shape[0]==0:
        return None
    boundary_tree = cKDTree(np.stack([boundary_W.real,boundary_W.imag],axis=1))
    stereo = stereo_w(ring_equi_deg,opposing_point)
    n = stereo.shape[0]
    query = np.stack([stereo.real,stereo.imag],axis=1)
    dist,idx = boundary_tree.query(query)
    exact = dist<tol
    exact_indices = np.where(exact)[0]
    if exact_indices.shape[0]==0:
        return None
    out = np.empty(n,dtype=complex)
    out[exact] = boundary_M[idx[exact]]
    for i in np.where(~exact)[0]:
        pos = np.searchsorted(exact_indices,i)
        right = exact_indices[pos%exact_indices.shape[0]]
        left = exact_indices[(pos-1)%exact_indices.shape[0]]
        if left==right:
            out[i] = out[left]
            continue
        span = (right-left)%n
        pos_in_span = (i-left)%n
        frac = pos_in_span/span
        out[i] = out[left]*(1-frac)+out[right]*frac
    return out

def _load_rotated_mesh(name):
    _build_optimized_mesh(name)
    grid_sample_w = np.load(f'pickle/interpolation_points/{name}_W.npy')/D
    mapped_sample_points = np.load(f'pickle/interpolation_points/{name}_optimized_M.npy')/D

    blocks_path = f'pickle/interpolation_points/{name}_blocks.npy'
    if os.path.isfile(blocks_path):
        cell_count,boundary_count = (int(v) for v in np.load(blocks_path))
    else:
        with open(f'pickle/matlab_saves/{name}_centers_file.txt','rb') as fh:
            cell_count = len(pickle.load(fh))
        boundary_count = grid_sample_w.shape[0]-13*cell_count

    is_mesh = np.ones(grid_sample_w.shape[0],dtype=bool)
    is_mesh[cell_count:cell_count+boundary_count] = False
    grid_sample_w = grid_sample_w[is_mesh]
    mapped_sample_points = mapped_sample_points[is_mesh]

    boundary_count = 0

    boundary_W = grid_sample_w[cell_count:cell_count + boundary_count]
    boundary_M = mapped_sample_points[cell_count:cell_count + boundary_count]

    rotation_phase,opposing_point,loaded_pickle = _compute_rotation_phase(name)
    mapped_sample_points = mapped_sample_points*np.exp(1j*rotation_phase)

    boundary_M = boundary_M*np.exp(1j*rotation_phase)

    return (grid_sample_w,mapped_sample_points,loaded_pickle,opposing_point,
        cell_count,boundary_count,boundary_W,boundary_M,rotation_phase)

def run_interpolation_algorithm(projection_name,vector_file_to_project=None,outline_out=None,projection_out=None,mesh_out=None,mesh_tree_out=None,mesh_tree_links_out=None,mesh_links_out=None,extra_outputs=True,graticule_out=None,graticule_step_deg=None,angle_outlier_deg=None):
    out_name = projection_name
    outline_out = outline_out or f'maps_projected/outline_{out_name}.shp'
    projection_out = projection_out or f'maps_projected/{out_name}.shp'
    mesh_out = mesh_out or f'maps_projected/mesh_{out_name}.shp'
    mesh_tree_out = mesh_tree_out or f'maps_projected/mesh_tree_{out_name}.shp'
    mesh_tree_links_out = mesh_tree_links_out or f'maps_projected/mesh_tree_links_{out_name}.shp'
    mesh_links_out = mesh_links_out or f'maps_projected/mesh_links_{out_name}.shp'

    (grid_sample_w,mapped_sample_points,loaded_pickle,opposing_point,
        cell_count,boundary_count,boundary_W,boundary_M,_rotation_phase) = _load_rotated_mesh(projection_name)

    tree_names,tree_depths = lattice_tree_names(projection_name,cell_count)
    mesh_names,mesh_depths,mesh_arms = (None,None,None) if tree_names is None else \
        mesh_point_names(tree_names,tree_depths,cell_count,mapped_sample_points.shape[0])

    if mesh_names is None:
        points_to_file(mesh_out,mapped_sample_points)
    else:
        named_points_to_file(mesh_out,mapped_sample_points,mesh_names,mesh_depths,mesh_arms)

    if tree_names is None:
        points_to_file(mesh_tree_out,mapped_sample_points[:cell_count])
    else:
        named_points_to_file(mesh_tree_out,mapped_sample_points[:cell_count],tree_names,tree_depths)

    tree_links = lattice_tree_links(projection_name,cell_count)
    if tree_links is not None:
        tree_links_to_file(mesh_tree_links_out,mapped_sample_points[:cell_count],tree_links)

    graph_links = lattice_mesh_links(projection_name,cell_count,mapped_sample_points.shape[0])
    if graph_links is not None:
        mesh_links_to_file(mesh_links_out,mapped_sample_points,graph_links)

    boundary_tree = cKDTree(np.stack([np.real(boundary_W),np.imag(boundary_W)],axis=1)) if boundary_W.shape[0]>0 else None
    BOUNDARY_EXACT_TOL = 1e-6

    node_xy = np.stack([np.real(grid_sample_w),np.imag(grid_sample_w)],axis=1)

    triangulation = lattice_triangulation(projection_name,grid_sample_w,cell_count)
    if triangulation is None:
        triangulation = Delaunay(node_xy)

    node_tree = cKDTree(node_xy)
    node_spacing = node_tree.query(node_xy,k=2)[0][:,1]

    _verts = triangulation.simplices
    _Wt = grid_sample_w[_verts]
    _Sys = np.empty((_verts.shape[0],3,3),dtype=complex)
    _Sys[:,:,0] = _Wt; _Sys[:,:,1] = np.conj(_Wt); _Sys[:,:,2] = 1.0
    _Sys_inv = np.zeros_like(_Sys)
    _nondegenerate = np.abs(np.linalg.det(_Sys))>1e-300
    _Sys_inv[_nondegenerate] = np.linalg.inv(_Sys[_nondegenerate])
    _A,_Bc,_Cc = _Sys_inv[:,0,:],_Sys_inv[:,1,:],_Sys_inv[:,2,:]
    _Mt = mapped_sample_points[_verts]
    _a = np.sum(_A*_Mt,axis=1); _b = np.sum(_Bc*_Mt,axis=1); _c = np.sum(_Cc*_Mt,axis=1)

    if isinstance(triangulation,ConstrainedTriangulation):
        real_rbf = triangulation.linear_interpolator(np.real(mapped_sample_points))
        imag_rbf = triangulation.linear_interpolator(np.imag(mapped_sample_points))
    else:
        real_rbf = sp.interpolate.LinearNDInterpolator(triangulation,np.real(mapped_sample_points))
        imag_rbf = sp.interpolate.LinearNDInterpolator(triangulation,np.imag(mapped_sample_points))

    _tri_centroids = node_xy[triangulation.simplices].mean(axis=1)
    _tri_tree = cKDTree(_tri_centroids)

    def nearest_triangle_affine(query):
        _,_tri = _tri_tree.query(query)
        _Wq = query[:,0]+1j*query[:,1]
        return _a[_tri]*_Wq+_b[_tri]*np.conj(_Wq)+_c[_tri]

    def sliver_reject(query,simplex):
        if LOCALITY<=0:
            return np.zeros(query.shape[0],dtype=bool)
        _,anchor_i = node_tree.query(query)
        corners = triangulation.simplices[np.maximum(simplex,0)]
        anchor = node_xy[anchor_i]
        reach = np.linalg.norm(node_xy[corners]-anchor[:,None,:],axis=2).max(axis=1)
        return reach>LOCALITY*node_spacing[anchor_i]

    LOCAL_LINEAR_K = 10
    GLOBAL_MEDIAN_SPACING = float(np.median(node_spacing))
    GLOBAL_M_SCALE = float(np.abs(mapped_sample_points-mapped_sample_points.mean()).std())
    LOCAL_LINEAR_RESID_TOL = 0.02

    def local_linear_fallback(query):
        dist,idx = node_tree.query(query,k=LOCAL_LINEAR_K)
        W_k = node_xy[idx]
        M_k = mapped_sample_points[idx]
        W0 = W_k.mean(axis=1)
        Wc = (W_k[:,:,0]-W0[:,0:1])+1j*(W_k[:,:,1]-W0[:,1:2])
        ones = np.ones_like(Wc)
        qc = (query[:,0]-W0[:,0])+1j*(query[:,1]-W0[:,1])
        est = np.zeros(query.shape[0],dtype=complex)
        untrustworthy = np.zeros(query.shape[0],dtype=bool)
        for i in range(query.shape[0]):
            A = np.stack([ones[i],Wc[i]],axis=1)
            sol,_,_,_ = np.linalg.lstsq(A,M_k[i],rcond=None)
            est[i] = sol[0]+sol[1]*qc[i]
            resid = np.abs(M_k[i]-A@sol).max()
            untrustworthy[i] = (resid>LOCAL_LINEAR_RESID_TOL*GLOBAL_M_SCALE) or (dist[i].max()>LOCALITY*GLOBAL_MEDIAN_SPACING)
        return est,untrustworthy

    def equi_to_stereo(equi):
        equi = equi*pi/180
        k = D/(1+sin(opposing_point.imag)*np.sin(np.imag(equi))+cos(opposing_point.imag)*np.cos(np.imag(equi))*np.cos(np.real(equi)-opposing_point.real))
        x = k*np.cos(np.imag(equi))*np.sin(np.real(equi)-opposing_point.real)
        y = k*(cos(opposing_point.imag)*np.sin(np.imag(equi))-sin(opposing_point.imag)*np.cos(np.imag(equi))*np.cos(np.real(equi)-opposing_point.real))

        x[np.isnan(x)] = 0
        y[np.isnan(y)] = 0
        x = np.minimum(x,10000000000000)
        x = np.maximum(x,-10000000000000)
        y = np.minimum(y,10000000000000)
        y = np.maximum(y,-10000000000000)

        return (x+y*1j)

    def reproj_boundary_ring(ring_equi):
        if boundary_tree is None:
            return None
        stereo = equi_to_stereo(ring_equi)
        n = stereo.shape[0]
        query = np.stack([np.real(stereo/D),np.imag(stereo/D)],axis=1)
        dist,idx = boundary_tree.query(query)
        exact = dist<BOUNDARY_EXACT_TOL
        exact_indices = np.where(exact)[0]
        if exact_indices.shape[0]==0:
            return None

        out = np.empty(n,dtype=complex)
        out[exact] = boundary_M[idx[exact]]
        for i in np.where(~exact)[0]:
            pos = np.searchsorted(exact_indices,i)
            right = exact_indices[pos%exact_indices.shape[0]]
            left = exact_indices[(pos-1)%exact_indices.shape[0]]
            if left==right:
                out[i] = out[left]
                continue
            span = (right-left)%n
            pos_in_span = (i-left)%n
            frac = pos_in_span/span
            out[i] = out[left]*(1-frac)+out[right]*frac
        return out

    def reproj(equi,stereoproj=True):
        if stereoproj:
            stereo = equi_to_stereo(equi)
        else:
            stereo = equi

        shape_ = np.shape(stereo)
        query = np.stack([np.real(stereo/D).ravel(),np.imag(stereo/D).ravel()],axis=1)

        out_real = np.empty(query.shape[0])
        out_imag = np.empty(query.shape[0])
        reject = np.zeros(query.shape[0],dtype=bool)

        exact = np.zeros(query.shape[0],dtype=bool)
        if boundary_tree is not None:
            dist,idx = boundary_tree.query(query)
            exact = dist<BOUNDARY_EXACT_TOL
            out_real[exact] = np.real(boundary_M[idx[exact]])
            out_imag[exact] = np.imag(boundary_M[idx[exact]])

        rest = ~exact
        if np.any(rest):
            rest_q = query[rest]
            rest_real = real_rbf(rest_q[:,0],rest_q[:,1])
            rest_imag = imag_rbf(rest_q[:,0],rest_q[:,1])
            rest_simplex = triangulation.find_simplex(rest_q)
            outside_hull = (rest_simplex<0)|np.isnan(rest_real)|np.isnan(rest_imag)
            sliver = np.zeros(rest_q.shape[0],dtype=bool)
            sliver[~outside_hull] = sliver_reject(rest_q[~outside_hull],rest_simplex[~outside_hull])

            rest_reject = outside_hull.copy()
            if np.any(outside_hull):
                affine_est = nearest_triangle_affine(rest_q[outside_hull])
                rest_real[outside_hull] = np.real(affine_est)
                rest_imag[outside_hull] = np.imag(affine_est)
            if np.any(sliver):
                fallback_est,still_bad = local_linear_fallback(rest_q[sliver])
                rest_real[sliver] = np.real(fallback_est)
                rest_imag[sliver] = np.imag(fallback_est)
                rest_reject[sliver] = still_bad
            out_real[rest] = rest_real
            out_imag[rest] = rest_imag
            reject[rest] = rest_reject

        reproj.rejected += int(reject.sum())
        reproj.total += reject.shape[0]
        reproj.last_reject = reject.reshape(shape_)
        return (out_real+1j*out_imag).reshape(shape_)
    reproj.rejected = 0
    reproj.total = 0
    reproj.last_reject = None

    def measure_distortion(filename,clip_poly_,progress_label=None):
        _,node_knn_idx = node_tree.query(node_xy,k=LOCAL_LINEAR_K)
        node_local_scale = np.zeros(node_xy.shape[0])
        for i in range(node_xy.shape[0]):
            W_k = node_xy[node_knn_idx[i]]
            M_k = mapped_sample_points[node_knn_idx[i]]
            W0 = node_xy[i]
            Wc = (W_k[:,0]-W0[0])+1j*(W_k[:,1]-W0[1])
            A = np.stack([np.ones_like(Wc),Wc],axis=1)
            sol,_,_,_ = np.linalg.lstsq(A,M_k,rcond=None)
            node_local_scale[i] = np.abs(sol[1])

        angle_errors = []
        length_errors = []
        angle_outliers = []
        with fiona.open(filename, layer=fiona.listlayers(filename)[0]) as layer:
            total_features = len(layer) if progress_label else 0
            start_time = time.time()
            for feature_i,feature in enumerate(layer):
                if progress_label:
                    progress_bar(progress_label,feature_i/max(total_features,1),start_time)
                geom = _clip_inverse(_clip_forward(shape(feature['geometry'])).intersection(clip_poly_))
                if geom.is_empty:
                    continue
                polys = geom.geoms if geom.geom_type in ('MultiPolygon','GeometryCollection') else [geom]
                for p in polys:
                    if p.geom_type!='Polygon' or p.is_empty:
                        continue
                    for ring in [p.exterior]+list(p.interiors):
                        coords = np.array(ring.coords)
                        n = coords.shape[0]
                        if n<4:
                            continue
                        equi = coords[:,0]+1j*coords[:,1]
                        orig = equi_to_stereo(equi)/D
                        proj = reproj(equi)

                        seg_lines = shapely.linestrings(np.stack([coords[:-1],coords[1:]],axis=1))
                        within_seg = shapely.intersects(clip_poly_,_clip_forward(seg_lines,repair=False))

                        m = n-1
                        idxs = np.arange(m)
                        prev_i = (idxs-1)%m
                        next_i = (idxs+1)%m

                        v,vp,vn = orig[:m],orig[prev_i],orig[next_i]
                        in1,in2 = vp-v,vn-v
                        valid_ang = (np.abs(in1)>0)&(np.abs(in2)>0)
                        pv,pvp,pvn = proj[:m],proj[prev_i],proj[next_i]
                        pin1,pin2 = pvp-pv,pvn-pv
                        valid_ang &= (np.abs(pin1)>0)&(np.abs(pin2)>0)

                        orig_ang = np.zeros(m); proj_ang = np.zeros(m)
                        orig_ang[valid_ang] = np.angle(in2[valid_ang]/in1[valid_ang])
                        proj_ang[valid_ang] = np.angle(pin2[valid_ang]/pin1[valid_ang])
                        ang_err = np.degrees(np.abs(np.angle(np.exp(1j*(proj_ang-orig_ang)))))

                        vertex_within = within_seg|within_seg[prev_i]
                        keep = valid_ang&vertex_within
                        angle_errors.extend(ang_err[keep].tolist())
                        if angle_outlier_deg is not None:
                            bad = keep&(ang_err>angle_outlier_deg)
                            for _i in np.where(bad)[0]:
                                angle_outliers.append((feature_i,float(np.real(equi[_i])),float(np.imag(equi[_i])),float(ang_err[_i])))

                        orig_vec,proj_vec = orig[1:]-orig[:-1],proj[1:]-proj[:-1]
                        orig_len = np.abs(orig_vec)
                        valid_len = orig_len>1e-15
                        mids = np.stack([(np.real(orig[:-1])+np.real(orig[1:]))/2,
                            (np.imag(orig[:-1])+np.imag(orig[1:]))/2],axis=1)
                        _,nearest = node_tree.query(mids)
                        expected_scale = node_local_scale[nearest]
                        valid_len &= expected_scale>0
                        rel_err = np.zeros(n-1)
                        rel_err[valid_len] = np.abs(np.abs(proj_vec[valid_len])/orig_len[valid_len]/expected_scale[valid_len]-1)
                        length_errors.extend(rel_err[within_seg&valid_len].tolist())
            if progress_label:
                progress_bar(progress_label,1.0,start_time,done=True)
        return np.array(angle_errors),np.array(length_errors),angle_outliers

    ##outline files
    boundary_ring_equi = loaded_pickle[0]*(180/pi)

    poly_ext = equi_to_stereo(boundary_ring_equi)
    poly_ext_out = np.stack([np.real(poly_ext),np.imag(poly_ext)],axis=1)
    multipolys = [[[poly_ext_out,[]]]]
    original_outline = Polygon(poly_ext_out)
    if extra_outputs:
        multi_to_file(f'maps_projected/original_outline_{out_name}.shp',multipolys)

    clip_poly = make_valid(Polygon(poly_ext_out))

    def _clip_forward(geom,repair=True):
        def _fwd(coords):
            z = equi_to_stereo(coords[:,0]+1j*coords[:,1])
            return np.stack([np.real(z),np.imag(z)],axis=1)
        out = shapely.transform(geom,_fwd)
        return out.buffer(0) if repair else out

    def _clip_inverse(geom):
        def _inv(coords):
            equi = stereo_to_equi_deg((coords[:,0]+1j*coords[:,1])/D,opposing_point)
            return np.stack([np.real(equi),np.imag(equi)],axis=1)
        return shapely.transform(geom,_inv)

    poly_ext = reproj_boundary_ring(boundary_ring_equi)
    if poly_ext is None:
        poly_ext = reproj(boundary_ring_equi)
    poly_ext_out = np.stack([np.real(poly_ext),np.imag(poly_ext)],axis=1)
    multipolys = [[[poly_ext_out,[]]]]
    multi_to_file(outline_out,multipolys)

    def _graticule_lines(step_deg):
        long_lines = np.arange(-180,180,step_deg)[:,None]+1j*np.arange(-90,90.1,0.1)[None,:]
        lat_lines = (np.arange(-180,180.1,0.1)[:,None]+1j*np.arange(-85,90,step_deg)[None,:]).T
        return list(long_lines)+list(lat_lines)

    if extra_outputs:
        line_list = []
        for gridline in _graticule_lines(5):
            line_comp = equi_to_stereo(gridline)
            line_split = np.stack([np.real(line_comp),np.imag(line_comp)],axis=1)
            line_list.append(line_split)
        lines_to_file(f'maps_projected/original_gridline_{out_name}.shp',line_list)

    if extra_outputs or graticule_step_deg is not None:
        step_deg = graticule_step_deg if graticule_step_deg is not None else 5
        out_path = graticule_out or f'maps_projected/gridline_{out_name}.shp'
        line_list = []
        gridlines = _graticule_lines(step_deg)
        _graticule_start = time.time()
        for _gridline_i,gridline in enumerate(gridlines):
            progress_bar('Projecting graticule',_gridline_i/max(len(gridlines),1),_graticule_start)
            line_comp = equi_to_stereo(gridline)
            line_split = np.stack([np.real(line_comp),np.imag(line_comp)],axis=1)
            intersected = LineString(line_split).intersection(make_valid(original_outline))
            pieces = intersected.geoms if intersected.geom_type=='MultiLineString' else [intersected]
            for piece in pieces:
                if piece.geom_type!='LineString' or len(piece.coords)==0:
                    continue
                reprojed = reproj(np.array(piece.coords)[:,0]+1j*np.array(piece.coords)[:,1],False)
                line_list.append(np.stack([np.real(reprojed),np.imag(reprojed)],axis=1))
        progress_bar('Projecting graticule',1.0,_graticule_start,done=True)
        lines_to_file(out_path,line_list)
        graticule_result_path = out_path
    else:
        graticule_result_path = None

    if vector_file_to_project is None:
        return {'outline':outline_out,'mesh':mesh_out,'mesh_tree':mesh_tree_out,'mesh_tree_links':mesh_tree_links_out,'mesh_links':mesh_links_out,'graticule':graticule_result_path}

    if extra_outputs:
        multipolys = read_in_file(vector_file_to_project,equi_to_stereo,clip_poly=clip_poly,clip_forward=_clip_forward,clip_inverse=_clip_inverse,progress_label='Projecting to stereographic (original)')
        multi_to_file(f'maps_projected/original_{out_name}.shp',multipolys)

    multipolys = read_in_file(vector_file_to_project,reproj,clip_poly=clip_poly,clip_forward=_clip_forward,clip_inverse=_clip_inverse,progress_label='Projecting to final position')
    multi_to_file(projection_out,multipolys)

    print(f'locality {LOCALITY}: {reproj.rejected} of {reproj.total} projected points fell back to nearest-point interpolation')

    _saved_rejected,_saved_total = reproj.rejected,reproj.total
    angle_errors,length_errors,angle_outliers = measure_distortion(vector_file_to_project,clip_poly,progress_label='Measuring distortion')
    reproj.rejected,reproj.total = _saved_rejected,_saved_total

    angle_error_table = _bucket_percentages(angle_errors,ANGLE_ERROR_BUCKETS_DEG)
    length_error_table = _bucket_percentages(length_errors,LENGTH_ERROR_BUCKETS)

    return {'outline':outline_out,'projection':projection_out,'mesh':mesh_out,'mesh_tree':mesh_tree_out,'mesh_tree_links':mesh_tree_links_out,'mesh_links':mesh_links_out,'rejected':reproj.rejected,'total':reproj.total,
        'angle_error_table':angle_error_table,'angle_error_count':int(angle_errors.shape[0]),
        'length_error_table':length_error_table,'length_error_count':int(length_errors.shape[0]),
        'angle_outliers':angle_outliers}

def boundary_ring_in_m(projection_name,grid_sample_w,mapped_sample_points,
        loaded_pickle,opposing_point,cell_count):
    ring_equi = np.real(loaded_pickle[0])*180/pi+1j*np.imag(loaded_pickle[0])*180/pi
    ring_w = stereo_w(ring_equi,opposing_point)

    tri_w = lattice_triangulation(projection_name,grid_sample_w,cell_count)
    if tri_w is None:
        tri_w = Delaunay(np.stack([grid_sample_w.real,grid_sample_w.imag],axis=1))
    if isinstance(tri_w,ConstrainedTriangulation):
        to_real = tri_w.linear_interpolator(mapped_sample_points.real)
        to_imag = tri_w.linear_interpolator(mapped_sample_points.imag)
    else:
        to_real = sp.interpolate.LinearNDInterpolator(tri_w,mapped_sample_points.real)
        to_imag = sp.interpolate.LinearNDInterpolator(tri_w,mapped_sample_points.imag)

    ring_m = (np.asarray(to_real(ring_w.real,ring_w.imag),dtype=float)
        +1j*np.asarray(to_imag(ring_w.real,ring_w.imag),dtype=float))
    off_mesh = ~np.isfinite(ring_m.real)|~np.isfinite(ring_m.imag)
    if np.any(off_mesh):
        ring_m[off_mesh] = map_through_nearest_triangle(tri_w,grid_sample_w,
            mapped_sample_points,ring_w[off_mesh])
    return ring_m

def neighbour_gap_on_globe(lon,lat,band=512):
    gap = np.zeros(lon.shape,dtype=np.float32)
    height = lon.shape[0]
    for start in range(0,height,band):
        stop = min(start+band,height)
        lo,la = lon[start:stop],lat[start:stop]
        across = great_circle_distance(lo[:,1:],la[:,1:],lo[:,:-1],la[:,:-1])
        np.maximum(gap[start:stop,1:],across,out=gap[start:stop,1:])
        np.maximum(gap[start:stop,:-1],across,out=gap[start:stop,:-1])
        del across
        if stop<height:
            down = great_circle_distance(lon[stop-1:stop+1][1:],lat[stop-1:stop+1][1:],
                lon[stop-1:stop+1][:-1],lat[stop-1:stop+1][:-1])
            np.maximum(gap[stop:stop+1],down,out=gap[stop:stop+1])
            np.maximum(gap[stop-1:stop],down,out=gap[stop-1:stop])
            del down
        if stop-start>1:
            inner = great_circle_distance(lo[1:],la[1:],lo[:-1],la[:-1])
            np.maximum(gap[start+1:stop],inner,out=gap[start+1:stop])
            np.maximum(gap[start:stop-1],inner,out=gap[start:stop-1])
            del inner
    return gap

def blank_folded_pixels(lon,lat,within_boundary,globe_radius,min_pixel_separation=4):
    height,width = lon.shape
    lat_rad,lon_rad = np.radians(lat),np.radians(lon)
    globe_xy = np.stack([np.cos(lat_rad)*np.cos(lon_rad),np.cos(lat_rad)*np.sin(lon_rad),
        np.sin(lat_rad)],axis=-1).reshape(-1,3).astype(np.float32)*np.float32(EARTH_RADIUS)
    del lat_rad,lon_rad
    chord = 2*EARTH_RADIUS*np.sin(globe_radius/(2*EARTH_RADIUS))

    to_boundary = ndimage.distance_transform_edt(~within_boundary).ravel().astype(np.float32)
    rows,cols = np.divmod(np.arange(height*width,dtype=np.int32),np.int32(width))
    inside_flat = within_boundary.ravel()

    order = np.lexsort((np.arange(height*width,dtype=np.int32),to_boundary)).astype(np.int32)
    alive = ~inside_flat.copy()
    blank = np.zeros(height*width,dtype=bool)
    separation_sq = min_pixel_separation*min_pixel_separation

    CHUNK = 4000
    position = 0
    while position<order.size and alive.any():
        live_idx = np.flatnonzero(alive).astype(np.int32)
        tree = cKDTree(np.asarray(globe_xy[live_idx],dtype=np.float64))
        live_at_build = live_idx.size
        while position<order.size:
            chunk = order[position:position+CHUNK]
            position += CHUNK
            found = tree.query_ball_point(np.asarray(globe_xy[chunk],dtype=np.float64),
                r=chord,workers=-1)
            for source,hits in zip(chunk,found):
                if not hits:
                    continue
                candidates = live_idx[hits]
                far_enough = ((rows[candidates].astype(np.int64)-int(rows[source]))**2
                    +(cols[candidates].astype(np.int64)-int(cols[source]))**2)>separation_sq
                candidates = candidates[far_enough&alive[candidates]
                    &(to_boundary[candidates]>to_boundary[source])]
                if candidates.size:
                    alive[candidates] = False
                    blank[candidates] = True
            del found
            if alive.sum()<live_at_build*0.5:
                break
    return blank.reshape(height,width)

def run_interpolation_algorithm_raster(projection_name,raster_path,out_tif=None,resolution=1200,margin_px=100):
    out_tif = out_tif or f'maps_projected/raster_{projection_name}.tif'
    import rasterio
    from rasterio.enums import ColorInterp

    (grid_sample_w,mapped_sample_points,loaded_pickle,opposing_point,
        cell_count,boundary_count,boundary_W,boundary_M,rotation_phase) = _load_rotated_mesh(projection_name)

    node_xy_m = np.stack([mapped_sample_points.real,mapped_sample_points.imag],axis=1)
    tri_m = lattice_triangulation(projection_name,mapped_sample_points,cell_count)
    if tri_m is None:
        tri_m = Delaunay(node_xy_m)
    if isinstance(tri_m,ConstrainedTriangulation):
        w_real_rbf = tri_m.linear_interpolator(grid_sample_w.real)
        w_imag_rbf = tri_m.linear_interpolator(grid_sample_w.imag)
    else:
        w_real_rbf = sp.interpolate.LinearNDInterpolator(tri_m,grid_sample_w.real)
        w_imag_rbf = sp.interpolate.LinearNDInterpolator(tri_m,grid_sample_w.imag)

    if isinstance(tri_m,ConstrainedTriangulation):
        sliver_tri = np.zeros(tri_m.simplices.shape[0],dtype=bool)
    else:
        node_spacing_m = cKDTree(node_xy_m).query(node_xy_m,k=2)[0][:,1]
        _corners = node_xy_m[tri_m.simplices]
        _edges = np.stack([np.linalg.norm(_corners[:,a]-_corners[:,b],axis=1) for a,b in ((0,1),(1,2),(2,0))],axis=1)
        sliver_tri = _edges.max(axis=1)>LOCALITY*node_spacing_m[tri_m.simplices].mean(axis=1)

    minx,maxx = float(node_xy_m[:,0].min()),float(node_xy_m[:,0].max())
    miny,maxy = float(node_xy_m[:,1].min()),float(node_xy_m[:,1].max())
    aspect = (maxx-minx)/(maxy-miny)
    nx,ny = (resolution,max(1,int(resolution/aspect))) if aspect>=1 else (max(1,int(resolution*aspect)),resolution)

    if margin_px:
        step_x = (maxx-minx)/max(nx-1,1)
        step_y = (maxy-miny)/max(ny-1,1)
        minx -= margin_px*step_x; maxx += margin_px*step_x
        miny -= margin_px*step_y; maxy += margin_px*step_y
        nx += 2*margin_px; ny += 2*margin_px

    xs = np.linspace(minx,maxx,nx)
    ys = np.linspace(maxy,miny,ny)
    gx,gy = np.meshgrid(xs,ys)

    simplex = tri_m.find_simplex(np.stack([gx.ravel(),gy.ravel()],axis=1)).reshape(gx.shape)
    on_mesh = simplex>=0
    rejected = on_mesh&sliver_tri[np.maximum(simplex,0)]
    inside = on_mesh&~rejected

    w_out_real = np.full(gx.shape,np.nan)
    w_out_imag = np.full(gx.shape,np.nan)
    w_out_real[inside] = w_real_rbf(gx[inside],gy[inside])
    w_out_imag[inside] = w_imag_rbf(gx[inside],gy[inside])

    outside = ~inside
    if np.any(outside):
        w_query = map_through_nearest_triangle(tri_m,mapped_sample_points,grid_sample_w,
            gx[outside]+1j*gy[outside],skip=sliver_tri)
        w_out_real[outside] = np.real(w_query)
        w_out_imag[outside] = np.imag(w_query)

    valid = np.isfinite(w_out_real)&np.isfinite(w_out_imag)
    extrapolated = valid&outside

    equi_out = stereo_to_equi_deg(w_out_real+1j*w_out_imag,opposing_point)
    lon_out,lat_out = np.real(equi_out),np.imag(equi_out)
    del equi_out,w_out_real,w_out_imag

    global max_neighbor_dist
    ring_m = boundary_ring_in_m(projection_name,grid_sample_w,mapped_sample_points,
        loaded_pickle,opposing_point,cell_count)
    m_boundary = make_valid(Polygon(np.stack([ring_m.real,ring_m.imag],axis=1)))
    within_boundary = shapely.contains_xy(m_boundary,gx,gy)

    on_mesh_fraction = float(inside.mean())
    extrapolated_fraction = float(extrapolated.mean())
    del gx,gy,simplex,on_mesh,rejected,inside,outside,extrapolated,node_xy_m,m_boundary,ring_m

    neighbour_gap = neighbour_gap_on_globe(lon_out,lat_out)
    measurable = within_boundary&np.isfinite(neighbour_gap)
    max_neighbor_dist = float(neighbour_gap[measurable].max()) if np.any(measurable) else None
    del neighbour_gap,measurable

    folded = blank_folded_pixels(lon_out,lat_out,within_boundary,2*max_neighbor_dist) \
        if max_neighbor_dist else np.zeros(valid.shape,dtype=bool)
    valid &= ~folded

    img = np.zeros((ny,nx,3),dtype=np.uint8)
    alpha = np.zeros((ny,nx),dtype=np.uint8)
    with rasterio.open(raster_path) as src:
        inv = ~src.transform
        vy,vx = np.where(valid)
        cols,rows = inv*(lon_out[valid],lat_out[valid])
        cols = np.clip(np.round(cols).astype(int),0,src.width-1)
        rows = np.clip(np.round(rows).astype(int),0,src.height-1)
        data = src.read()
        img[vy,vx] = data[:,rows,cols].T
        alpha[vy,vx] = 255

    px = (maxx-minx)/(nx-1) if nx>1 else (maxx-minx)
    py = (maxy-miny)/(ny-1) if ny>1 else (maxy-miny)
    transform = rasterio.transform.from_origin(minx-px/2,maxy+py/2,px,py)
    with rasterio.open(out_tif,'w',driver='GTiff',width=nx,height=ny,count=4,dtype='uint8',
            transform=transform,crs=None,photometric='RGB',alpha='unspecified',
            compress='deflate',tiled=True) as dst:
        for band in range(3):
            dst.write(img[:,:,band],band+1)
        dst.write(alpha,4)
        dst.colorinterp = [ColorInterp.red,ColorInterp.green,ColorInterp.blue,ColorInterp.alpha]

    return {'raster':out_tif,'width':nx,'height':ny,'valid_fraction':float(valid.mean()),
        'on_mesh_fraction':on_mesh_fraction,
        'extrapolated_fraction':extrapolated_fraction,
        'within_boundary_fraction':float(within_boundary.mean()),
        'blanked_fraction':float(folded.mean()),
        'max_neighbor_dist':max_neighbor_dist}

AREA_DISTORTION_BUCKETS = [0,0.99,0.995,0.998,0.999,1.0,float('inf')]

def _area_distortion_label(ratio):
    if not np.isfinite(ratio):
        return 'degenerate'
    for lo,hi in zip(AREA_DISTORTION_BUCKETS,AREA_DISTORTION_BUCKETS[1:]):
        if hi==float('inf'):
            if ratio>=lo:
                return f'>= {lo:g}'
        elif lo<=ratio<hi:
            return f'{lo:g}-{hi:g}'
    return f'< {AREA_DISTORTION_BUCKETS[0]:g}'

def run_interpolation_algorithm_distortion(projection_name,out_shp=None):
    out_shp = out_shp or f'maps_projected/distortion_{projection_name}.shp'

    (grid_sample_w,mapped_sample_points,loaded_pickle,opposing_point,
        cell_count,boundary_count,boundary_W,boundary_M,rotation_phase) = _load_rotated_mesh(projection_name)

    bw_xy = np.stack([boundary_W.real,boundary_W.imag],axis=1)
    _,nearest = cKDTree(bw_xy).query(bw_xy,k=2)
    nearest = nearest[:,1]
    equi_b = stereo_to_equi_deg(boundary_W,opposing_point)
    lon_b,lat_b = np.radians(equi_b.real),np.radians(equi_b.imag)
    vb = np.stack([np.cos(lat_b)*np.cos(lon_b),np.cos(lat_b)*np.sin(lon_b),np.sin(lat_b)],axis=1)
    sph_dist = np.arccos(np.clip(np.sum(vb*vb[nearest],axis=1),-1,1))
    m_dist = np.abs(boundary_M-boundary_M[nearest])
    valid_b = sph_dist>1e-9
    area_scale_ref = float(np.median((m_dist[valid_b]/sph_dist[valid_b])**2))

    node_xy = np.stack([grid_sample_w.real,grid_sample_w.imag],axis=1)
    triangulation = lattice_triangulation(projection_name,grid_sample_w,cell_count)
    if triangulation is None:
        triangulation = Delaunay(node_xy)
    verts = triangulation.simplices
    T = verts.shape[0]
    Wt = grid_sample_w[verts]
    Mt = mapped_sample_points[verts]

    equi_deg = stereo_to_equi_deg(Wt,opposing_point)
    lon,lat = np.radians(equi_deg.real),np.radians(equi_deg.imag)
    v = np.stack([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)],axis=-1)
    v0,v1,v2 = v[:,0],v[:,1],v[:,2]
    triple = np.sum(v0*np.cross(v1,v2),axis=-1)
    denom = 1+np.sum(v0*v1,axis=-1)+np.sum(v1*v2,axis=-1)+np.sum(v2*v0,axis=-1)
    spherical_area = 2*np.abs(np.arctan2(triple,denom))

    Mx,My = Mt.real,Mt.imag
    area_m = 0.5*np.abs((Mx[:,1]-Mx[:,0])*(My[:,2]-My[:,0])-(Mx[:,2]-Mx[:,0])*(My[:,1]-My[:,0]))

    Sys = np.empty((T,3,3),dtype=complex)
    Sys[:,:,0] = Wt
    Sys[:,:,1] = np.conj(Wt)
    Sys[:,:,2] = 1.0
    Sys_inv = np.zeros_like(Sys)
    nondegenerate = np.abs(np.linalg.det(Sys))>1e-300
    Sys_inv[nondegenerate] = np.linalg.inv(Sys[nondegenerate])
    A,Bc = Sys_inv[:,0,:],Sys_inv[:,1,:]
    a = np.sum(A*Mt,axis=1)
    b = np.sum(Bc*Mt,axis=1)
    s = (a*np.conj(a)).real-(b*np.conj(b)).real
    coeff_mag = np.maximum(np.abs(A).max(axis=1),np.abs(Bc).max(axis=1))

    valid = (s>0)&np.isfinite(s)&(coeff_mag<1e5)&(spherical_area>1e-15)&(area_m>0)
    ratio = np.full(T,np.nan)
    ratio[valid] = (area_m[valid]/spherical_area[valid])/area_scale_ref

    labels = [_area_distortion_label(r) for r in ratio]
    polys_xy = [np.stack([Mt[i].real,Mt[i].imag],axis=1) for i in range(T)]
    polygons_with_attrs_to_file(out_shp,polys_xy,{
        'bucket':labels,
        'ratio':[float(r) if np.isfinite(r) else -1.0 for r in ratio],
    })
    return {'distortion':out_shp,'triangle_count':T}

def write_shifted_mesh(projection_name,N,w_pos,delta_out,rotation_phase,out_suffix):
    delta_raw = delta_out*D*np.exp(-1j*rotation_phase)

    M = np.load(f'pickle/interpolation_points/{projection_name}_M.npy').astype(complex)
    W = np.load(f'pickle/interpolation_points/{projection_name}_W.npy')
    blocks_path = f'pickle/interpolation_points/{projection_name}_blocks.npy'
    if os.path.isfile(blocks_path):
        cell_count,boundary_count = (int(v) for v in np.load(blocks_path))
        if cell_count!=N:
            raise ValueError(f'{projection_name}: blocks.npy cell_count ({cell_count}) != centers_file.txt length ({N})')
    else:
        cell_count = N
        boundary_count = M.shape[0]-13*cell_count

    M[:cell_count] += delta_raw

    if boundary_count>0:
        boundary_W = W[cell_count:cell_count+boundary_count]
        cell_tree = cKDTree(np.stack([w_pos.real,w_pos.imag],axis=1))
        _,owner = cell_tree.query(np.stack([boundary_W.real,boundary_W.imag],axis=1))
        M[cell_count:cell_count+boundary_count] += delta_raw[owner]

    extra_start = cell_count+boundary_count
    n_extra = M.shape[0]-extra_start
    if n_extra%cell_count!=0:
        raise ValueError(f'{projection_name}: {n_extra} leftover M.npy points not evenly divisible by {cell_count} cells')
    per_cell = n_extra//cell_count
    M[extra_start:] = (M[extra_start:].reshape(cell_count,per_cell)+delta_raw[:,None]).reshape(-1)

    out_path = f'pickle/interpolation_points/{projection_name}{out_suffix}_M.npy'
    np.save(out_path,M)
    return {'relaxed_M':out_path,'cell_count':cell_count,'boundary_count':boundary_count,
        'detail_points_per_cell':per_cell}

if __name__=='__main__':
    run_interpolation_algorithm(sys.argv[1],None if len(sys.argv)<2 else sys.argv[2])
