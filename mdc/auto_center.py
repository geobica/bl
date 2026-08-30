import numpy as np
from math import pi

def wrap_lon(lon):
    return (lon+180)%360-180

def antipode(lon,lat):
    return wrap_lon(lon+180),-lat

def _rings_of(boundary):
    if hasattr(boundary,'geom_type'):
        polys = [boundary] if boundary.geom_type=='Polygon' else list(boundary.geoms)
        rings = []
        for poly in polys:
            if poly.geom_type!='Polygon':
                continue
            rings.append(np.asarray(poly.exterior.coords,dtype=float))
            rings.extend(np.asarray(hole.coords,dtype=float) for hole in poly.interiors)
        return rings
    return [np.asarray(ring,dtype=float) for ring in boundary]

def _densify(ring,max_step_deg):
    lon,lat = ring[:,0].copy(),ring[:,1]
    lon[1:] = lon[0]+np.cumsum((np.diff(lon)+180)%360-180)
    out = [ring[:1]]
    for i in range(len(ring)-1):
        steps = int(np.ceil(max(abs(lon[i+1]-lon[i]),abs(lat[i+1]-lat[i]))/max_step_deg))
        frac = np.linspace(0,1,max(steps,1)+1)[1:,None]
        out.append(np.column_stack([lon[i],lat[i]])+frac*np.array([[lon[i+1]-lon[i],lat[i+1]-lat[i]]]))
    return np.vstack(out)

def _unit_vectors(lonlat_deg):
    lon,lat = np.radians(lonlat_deg[:,0]),np.radians(lonlat_deg[:,1])
    return np.column_stack([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)])

def _turning(ring_rad,lon,lat):
    dlon = ring_rad[:,0]-np.radians(lon)
    rlat = ring_rad[:,1]
    plat = np.radians(lat)
    bearing = np.arctan2(np.sin(dlon)*np.cos(rlat),
        np.cos(plat)*np.sin(rlat)-np.sin(plat)*np.cos(rlat)*np.cos(dlon))
    turn = (np.diff(bearing,append=bearing[:1])+pi)%(2*pi)-pi
    return turn.sum()

class _Boundary:
    def __init__(self,boundary,max_step_deg=0.5):
        self.rings = [_densify(ring,max_step_deg) for ring in _rings_of(boundary)]
        if not self.rings:
            raise ValueError('the boundary has no rings')
        self.rings_rad = [np.radians(ring) for ring in self.rings]
        self.vectors = np.vstack([_unit_vectors(ring) for ring in self.rings])

    def contains(self,lon,lat):
        return sum(_turning(ring,lon,lat) for ring in self.rings_rad)>pi

    def distance(self,lon,lat):
        v = _unit_vectors(np.array([[lon,lat]]))[0]
        chord = np.sqrt(np.maximum(((self.vectors-v)**2).sum(axis=1).min(),0.0))
        return np.degrees(2*np.arcsin(min(chord/2,1.0)))

def _score(edge,lon,lat):
    alon,alat = antipode(lon,lat)
    if not edge.contains(lon,lat):
        return None
    if edge.contains(alon,alat):
        return None
    return min(edge.distance(lon,lat),edge.distance(alon,alat))

def find_center(boundary,grid_n=40):
    from scipy.optimize import minimize
    edge = _Boundary(boundary)

    best = None
    best_score = -1.0
    for x in np.linspace(-180,180,grid_n,endpoint=False):
        for y in np.linspace(-90,90,grid_n):
            s = _score(edge,x,y)
            if s is not None and s>best_score:
                best_score = s
                best = (x,y)

    if best is None:
        raise ValueError('No point that is within the boundary and which has an antipode outside the boundary could be found.')

    def neg_score(pt):
        s = _score(edge,wrap_lon(pt[0]),pt[1])
        return -s if s is not None else 1e9

    refined = minimize(neg_score,best,method='Nelder-Mead',
        options={'xatol':1e-3,'fatol':1e-3,'maxiter':500})
    if refined.success:
        lon,lat = wrap_lon(float(refined.x[0])),float(refined.x[1])
        s = _score(edge,lon,lat)
        if s is not None and s>best_score:
            best_score = s
            best = (lon,lat)

    return float(best[0]),float(best[1]),float(best_score)
