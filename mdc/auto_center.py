import numpy as np
from shapely.geometry import Point
from scipy.optimize import minimize

def wrap_lon(lon):
    return (lon+180)%360-180

def antipode(lon,lat):
    return wrap_lon(lon+180),-lat

def _score(poly,boundary,lon,lat):
    alon,alat = antipode(lon,lat)
    p = Point(lon,lat)
    ap = Point(alon,alat)
    if not poly.contains(p):
        return None
    if poly.contains(ap):
        return None
    return min(boundary.distance(p),boundary.distance(ap))

def find_center(poly,grid_n=40):
    boundary = poly.boundary
    minx,miny,maxx,maxy = poly.bounds
    xs = np.linspace(minx,maxx,grid_n)
    ys = np.linspace(miny,maxy,grid_n)

    best = None
    best_score = -1.0
    for x in xs:
        for y in ys:
            s = _score(poly,boundary,x,y)
            if s is not None and s>best_score:
                best_score = s
                best = (x,y)

    if best is None:
        raise ValueError('No point that is within the boundary and which has an antipode outside the boundary could be found.')

    def neg_score(pt):
        s = _score(poly,boundary,pt[0],pt[1])
        return -s if s is not None else 1e9

    refined = minimize(neg_score,best,method='Nelder-Mead',
        options={'xatol':1e-3,'fatol':1e-3,'maxiter':500})
    if refined.success:
        s = _score(poly,boundary,refined.x[0],refined.x[1])
        if s is not None and s>best_score:
            best_score = s
            best = (float(refined.x[0]),float(refined.x[1]))

    return best[0],best[1],best_score
