import numpy as np
from math import e,pi,sqrt,sin,cos,tan,atan
import os
import warnings
warnings.filterwarnings('ignore', message='Unable to import Axes3D')
import matplotlib.pyplot as plt
from scipy import io
import scipy.special
import copy
from shapely.geometry import Polygon, MultiPolygon, Point
from shapely.prepared import prep
from concurrent.futures import ThreadPoolExecutor
import pickle
from scipy.optimize import minimize
import time
import sys
from utils import *

import contextlib

_DISPLAY_VARS = ('DISPLAY', 'WAYLAND_DISPLAY', 'XDG_SESSION_TYPE', 'DBUS_SESSION_BUS_ADDRESS')

@contextlib.contextmanager
def _spawning_octave():
    saved = {v: os.environ.pop(v, None) for v in _DISPLAY_VARS}
    saved['QT_QPA_PLATFORM'] = os.environ.pop('QT_QPA_PLATFORM', None)
    os.environ.setdefault('NO_AT_BRIDGE', '1')
    os.environ.setdefault('QT_LOGGING_RULES', '*.warning=false')
    os.environ['QT_QPA_PLATFORM'] = 'offscreen'
    try:
        yield
    finally:
        for v, val in saved.items():
            if val is None:
                os.environ.pop(v, None)
            else:
                os.environ[v] = val

oct2py = None
_oc_instance = None

def _oc():
    global oct2py,_oc_instance
    if _oc_instance is None:
        with _spawning_octave():
            import oct2py as _oct2py_module
            oct2py = _oct2py_module
            _oc_instance = oct2py.Oct2Py()
    return _oc_instance

D = 12756274

FOURIER_POWER = 8
BUFFER_DIVISOR = 240
WORKERS = 1
LEGACY = False
ETA_REFERENCE_CELLS = None
beta_raw = None

def f_W_ring(w,z,beta,power=8):
    theta = np.linspace(0,2*pi,2**power)
    eth_k = np.remainder(np.imag(np.log(z)),2*pi)
    eth_k[eth_k>2*pi-0.0000000000001]=0
    sums_d2 = np.prod((2-2*np.cos(theta[:,None]-eth_k[None,:]))**(beta_raw[None,:]/2),axis=1)
    prev_eth = np.zeros_like(theta)
    w_pos = np.zeros_like(theta+0j)
    for i in range(z.shape[0]):
        if eth_k[i]<eth_k[(i+1)%z.shape[0]]:
            prev_eth[np.logical_and(theta>=eth_k[i],theta<eth_k[(i+1)%z.shape[0]])] = i
    for i in range(z.shape[0]):
        wheres = np.where(prev_eth==i)[0]
        if wheres.shape[0]==1:
            total_len = eth_k[(i+1)%z.shape[0]]-eth_k[i]
            if total_len<=theta[1]:
                w_pos[wheres[0]] = w[i]
            else:
                lambd = np.nan_to_num((theta[wheres[0]]-eth_k[i])/total_len)
                w_pos[wheres[0]] = w[i]*(1-lambd)+w[(i+1)%z.shape[0]]*lambd
        elif wheres.shape[0]>1:
            test = np.prod((2-2*np.cos(eth_k[i]-eth_k[:i]))**(beta_raw[:i]/2))*np.prod((2-2*np.cos(eth_k[i]-eth_k[i+1:]))**(beta_raw[i+1:]/2))
            test = np.nan_to_num(test)
            test3 = np.prod((2-2*np.cos(eth_k[(i+1)%z.shape[0]]-eth_k[:(i+1)%z.shape[0]]))**(beta_raw[:(i+1)%z.shape[0]]/2))*np.prod((2-2*np.cos(eth_k[(i+1)%z.shape[0]]-eth_k[(i+1)%z.shape[0]+1:]))**(beta_raw[(i+1)%z.shape[0]+1:]/2))
            test3 = np.nan_to_num(test3)

            #integral from eth_k[i] to theta[wheres[0]]:
            b = beta_raw[i]/2
            v = eth_k[i]
            br = beta_raw[(i+1)%z.shape[0]]/2
            vr = eth_k[(i+1)%z.shape[0]]

            integral_left = lambda x:-(np.sqrt(2)*np.sqrt(np.cos(v - x) + 1)*np.tan((v - x)/2)*(2 - 2*np.cos(v - x))**b*scipy.special.hyp2f1(1/2, b + 1/2, b + 3/2, np.sin((v - x)/2)**2))/(2*b + 1)
            integral_right = lambda x:-(np.sqrt(2)*np.sqrt(np.cos(vr - x) + 1)*np.tan((vr - x)/2)*(2 - 2*np.cos(vr - x))**br*scipy.special.hyp2f1(1/2, br + 1/2, br + 3/2, np.sin((vr - x)/2)**2))/(2*br + 1)

            integral_dtheta = np.zeros_like(theta[wheres]+0j)
            int_left_eval = integral_left(theta[wheres[0]])
            int_right_eval = integral_right(theta[wheres[-1]])
            if not np.any(np.isinf(int_left_eval)):
                integral_dtheta += test*int_left_eval
            sums_here = sums_d2[wheres]
            if np.any(np.isinf(sums_here[0])):
                sums_here[0] = sums_here[1]

            integral_dtheta[1:]+=np.cumsum((sums_here[:-1]+sums_here[1:])/2)*theta[1]

            total_len = integral_dtheta[-1]
            if not np.any(np.isinf(int_right_eval)):
                total_len-=test3*integral_right(theta[wheres[-1]]) #this is correct
            if total_len<=0:
                w_pos[wheres] = w[i]
            else:
                lambd = np.nan_to_num(integral_dtheta/total_len)
                w_pos[wheres] = w[i]*(1-lambd)+w[(i+1)%z.shape[0]]*(lambd)
                if np.any(np.isnan(w_pos[wheres])):
                    w_pos[wheres] = w[i]

    return w_pos

def prepare_coefficients(w,z,beta,power=8):
    sample_density = 2**power
    f_M_z_in_pre = np.linspace(0,2*pi,sample_density,endpoint=False)

    f_M_z_loopup_table = f_W_ring(w,z,beta,power=power)

    t1_ln_tlu = np.log(np.absolute(f_M_z_loopup_table/D)**2+1)

    fft_get = -np.fft.fft(t1_ln_tlu)/2**power

    return fft_get

_gl_cache = {}

def gauss_legendre_01(order):
    if order not in _gl_cache:
        nodes,weights = np.polynomial.legendre.leggauss(order)
        _gl_cache[order] = ((nodes+1)/2,weights/2)
    return _gl_cache[order]

def f_prime_W_grid(zeta,z,beta,p):
    out = np.empty(zeta.shape,dtype=complex)
    block = max(1,int(2000000/max(1,z.shape[0])))
    flat = zeta.reshape(-1)
    res = out.reshape(-1)
    for s in range(0,flat.shape[0],block):
        chunk = flat[s:s+block]
        res[s:s+block] = np.exp(np.sum(beta[:,None]*np.log(1-chunk[None,:]/z[:,None]),axis=0)
            + np.polynomial.polynomial.polyval(chunk,p))
    return out

def sc_map_continuous(zp,z,beta,p):
    zp = np.asarray(zp)
    flat = zp.reshape(-1)
    wp = np.zeros(flat.shape,dtype=complex)
    radius = np.absolute(flat)
    orders = np.select([radius<0.5,radius<0.8,radius<0.95],[16,32,64],96)
    for order in (16,32,64,96):
        sel = np.where((radius>0) & (orders==order))[0]
        if sel.shape[0]==0:
            continue
        nodes,weights = gauss_legendre_01(order)
        here = flat[sel]
        zeta = here[:,None]*nodes[None,:]
        wp[sel] = here*np.sum(f_prime_W_grid(zeta,z,beta,p)*weights[None,:],axis=1)
    wp = wp.reshape(zp.shape)
    return wp

def tlu_M_over_tlu_W(zp,z,beta,p):
    k = np.arange(p.shape[0])
    return np.absolute(np.exp(np.sum(p[:,None]*zp[None,:]**k[:,None],axis=0)))

def unpack_center(new_centers,j,shared):
    record = ['']
    for i in range(1,10):
        value = new_centers[j,i][:,0]
        if i>=7:
            ref = shared[i-7]
            if value.shape==ref.shape and np.array_equal(value,ref):
                record.append(ref)
                continue
        record.append(copy.copy(value))
    record.append('')
    for i in range(11,13):
        record.append(copy.copy(new_centers[j,i][:,0]))
    return record

def slim_record(record):
    slim = list(record)
    slim[7] = ''
    slim[8] = ''
    slim[9] = ''
    return slim

level_cell_count = {}
worker_sessions = []

def save_level(path,centers,current_level,index_within_last,r):
    level_cell_count[path] = len(current_level)
    compact = [slim_record(centers[i]) for i in current_level]
    scipy.io.savemat(path, mdict={'centers': compact,'r': r,
        'current_level': list(range(len(current_level))),'index_within_last': index_within_last})

def run_layer_b(preceeder,level_i):
    infile_sc = f'pickle/matlab_saves/{preceeder}sc.mat'
    infile_level = f'pickle/matlab_saves/{preceeder}level_{level_i}.mat'
    infile_prev = f'pickle/matlab_saves/{preceeder}level_{level_i-1}_computed.mat'
    outfile = f'pickle/matlab_saves/{preceeder}level_{level_i}_computed.mat'
    cells = 2*level_cell_count.get(infile_level,0)
    if WORKERS<=1 or LEGACY or cells<8*WORKERS:
        _oc().feval('sc_toolbox_octave/layer_b.m',infile_sc,infile_level,infile_prev,outfile)
        return
    _oc()
    prefix = f'pickle/matlab_saves/{preceeder}level_{level_i}_'
    with _spawning_octave():
        while len(worker_sessions)<WORKERS:
            worker_sessions.append(oct2py.Oct2Py())
    def chunk(worker):
        worker_sessions[worker-1].feval('sc_toolbox_octave/layer_b_chunk.m',infile_sc,infile_level,
            infile_prev,float(worker),float(WORKERS),prefix)
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        list(pool.map(chunk,range(1,WORKERS+1)))
    _oc().feval('sc_toolbox_octave/layer_b_merge.m',prefix,float(WORKERS),outfile)

def explore_polygon(poly_coords,preceeder,inverted_buffer_metric=False,crpfun_workers=0):
    poly_string = 'A = ['
    x = np.real(poly_coords)
    y = np.imag(poly_coords)
    for i in range(x.shape[0]):
        poly_string = poly_string+str(x[i])
        if y[i]>=0:
            poly_string = poly_string+'+'+str(y[i])+'i'
        else:
            poly_string = poly_string+str(y[i])+'i'
        if i<x.shape[0]-1:
            poly_string = poly_string+','
        else:
            poly_string = poly_string+'];'
    sides = 7
    r = sqrt(-5 + sqrt(3*sides)*cos(1/3*atan(1/(3*sqrt(3)))) + 3*sqrt(sides)*sin(1/3*atan(1/(3*sqrt(3)))))

    if not os.path.isfile(f'pickle/matlab_saves/{preceeder}sc.mat'):
        _oc().eval(poly_string)
        _oc().eval(f'A_D = {r};')
        _oc().eval(f'outfile = strcat(\'pickle/matlab_saves/{preceeder}sc.mat\');')
        print(f'sc.mat: solving crdiskmap for {poly_coords.shape[0]} boundary vertices...')
        _sc_start = time.time()
        _oc().feval('sc_toolbox_octave/layer_0.m',_oc().pull('A'),_oc().pull('A_D'),_oc().pull('outfile'),crpfun_workers)
        print(f'sc.mat: done in {time.time()-_sc_start:.1f}s')
    #data = np.loadtxt(f'pickle/matlab_saves/{preceeder}sc.mat')
    loaded = io.loadmat(f'pickle/matlab_saves/{preceeder}sc', verify_compressed_data_integrity=False)
    w_raw = loaded['w'][:,0]
    beta_raw = loaded['beta'][:,0]
    z_raw = loaded['z'][:,0]

    centers = []
    centers.append(['',np.array([0j]),loaded['A_W'][:,0],0,loaded['A_D'][:,0],loaded['R_D'][:,0],loaded['L_D'][:,0],w_raw,z_raw,beta_raw,'',loaded['R_W'][:,0],loaded['L_W'][:,0]])
    save_level(f'pickle/matlab_saves/{preceeder}level_0.mat',centers,[0],[],r)
    if not os.path.isfile(f'pickle/matlab_saves/{preceeder}level_0_computed.mat'):
        _oc().eval(f'infile_sc = \'pickle/matlab_saves/{preceeder}sc.mat\';')
        _oc().eval(f'infile_level = \'pickle/matlab_saves/{preceeder}level_0.mat\';')
        _oc().eval(f'outfile_level_computed = strcat(\'pickle/matlab_saves/{preceeder}level_0_computed.mat\');')
        _oc().feval('sc_toolbox_octave/layer_a.m',_oc().pull('infile_sc'),_oc().pull('infile_level'),_oc().pull('outfile_level_computed'))

    loaded = io.loadmat(f'pickle/matlab_saves/{preceeder}level_0_computed', verify_compressed_data_integrity=False)
    new_centers = loaded['new_centers']
    shared = (w_raw,z_raw,beta_raw)
    current_level = []
    index_within_last = []
    for j in range(len(new_centers)):
        current_level.append(len(centers))
        centers.append(unpack_center(new_centers,j,shared))
        index_within_last.append(j)

    save_level(f'pickle/matlab_saves/{preceeder}level_1.mat',centers,current_level,index_within_last,r)

    w_polygon = Polygon(np.stack([np.real(poly_coords),np.imag(poly_coords)],axis=1))
    w_prepared = prep(w_polygon)
    w_edge = w_polygon.exterior
    w_placed_arr = np.array([complex(np.ravel(center[1])[0]) for center in centers])

    done_exploring = False
    level_i = 0
    _explore_start = time.time()
    _level_start = _explore_start
    while not done_exploring:
        level_i += 1
        run_layer_b(preceeder,level_i)

        loaded = io.loadmat(f'pickle/matlab_saves/{preceeder}level_{level_i}_computed', verify_compressed_data_integrity=False)
        new_centers = loaded['new_centers']
        current_level = []
        index_within_last = []
        for j in range(len(new_centers)):
            new_center_to_add = unpack_center(new_centers,j,shared)
            w_cen = complex(np.ravel(new_center_to_add[1])[0])
            min_val = np.amin(np.absolute(w_placed_arr-w_cen))
            w_point = Point(w_cen.real,w_cen.imag)
            if inverted_buffer_metric:
                if w_cen == 0:
                    in_range = True
                else:
                    buffering = np.absolute(1/w_cen)/80
                    inv_point = Point(np.real(1/w_cen),np.imag(1/w_cen))
                    inv_polygon = Polygon(np.stack([np.real(1/w_raw),np.imag(1/w_raw)],axis=1)).buffer(buffering)
                    in_range = inv_polygon.contains(inv_point)
            else:
                buffering = float((np.absolute(w_cen/D)**2+1)*D/BUFFER_DIVISOR)
                in_range = w_prepared.contains(w_point) and w_edge.distance(w_point)>buffering
            reach = np.absolute(new_center_to_add[11]-new_center_to_add[1])
            if min_val>0.01*reach and in_range and np.absolute(new_center_to_add[2]-new_center_to_add[1])<3*reach:
                current_level.append(len(centers))
                centers.append(new_center_to_add)
                w_placed_arr = np.append(w_placed_arr,w_cen)
                index_within_last.append(j)

        _now = time.time()
        _level_secs = _now - _level_start
        _total_secs = _now - _explore_start
        _level_start = _now
        _added = len(current_level)
        _total = len(centers)
        _pace = _total_secs/_total if _total else float('nan')
        eta_str = ''
        if ETA_REFERENCE_CELLS:
            _remaining = max(0,ETA_REFERENCE_CELLS-_total)
            _eta_secs = _remaining*_pace
            eta_str = f' | ETA ~{_eta_secs/60:.1f}min to reach {ETA_REFERENCE_CELLS} cells'
        print(f'level {level_i}: +{_added} cells (total {_total}) | '
            f'{_level_secs:.1f}s this level, {_total_secs/60:.1f}min elapsed, '
            f'{_pace:.3f}s/cell avg{eta_str}')
        with open(f'pickle/matlab_saves/{preceeder}level_{level_i}_centers.txt', 'wb') as fh:
            pickle.dump(centers, fh)
        for old_level in range(1,level_i-1):
            stale = f'pickle/matlab_saves/{preceeder}level_{old_level}_centers.txt'
            if os.path.isfile(stale):
                os.remove(stale)
        with open(f'pickle/matlab_saves/{preceeder}level_{level_i}_current_level.txt', 'wb') as fh:
            pickle.dump(current_level, fh)
        with open(f'pickle/matlab_saves/{preceeder}level_{level_i}_index_within_last.txt', 'wb') as fh:
            pickle.dump(index_within_last, fh)

        if len(current_level)==0:
            done_exploring = True
        else:
            save_level(f'pickle/matlab_saves/{preceeder}level_{level_i+1}.mat',centers,current_level,index_within_last,r)

    with open(f'pickle/matlab_saves/{preceeder}centers_file.txt', 'wb') as fh:
         pickle.dump(centers, fh)

    return w_raw,beta_raw,z_raw,centers

def M_centers_from_W(centers_W,w_raw,beta_raw,z_raw,preceeder):
    line_placements = []
    line_posers = []
    line_placements_w = []

    line_placements_raw = []
    line_placements_raw_w = []

    ii=0
    w_center_array = np.hstack(np.array([cent[1] for cent in centers_W]))
    w_anchor_array = np.hstack(np.array([cent[2] for cent in centers_W]))
    w_right_array = np.hstack(np.array([cent[11] for cent in centers_W]))
    power = FOURIER_POWER

    resize_factor = []
    previous_index = []
    associates = []
    tlu_M_list = []
    extra_tlu_M_list = []

    extra_tlu_M_over_tlu_W_list = []

    _start_time = time.time()
    for center in centers_W:
        progress_bar(f'Projected {ii}/{len(centers_W)} centers.',ii/max(len(centers_W),1),_start_time)
        previous_index.append(0)

        associates.append(np.absolute(center[8]-np.roll(center[8],1))/np.sum(np.absolute(center[8]-np.roll(center[8],1))))
        p = prepare_coefficients(center[7],center[8],center[9],power)
        detail_w = np.ndarray.flatten(np.array(center[7])[:,None]*np.linspace(0,1,10,endpoint=False)[None,:]+np.roll(np.array(center[7]),1)[:,None]*(1-np.linspace(0,1,10,endpoint=False))[None,:])
        detail_w/=D
        a = center[1]/D
        detail_w = (detail_w-a)/(1+detail_w*np.conj(a))
        stereo_dist = np.amin(np.absolute(detail_w))

        exp_p_0 = sc_map_continuous(np.array([center[4]+0j]),center[8],center[9],np.ones((1))*p[0])[0]/(center[2]-center[1])

        extra_D = np.array([center[4],center[5],center[6],center[4],center[5],center[6],center[4],center[5],center[6],center[4],center[5],center[6]])[:,0]
        extra_D[3:6] *= 1/3
        extra_D[6:9] *= e**(pi*1j/3)/4
        extra_D[9:12] *= e**(pi*1j/3)*3/4
        M_vals = sc_map_continuous(extra_D,center[8],center[9],p*2)
        W_vals = sc_map_continuous(extra_D,center[8],center[9],np.ones((1)))

        tlu_W = np.absolute(center[1]/D)**2+1
        tlu_M = np.absolute(tlu_W/np.exp(np.absolute(p[0])))
        exp_p_0 /= (1/(np.absolute(stereo_dist)**2+1)/tlu_M)
        tlu_M = 1/(np.absolute(stereo_dist)**2+1)

        extra_tlu_M_over_tlu_W = np.hstack([p[0]*2,*list(tlu_M_over_tlu_W(extra_D,center[8],center[9],p*2))])
        extra_tlu_M_over_tlu_W *= tlu_W/tlu_M

        M_vals = np.nan_to_num(np.hstack(M_vals))

        ii+= 1
        if len(line_placements)==0:
            line_placements.append([0,*list(M_vals)])
            line_placements_raw.append([0,*list(M_vals)])

            line_posers.append(1+0j)
            resize_factor.append(exp_p_0)
        else:
            # find center[2] in centers[:,1]
            anchor_i = np.argmin(np.absolute(w_center_array-center[2]))
            previous_index[len(previous_index)-1] = anchor_i
            # so this one's 0 corresponds to
            this_cent = line_placements[anchor_i][int(1+center[3][0])]

            this_anchor = line_placements[anchor_i][0]
            resize_factor.append(exp_p_0)

            line_placements_raw.append([0,*list(M_vals)])
            line_posers.append((this_anchor-this_cent)/M_vals[0])
            M_vals = this_cent+M_vals*(this_anchor-this_cent)/M_vals[0]

            line_placements.append([this_cent,*list(M_vals)])
        tlu_M_list.append(tlu_M)
        extra_tlu_M_over_tlu_W_list.append(extra_tlu_M_over_tlu_W)

        line_placements_w.append([0,*list(W_vals)])
        line_placements_raw_w.append([0,*list(W_vals)])
    progress_bar('',1.0,_start_time,done=True)

    with open(f'pickle/matlab_saves/{preceeder}line_placements.txt', 'wb') as fh:
         pickle.dump(line_placements, fh)
    with open(f'pickle/matlab_saves/{preceeder}line_posers.txt', 'wb') as fh:
         pickle.dump(line_posers, fh)
    with open(f'pickle/matlab_saves/{preceeder}line_placements_raw.txt', 'wb') as fh:
         pickle.dump(line_placements_raw, fh)
    with open(f'pickle/matlab_saves/{preceeder}resize_factor.txt', 'wb') as fh:
         pickle.dump(resize_factor, fh)
    rescaled_raw = np.absolute(np.hstack(resize_factor))

    line_posers = line_posers/np.absolute(line_posers)/rescaled_raw

    line_placements = []
    line_placements_w = []
    line_placements.append(list(np.hstack(line_placements_raw[0])*line_posers[0]))
    centers_M = []

    center = centers_W[0]
    detail_w = np.ndarray.flatten(np.array(center[7])[:,None]*np.linspace(0,1,10,endpoint=False)[None,:]+np.roll(np.array(center[7]),1)[:,None]*(1-np.linspace(0,1,10,endpoint=False))[None,:])
    detail_w/=D

    for i in range(0,len(line_placements_raw)):
        line_placements_w.append(w_center_array[i]+(w_anchor_array[i]-w_center_array[i])*np.hstack(line_placements_raw_w[i])/np.hstack(line_placements_raw_w[i])[1])
        extra_tlu_M_list.append(extra_tlu_M_over_tlu_W_list[i]*(np.absolute(line_placements_w[i]/D)**2+1))
        for ii in range(len(line_placements_w[i])):
            a = line_placements_w[i][ii]/D
            detail_w_mod = (detail_w-a)/(1+detail_w*np.conj(a))
            #detail_w*=D
            stereo_dist = np.amin(np.absolute(detail_w_mod))
            extra_tlu_M_list[i][ii] = 1/(np.absolute(stereo_dist)**2+1)

        if i > 0:
            if np.amin(np.absolute(w_center_array[:i]-w_right_array[i]))<0.01*np.absolute(w_center_array[i]-w_right_array[i]):
                print(np.amin(np.absolute(w_center_array[:i]-w_right_array[i])))
                line_placements.append((np.hstack(line_placements[np.argmin(np.absolute(w_center_array[:i]-w_right_array[i]))])[0]+line_posers[i]*(np.hstack(line_placements_raw[i])-np.hstack(line_placements_raw[i])[2])+\
                    np.hstack(line_placements[previous_index[i]])[0]+line_posers[i]*(np.hstack(line_placements_raw[i])-np.hstack(line_placements_raw[i])[1]))/2)
            else:
                line_placements.append(np.hstack(line_placements[previous_index[i]])[0]+line_posers[i]*(np.hstack(line_placements_raw[i])-np.hstack(line_placements_raw[i])[1]))
        centers_M.append(['',line_placements[len(line_placements)-1][0],line_placements[len(line_placements)-1][1],line_posers[i],0,0,0,0,0,0,0,line_placements[len(line_placements)-1][2],line_placements[len(line_placements)-1][3]])

    return centers_M,associates,tlu_M_list,line_placements,line_placements_w,extra_tlu_M_list#,line_placements_W_inv,line_placements_M_inv

'''
Create correspondence between:
string representation of heptagon vertex
position in W
position of its anchor in W (backward in tree), already known, store as A_W
whether this corresponds to the A (0), R (1), or L (2) of its anchor
position of its anchor in D calculated with A_D = evalinv(f,A_W)
position of its right branch in D
position of its left branch in D
w, z, and beta
'''

def run_pipeline(name):
    global beta_raw
    boundary_stride = 10

    from_start = True

    preceeder = f'{name}_'

    pickle_off = open(f"pickle/bubble_wrap/{preceeder}full_wrapping.txt", "rb")
    loaded_pickle = pickle.load(pickle_off)
    full_wrapping = np.array(loaded_pickle[0])[::-boundary_stride]

    nemo_point = loaded_pickle[1]-pi

    opposing_point = np.conj(nemo_point)
    k = D/(1+sin(opposing_point.imag)*np.sin(np.imag(full_wrapping))+cos(opposing_point.imag)*np.cos(np.imag(full_wrapping))*np.cos(np.real(full_wrapping)-opposing_point.real))
    x = k*np.cos(np.imag(full_wrapping))*np.sin(np.real(full_wrapping)-opposing_point.real)
    y = k*(cos(opposing_point.imag)*np.sin(np.imag(full_wrapping))-sin(opposing_point.imag)*np.cos(np.imag(full_wrapping))*np.cos(np.real(full_wrapping)-opposing_point.real))

    if from_start:
        poly_use = (x+y*1j)

        w_raw,beta_raw,z_raw,centers = explore_polygon(poly_use,preceeder,inverted_buffer_metric=False)

        centers_M,associates,tlu_M_list,line_placements_m,line_placements_w,extra_tlu_M_list = M_centers_from_W(centers,w_raw,beta_raw,z_raw,preceeder)
        with open(f'pickle/matlab_saves/{preceeder}centers_M.txt', 'wb') as fh:
             pickle.dump(centers_M, fh)
        with open(f'pickle/matlab_saves/{preceeder}associates.txt', 'wb') as fh:
             pickle.dump(associates, fh)
        with open(f'pickle/matlab_saves/{preceeder}tlu_M_list.txt', 'wb') as fh:
             pickle.dump(tlu_M_list, fh)
        with open(f'pickle/matlab_saves/{preceeder}line_placements_m.txt', 'wb') as fh:
             pickle.dump(line_placements_m, fh)
        with open(f'pickle/matlab_saves/{preceeder}line_placements_w.txt', 'wb') as fh:
             pickle.dump(line_placements_w, fh)
        with open(f'pickle/matlab_saves/{preceeder}extra_tlu_M_list.txt', 'wb') as fh:
             pickle.dump(extra_tlu_M_list, fh)
    else:
        poly_use = (x+y*1j)

        w_raw,beta_raw,z_raw,centers = explore_polygon(poly_use,preceeder,inverted_buffer_metric=False)

        pickle_off = open(f"pickle/matlab_saves/{preceeder}centers_M.txt", "rb")
        centers_M = pickle.load(pickle_off)
        pickle_off = open(f"pickle/matlab_saves/{preceeder}associates.txt", "rb")
        associates = pickle.load(pickle_off)
        pickle_off = open(f"pickle/matlab_saves/{preceeder}tlu_M_list.txt", "rb")
        tlu_M_list = pickle.load(pickle_off)
        pickle_off = open(f"pickle/matlab_saves/{preceeder}line_placements_m.txt", "rb")
        line_placements_m = pickle.load(pickle_off)
        pickle_off = open(f"pickle/matlab_saves/{preceeder}line_placements_w.txt", "rb")
        line_placements_w = pickle.load(pickle_off)
        pickle_off = open(f"pickle/matlab_saves/{preceeder}extra_tlu_M_list.txt", "rb")
        extra_tlu_M_list = pickle.load(pickle_off)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    power = FOURIER_POWER
    w_center_array = np.hstack(np.array([cent[1] for cent in centers]))
    w_anchor_array = np.hstack(np.array([cent[2] for cent in centers]))
    m_center_array = np.hstack(np.array([cent[1] for cent in centers_M]))
    m_anchor_array = np.hstack(np.array([cent[2] for cent in centers_M]))
    m_posers_array = np.hstack(np.array([cent[3] for cent in centers_M]))

    nearest_to_W_vert = np.argmin(np.absolute(w_center_array[:,None]-w_raw[None,:]),axis=0)

    M_vertices_est = []
    ii = 0
    for i in range(nearest_to_W_vert.shape[0]):
        if ii<poly_use.shape[0] and np.absolute(w_raw[i]-poly_use[ii])<0.00000001:
            center = centers[nearest_to_W_vert[i]]
            p = prepare_coefficients(center[7],center[8],center[9],power)
            M_vals = sc_map_continuous(np.array([center[8][i]]),center[8],center[9],p)
            M_vertices_est.append(m_center_array[nearest_to_W_vert[i]]+m_posers_array[nearest_to_W_vert[i]]*M_vals[0])
            ii+=1

    if len(M_vertices_est)!=poly_use.shape[0]:
        print(f'warning: {len(M_vertices_est)} of {poly_use.shape[0]} boundary vertices matched a solved prevertex')

    np.save(f'pickle/interpolation_points/{preceeder}blocks.npy',np.array([w_center_array.shape[0],len(M_vertices_est)]))
    np.save(f'pickle/interpolation_points/{preceeder}W.npy',np.hstack([w_center_array,poly_use[:len(M_vertices_est)],np.ndarray.flatten(np.array(line_placements_w)[:,4:])]))
    np.save(f'pickle/interpolation_points/{preceeder}M.npy',np.hstack([m_center_array,M_vertices_est,np.ndarray.flatten(np.array(line_placements_m)[:,4:])]))
    np.save(f'pickle/interpolation_points/{preceeder}tlu_M.npy',np.hstack([np.array(extra_tlu_M_list)[:,0],np.ones_like(M_vertices_est),np.ndarray.flatten(np.array(extra_tlu_M_list)[:,4:])]))

if __name__=='__main__':
    if len(sys.argv)>1:
        run_pipeline(sys.argv[1])
    else:
        import interactive_menu
        interactive_menu.run(run_pipeline)