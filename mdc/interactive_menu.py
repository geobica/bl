import os
import glob
import pickle
import time
import numpy as np
import fiona
from shapely.geometry import shape, Polygon, mapping
from shapely.ops import unary_union
import shapely.wkb

import menu_ui
import auto_center
import bubble_wrap
import interpolate
from utils import *

VECTOR_EXTS = ('gpkg','shp')
RASTER_EXTS = ('tif','tiff')
POLYGON_CACHE_DIR = 'pickle/vector_polygon_cache'

def _basenames(directory,exts=VECTOR_EXTS):
    names = set()
    for ext in exts:
        for path in glob.glob(f'{directory}/*.{ext}'):
            names.add(os.path.basename(path)[:-(len(ext)+1)])
    return sorted(names)

def _find_file(directory,name,exts):
    for ext in exts:
        path = f'{directory}/{name}.{ext}'
        if os.path.isfile(path):
            return path
    raise FileNotFoundError(f'{directory}/{name}.({"|".join(exts)})')

def _find_vector_file(directory,name):
    return _find_file(directory,name,VECTOR_EXTS)

def _find_raster_file(directory,name):
    return _find_file(directory,name,RASTER_EXTS)

def _load_polygon(path,on_feature=None):
    with fiona.open(path, layer=fiona.listlayers(path)[0]) as layer:
        geoms = []
        for feature in layer:
            geoms.append(shape(feature['geometry']).buffer(0))
            if on_feature is not None:
                on_feature()
    return unary_union(geoms)

def _load_polygon_cached(path,on_feature=None):
    os.makedirs(POLYGON_CACHE_DIR,exist_ok=True)
    st = os.stat(path)
    cache_path = f'{POLYGON_CACHE_DIR}/{os.path.basename(path)}.pkl'
    if os.path.isfile(cache_path):
        try:
            with open(cache_path,'rb') as fh:
                cached = pickle.load(fh)
        except (pickle.UnpicklingError,EOFError):
            cached = None
        if cached and cached.get('mtime')==st.st_mtime and cached.get('size')==st.st_size:
            if on_feature is not None:
                for _ in range(cached.get('feature_count',0)):
                    on_feature()
            return shapely.wkb.loads(cached['wkb'])

    with fiona.open(path, layer=fiona.listlayers(path)[0]) as layer:
        feature_count = len(layer)
    poly = _load_polygon(path,on_feature=on_feature)
    with open(cache_path,'wb') as fh:
        pickle.dump({'mtime':st.st_mtime,'size':st.st_size,'feature_count':feature_count,
            'wkb':shapely.wkb.dumps(poly)},fh)
    return poly

def _wrapping_to_lonlat_polygon(points_complex):
    pts = np.asarray(points_complex)
    lon = np.degrees(np.real(pts))
    lat = np.degrees(np.imag(pts))
    lon = (lon+180)%360-180
    return Polygon(np.column_stack([lon,lat]))

def _boundary_lonlat_polygon(name):
    with open(f'pickle/bubble_wrap/{name}_full_wrapping.txt','rb') as fh:
        pts,nemo = pickle.load(fh)
    return _wrapping_to_lonlat_polygon(pts)

def _save_boundary_polygon(name,poly):
    schema = {'geometry':'Polygon'}
    parts = list(poly.geoms) if poly.geom_type=='MultiPolygon' else [poly]
    with fiona.open(f'boundary_polygons/{name}.gpkg','w','GPKG',schema) as c:
        for part in parts:
            c.write({'geometry': mapping(Polygon(part.exterior))})

def _precomputed_projection_items():
    names = set()
    for suffix in ('_tlu_M_list.txt','_level_0.mat'):
        for path in glob.glob(f'pickle/matlab_saves/*{suffix}'):
            names.add(os.path.basename(path)[:-len(suffix)])
    names = {n for n in names if os.path.isfile(f'pickle/bubble_wrap/{n}_full_wrapping.txt')}
    items = []
    for name in sorted(names):
        done = os.path.isfile(f'pickle/matlab_saves/{name}_tlu_M_list.txt')
        items.append(menu_ui.MenuItem(name,('precomputed',name),
            color='light_gray' if done else 'dark_gray',
            annotation=None if done else 'in progress'))
    return items

def _boundary_polygon_items():
    items = [menu_ui.MenuItem(name,('boundary',name)) for name in _basenames('boundary_polygons')]
    items.append(menu_ui.MenuItem('New Boundary',('new_boundary',None)))
    return items

def _main_menu():
    sections = [
        menu_ui.MenuSection('A precomputed projection',_precomputed_projection_items()),
        menu_ui.MenuSection('Create a new projection from a boundary',_boundary_polygon_items()),
        menu_ui.MenuSection('Cancel',items=None),
    ]
    return menu_ui.run_menu(sections,'Which projection do you want to use?',enable_overwrite_toggle=True)

def _new_boundary_flow():
    names = _basenames('input_sample')
    if not names:
        print('No shapefiles found in input_sample/.')
        return None
    sections = [
        menu_ui.MenuSection('input_sample/',[menu_ui.MenuItem(n,n) for n in names]),
        menu_ui.MenuSection('Cancel',items=None),
    ]
    name,_overwrite = menu_ui.run_menu(sections,'Select a shapefile to run bubble_wrap on:')
    if name is None:
        return None

    points,nemo_point = bubble_wrap.run_bubble_wrap(name,input_path='input_sample',pickle_path='pickle/bubble_wrap')

    poly = _wrapping_to_lonlat_polygon(points)
    _save_boundary_polygon(name,poly)
    return name

def _boundary_rings(path):
    rings = []
    with fiona.open(path, layer=fiona.listlayers(path)[0]) as layer:
        for feature in layer:
            geom = shape(feature['geometry'])
            for poly in (geom.geoms if geom.geom_type == 'MultiPolygon' else [geom]):
                if poly.geom_type == 'Polygon':
                    rings.append(np.asarray(poly.exterior.coords, dtype=float))
                    rings.extend(np.asarray(hole.coords, dtype=float) for hole in poly.interiors)
    return rings

def _ensure_auto_center_wrapping(name):
    rings = _boundary_rings(_find_vector_file('boundary_polygons', name))
    lon, lat, score = auto_center.find_center(rings)
    print(f'Picked the center for {name} as: {abs(lat):.0f}°{"NS"[lat<0]}, {abs(lon):.0f}°{"EW"[lon<0]}')
    bubble_wrap.run_bubble_wrap(name,lon,lat,input_path='boundary_polygons',pickle_path='pickle/bubble_wrap')

def _overlap_color(boundary_poly,map_poly):
    if not boundary_poly.intersects(map_poly):
        return 'red'
    if boundary_poly.contains(map_poly):
        return 'light_gray'
    return 'orange'

def _print_result_path(label,path):
    print(f'{label} can be found at:')
    print(f' - {path}')

def _print_error_table(title,total_count,rows,unit=''):
    print(f'{title} ({total_count} measured):')
    labels = [f'{label} {unit}'.rstrip() for label,_,_ in rows]
    width = max((len(l) for l in labels),default=0)
    for label,(_,count,pct) in zip(labels,rows):
        print(f'  {label.ljust(width)}  {pct:6.2f}%  ({count})')

RASTER_MARGIN_CHOICES = (0,100,250,500,1000)
RASTER_RESOLUTION_CHOICES = (600,1200,2400,6000,12000)

def _raster_options():
    margin,_ = menu_ui.run_menu(
        [menu_ui.MenuSection('Add buffer on each side',
            [menu_ui.MenuItem(f'{m} px' if m else 'none',m) for m in RASTER_MARGIN_CHOICES]),
         menu_ui.MenuSection('Cancel',items=None)],
        'How many pixels of buffer should be on each side of the projected shapes?')
    if margin is None:
        return None,None
    resolution,_ = menu_ui.run_menu(
        [menu_ui.MenuSection('Pixels along the longer side',
            [menu_ui.MenuItem(f'{r} px'+(' (default)' if r==1200 else ''),r)
                for r in RASTER_RESOLUTION_CHOICES]),
         menu_ui.MenuSection('Cancel',items=None)],
        'What resolution should the raster be?')
    if resolution is None:
        return None,None
    return margin,resolution

def _project_flow(name):
    map_names = _basenames('maps_to_project')
    raster_names = _basenames('rasters_to_project',RASTER_EXTS)
    if not map_names and not raster_names:
        print(f'No shapefiles found in maps_to_project/ with extensions in {VECTOR_EXTS} and no rasters found in rasters_to_project/ with extensions in {RASTER_EXTS}.')
        return

    boundary_poly = _boundary_lonlat_polygon(name)

    map_paths = {mn:_find_vector_file('maps_to_project',mn) for mn in map_names}
    total_features = 0
    for path in map_paths.values():
        with fiona.open(path, layer=fiona.listlayers(path)[0]) as layer:
            total_features += len(layer)
    total_features = max(total_features,1)

    start_time = time.time()
    done_features = [0]
    def on_feature():
        done_features[0] += 1
        progress_bar(f'Loading {name} projection',done_features[0]/total_features,start_time)

    items = []
    for map_name in map_names:
        map_poly = _load_polygon_cached(map_paths[map_name],on_feature=on_feature)
        items.append(menu_ui.MenuItem(map_name,('vector',map_name),color=_overlap_color(boundary_poly,map_poly)))
    progress_bar('',1.0,start_time,done=True)

    raster_items = [menu_ui.MenuItem(rn,('raster',rn)) for rn in raster_names]

    other_items = [
        menu_ui.MenuItem('Graticule',('other','graticule')),
        menu_ui.MenuItem('Distortion',('other','distortion')),
    ]

    sections = [menu_ui.MenuSection('maps_to_project/',items)] if items else []
    if raster_items:
        sections.append(menu_ui.MenuSection('rasters_to_project/',raster_items))
    sections.append(menu_ui.MenuSection('Other',other_items))
    sections.append(menu_ui.MenuSection('Cancel',items=None))
    selected,_overwrite = menu_ui.run_menu(sections,'Which shapefile or raster do you want to project?')
    if selected is None:
        return
    kind,selected_name = selected

    if kind=='raster':
        raster_file = _find_raster_file('rasters_to_project',selected_name)
        margin,resolution = _raster_options()
        if margin is None:
            return
        raster_out = f'maps_projected/raster_{name}_{selected_name}.tif'
        result = interpolate.run_interpolation_algorithm_raster(name,raster_file,out_tif=raster_out,
            resolution=resolution,margin_px=margin)
        _print_result_path(f'The projection of {selected_name} within {name}',result['raster'])
        print(f'({result["width"]}x{result["height"]}, {100*result["valid_fraction"]:.1f}% of pixels filled; '
            f'{100*result["on_mesh_fraction"]:.1f}% sit on the mesh, the other '
            f'{100*result["extrapolated_fraction"]:.1f}% are extrapolated from the nearest triangle)')
        print(f'georeferenced to the same coordinates as boundary_{name}.shp and mesh_{name}.shp')
        return

    if kind=='other' and selected_name=='graticule':
        graticule_out = f'maps_projected/graticule_{name}.shp'
        result = interpolate.run_interpolation_algorithm(name,vector_file_to_project=None,outline_out=f'maps_projected/boundary_{name}.shp',
            mesh_out=f'maps_projected/mesh_{name}.shp',extra_outputs=False,
            graticule_out=graticule_out,graticule_step_deg=1.0)
        _print_result_path(f'The projection of a 1° graticule within {name}',result['graticule'])
        return

    if kind=='other' and selected_name=='distortion':
        distortion_out = f'maps_projected/distortion_{name}.shp'
        result = interpolate.run_interpolation_algorithm_distortion(name,out_shp=distortion_out)
        _print_result_path(f'The triangulation with an area distortion field containing {result["triangle_count"]} triangles for {name}',result['distortion'])
        print("(fields: 'bucket' -- a label like '0.99-0.995', 'ratio' -- the raw ratio, 1.0 at the boundary)")
        return

    vector_file = _find_vector_file('maps_to_project',selected_name)
    boundary_out = f'maps_projected/boundary_{name}.shp'
    projection_out = f'maps_projected/projection_{name}_{selected_name}.shp'
    mesh_out = f'maps_projected/mesh_{name}.shp'
    mesh_tree_out = f'maps_projected/mesh_tree_{name}.shp'
    mesh_tree_links_out = f'maps_projected/mesh_tree_links_{name}.shp'
    mesh_links_out = f'maps_projected/mesh_links_{name}.shp'
    result = interpolate.run_interpolation_algorithm(name,vector_file,outline_out=boundary_out,projection_out=projection_out,mesh_out=mesh_out,mesh_tree_out=mesh_tree_out,mesh_tree_links_out=mesh_tree_links_out,mesh_links_out=mesh_links_out,extra_outputs=False)

    _print_error_table(f'Angle error at vertices within {name}',result['angle_error_count'],
        result['angle_error_table'],unit='deg')
    _print_error_table(f'Length distortion of segments within {name}',result['length_error_count'],
        result['length_error_table'],unit='')
    _print_result_path('The projected points used to compute the triangulation',mesh_out)
    _print_result_path('Just the lattice tree those points were grown from',mesh_tree_out)
    _print_result_path("That tree's edges, plus the neighbours find_RL recovers",mesh_tree_links_out)
    _print_result_path('The same lattice drawn through the surrounding samples',mesh_links_out)
    _print_result_path(f'The projection of {selected_name} within {name}',projection_out)
    print('with the projected boundary at:')
    print(f' - {boundary_out}')

def run(run_pipeline_fn):
    selection,overwrite = _main_menu()
    if selection is None:
        return

    kind,name = selection
    if kind=='new_boundary':
        name = _new_boundary_flow()
        if name is None:
            return
        kind = 'boundary'

    if kind=='boundary' and (overwrite or not os.path.isfile(f'pickle/bubble_wrap/{name}_full_wrapping.txt')):
        _ensure_auto_center_wrapping(name)

    if overwrite or not os.path.isfile(f'pickle/matlab_saves/{name}_tlu_M_list.txt'):
        run_pipeline_fn(name)

    _project_flow(name)
