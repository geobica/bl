import copy
import math
import os
import pickle as pkl
from datetime import datetime
from os import pread
from typing import List, Dict

import numpy as np
import shapely
from matplotlib import pyplot as plt
from shapely import (
    ops,
    wkb,
    MultiPolygon,
    Polygon,
    LineString,
    MultiLineString,
    GeometryCollection,
)
import geopandas as gpd

from utils import get_closest_in_direction, equi_to_stereo, stereo_to_equi


def ordinal(n):
    return str(n) + (
        "th" if 4 <= n % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    )


merge_the_states = False

numstr = "100"
states_merged_str = ""
if merge_the_states:
    states_merged_str = "_MERGED"


def determine_boundary_simplification(multiboundary, for_tolerance=0.0001):
    simplified_boundaries = shapely.simplify(
        multiboundary, tolerance=for_tolerance, preserve_topology=True
    )  # MultiLineString(new_multiboundary)
    if type(simplified_boundaries) is LineString:
        if not simplified_boundaries.is_empty:
            simplified_boundaries = MultiLineString([simplified_boundaries])
        else:
            return simplified_boundaries

    cannot_be_destroyed = []
    for geom_i, geom in enumerate(simplified_boundaries.geoms):
        geoxy = np.array(geom.coords.xy).T
        if geoxy[0, 0] == geoxy[-1, 0] and geoxy[0, 1] == geoxy[-1, 1]:
            cannot_be_destroyed.append(False)
        else:
            cannot_be_destroyed.append(True)

    simplified_areas = []
    area_ids = []
    if len(simplified_boundaries.geoms) <= 1:
        return simplified_boundaries
    which_poly_points_are_from = np.hstack(
        [
            np.ones_like(geom.coords.xy) * geom_i
            for geom_i, geom in enumerate(simplified_boundaries.geoms)
        ]
    ).T[:, 0]
    all_points = np.hstack(
        [np.array(geom.coords.xy) for geom in simplified_boundaries.geoms]
    ).T
    for geom_i, geom in enumerate(simplified_boundaries.geoms):
        simplified_areas.append(geom.length)
        area_ids.append(geom_i)
    sorted_areas, sorted_ids = zip(*sorted(zip(simplified_areas, area_ids)))

    for_destruction = []
    undestroyed_points_for_comparison = which_poly_points_are_from != -1
    for sort_i, sort_area in enumerate(sorted_areas):
        geom_i = sorted_ids[sort_i]
        if not cannot_be_destroyed[geom_i]:
            to_compare_to = np.logical_and(
                undestroyed_points_for_comparison, which_poly_points_are_from != geom_i
            )
            if np.sum(to_compare_to) == 0:
                cannot_be_destroyed[geom_i] = True
                # there's no points to compare to
                break
            this_poly_xy = np.array(simplified_boundaries.geoms[geom_i].coords.xy).T
            for point_xy_i in range(this_poly_xy.shape[0]):
                if (
                        np.min(
                            np.linalg.norm(
                                (all_points[to_compare_to])
                                - (this_poly_xy[point_xy_i])[None, :],
                                axis=-1,
                            )
                        )
                        > for_tolerance
                ):
                    cannot_be_destroyed[geom_i] = True
                    break
        if not cannot_be_destroyed[geom_i]:
            for_destruction.append(geom_i)
            undestroyed_points_for_comparison[which_poly_points_are_from == geom_i] = (
                False
            )

    undestroyed_multiboundary = []
    for geom_i, geom in enumerate(simplified_boundaries.geoms):
        if geom_i not in for_destruction:
            undestroyed_multiboundary.append(geom)
    return MultiLineString(undestroyed_multiboundary)


class Mouseoverable:
    def __init__(self, id: str, boundaries: List[Polygon], attributes: Dict):
        self.id = id
        self.boundaries = []
        for boundary in boundaries:
            if type(boundary) is Polygon:
                self.boundaries.append(LineString(boundary.exterior))
            elif type(boundary) is MultiPolygon:
                for poly in list(boundary.geoms):
                    self.boundaries.append(LineString(poly.exterior))
        self.attributes = attributes


def to_js_strings(in_direction_all_in_direction):
    in_direction = in_direction_all_in_direction[0]
    all_in_direction = in_direction_all_in_direction[1]
    connection_grid = in_direction_all_in_direction[2]
    all_data_string = in_direction_all_in_direction[3] + ";"
    centroid_arr = in_direction_all_in_direction[4]
    in_direction_string = (
            "["
            + ",".join(
        [f'[{",".join([str(int(el_el)) for el_el in el])}]' for el in in_direction]
    )
            + "];"
    )
    all_in_direction_string = (
            "["
            + ",".join(
        [
            "["
            + ",".join(
                [
                    "[" + ",".join([str(int(el_el_el)) for el_el_el in el_el]) + "]"
                    for el_el in el
                ]
            )
            + "]"
            for el in all_in_direction
        ]
    )
            + "];"
    )
    connection_grid_string = (
            "["
            + ",".join(
        [
            "["
            + ",".join(
                [
                    "[" + ",".join([str(int(el_el_el)) for el_el_el in el_el]) + "]"
                    for el_el in el
                ]
            )
            + "]"
            for el in connection_grid
        ]
    )
            + "];"
    )
    print([list(el) for el in centroid_arr])
    centroid_arr_string = (
            "["
            + ",".join(
        ["[" + ",".join([str(el_el) for el_el in el]) + "]" for el in centroid_arr]
    )
            + "];"
    )
    print("centroid_arr_string", centroid_arr_string)
    return (
        in_direction_string,
        all_in_direction_string,
        connection_grid_string,
        all_data_string,
        centroid_arr_string,
    )


cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]


shapefile_path = "shapefiles/quadria_districts_with_geom_info.shp"
gdf = gpd.read_file(shapefile_path)

district_data = []

USE_ALL_STATES = True

features_to_use = [4, 21, 30, 34, 35, 51]  # northeastern states
ii = 0
district_centroids = []
total_bounds = [None, None, None, None]
for i, feature in gdf.iterrows():
    if i not in features_to_use and not USE_ALL_STATES:
        continue
    if feature.geometry is None:
        continue
    # print(feature["STATENAME"]+" "+ordinal(int(feature["DISTRICT"])) if int(feature["DISTRICT"])>0 else feature["STATENAME"])
    # district_data.append(
    #     Mouseoverable(id=ii, boundaries=[feature.geometry], attributes={"NAME": feature["STATENAME"]+" "+ordinal(int(feature["DISTRICT"])) if int(feature["DISTRICT"])>0 else feature["STATENAME"]})
    # )
    # print(feature["STATENAME"]+" "+ordinal(int(feature["DISTRICT"])) if int(feature["DISTRICT"])>0 else feature["STATENAME"])
    try:
        name_attribute = feature["STATENAME"]+" "+ordinal(int(feature["DISTRICT"])) if int(feature["DISTRICT"])>0 else feature["STATENAME"]
    except:
        if "NAME" in feature:
            name_attribute = feature["NAME"]
        else:
            name_attribute = f"{ii}"
    district_data.append(
        Mouseoverable(
            id=ii, boundaries=[feature.geometry], attributes={"NAME": name_attribute}
        )
    )
    ii += 1
    print(feature.geometry.bounds)
    bounds = feature.geometry.bounds
    for bound_i in range(4):
        if total_bounds[bound_i] is None:
            total_bounds[bound_i] = bounds[bound_i]
    for bound_i in range(4):
        if bound_i < 2:
            total_bounds[bound_i] = min(total_bounds[bound_i], bounds[bound_i])
        else:
            total_bounds[bound_i] = max(total_bounds[bound_i], bounds[bound_i])
    # district_centroids.append(np.array(feature.geometry.bounds))
# print("district_centroids centroid",np.mean(district_centroids,axis=0)[:,0])
total_bounds = np.array(total_bounds)
print("total_bounds", total_bounds)
# just east coast
stereo_center_of_lines = np.array([-78, 38]) * math.pi / 180
# center is kansas
stereo_center_of_lines = np.array([-98, 38]) * math.pi / 180
stereo_center_of_lines = (total_bounds[:2] + total_bounds[2:]) / 2 * math.pi / 180


def xxy_to_xxy_stereo(xxy, stereo_center):
    xxy_equi = np.array(xxy).T * math.pi / 180
    segment_dists = np.linalg.norm(xxy_equi - np.roll(xxy_equi, 1, axis=0), axis=1)
    while np.max(segment_dists) > 0.01:
        where_break = np.argmax(segment_dists)
        xxy_equi = np.vstack(
            [
                xxy_equi[:where_break],
                (xxy_equi[where_break - 1] + xxy_equi[where_break]) / 2,
                xxy_equi[where_break:],
            ]
        )
        segment_dists = np.linalg.norm(xxy_equi - np.roll(xxy_equi, 1, axis=0), axis=1)
    xxy_stereo = equi_to_stereo(xxy_equi, stereo_center)
    return xxy_stereo


def stereo_bounds(gdf_, stereo_center):
    total_bounds = [0, 0, 0, 0]
    for i, feature in gdf_.iterrows():
        if type(feature.geometry) is Polygon:
            geom = feature.geometry
            xxy_stereo = xxy_to_xxy_stereo(geom.exterior.xy, stereo_center)
            # plt.plot(xxy_stereo[:,0], xxy_stereo[:,1])
            for bound_i in range(2):
                total_bounds[bound_i] = min(
                    total_bounds[bound_i], np.min(xxy_stereo[:, bound_i])
                )
                total_bounds[2 + bound_i] = max(
                    total_bounds[2 + bound_i], np.max(xxy_stereo[:, bound_i])
                )
        elif type(feature.geometry) is MultiPolygon:
            for geom in feature.geometry.geoms:
                xxy_stereo = xxy_to_xxy_stereo(geom.exterior.xy, stereo_center)
                # plt.plot(xxy_stereo[:,0], xxy_stereo[:,1])
                for bound_i in range(2):
                    total_bounds[bound_i] = min(
                        total_bounds[bound_i], np.min(xxy_stereo[:, bound_i])
                    )
                    total_bounds[2 + bound_i] = max(
                        total_bounds[2 + bound_i], np.max(xxy_stereo[:, bound_i])
                    )
    return np.array(total_bounds)


# plt.show()
total_bounds = stereo_bounds(gdf, stereo_center_of_lines)
new_stereo_center = (total_bounds[:2] + total_bounds[2:]) / 2
print(
    "rsddsr",
    stereo_to_equi(new_stereo_center[None, :], stereo_center_of_lines),
    stereo_center_of_lines,
)
stereo_center_of_lines = stereo_to_equi(
    new_stereo_center[None, :], stereo_center_of_lines
)[0]
print("stereo_center_of_lines", stereo_center_of_lines)
# for use in deciding SVG bounds
total_bounds = 1000 * stereo_bounds(gdf, stereo_center_of_lines)
ultimate_bounds = [
    total_bounds[0],
    (total_bounds[1] + total_bounds[3]) / 2
    - 600 / 600 * (total_bounds[2] - total_bounds[0]) / 2,
    total_bounds[2],
    (total_bounds[1] + total_bounds[3]) / 2
    + 600 / 600 * (total_bounds[2] - total_bounds[0]) / 2,
]
print("new_ total_bounds", total_bounds, ultimate_bounds)
# get in_direction data from utils.py
in_direction, all_in_direction, connection_grid, all_data_string, centroid_arr = (
    get_closest_in_direction(shapefile_path)
)
in_direction_dict = {}
# ["MapSVG_states_trivial_1", "MapSVG_states_trivial_2"]
# in_direction_dict["MapSVG_states_trivial_1"] = to_js_strings(
#     get_closest_in_direction(shapefile_path)
# )
in_direction_dict["MapSVG_states_trivial_2"] = to_js_strings(
    get_closest_in_direction(
        shapefile_path,
        secondary_movement=True,
    )
)
# in_direction_dict["MapSVG_states_trivial_3"] = to_js_strings(
#     get_closest_in_direction(shapefile_path, secondary_movement=True, remove_loops=True)
# )

print(in_direction_dict)
# indirection_string,all_in_direction_string = to_js_strings(in_direction, all_in_direction)

svg_name = list(in_direction_dict.keys())[0]

with open("svg_test.html", "w") as svg_f:
    svg_f.write("<!DOCTYPE html>")
    # just use the first svg_name for the style, though I don't remember why it's in here
    svg_f.write(
        f"<style>body {{height: 100%;}} .MapSVG{{border: 2px solid black;}}.MapSVG:has(>#{svg_name}_obj_{district_data[-1].id}:hover) ~ .MapSVG>#{svg_name}_obj_{district_data[-1].id}{{fill: green;}}</style>"
    )
    for svg_name_i, svg_name in enumerate(list(in_direction_dict.keys())):
        svg_f.write(
            f'\n<svg class="MapSVG" id="{svg_name}" style="position: absolute;left: 0px;top: {0 * svg_name_i}px" width="800" height="800" viewBox="{ultimate_bounds[0]} {-ultimate_bounds[3]} {ultimate_bounds[2]} {-ultimate_bounds[1]}">'
        )
        svg_f.write('\n<g id="map_group">')
        svg_f.write(f"\n<defs>")
        total_number_of_points = 0
        for i in range(len(district_data)):
            print(district_data[i])
            print(
                f"boundaries: {i}/{len(district_data)} ({district_data[i].boundaries})"
            )
            welded = None
            for bound in district_data[i].boundaries:
                if welded is None:
                    welded = bound
                else:
                    if type(bound) is LineString or type(bound) is MultiLineString:
                        if not bound.is_empty:
                            welded = welded.union(bound)
                    else:
                        if type(bound) is GeometryCollection:
                            for geom in bound.geoms:
                                if (
                                        type(geom) is LineString
                                        or type(geom) is MultiLineString
                                ):
                                    if not geom.is_empty:
                                        welded = welded.union(geom)
            svg_f.write(f'\n<g id="{svg_name}_{district_data[i].id}">')

            M_strings = []
            if type(welded) is MultiLineString:
                merged_welded = ops.linemerge(welded)
                if type(merged_welded) is GeometryCollection:
                    new_merged_welded_list = []
                    for geom in merged_welded.geoms:
                        if type(geom) is LineString:
                            new_merged_welded_list.append(geom)
                        elif type(geom) is MultiLineString:
                            for poly in geom.geoms:
                                new_merged_welded_list.append(poly)
                    merged_welded = MultiLineString(new_merged_welded_list)

                for geom in shapely.polygonize(
                        [
                            LineString([*[coord for coord in np.array(geom_i.coords.xy).T]])
                            for geom_i in list(
                            (
                                    merged_welded
                                    if (type(merged_welded) is MultiLineString)
                                    else (
                                            MultiLineString([merged_welded])
                                            if (type(merged_welded) is LineString)
                                            else merged_welded
                                    )
                            ).geoms
                        )
                        ]
                ).geoms:
                    xxy = geom.exterior.coords.xy
                    xxy_stereo = xxy_to_xxy_stereo(xxy, stereo_center_of_lines)

                    joined_L_parts = " L ".join(
                        [
                            f"{coord[0] * 1000} {-coord[1] * 1000}"
                            for coord in xxy_stereo
                        ]
                    )
                    M_strings.append(f"M {joined_L_parts} Z")
                    total_number_of_points += xxy_stereo.shape[0]
            else:
                geom = welded
                xxy = geom.coords.xy
                xxy_stereo = xxy_to_xxy_stereo(xxy, stereo_center_of_lines)
                joined_L_parts = " L ".join(
                    [f"{coord[0] * 1000} {-coord[1] * 1000}" for coord in xxy_stereo]
                )
                M_strings.append(f"M {joined_L_parts} Z")
                total_number_of_points += xxy_stereo.shape[0]
                plt.plot(xxy_stereo[:, 0], xxy_stereo[:, 1], c=cycle[i % len(cycle)])
            joined_M_strings = " ".join(M_strings)
            path_svg = f'\n<path id="path1" fill-rule="evenodd" d="{joined_M_strings}"/>'
            # add centroid center
            centroid_projected = (
                    equi_to_stereo(
                        np.array(welded.centroid.xy)[None, :] * math.pi / 180,
                        stereo_center_of_lines,
                    )[0]
                    * 1000
            )
            centroid_projected[1] *= -1
            print(centroid_projected)
            path_svg = (
                    path_svg
                    + f'\n<circle class="{svg_name}_district_circle" id="{svg_name}_circle_{district_data[i].id}" district_id="{district_data[i].id}" r="5" cx="{centroid_projected[0]}" cy="{centroid_projected[1]}" fill="red" />'
            )
            svg_f.write(path_svg)
            svg_f.write(f"\n</g>")
            svg_f.write(
                f'\n<mask id="{svg_name}_inner_{district_data[i].id}"><rect x="-1000" y="-1000" width="2000" height="2000" fill="black"/><use href = "#{svg_name}_{district_data[i].id}" fill="white"/> </mask>'
            )
        for color in ["black", "white", "blue", "red"]:
            svg_f.write(
                f"""
                <marker 
                  id='arrowhead_{color}' 
                  orient="auto" 
                  markerWidth='5' 
                  markerHeight='4' 
                  refX='3' 
                  refY='2'
                >
                  <path d='M0,0 V4 L4,2 Z' fill="{color}" />
                </marker>
            """
            )
        svg_f.write(f"</defs>")
        for i in range(len(district_data)):
            svg_f.write(
                f'\n<use id="{svg_name}_obj_{district_data[i].id}" href="#{svg_name}_{district_data[i].id}" stroke-width="0.1" stroke="red" fill="none" mask="url(#{svg_name}_inner_{district_data[i].id})" stroke-linejoin="round"/>'
            )
        for i in range(len(district_data)):
            svg_f.write(
                f'\n<use id="{svg_name}_g_{district_data[i].id}" href = "#{svg_name}_{district_data[i].id}" fill="white" fill-opacity="0.0" stroke="none"/>'
            )

        cardinal_directions = ["north", "east", "south", "west"]
        cardinal_colors = ["black", "blue", "red", "white"]
        for i in range(len(district_data)):
            for cardinal_i in range(len(cardinal_directions)):
                svg_f.write(
                    f'\n<line id="{svg_name}_{district_data[i].id}_{cardinal_directions[cardinal_i]}_miniarrowline" x1="0" y1="0" x2="0" y2="0" style="pointer-events: none;stroke:{cardinal_colors[cardinal_i]};stroke-width:2" marker-end="url(#arrowhead_{cardinal_colors[cardinal_i]})" />'
                )
        for cardinal_i in range(len(cardinal_directions)):
            svg_f.write(
                f'\n<line id="{svg_name}_{cardinal_directions[cardinal_i]}_arrowline" x1="0" y1="0" x2="0" y2="0" style="pointer-events: none;stroke:{cardinal_colors[cardinal_i]};stroke-width:2" marker-end="url(#arrowhead_{cardinal_colors[cardinal_i]})" />'
            )
        svg_f.write("</svg>")

        #### MOUSEOVER
        svg_f.write(
            f'\n<svg class="MapSVG_mouseover" id="{svg_name}_mouseover_area" style="position: absolute;left: 0px;top: {0 * svg_name_i}px;z-index:2;pointer-events: none;user-select: none;-webkit-user-drag: none;" width="800" height="800" viewBox="0 0 800 800">'
        )
        svg_f.write(
            f'\n<text id="{svg_name}_label_text" x="{50}" y="50" fill="black" fill-opacity="1.0" stroke="black">None</text>'
        )
        # mouseover rect to display information
        svg_f.write(
            f'\n<rect id="{svg_name}_mouseover_rect" x="0" y="0" width="10" height="10" rx="10" ry="10" fill="black" fill-opacity="0.0" stroke="none" style="z-index:2;"/>'
        )
        svg_f.write(
            f'\n<text id="{svg_name}_mouseover_text" x="0" y="0" fill="white" fill-opacity="0.0" stroke="none" style="z-index:3;" dominant-baseline="hanging">None</text>'
        )
        svg_f.write("\n</svg>")

    js = """
    var hovered_district = {};
    var clicked_district = {};
    var nav_arr_dict = {};
    var map_svg_list = document.getElementsByClassName("MapSVG");
    for(var map_i=0;map_i<map_svg_list.length;map_i++){
        let map_id = map_svg_list[map_i].id;
        hovered_district[map_id] = "";
        clicked_district[map_id] = "";
        nav_arr_dict[map_id] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    }
    var in_direction_dict = {};
    var all_in_direction_dict = {};
    var connection_grid_dict = {};
    var centroid_arr_dict = {};
    var all_data_string_dict = {};
    var mouseoverable_list = [];
    var mouseoverable_dict = [];
"""
    for key in in_direction_dict:
        js = js + (
                """
    
            in_direction_dict[\""""
                + key
                + """\"] = """
                + in_direction_dict[key][0]
                + """;
        all_in_direction_dict[\""""
                + key
                + """\"] = """
                + in_direction_dict[key][1].replace(
            "]],[[",
            "]],\n"
            + " "
            * (len('        all_in_direction_dict["') + len(key) + len('"] = ['))
            + "[[",
        )
                + """;
        connection_grid_dict[\""""
                + key
                + """\"] = """
                + in_direction_dict[key][2]
                + """;
        all_data_string_dict[\""""
                + key
                + """\"] = """
                + in_direction_dict[key][3]
                + """;
        centroid_arr_dict[\""""
                + key
                + """\"] = """
                + in_direction_dict[key][4]
                + """;
            """
        )

    for i in range(len(district_data)):
        attribute_string = ",".join(
            [
                f'{key}:"{district_data[i].attributes[key]}"'
                for key in district_data[i].attributes
            ]
        )
        print(attribute_string)
        js += f'mouseoverable_list.push({{id:"{district_data[i].id}", attributes:{{{attribute_string}}}}});\n'
        js += f'mouseoverable_dict["{district_data[i].id}"] = {{id:"{district_data[i].id}", attributes:{{{attribute_string}}}}}\n'
    js += (
            """
    const mouseover_buffer = 12;
    const ARROWHEAD_PIXEL_GAP = 10;
    const STANDARD_FILL = false;
    function rgbToHex(r, g, b) {
        return "#" + ((1 << 24) + (r << 16) + (g << 8) + b).toString(16).slice(1).toUpperCase();
    }
    function reset_hovered(hovered_district_id,map_id){
        let district_name = "";
        if(hovered_district_id in mouseoverable_dict){
            district_name = mouseoverable_dict[hovered_district_id].attributes.NAME;
        }
        document.getElementById(map_id+'_label_text').innerHTML=''+district_name;
        document.getElementById(map_id+"_mouseover_text").innerHTML = ''+district_name;
        if(hovered_district_id==''){
            document.getElementById(map_id+"_mouseover_rect").setAttribute('fill-opacity', "0.0");
            document.getElementById(map_id+"_mouseover_text").setAttribute('fill-opacity', "0.0");
        }else{
            document.getElementById(map_id+"_mouseover_rect").setAttribute('fill-opacity', "0.1");
            //document.getElementById(map_id+"_mouseover_rect").setAttribute('width', "100");
            //document.getElementById(map_id+"_mouseover_rect").setAttribute('height', "40");
            document.getElementById(map_id+"_mouseover_rect").setAttribute('width', ""+(mouseover_buffer*2+document.getElementById(map_id+"_mouseover_text").getBBox().width));
            document.getElementById(map_id+"_mouseover_rect").setAttribute('height', ""+(mouseover_buffer*2-5+document.getElementById(map_id+"_mouseover_text").getBBox().height));
            document.getElementById(map_id+"_mouseover_text").setAttribute('fill-opacity', "1.0");
        }
    }
    for(var map_i=0;map_i<map_svg_list.length;map_i++){
        let map_id = map_svg_list[map_i].id;
        for(var i=0;i<mouseoverable_list.length;i+=1){
            var current_member_code = mouseoverable_list[i].id;
            var fill_color = "grey";
            if(current_member_code!=''){
                fill_color = "grey";
            }
            if(STANDARD_FILL){
                document.getElementById(map_id+'_obj_'+mouseoverable_list[i].id).style.fill = fill_color;
            }
            document.getElementById(map_id+'_g_'+mouseoverable_list[i].id).addEventListener('mouseenter', function(e) {
                var current_member_code = e.target.id.substring(3+map_id.length);
                var hovered_district_was = ''+hovered_district[map_id];
                hovered_district[map_id] = current_member_code;
    
                var fill_color = "orange";
                if(STANDARD_FILL){
                    if(hovered_district_was!=''){
                        if(hovered_district_was==clicked_district[map_id]){
                            document.getElementById(map_id+'_obj_'+hovered_district_was).style.fill = "green";
                        }else{
                            document.getElementById(map_id+'_obj_'+hovered_district_was).style.fill = "grey";
                        }
                    }
                    document.getElementById(map_id+'_obj_'+current_member_code).style.fill = fill_color;
                }
                reset_hovered(hovered_district[map_id],map_id);
            });
            document.getElementById(map_id+'_g_'+mouseoverable_list[i].id).addEventListener('mouseleave', function(e) {
                var current_member_code = e.target.id.substring(3+map_id.length);
                if(STANDARD_FILL){
                    if(current_member_code==clicked_district[map_id]){
                        document.getElementById(map_id+'_obj_'+current_member_code).style.fill = "green";
                    }else{
                        document.getElementById(map_id+'_obj_'+current_member_code).style.fill = "grey";
                    }
                }
                if(current_member_code==hovered_district[map_id]){
                    hovered_district[map_id] = '';
                }
                reset_hovered(hovered_district[map_id],map_id);
            });
        }
    }
    let view_box_bounds_by_id = {};
    let ultimate_bounds_by_id = {};
    
    var zooms_initialized = false;
    
    function initialize_view_box_dicts(){
        var map_svg_list = document.getElementsByClassName("MapSVG");
        for(var map_i=0;map_i<map_svg_list.length;map_i++){
            let map_id = map_svg_list[map_i].id;
            if (!(map_id in view_box_bounds_by_id)) {
                var rect = document.getElementById(map_id).getBoundingClientRect();
                let box_width = """
            + str((ultimate_bounds[2] - ultimate_bounds[0]))
            + """;
            let box_height = box_width*(rect.bottom-rect.top)/(rect.right-rect.left)
            let center_pos = ["""
            + str((ultimate_bounds[0] + ultimate_bounds[2]) / 2)
            + ""","""
            + str(-(ultimate_bounds[1] + ultimate_bounds[3]) / 2)
            + """];
            view_box_bounds_by_id[map_id] = [center_pos[0]-box_width/2,center_pos[1]-box_height/2,box_width,box_height];
            ultimate_bounds_by_id[map_id] = [center_pos[0]-box_width/2,center_pos[1]-box_height/2,box_width,box_height];
        }
    }
}

function initialize_zooms(){
    var map_svg_list = document.getElementsByClassName("MapSVG");
    initialize_view_box_dicts();
    select_district();
    for(var map_i=0;map_i<map_svg_list.length;map_i++){
        let map_id = map_svg_list[map_i].id;
        // REDRAW ARROWS
        var cardinal_directions = ["east","north","west","south"];
        var district_circles = document.getElementsByClassName(map_id+"_district_circle");
        for(var theta_i=0;theta_i<cardinal_directions.length;theta_i++){
            document.getElementById(map_id+"_"+cardinal_directions[theta_i]+"_arrowline").style["stroke-width"] = 12*view_box_bounds_by_id[map_id][2]/1200+"px";
            if(document.getElementById(map_id+"_"+cardinal_directions[theta_i]+"_arrowline").getAttribute('stroke-dasharray') !== null){
                var dasharray_list = document.getElementById(map_id+"_"+cardinal_directions[theta_i]+"_arrowline").getAttribute('stroke-dasharray').split(" ");
                var dasharray_max = parseFloat(document.getElementById(map_id+"_"+cardinal_directions[theta_i]+"_arrowline").getAttribute('stroke-dasharray-max'));
                for(var dash_i=0;dash_i<dasharray_list.length;dash_i++){
                    dasharray_list[dash_i] = parseFloat(dasharray_list[dash_i]);
                }
                var max_in_list = Math.max(...dasharray_list);
                for(var dash_i=0;dash_i<dasharray_list.length;dash_i++){
                    dasharray_list[dash_i] = dasharray_list[dash_i]/max_in_list*dasharray_max*6*view_box_bounds_by_id[map_id][2]/1200;
                }
                if(dasharray_list.join(" ")!="NaN"){
                    document.getElementById(map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute("stroke-dasharray",dasharray_list.join(" "));
                }
            }





            for(var district_i=0;district_i<district_circles.length;district_i++){
                var district_id = district_circles[district_i].getAttribute("district_id");
                document.getElementById(map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").style["stroke-width"] = 4*view_box_bounds_by_id[map_id][2]/1200+"px";
                if(document.getElementById(map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").getAttribute('stroke-dasharray') !== null){
                    var dasharray_list = document.getElementById(map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").getAttribute('stroke-dasharray').split(" ");
                    var dasharray_max = parseFloat(document.getElementById(map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").getAttribute('stroke-dasharray-max'));
                    for(var dash_i=0;dash_i<dasharray_list.length;dash_i++){
                        dasharray_list[dash_i] = parseFloat(dasharray_list[dash_i]);
                    }
                    var max_in_list = Math.max(...dasharray_list);
                    for(var dash_i=0;dash_i<dasharray_list.length;dash_i++){
                        dasharray_list[dash_i] = dasharray_list[dash_i]/max_in_list*dasharray_max*6*view_box_bounds_by_id[map_id][2]/1200;
                    }
                    if(dasharray_list.join(" ")!="NaN"){
                        document.getElementById(map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute("stroke-dasharray",dasharray_list.join(" "));
                    }
                }
            }
        }

        var district_circles = document.getElementsByClassName(map_id+"_district_circle");
        for(var i=0;i<district_circles.length;i++){
            district_circles[i].setAttribute("r",20*view_box_bounds_by_id[map_id][2]/1200);
        }
    }
    //arrow endpoint should be slightly before target end, by a few pixels
    zooms_initialized = true;
}

var map_svg_list = document.getElementsByClassName("MapSVG");
for(var i=0;i<map_svg_list.length;i++){
    map_svg_list[i].addEventListener('wheel', preventScroll, {passive: false});
}
window.addEventListener("keydown", function(e) {
    // space and arrow keys
    if([37, 38, 39, 40].indexOf(e.keyCode) > -1) {
        e.preventDefault();
    }
}, false);
function preventScroll(e){
    e.preventDefault();

    return false;
}

var current_map_id = "MapSVG_states_trivial_2";
window.addEventListener("wheel", event => {
    var map_svg_list = document.getElementsByClassName("MapSVG");
    initialize_view_box_dicts();
    for(var i=0;i<map_svg_list.length;i++){
        let map_id = map_svg_list[i].id;
        let view_box_bounds = view_box_bounds_by_id[map_id];
        let ultimate_bounds = ultimate_bounds_by_id[map_id];
        var rect = document.getElementById(map_id).getBoundingClientRect();
        if(event.clientX>rect.left&&event.clientX<rect.right&&event.clientY>rect.top&&event.clientY<rect.bottom){
            event.preventDefault();
            event.stopPropagation();
            current_map_id = map_id;
            var mouse_x = (event.clientX - rect.left)/(rect.right-rect.left);
            var mouse_y = (event.clientY - rect.top)/(rect.bottom-rect.top);
            var mouse_x_in_svg = view_box_bounds[0]+(view_box_bounds[2])*mouse_x;
            var mouse_y_in_svg = view_box_bounds[1]+(view_box_bounds[3])*mouse_y;
            var product = 1.002**event.deltaY;
            view_box_bounds = [mouse_x_in_svg+(view_box_bounds[0]-mouse_x_in_svg)*product,mouse_y_in_svg+(view_box_bounds[1]-mouse_y_in_svg)*product,view_box_bounds[2]*product,view_box_bounds[3]*product]
            view_box_bounds[0] = Math.max(ultimate_bounds[0],view_box_bounds[0]);
            view_box_bounds[1] = Math.max(ultimate_bounds[1],view_box_bounds[1]);
            view_box_bounds[2] = Math.min(view_box_bounds[2],ultimate_bounds[2]);
            view_box_bounds[3] = Math.min(view_box_bounds[3],ultimate_bounds[3]);
            if(view_box_bounds[2]+view_box_bounds[0]>ultimate_bounds[2]+ultimate_bounds[0]){
                view_box_bounds[0] = Math.min(ultimate_bounds[2]+ultimate_bounds[0],view_box_bounds[2]+view_box_bounds[0])-view_box_bounds[2];
            }
            if(view_box_bounds[3]+view_box_bounds[1]>ultimate_bounds[3]+ultimate_bounds[1]){
                view_box_bounds[1] = Math.min(ultimate_bounds[3]+ultimate_bounds[1],view_box_bounds[3]+view_box_bounds[1])-view_box_bounds[3];
            }
            //view_box_bounds[3] = Math.min(ultimate_bounds[3],view_box_bounds[3]);
            document.getElementById(map_id).setAttribute("viewBox",""+view_box_bounds[0]+" "+view_box_bounds[1]+" "+view_box_bounds[2]+" "+view_box_bounds[3]);

            view_box_bounds_by_id[map_id] = view_box_bounds;
        }
    }
    initialize_zooms();

});

addEventListener("mousemove", (event) => {
    for(var i=0;i<map_svg_list.length;i++){
        let map_id = map_svg_list[i].id;
        var rect = document.getElementById(map_id).getBoundingClientRect();
        if(event.clientX>rect.left&&event.clientX<rect.right&&event.clientY>rect.top&&event.clientY<rect.bottom){
            current_map_id = map_id;
        }


        var rect = document.getElementById(map_id).getBoundingClientRect();
        var mouse_x = event.clientX - rect.left;
        var mouse_y = event.clientY - rect.top;
        document.getElementById(map_id+"_mouseover_rect").setAttribute('x', ""+mouse_x);
        document.getElementById(map_id+"_mouseover_rect").setAttribute('y', ""+mouse_y);
        document.getElementById(map_id+"_mouseover_text").setAttribute('x', ""+(mouseover_buffer+mouse_x));
        document.getElementById(map_id+"_mouseover_text").setAttribute('y', ""+(mouseover_buffer+mouse_y));
    }
});

function color_map(map_id,new_position){

    var district_circles = document.getElementsByClassName(map_id+"_district_circle");
    for(var district_i=0;district_i<district_circles.length;district_i++){
        var district_id = district_circles[district_i].getAttribute("district_id").toString();
        var district_obj = document.getElementById(map_id+'_obj_'+district_id);
        var district_circle = district_circles[district_i];
        if(district_id==new_position.toString()){
            district_obj.style.fill = "#aa0000";
        }else{
            district_obj.style.fill = "#3333ff";
        }
        district_circle.style.fill = district_obj.style.fill;
    }
}

var already_passed = [];
var scores_storage = [];
var score_stacking = 1;
var prev = [];
function get_in_direction(map_id,state_i,theta_i,simple=true,simulate=false) {
    var in_both = [];
    var new_position = -1;

    var cx = parseFloat(document.getElementById(map_id+"_circle_"+state_i).getAttribute("cx"));
    var cy = parseFloat(document.getElementById(map_id+"_circle_"+state_i).getAttribute("cy"));

    centroid_arr_dict[map_id];
    connection_grid_dict[map_id];
    var district_circles = document.getElementsByClassName(map_id+"_district_circle");
    var cleaned_scores = [];
    var score_storage_was = scores_storage.slice();
    var score_stacking_was = score_stacking+0;
    if(!simple){
        var scores = [];
        var best_score = 45;
        var best_score_i = state_i;
        for(var district_i=0;district_i<district_circles.length;district_i++){
            var district_id = parseInt(district_circles[district_i].getAttribute("district_id"));
            var diff_x = centroid_arr_dict[map_id][district_id][0]-centroid_arr_dict[map_id][state_i][0];
            var diff_y = centroid_arr_dict[map_id][district_id][1]-centroid_arr_dict[map_id][state_i][1];
            var angle = 180/Math.PI*((Math.atan2(diff_y,diff_x)+Math.PI*4+Math.PI/4+Math.PI/2*(2-theta_i))%(Math.PI*2));
            angle = Math.max(0,90-Math.abs(angle-180));
            var score = connection_grid_dict[map_id][state_i][district_id][theta_i]+connection_grid_dict[map_id][district_id][state_i][(theta_i+2)%4]/2;
            cleaned_scores.push(score+0);
            if(scores_storage.length<=district_id){
                scores_storage.push(score);
            }else{
                score += scores_storage[district_id]*2/3;
                scores_storage[district_id] = score;
            }
            var multiplier = 1;
            if(prev[0]==district_id&&(prev[1]+theta_i)%2==1){
                multiplier /= 5;
            }
            if(state_i==district_id){
                multiplier = 0;
            }
            multiplier /= score_stacking;
            if(connection_grid_dict[map_id][state_i][district_id][theta_i]+connection_grid_dict[map_id][district_id][state_i][(theta_i+2)%4]/2==0){
                multiplier = 0;
            }
            var district_obj = document.getElementById(map_id+'_obj_'+district_id);
            var district_circle = district_circles[district_i];
            //district_obj.style.fill = rgbToHex(Math.floor((multiplier*(angle+6*score))/2/360*255),Math.floor(angle*255/360),255-Math.floor((multiplier*(angle+6*score))/2/360*255));
            scores.push(multiplier*(angle+8*score));
            if(multiplier*(angle+8*score)>best_score){
                best_score = multiplier*(angle+8*score);
                best_score_i = district_id;
            }
        }
        score_stacking = score_stacking*2/3+1;
        new_position = best_score_i;
        if(new_position!=state_i){
            for(var district_i=0;district_i<district_circles.length;district_i++){
                var district_id = parseInt(district_circles[district_i].getAttribute("district_id"));
                var district_obj = document.getElementById(map_id+'_obj_'+district_id);
                district_obj.style.fill = rgbToHex(Math.floor(scores[district_i]/best_score*255*12/15+255*3/15),Math.floor(scores[district_i]/best_score*255*12/15+255*3/15),255-Math.floor(scores[district_i]/best_score*255));
                if(district_id==new_position.toString()){
                    district_obj.style.fill = "#aa0000";
                }
                district_circles[district_i].style.fill = district_obj.style.fill;
            }
        }
    }else if(simulate){
        var scores = [];
        var best_score = 45;
        var best_score_i = state_i;
        var sim_score_stacking = score_stacking+0;
        for(var district_i=0;district_i<district_circles.length;district_i++){
            var district_id = parseInt(district_circles[district_i].getAttribute("district_id"));
            var diff_x = centroid_arr_dict[map_id][district_id][0]-centroid_arr_dict[map_id][state_i][0];
            var diff_y = centroid_arr_dict[map_id][district_id][1]-centroid_arr_dict[map_id][state_i][1];
            var angle = 180/Math.PI*((Math.atan2(diff_y,diff_x)+Math.PI*4+Math.PI/4+Math.PI/2*(2-theta_i))%(Math.PI*2));
            angle = Math.max(0,90-Math.abs(angle-180));
            var score = connection_grid_dict[map_id][state_i][district_id][theta_i]+connection_grid_dict[map_id][district_id][state_i][(theta_i+2)%4]/2;

            if(scores_storage.length<=district_id){
            }else{
                score += scores_storage[district_id]*2/3;
            }
            var multiplier = 1;
            if(prev[0]==district_id&&(prev[1]+theta_i)%2==1){
                multiplier /= 5;
            }
            if(state_i==district_id){
                multiplier = 0;
            }
            multiplier /= sim_score_stacking;
            if(connection_grid_dict[map_id][state_i][district_id][theta_i]+connection_grid_dict[map_id][district_id][state_i][(theta_i+2)%4]/2==0){
                multiplier = 0;
            }
            var district_obj = document.getElementById(map_id+'_obj_'+district_id);
            var district_circle = district_circles[district_i];
            scores.push(multiplier*(angle+8*score));
            if(multiplier*(angle+8*score)>best_score){
                best_score = multiplier*(angle+8*score);
                best_score_i = district_id;
            }
        }
        sim_score_stacking = sim_score_stacking*2/3+1;
        new_position = best_score_i;

    }

    if(!simple){
    }else if(!simulate){
        if(state_i==clicked_district[map_id]){
        }else{
        }
    }
    if(!simple){
        //color_map(map_id,new_position);
        if((prev[1]+theta_i)%2==1){
            score_stacking = 1;
            scores_storage = cleaned_scores;
        }
        scores_storage[state_i] = 0;
        if(new_position==state_i){
            scores_storage = score_storage_was.slice();
            score_stacking = score_stacking_was+0;
        }else{
            prev = [state_i,theta_i+0];
        }
    }


    return new_position;
}


document.onkeydown = function(e) {
    var theta_i = 0;
    switch (e.keyCode) {
        case 37:
            theta_i = 2;
            break;
        case 38:
            theta_i = 1;
            break;
        case 39:
            theta_i = 0;
            break;
        case 40:
            theta_i = 3;
            break;
    }
    var new_position = get_in_direction(current_map_id,clicked_district[current_map_id],theta_i,false);

    if(new_position!=-1){
        if(STANDARD_FILL){
            document.getElementById(current_map_id+'_obj_'+clicked_district[current_map_id]).style.fill = "grey";
            document.getElementById(current_map_id+'_obj_'+new_position).style.fill = "green";
        }
        //if(!in_direction[new_position].includes(clicked_district[current_map_id])){
        //    in_direction[new_position][(theta_i+2)%4] = clicked_district[current_map_id];
        //}
        clicked_district[current_map_id] = new_position;
    }
    initialize_zooms();
};
function select_district(){
    initialize_view_box_dicts();
    var cardinal_directions = ["east","north","west","south"];
    var count_dict = {};
    var count_dict_by_district_id = {};
    var viewbox_size = (view_box_bounds_by_id[current_map_id][2])/
                       (ultimate_bounds_by_id[current_map_id][2]);
    for(var theta_i=0;theta_i<cardinal_directions.length;theta_i++){
        if(document.getElementById(current_map_id+"_circle_"+clicked_district[current_map_id])){
            document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('x1', document.getElementById(current_map_id+"_circle_"+clicked_district[current_map_id]).getAttribute('cx'));
            document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('y1', document.getElementById(current_map_id+"_circle_"+clicked_district[current_map_id]).getAttribute('cy'));
            var state_north_of_here = get_in_direction(current_map_id,clicked_district[current_map_id],theta_i,simple=true,simulate=true);
            if(state_north_of_here<0){
                state_north_of_here = clicked_district[current_map_id];
            }
            var target_x1 = parseFloat(document.getElementById(current_map_id+"_circle_"+clicked_district[current_map_id]).getAttribute('cx'));
            var target_y1 = parseFloat(document.getElementById(current_map_id+"_circle_"+clicked_district[current_map_id]).getAttribute('cy'));
            var target_x2 = parseFloat(document.getElementById(current_map_id+"_circle_"+state_north_of_here).getAttribute('cx'));
            var target_y2 = parseFloat(document.getElementById(current_map_id+"_circle_"+state_north_of_here).getAttribute('cy'));

            var full_line_length = Math.sqrt((target_x2-target_x1)**2 + (target_y2-target_y1)**2);
            var display_x2 = target_x2;
            var display_y2 = target_y2;
            if(full_line_length>0){
                display_x2 = target_x1+(target_x2-target_x1)*Math.max(full_line_length/2,full_line_length-(ARROWHEAD_PIXEL_GAP*viewbox_size))/full_line_length;
                display_y2 = target_y1+(target_y2-target_y1)*Math.max(full_line_length/2,full_line_length-(ARROWHEAD_PIXEL_GAP*viewbox_size))/full_line_length;
                document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").style.display = "block";
                document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('target_x2', target_x2);
                document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('target_y2', target_y2);
                document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('x2', display_x2.toString());
                document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('y2', display_y2.toString());
            }else{
                document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").style.display = "none";
            }

            if(state_north_of_here in count_dict){
                count_dict[state_north_of_here] = count_dict[state_north_of_here]+1;
            }else{
                count_dict[state_north_of_here] = 1;
            }


        }
        //MINILINES
        //there's one for each district and cardinal direction pair
        var district_circles = document.getElementsByClassName(current_map_id+"_district_circle");
        for(var district_i=0;district_i<district_circles.length;district_i++){
            var district_id = district_circles[district_i].getAttribute("district_id");

            var state_north_of_here = get_in_direction(current_map_id,district_id,theta_i,simple=true,simulate=true);
            if(state_north_of_here<0){
                state_north_of_here = district_id;
            }
            var target_x1 = parseFloat(document.getElementById(current_map_id+"_circle_"+district_id).getAttribute('cx'));
            var target_y1 = parseFloat(document.getElementById(current_map_id+"_circle_"+district_id).getAttribute('cy'));
            var target_x2 = parseFloat(document.getElementById(current_map_id+"_circle_"+state_north_of_here).getAttribute('cx'));
            var target_y2 = parseFloat(document.getElementById(current_map_id+"_circle_"+state_north_of_here).getAttribute('cy'));

            var full_line_length = Math.sqrt((target_x2-target_x1)**2 + (target_y2-target_y1)**2);

            var display_x_offset = (target_y2-target_y1)/full_line_length;
            var display_y_offset = -(target_x2-target_x1)/full_line_length;
            var display_x1 = target_x1+display_x_offset*ARROWHEAD_PIXEL_GAP/3*viewbox_size;
            var display_y1 = target_y1+display_y_offset*ARROWHEAD_PIXEL_GAP/3*viewbox_size;


            document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline")
                .setAttribute('x1', 
                display_x1);
            document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline")
                .setAttribute('y1', 
                display_y1);

            var display_x2 = target_x2;
            var display_y2 = target_y2;
            if(false&&full_line_length>0){
                display_x2 = target_x1+(target_x2-target_x1)*Math.max(full_line_length/2,full_line_length-(ARROWHEAD_PIXEL_GAP*2*viewbox_size))/full_line_length+display_x_offset*ARROWHEAD_PIXEL_GAP/3*viewbox_size;
                display_y2 = target_y1+(target_y2-target_y1)*Math.max(full_line_length/2,full_line_length-(ARROWHEAD_PIXEL_GAP*2*viewbox_size))/full_line_length+display_y_offset*ARROWHEAD_PIXEL_GAP/3*viewbox_size;
                document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").style.display = "block";
                document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('target_x2', target_x2);
                document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('target_y2', target_y2);
                document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('x2', display_x2.toString());
                document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('y2', display_y2.toString());
            }else{
                document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").style.display = "none";
            }

            if(!(district_id in count_dict_by_district_id)){
                count_dict_by_district_id[district_id] = {};
            }
            if(state_north_of_here in count_dict_by_district_id[district_id]){
                count_dict_by_district_id[district_id][state_north_of_here] = count_dict_by_district_id[district_id][state_north_of_here]+1;
            }else{
                count_dict_by_district_id[district_id][state_north_of_here] = 1;
            }


        }
    }

    var count_dict_2_by_district_id = {};
    var district_circles = document.getElementsByClassName(current_map_id+"_district_circle");
    for(var district_i=0;district_i<district_circles.length;district_i++){
        var district_id = district_circles[district_i].getAttribute("district_id");
        count_dict_2_by_district_id[district_id] = {};
        for(var theta_i=0;theta_i<cardinal_directions.length;theta_i++){
            var state_north_of_here = get_in_direction(current_map_id,district_id,theta_i,simple=true,simulate=true);
            if(state_north_of_here<0){
                state_north_of_here = district_id;
            }
            if(state_north_of_here in count_dict_2_by_district_id[district_id]){
                count_dict_2_by_district_id[district_id][state_north_of_here] = count_dict_2_by_district_id[district_id][state_north_of_here]+1;
            }else{
                count_dict_2_by_district_id[district_id][state_north_of_here] = 1;
            }
            if(count_dict_by_district_id[district_id][state_north_of_here]==1){
                document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('stroke-dasharray',"");
            }
            if(count_dict_by_district_id[district_id][state_north_of_here]==2){
                if(count_dict_2_by_district_id[district_id][state_north_of_here]==1){
                    document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('stroke-dasharray',"4");
                    document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('stroke-dasharray-max',"4");
                }else{
                    document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('stroke-dasharray',"0 4 0");
                    document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('stroke-dasharray-max',"4");
                }
            }
            if(count_dict_by_district_id[district_id][state_north_of_here]==3){
                if(count_dict_2_by_district_id[district_id][state_north_of_here]==1){
                    document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('stroke-dasharray',"3 6");
                    document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('stroke-dasharray-max',"6");
                }else if(count_dict_2_by_district_id[district_id][state_north_of_here]==2){
                    document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('stroke-dasharray',"0 3 3 3");
                    document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('stroke-dasharray-max',"3");
                }else{
                    document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('stroke-dasharray',"0 6 3 0");
                    document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('stroke-dasharray-max',"6");
                }
            }
            if(count_dict_by_district_id[district_id][state_north_of_here]==4){
                if(count_dict_2_by_district_id[district_id][state_north_of_here]==1){
                    document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('stroke-dasharray',"2.5 7.5");
                    document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('stroke-dasharray-max',"7.5");
                }else if(count_dict_2_by_district_id[district_id][state_north_of_here]==2){
                    document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('stroke-dasharray',"0 2.5 2.5 5");
                    document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('stroke-dasharray-max',"5");
                }else if(count_dict_2_by_district_id[district_id][state_north_of_here]==3){
                    document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('stroke-dasharray',"0 5 2.5 2.5");
                    document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('stroke-dasharray-max',"5");
                }else{
                    document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('stroke-dasharray',"0 7.5 2.5 0");
                    document.getElementById(current_map_id+"_"+district_id+"_"+cardinal_directions[theta_i]+"_miniarrowline").setAttribute('stroke-dasharray-max',"7.5");
                }
            }
        }
    }

    var count_dict_2 = {};
    for(var theta_i=0;theta_i<cardinal_directions.length;theta_i++){
        if(document.getElementById(current_map_id+"_circle_"+clicked_district[current_map_id])){
            var state_north_of_here = get_in_direction(current_map_id,clicked_district[current_map_id],theta_i,simple=true,simulate=true);
            if(state_north_of_here<0){
                state_north_of_here = clicked_district[current_map_id];
            }
            if(state_north_of_here in count_dict_2){
                count_dict_2[state_north_of_here] = count_dict_2[state_north_of_here]+1;
            }else{
                count_dict_2[state_north_of_here] = 1;
            }
            if(count_dict[state_north_of_here]==1){
                document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('stroke-dasharray',"");
            }
            if(count_dict[state_north_of_here]==2){
                if(count_dict_2[state_north_of_here]==1){
                    document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('stroke-dasharray',"4");
                    document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('stroke-dasharray-max',"4");
                }else{
                    document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('stroke-dasharray',"0 4 0");
                    document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('stroke-dasharray-max',"4");
                }
            }
            if(count_dict[state_north_of_here]==3){
                if(count_dict_2[state_north_of_here]==1){
                    document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('stroke-dasharray',"3 6");
                    document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('stroke-dasharray-max',"6");
                }else if(count_dict_2[state_north_of_here]==2){
                    document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('stroke-dasharray',"0 3 3 3");
                    document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('stroke-dasharray-max',"3");
                }else{
                    document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('stroke-dasharray',"0 6 3 0");
                    document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('stroke-dasharray-max',"6");
                }
            }
            if(count_dict[state_north_of_here]==4){
                if(count_dict_2[state_north_of_here]==1){
                    document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('stroke-dasharray',"2.5 7.5");
                    document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('stroke-dasharray-max',"7.5");
                }else if(count_dict_2[state_north_of_here]==2){
                    document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('stroke-dasharray',"0 2.5 2.5 5");
                    document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('stroke-dasharray-max',"5");
                }else if(count_dict_2[state_north_of_here]==3){
                    document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('stroke-dasharray',"0 5 2.5 2.5");
                    document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('stroke-dasharray-max',"5");
                }else{
                    document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('stroke-dasharray',"0 7.5 2.5 0");
                    document.getElementById(current_map_id+"_"+cardinal_directions[theta_i]+"_arrowline").setAttribute('stroke-dasharray-max',"7.5");
                }
            }
        }
    }
}
addEventListener("mousedown", (event) => {
    scores_storage = [];
    prev = [];
    already_passed = [];
    score_stacking = 1;
    var rect = document.getElementById(current_map_id).getBoundingClientRect();
    var mouse_x = event.clientX - rect.left;
    var mouse_y = event.clientY - rect.top;
    if(STANDARD_FILL){
        if(hovered_district[current_map_id]!=clicked_district[current_map_id]&&clicked_district[current_map_id]!=""){
            document.getElementById(current_map_id+'_obj_'+clicked_district[current_map_id]).style.fill = "grey";
        }
    }
    clicked_district[current_map_id] = hovered_district[current_map_id];
    document.getElementById(current_map_id+"_mouseover_rect").setAttribute('x', ""+mouse_x);
    document.getElementById(current_map_id+"_mouseover_rect").setAttribute('y', ""+mouse_y);
    document.getElementById(current_map_id+"_mouseover_text").setAttribute('x', ""+(mouseover_buffer+mouse_x));
    document.getElementById(current_map_id+"_mouseover_text").setAttribute('y', ""+(mouseover_buffer+mouse_y));
    initialize_zooms();
    color_map(current_map_id,clicked_district[current_map_id]);
    //document.getElementById("lean_mouseover_text").setAttribute('y', ""+(20+mouse_y))
});

initialize_zooms();
var map_svg_list = document.getElementsByClassName("MapSVG");
for(var i=0;i<map_svg_list.length;i++){
    var map_id = map_svg_list[i].id;
    var district_circles = document.getElementsByClassName(map_id+"_district_circle");
    for(var district_i=0;district_i<district_circles.length;district_i++){
        var district_id = district_circles[district_i].getAttribute("district_id").toString();
        var district_obj = document.getElementById(map_id+'_obj_'+district_id);
        var district_circle = district_circles[district_i];
        if(false){
            district_obj.style.fill = "red";
        }else{
            district_obj.style.fill = "grey";
        }
        district_circle.style.fill = district_obj.style.fill;
    }
}
"""
    )
    svg_f.write(
        f'<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.7.1/jquery.min.js"></script><script>{js}</script>'
    )
