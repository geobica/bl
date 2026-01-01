import itertools
import os
import pickle

import geopandas as gpd
import matplotlib.pyplot as plt
import random
import numpy as np
import shapely
from shapely import MultiPolygon
from shapely.geometry import Point, LineString, Polygon
import math

from shapely.geometry.multilinestring import MultiLineString
import copy
import matplotlib.colors as mcolors


def equi_to_stereo(equi_arr, stereo_center_as_equi_coord):
    lambd_0 = stereo_center_as_equi_coord[0]
    phi_1 = stereo_center_as_equi_coord[1]
    lambd = equi_arr[:, 0]
    phi = equi_arr[:, 1]
    R = 1
    k = (2 * R) / (
        1
        + np.sin(phi_1) * np.sin(phi)
        + np.cos(phi_1) * np.cos(phi) * np.cos(lambd - lambd_0)
    )
    stereo_x = k * np.cos(phi) * np.sin(lambd - lambd_0)
    stereo_y = k * (
        np.cos(phi_1) * np.sin(phi)
        - np.sin(phi_1) * np.cos(phi) * np.cos(lambd - lambd_0)
    )
    stereo_arr = np.vstack([stereo_x, stereo_y]).T
    return stereo_arr


def stereo_to_equi(stereo_arr, stereo_center_as_equi_coord):
    """
    Inverse of equi_to_stereo.

    Parameters
    ----------
    stereo_arr : (N, 2) array
        Stereographic coordinates [x, y]
    stereo_center_as_equi_coord : (2,) array
        [lambda_0, phi_1] in radians

    Returns
    -------
    equi_arr : (N, 2) array
        Equirectangular coordinates [lambda, phi] in radians
    """
    x = stereo_arr[:, 0]
    y = stereo_arr[:, 1]

    lambd_0 = stereo_center_as_equi_coord[0]
    phi_1 = stereo_center_as_equi_coord[1]

    R = 1.0

    rho = np.sqrt(x**2 + y**2)
    c = 2 * np.arctan2(rho, 2 * R)

    sin_c = np.sin(c)
    cos_c = np.cos(c)

    # Latitude
    phi = np.arcsin(
        cos_c * np.sin(phi_1) + (y * sin_c * np.cos(phi_1)) / np.where(rho == 0, 1, rho)
    )

    # Longitude
    lambd = lambd_0 + np.arctan2(
        x * sin_c, rho * np.cos(phi_1) * cos_c - y * np.sin(phi_1) * sin_c
    )

    equi_arr = np.vstack([lambd, phi]).T
    return equi_arr


def connection_grid_first_collision_method(
    gdf, density_used=40, return_rays=False, stereo_center=None, random_points=False
):
    # center is kansas
    if stereo_center is None:
        stereo_center = np.array([-98, 38]) * math.pi / 180
    theta_list = [0, math.pi / 2, math.pi, math.pi * 3 / 2]
    connection_grid = np.zeros([gdf.shape[0], gdf.shape[0], 4])
    centroid_list = []
    connection_grid_i = 0
    ray_list = []

    projected_feature_geoms = []
    for i, feature in gdf.iterrows():
        if feature.geometry is None:
            continue
        projected_feature_geoms.append(
            project_feature_geometry(feature, stereo_center=stereo_center)
        )

    for i, feature in gdf.iterrows():
        if feature.geometry is None:
            continue
        ray_list_i = []
        print(f"Collision Grid [First Collision Method]: {i}/{gdf.shape[0]}")
        geom_projected = project_feature_geometry(feature, stereo_center=stereo_center)
        minx, miny, maxx, maxy = geom_projected.bounds
        centroid_list.append(geom_projected.centroid.xy)
        # ax.plot([minx,minx,maxx,maxx], [miny,maxy,maxy,miny],)
        xx, yy = np.meshgrid(
            np.linspace(minx, maxx, num=density_used, endpoint=True)[1:-1],
            np.linspace(miny, maxy, num=density_used, endpoint=True)[1:-1],
        )
        if random_points:
            xx = minx + (maxx - minx) * np.random.rand(density_used, density_used)
            yy = miny + (maxy - miny) * np.random.rand(density_used, density_used)
        point_coords = np.array(
            list(zip(np.ndarray.flatten(xx), np.ndarray.flatten(yy)))
        )
        are_points_contained = [geom_projected.contains(Point(i)) for i in point_coords]
        for theta_i, theta in enumerate(theta_list):
            ray_list_theta = []
            rad = 100 * max(maxx - minx, maxy - miny) / 2
            thetaward_lines = [
                LineString([el, el + rad * np.array([np.cos(theta), np.sin(theta)])])
                for el in point_coords[are_points_contained]
            ]

            collision_counter = []
            for line in thetaward_lines:
                possible_rays = []
                distances = []
                for j, feature_j in enumerate(projected_feature_geoms):
                    try:
                        distances.append(
                            line.length
                            if i == j
                            else shapely.ops.split(line, feature_j).geoms[0].length
                        )
                        # print("shapely.ops.split(line, feature_j).geoms[0]",shapely.ops.split(line, feature_j).geoms[0],shapely.ops.split(line, feature_j).geoms[0].xy)
                        possible_rays.append(
                            shapely.ops.split(line, feature_j).geoms[0].xy
                        )
                    except:
                        distances.append(line.length)
                        possible_rays.append(line.xy)
                if np.max(distances) == np.min(distances):
                    # no collision
                    first_collision = -1
                    ray_list_theta.append([line.xy[0], line.xy[1], i])
                    collision_counter.append(i)
                else:
                    first_collision = np.argmin(distances)
                    ray_list_theta.append(
                        [
                            possible_rays[first_collision][0],
                            possible_rays[first_collision][1],
                            first_collision,
                        ]
                    )
                    collision_counter.append(first_collision)
            ray_list_i.append(ray_list_theta)
            if len(collision_counter):
                uniq_values, uniq_counts = np.unique(
                    collision_counter, return_counts=True
                )

                most_common_first_collision = uniq_values[np.argmax(uniq_counts)]

                # from Feature i in direction theta_i the first encountersed most of the time is most_common_first_collision
                for value_i in range(len(uniq_values)):
                    connection_grid[
                        connection_grid_i, uniq_values[value_i], theta_i
                    ] = uniq_counts[value_i]
        ray_list.append(ray_list_i)
        connection_grid_i += 1
    if return_rays:
        return connection_grid, np.array(centroid_list)[:, :, 0], ray_list
    return connection_grid, np.array(centroid_list)[:, :, 0], None


def project_feature_geometry(feature, stereo_center=None):
    # center is kansas
    if stereo_center is None:
        stereo_center = np.array([-98, 38]) * math.pi / 180
    if type(feature.geometry) is MultiPolygon:
        projected_polygons = []
        for polygon in feature.geometry.geoms:
            projected_polygon = Polygon(
                equi_to_stereo(
                    np.array(polygon.exterior.coords.xy).T * math.pi / 180,
                    stereo_center,
                ),
                [
                    equi_to_stereo(
                        np.array(interior.coords.xy).T * math.pi / 180, stereo_center
                    )
                    for interior in polygon.interiors
                ],
            )
            projected_polygons.append(projected_polygon)
        projected_multipolygon = MultiPolygon(projected_polygons)
        return projected_multipolygon
    elif type(feature.geometry) is Polygon:
        polygon = feature.geometry
        projected_polygon = Polygon(
            equi_to_stereo(
                np.array(polygon.exterior.coords.xy).T * math.pi / 180, stereo_center
            ),
            [
                equi_to_stereo(
                    np.array(interior.coords.xy).T * math.pi / 180, stereo_center
                )
                for interior in polygon.interiors
            ],
        )
        return projected_polygon


def RayPlurality(shapefile_path, plot=False, stereo_center=None):
    if stereo_center is None:
        stereo_center = np.array([0, 0]) * math.pi / 180
    gdf = gpd.read_file(shapefile_path)
    names = []
    for i, feature in gdf.iterrows():
        if feature.geometry is None:
            continue
        try:
            name_attribute = (
                feature["STATENAME"] + " " + ordinal(int(feature["DISTRICT"]))
                if int(feature["DISTRICT"]) > 0
                else feature["STATENAME"]
            )
        except:
            if "CODE" in feature:
                name_attribute = str(feature["CODE"])
            elif "NAME" in feature:
                name_attribute = str(feature["NAME"])
            elif "Name" in feature:
                name_attribute = str(feature["Name"])
            else:
                name_attribute = f"{i}"
        if name_attribute.startswith("District "):
            names.append(" ".join(name_attribute.split(" ")[1:]))
        else:
            names.append(name_attribute)

    pad_amount = max([len(el) for el in names])

    connection_grid, centroid_arr, ray_list = connection_grid_first_collision_method(
        gdf,
        density_used=10,
        return_rays=plot,
        stereo_center=stereo_center,
        random_points=False,
    )

    direction_table = np.zeros((connection_grid.shape[0], 4))
    directions_in_conn_grid_order = ["r", "u", "l", "d"]
    directions_in_wanted_order = ["u", "d", "l", "r"]
    # directions_in_wanted_order = ['r','u','l','d']
    direction_dict = {}
    for el in directions_in_wanted_order:
        direction_dict[el] = directions_in_conn_grid_order.index(el)

    direction_strings = [
        " " * pad_amount + "  " + (" " * (pad_amount)).join(directions_in_wanted_order)
    ]
    for i in range(connection_grid.shape[0]):
        line_string = f"{' '*(pad_amount-len(names[i])) + names[i]}:"
        for theta_i_i in range(4):
            theta_i = direction_dict[directions_in_wanted_order[theta_i_i]]
            if (
                connection_grid[i, i, theta_i]
                < np.sum(connection_grid[i, :, theta_i]) / 2
            ):
                connection_grid[i, i, theta_i] = 0
                direction_table[i, theta_i] = np.argmax(connection_grid[i, :, theta_i])
                line_string = f"{line_string} {' '*(pad_amount-len(names[int(direction_table[i,theta_i])])) + names[int(direction_table[i,theta_i])]}"
            else:
                direction_table[i, theta_i] = -1
                line_string = f"{line_string} {'-'*pad_amount}"
        direction_strings.append(line_string)
    direction_strings[1:] = sorted(direction_strings[1:])
    print(f"RayPlurality [{shapefile_path}]")
    print("\n".join(direction_strings))
    tab_color_list = list(mcolors.TABLEAU_COLORS.keys())
    if plot:
        ax = plt.gca()
        for i, feature in gdf.iterrows():
            if feature.geometry is None:
                continue
            feature.geometry = project_feature_geometry(
                feature, stereo_center=stereo_center
            )
            feature_gdf = gpd.GeoDataFrame([feature], crs=gdf.crs)
            as_hsv = mcolors.rgb_to_hsv(
                mcolors.to_rgb(tab_color_list[i % len(tab_color_list)])
            )
            whitening = 4
            as_hsv[2] = (1 * (whitening - 1) + as_hsv[2]) / whitening
            as_hsv[1] = (0 + as_hsv[1]) / whitening
            feature_gdf.plot(ax=ax, color=mcolors.hsv_to_rgb(as_hsv), edgecolor="black")
        for theta_i in range(4):
            for ray in ray_list[0][theta_i]:
                plt.plot(
                    [ray[0][0], ray[0][1]],
                    [ray[1][0], ray[1][1]],
                    c=tab_color_list[ray[2] % len(tab_color_list)],
                    alpha=(1.0 if theta_i == 0 else 0.3),
                    linewidth=(2.0 if theta_i == 0 else 0.5),
                )
        plt.xlim((-5 * math.pi / 180, 5 * math.pi / 180))
        plt.ylim((-5 * math.pi / 180, 5 * math.pi / 180))
        plt.show()
    return np.array(direction_table).astype(int), connection_grid


def get_min_presses(direction_table, return_paths=False):
    optimal_path = [
        [None for i in range(direction_table.shape[0])]
        for j in range(direction_table.shape[0])
    ]
    min_presses = np.zeros((direction_table.shape[0], direction_table.shape[0])) + 1000
    for i in range(direction_table.shape[0]):
        for theta_i in range(4):
            if direction_table[i, theta_i] >= 0:
                if 1 < min_presses[i, direction_table[i, theta_i]]:
                    optimal_path[i][direction_table[i, theta_i]] = [theta_i]
                min_presses[i, direction_table[i, theta_i]] = min(
                    min_presses[i, direction_table[i, theta_i]], 1
                )
        min_presses[i, i] = 0
        optimal_path[i][i] = []
    total = np.sum(min_presses)
    finished = False
    while not finished:
        for j in range(direction_table.shape[0]):
            for i in range(direction_table.shape[0]):
                for theta_i in range(4):
                    # j to i
                    if direction_table[i, theta_i] >= 0:
                        if min_presses[j, i] >= 0:
                            if (
                                min_presses[j, i] + 1
                                < min_presses[j, direction_table[i, theta_i]]
                            ):
                                optimal_path[j][
                                    direction_table[i, theta_i]
                                ] = copy.copy(optimal_path[j][i])
                                optimal_path[j][direction_table[i, theta_i]].append(
                                    theta_i
                                )
                            min_presses[j, direction_table[i, theta_i]] = min(
                                min_presses[j, direction_table[i, theta_i]],
                                min_presses[j, i] + 1,
                            )
        if np.sum(min_presses) == total:
            finished = True
        else:
            total = np.sum(min_presses)

    if return_paths:
        return min_presses, optimal_path
    return min_presses


def loss(direction_table, connection_grid):
    total = 0
    for i in range(direction_table.shape[0]):
        for theta_i in range(4):
            total += connection_grid[i, direction_table[i, theta_i], theta_i]
    return total


def ImprovedRayPlurality(shapefile_path, plot=False, stereo_center=None):
    direction_table, connection_grid = RayPlurality(
        shapefile_path, plot=plot, stereo_center=stereo_center
    )  # ,stereo_center=np.array([-98, 38]) * math.pi / 180)

    # potential modifications
    finished = False
    while not finished:
        best_improvement = None
        best_improvement_score = None
        for i in range(direction_table.shape[0]):
            for theta_i in range(4):
                for j in range(direction_table.shape[0]):
                    if connection_grid[i, j, theta_i] > 0:
                        new_direction_table = np.array(direction_table)
                        new_direction_table[i, theta_i] = j
                        if np.max(get_min_presses(new_direction_table)) < np.max(
                            get_min_presses(direction_table)
                        ):
                            if best_improvement is None:
                                best_improvement = [i, theta_i, j]
                                best_improvement_score = [
                                    np.max(get_min_presses(new_direction_table)),
                                    loss(new_direction_table, connection_grid),
                                ]
                            elif (
                                loss(new_direction_table, connection_grid)
                                > best_improvement_score[1]
                            ):
                                best_improvement = [i, theta_i, j]
                                best_improvement_score = [
                                    np.max(get_min_presses(new_direction_table)),
                                    loss(new_direction_table, connection_grid),
                                ]

        if best_improvement is not None:
            direction_table[
                best_improvement[0], best_improvement[1]
            ] = best_improvement[2]
        else:
            finished = True

    return np.array(direction_table).astype(int), connection_grid



def get_closest_in_direction(file_path, secondary_movement=False, remove_loops=False):
    # Load shapefile using geopandas
    gdf = gpd.read_file(file_path)

    print("oss", os.path.split(file_path)[-1].split(".")[0])
    file_identifier = os.path.split(file_path)[-1].split(".")[0]

    if not (
        os.path.exists(f"data/connection_grid_{file_identifier}.npy")
        and os.path.exists(f"data/centroid_arr_{file_identifier}.npy")
    ):
        connection_grid, centroid_arr, _ = connection_grid_first_collision_method(gdf)
        np.save(f"data/connection_grid_{file_identifier}.npy", connection_grid)
        np.save(f"data/centroid_arr_{file_identifier}.npy", centroid_arr)

    connection_grid = np.load(f"data/connection_grid_{file_identifier}.npy")
    centroid_arr = np.load(f"data/centroid_arr_{file_identifier}.npy")
    # if it's overwhelmingly obvious where a direction should go,
    # like Nigeria goes east to Cameroon,
    # and the reverse of an alternative, like Chad to Nigeria, isn't especially strong,
    # eliminate the alternative (Nigeria to Chad)
    for district_i in range(connection_grid.shape[0]):
        for district_j in range(connection_grid.shape[0]):
            for theta_i in range(4):
                if (
                    np.max(connection_grid[district_i, :, theta_i])
                    / np.sum(connection_grid[district_i, :, theta_i])
                    > 0.75
                ):
                    if np.argmax(connection_grid[district_i, :, theta_i]) != district_j:
                        if (
                            connection_grid[district_j, district_i, (theta_i + 2) % 4]
                            < np.max(connection_grid[district_j, :, (theta_i + 2) % 4])
                            / 2
                        ):
                            connection_grid[district_i, district_j, theta_i] = min(
                                1, connection_grid[district_i, district_j, theta_i]
                            )

    print("max", np.max(connection_grid))

    stereo_center_of_lines = np.array([0.05458667, -0.06561405])
    centroid_arr = equi_to_stereo(centroid_arr * math.pi / 180, stereo_center_of_lines)

    theta_signs = [1, 1, -1, -1]
    # print("cent", centroid_arr[17])
    # print(centroid_arr[29])
    # gsn
    # plt.plot(connection_grid[22,34])
    # plt.plot(connection_grid[18,:,0])
    # plt.plot(connection_grid[27,:,1])
    # plt.plot(connection_grid[2,:,0])
    # plt.show()
    # plt.plot(connection_grid[2,18,:])
    # plt.plot(connection_grid[2,27,:])
    # plt.show()
    # plt.plot(connection_grid[53,24,:])
    # plt.plot(connection_grid[53,43,:])
    # plt.show()
    # plt.plot(np.sum(np.sum(connection_grid,axis=1),axis=1))
    # plt.show()
    # print("centroid_arr.shape")
    # print(centroid_arr.shape)
    # where_filled = np.where(np.sum(np.sum(connection_grid,axis=1),axis=1)>0)[0]
    # connection_grid = connection_grid[where_filled]
    # connection_grid = connection_grid[:,where_filled]
    # np.save(f"data/connection_grid_{file_identifier}.npy", connection_grid)
    # centroid_arr = centroid_arr[where_filled]
    for district_i in range(connection_grid.shape[0]):
        for district_j in range(connection_grid.shape[0]):
            for theta_i in range(4):
                if (
                    (
                        centroid_arr[district_i, theta_i % 2] * theta_signs[theta_i]
                        > centroid_arr[district_j, theta_i % 2] * theta_signs[theta_i]
                    )
                    and connection_grid[district_i, district_j, theta_i] < 1600 / 5
                    and connection_grid[district_i, district_j, theta_i]
                    < np.max(connection_grid[district_i, :, :]) / 2
                ):
                    connection_grid[district_i, district_j, theta_i] = 0

                # if the angle is more than 45° from ideal
                major_offset = abs(
                    centroid_arr[district_i, theta_i % 2]
                    - centroid_arr[district_j, theta_i % 2]
                )
                minor_offset = abs(
                    centroid_arr[district_i, (theta_i + 1) % 2]
                    - centroid_arr[district_j, (theta_i + 1) % 2]
                )
                if (
                    minor_offset > major_offset
                    and connection_grid[district_i, district_j, theta_i] < 1600 / 5
                    and connection_grid[district_i, district_j, theta_i]
                    < np.max(connection_grid[district_i, :, :]) / 2
                ):
                    connection_grid[district_i, district_j, theta_i] = 0

    for district_i in range(connection_grid.shape[0]):
        connection_grid[district_i, district_i, :] = 0

    # it's necessary that if A goes to B in opposite directions that the lesser one be zeroed
    # excluding stuff where it's a massive amount, like The Gambia going south
    connection_grid[
        np.logical_and(
            connection_grid - np.roll(connection_grid, 2, axis=2) < 0,
            connection_grid < 1600 / 3,
        )
    ] = 0

    # needs to be at least 20% of the connections to count
    # connection_grid[connection_grid<0.2*np.sum(connection_grid,axis=1)[:,None,:]] = 0

    # so large A to tiny B still has at least something for A to B
    connection_grid[
        np.logical_and(
            np.logical_and(
                connection_grid == 0,
                np.transpose(np.roll(connection_grid, 2, axis=2), (1, 0, 2))
                > np.max(
                    np.transpose(np.roll(connection_grid, 2, axis=2), (1, 0, 2)), axis=0
                )[None, :, :]
                / 2,
            ),
            np.transpose(np.roll(connection_grid, 2, axis=2), (1, 0, 2))
            >= np.max(
                np.transpose(np.roll(connection_grid, 2, axis=2), (1, 0, 2)), axis=2
            )[:, :, None],
        )
    ] = 1  # 0.5

    if remove_loops:
        for district_i in range(connection_grid.shape[0]):
            connection_grid[district_i, district_i, :] = 0

        # it's necessary that if A goes to B in opposite directions that the lesser one be zeroed
        connection_grid[connection_grid - np.roll(connection_grid, 2, axis=2) < 0] = 0

        # REMOVE LOOPS
        finished_finding_loops = False
        while not finished_finding_loops:
            loops_found = 0
            for theta_i in range(4):
                for district_i in range(connection_grid.shape[0]):
                    # district_i is the root, with height 0, all else rises from that

                    quickest_path_to_each = [
                        None for i in range(connection_grid.shape[0])
                    ]
                    quickest_path_to_each[district_i] = [district_i]
                    finished = False
                    while not finished:
                        num_added = 0
                        for from_i in range(connection_grid.shape[0]):
                            # print(from_i,quickest_path_to_each[from_i])
                            if quickest_path_to_each[from_i] is not None:
                                next_list = np.where(
                                    connection_grid[from_i, :, theta_i] > 0
                                )[0]
                                for next_i in next_list:
                                    if next_i == district_i:
                                        # LOOP FOUND!!!
                                        loops_found += 1
                                        loop = quickest_path_to_each[from_i]
                                        if 15 in loop:
                                            print(
                                                "Loop fuond!!!",
                                                quickest_path_to_each[from_i],
                                            )
                                        loop_connection_strengths = []
                                        for loop_i in range(
                                            len(quickest_path_to_each[from_i])
                                        ):
                                            loop_connection_strengths.append(
                                                connection_grid[
                                                    loop[loop_i],
                                                    loop[(loop_i + 1) % len(loop)],
                                                    theta_i,
                                                ]
                                            )

                                        weakest_i = np.argmin(loop_connection_strengths)
                                        if 15 in loop:
                                            print(
                                                loop_connection_strengths,
                                                loop[weakest_i],
                                                loop[(weakest_i + 1) % len(loop)],
                                                theta_i,
                                            )
                                        # break the weakest link in the loop
                                        connection_grid[
                                            loop[weakest_i],
                                            loop[(weakest_i + 1) % len(loop)],
                                            theta_i,
                                        ] = 0
                                        finished = True
                                        break
                                    if quickest_path_to_each[next_i] is None:
                                        quickest_path_to_each[next_i] = list(
                                            quickest_path_to_each[from_i]
                                        )
                                        quickest_path_to_each[next_i].append(next_i)
                                        num_added += 1
                            if finished:
                                break
                        if num_added == 0:
                            finished = True
                    # print(quickest_path_to_each)
            if loops_found == 0:
                finished_finding_loops = True

            # heights = np.zeros((connection_grid.shape[0]))-1
            # heights[district_i] = 0
            # current_height = 0
            # new_heights = np.array(heights)
            # for next_i in range(connection_grid.shape[0]):
            #     where_branches_to_next_i = np.logical_and(connection_grid[:, next_i, theta_i]>0,heights==current_height)
            #     if np.sum(where_branches_to_next_i) > 0:
            #         if next_i == district_i:
            #             # LOOP FOUND!!!
            #
            #         new_heights[next_i] = current_height+1

    closest_in_direction = np.zeros((connection_grid.shape[0], 4)).astype(int) - 1
    closest_in_direction[np.max(connection_grid, axis=1) > 0] = np.argmax(
        connection_grid, axis=1
    )[np.max(connection_grid, axis=1) > 0].astype(int)

    directly_adjacent_by_starting_from = []
    for starting_from in range(connection_grid.shape[0]):
        directly_adjacent = []
        for district_i in range(connection_grid.shape[0]):
            if np.max(connection_grid[starting_from, district_i, :]) > 0:
                directly_adjacent.append(district_i)
        directly_adjacent_by_starting_from.append(directly_adjacent)

    def evaluate_solution(solution_arr, starting_from, ignore_doubles=False):
        directly_adjacent = directly_adjacent_by_starting_from[starting_from]
        solution_arr = np.array(solution_arr)
        # solution_arr contains 20 values, which are the immediate next ones in directions 0,1,2,3
        # then the 16 secondaries

        adjacents_included = 0
        for i in range(len(directly_adjacent)):
            if directly_adjacent[i] in solution_arr:
                adjacents_included += 1
        if not ignore_doubles:
            double_adjacents = (
                []
            )  # list(set(list(itertools.chain.from_iterable([directly_adjacent_by_starting_from[solution_arr[i]] for i in range(4)]))).difference(set(directly_adjacent)))
            # print(double_adjacents)
            for i in range(4):
                for el in directly_adjacent_by_starting_from[solution_arr[i]]:
                    if el not in double_adjacents and el not in directly_adjacent:
                        double_adjacents.append(el)
                        adjacents_included += 0.5
            # for i in range(len(double_adjacents)):
            #     if double_adjacents[i] in solution_arr:
            #         adjacents_included += 0.5
        bonus = 0.4 * np.sum(
            [
                [
                    connection_grid[starting_from, solution_arr[4 + 4 * i + j], i]
                    * (1 if (j - i) % 2 == 1 else 0)
                    for j in range(4)
                ]
                for i in range(4)
            ]
        )

        # so going back should be where you were
        for i in range(4):
            if solution_arr[4 + 4 * i + (i + 2) % 4] == starting_from:
                bonus += 500

        if not ignore_doubles:
            for i in range(4):
                if solution_arr[i] == starting_from:
                    if np.any(
                        solution_arr[:4] != solution_arr[4 + 4 * i : 4 + 4 * i + 4]
                    ):
                        bonus -= 100000

        # to avoid down-left etc going back to the same thing
        for i in [5, 7, 8, 10, 13, 15, 16, 18]:
            if (
                solution_arr[i] == starting_from
                and solution_arr[int((i - 4) / 4)] != starting_from
            ):
                bonus -= 100000

        return bonus + (
            100000 * adjacents_included
            + np.sum(
                [
                    connection_grid[starting_from, solution_arr[i], i]
                    - (
                        1000
                        if connection_grid[starting_from, solution_arr[i], i] == 0
                        else 0
                    )
                    for i in range(4)
                ]
            )
            + np.sum(
                [
                    [
                        connection_grid[solution_arr[i], solution_arr[4 + 4 * i + j], j]
                        - (
                            1000
                            if connection_grid[
                                solution_arr[i], solution_arr[4 + 4 * i + j], j
                            ]
                            == 0
                            else 0
                        )
                        for i in range(4)
                    ]
                    for j in range(4)
                ]
            )
        )

    def find_best_for_possible_start(
        starting_from, possible_start, possible_value_dicts
    ):
        default_solution = possible_start[:]

        for theta_i in range(4):
            for theta_j in range(4):
                default_solution.append(
                    np.argmax(connection_grid[default_solution[theta_i], :, theta_j])
                )
        best_switches = np.zeros((20, 20))
        best_switch_val = -1000000000000000000
        best_switch_solution = None
        for i in range(16):
            for j in range(16):
                for k in range(16):
                    if i == j or i == k or k == j:
                        continue
                    best_switch = -1000000000000000000
                    for p_i in range(
                        len(possible_value_dicts[4 + i][default_solution[int(i / 4)]])
                    ):
                        for p_j in range(
                            len(
                                possible_value_dicts[4 + j][
                                    default_solution[int(j / 4)]
                                ]
                            )
                        ):
                            for p_k in range(
                                len(
                                    possible_value_dicts[4 + k][
                                        default_solution[int(k / 4)]
                                    ]
                                )
                            ):
                                new_solution = default_solution[:]
                                new_solution[4 + i] = possible_value_dicts[4 + i][
                                    default_solution[int(i / 4)]
                                ][p_i]
                                new_solution[4 + j] = possible_value_dicts[4 + j][
                                    default_solution[int(j / 4)]
                                ][p_j]
                                new_solution[4 + k] = possible_value_dicts[4 + k][
                                    default_solution[int(k / 4)]
                                ][p_k]
                                for ii in range(4):
                                    if new_solution[ii] == starting_from:
                                        new_solution[4 + 4 * ii : 4 + 4 * ii + 4] = (
                                            new_solution[:4]
                                        )
                                # to avoid down-left etc going back to the same thing
                                avoid_this_one = 0
                                for ii in [5, 7, 8, 10, 13, 15, 16, 18]:
                                    if (
                                        new_solution[ii] == starting_from
                                        and new_solution[int((ii - 4) / 4)]
                                        != starting_from
                                    ):
                                        # this can either replace the leaf, or the whole branch
                                        branch_strength = connection_grid[
                                            starting_from,
                                            new_solution[ii],
                                            int((ii - 4) / 4),
                                        ]
                                        leaf_strength = connection_grid[
                                            new_solution[int((ii - 4) / 4)],
                                            new_solution[ii],
                                            ii % 4,
                                        ]
                                        if branch_strength > leaf_strength:
                                            potential_alternative_scores = np.array(
                                                connection_grid[
                                                    new_solution[int((ii - 4) / 4)],
                                                    :,
                                                    ii % 4,
                                                ]
                                            )
                                            potential_alternative_scores[
                                                new_solution[ii]
                                            ] = 0
                                            if (
                                                np.max(potential_alternative_scores)
                                                == 0
                                            ):
                                                if (
                                                    leaf_strength
                                                    < np.sqrt(1600 / 5) + 10
                                                ):
                                                    new_solution[ii] = new_solution[
                                                        int((ii - 4) / 4)
                                                    ]
                                            else:
                                                new_solution[ii] = np.argmax(
                                                    potential_alternative_scores
                                                )
                                        else:
                                            avoid_this_one = -100000
                                            # potential_alternative_scores = np.array(connection_grid[starting_from,:,int((ii-4) / 4)])
                                            # potential_alternative_scores[new_solution[int((ii-4) / 4)]] = 0
                                            # if np.max(potential_alternative_scores)==0:
                                            #     if branch_strength<np.sqrt(1600/5)+10:
                                            #         new_solution[int((ii-4) / 4)] = starting_from
                                            #         new_solution[4+4*int((ii - 4) / 4):4+4*int((ii - 4) / 4)+4] = new_solution[:4]
                                            # else:
                                            #     new_solution[int((ii-4) / 4)] = np.argmax(potential_alternative_scores)
                                val = avoid_this_one + evaluate_solution(
                                    new_solution,
                                    starting_from,
                                )
                                best_switch = max(best_switch, val)
                                # print("val",val,best_switch_val,val > best_switch_val)
                                if val > best_switch_val:
                                    best_switch_val = val
                                    best_switch_solution = new_solution[:]

                    best_switches[i, j] = best_switch
        for t in range(2):
            default_solution = best_switch_solution[:]
            best_switches = np.zeros((20, 20))
            best_switch_val = -1000000000000000000
            best_switch_solution = None
            for i in range(16):
                for j in range(16):
                    for k in range(16):
                        if i == j or i == k or k == j:
                            continue
                        best_switch = -1000000000000000000
                        # print(possible_value_dicts[4 + i],default_solution,default_solution[int(i/4)],i)
                        for p_i in range(
                            len(
                                possible_value_dicts[4 + i][
                                    default_solution[int(i / 4)]
                                ]
                            )
                        ):
                            for p_j in range(
                                len(
                                    possible_value_dicts[4 + j][
                                        default_solution[int(j / 4)]
                                    ]
                                )
                            ):
                                for p_k in range(
                                    len(
                                        possible_value_dicts[4 + k][
                                            default_solution[int(k / 4)]
                                        ]
                                    )
                                ):
                                    new_solution = default_solution[:]
                                    new_solution[4 + i] = possible_value_dicts[4 + i][
                                        default_solution[int(i / 4)]
                                    ][p_i]
                                    new_solution[4 + j] = possible_value_dicts[4 + j][
                                        default_solution[int(j / 4)]
                                    ][p_j]
                                    new_solution[4 + k] = possible_value_dicts[4 + k][
                                        default_solution[int(k / 4)]
                                    ][p_k]
                                    for ii in range(4):
                                        if new_solution[ii] == starting_from:
                                            new_solution[
                                                4 + 4 * ii : 4 + 4 * ii + 4
                                            ] = new_solution[:4]
                                    # to avoid down-left etc going back to the same thing
                                    avoid_this_one = 0
                                    for ii in [5, 7, 8, 10, 13, 15, 16, 18]:
                                        if (
                                            new_solution[ii] == starting_from
                                            and new_solution[int((ii - 4) / 4)]
                                            != starting_from
                                        ):
                                            # this can either replace the leaf, or the whole branch
                                            branch_strength = connection_grid[
                                                starting_from,
                                                new_solution[ii],
                                                int((ii - 4) / 4),
                                            ]
                                            leaf_strength = connection_grid[
                                                new_solution[int((ii - 4) / 4)],
                                                new_solution[ii],
                                                ii % 4,
                                            ]
                                            if branch_strength > leaf_strength:
                                                potential_alternative_scores = np.array(
                                                    connection_grid[
                                                        new_solution[int((ii - 4) / 4)],
                                                        :,
                                                        ii % 4,
                                                    ]
                                                )
                                                potential_alternative_scores[
                                                    new_solution[ii]
                                                ] = 0
                                                if (
                                                    np.max(potential_alternative_scores)
                                                    == 0
                                                ):
                                                    if (
                                                        leaf_strength
                                                        < np.sqrt(1600 / 5) + 10
                                                    ):
                                                        new_solution[ii] = new_solution[
                                                            int((ii - 4) / 4)
                                                        ]
                                                else:
                                                    new_solution[ii] = np.argmax(
                                                        potential_alternative_scores
                                                    )
                                            else:
                                                avoid_this_one = -100000
                                                # potential_alternative_scores = np.array(connection_grid[starting_from,:,int((ii-4) / 4)])
                                                # potential_alternative_scores[new_solution[int((ii-4) / 4)]] = 0
                                                # if np.max(potential_alternative_scores)==0:
                                                #     if branch_strength<np.sqrt(1600/5)+10:
                                                #         new_solution[int((ii-4) / 4)] = starting_from
                                                #         new_solution[4+4*int((ii - 4) / 4):4+4*int((ii - 4) / 4)+4] = new_solution[:4]
                                                # else:
                                                #     new_solution[int((ii-4) / 4)] = np.argmax(potential_alternative_scores)

                                    val = avoid_this_one + evaluate_solution(
                                        new_solution,
                                        starting_from,
                                    )
                                    best_switch = max(best_switch, val)
                                    if val > best_switch_val:
                                        best_switch_val = val
                                        best_switch_solution = new_solution[:]

                        best_switches[i, j] = best_switch
        best_switch_solution = [int(el) for el in best_switch_solution]
        return best_switch_solution

    # connection_grid = np.sqrt(connection_grid)
    # connection_grid[connection_grid > 0] += 10

    # so you don't get stuck in the Seychelles
    for district_i in range(connection_grid.shape[0]):
        places_from_there = connection_grid[district_i, :, :]
        places_from_there[district_i, :] = 0
        if np.max(connection_grid[district_i, :, :]) == 0:
            places_that_get_there = np.where(connection_grid[:, district_i, :] > 0)
            print("places_that_get_there", places_that_get_there)
            for i in range(len(places_that_get_there[0])):
                if places_that_get_there[0][i] != district_i:
                    connection_grid[
                        district_i,
                        places_that_get_there[0][i],
                        (places_that_get_there[1][i] + 2) % 4,
                    ] = 1

    return (
        closest_in_direction,
        [[[el_el] for el_el in el] for el in closest_in_direction],
        connection_grid,
        str({}),
        centroid_arr,
    )

    possible_value_dicts_by_starting_from = []
    if not os.path.exists(f"data/all_data_{file_identifier}.pkl"):
        all_data = {}
        for starting_from in range(connection_grid.shape[0]):
            directly_adjacent = []
            for district_i in range(connection_grid.shape[0]):
                if np.max(connection_grid[starting_from, district_i, :]) > 0:
                    directly_adjacent.append(district_i)
            print("Starting from ", starting_from)
            all_data[starting_from] = {"best_possible": None, "by_starts": {}}
            directly_adjacent = []
            for district_i in range(connection_grid.shape[0]):
                if np.max(connection_grid[starting_from, district_i, :]) > 0:
                    directly_adjacent.append(district_i)
            # so if there's nowhere else to go, the district connects to itself in that direction
            for district_i in range(connection_grid.shape[0]):
                for theta_i in range(4):
                    if np.max(connection_grid[district_i, :, theta_i]) == 0:
                        connection_grid[district_i, district_i, theta_i] = 1

            possible_values = [[] for i in range(20)]
            possible_value_dicts = [{} for i in range(20)]
            for theta_i in range(4):
                for j in range(connection_grid.shape[1]):
                    if connection_grid[starting_from, j, theta_i] > 0:
                        possible_values[theta_i].append(j)
                        possible_value_dicts[theta_i][j] = [j]
                        for theta_j in range(4):
                            possible_value_dicts[4 + 4 * theta_i + theta_j][j] = []
            for theta_i in range(4):
                for theta_j in range(4):
                    for possible_i in possible_values[theta_i]:
                        for j in range(connection_grid.shape[1]):
                            if connection_grid[possible_i, j, theta_j] > 0:
                                if j not in possible_values[4 + 4 * theta_i + theta_j]:
                                    possible_values[4 + 4 * theta_i + theta_j].append(j)
                                    possible_value_dicts[4 + 4 * theta_i + theta_j][
                                        possible_i
                                    ].append(j)
            possible_value_dicts_by_starting_from.append(possible_value_dicts)

            coord_4d = [0, 0, 0, 0]
            possible_starts = []
            val_lens = [len(el) for el in possible_values[:4]]
            for i in range(np.prod(val_lens)):
                possible_starts.append(
                    [possible_values[el_i][coord_4d[el_i]] for el_i in range(4)]
                )
                coord_4d[-1] += 1
                for v_i in range(4):
                    if coord_4d[3 - v_i] == val_lens[3 - v_i]:
                        coord_4d[3 - v_i] = 0
                        coord_4d[3 - v_i - 1] += 1
            print("Possible starts:", len(possible_starts))
            best_possible_start_val = 0
            best_possible_start = None
            possible_start_val_list = []
            for possible_start in possible_starts:

                best_switch_solution = find_best_for_possible_start(
                    starting_from, possible_start, possible_value_dicts
                )
                all_data[starting_from]["by_starts"][
                    str(possible_start)
                    .replace("[", "")
                    .replace("]", "")
                    .replace(" ", "")
                ] = best_switch_solution[:]
                possible_start_val_list.append(
                    evaluate_solution(best_switch_solution, starting_from)
                )
                if (
                    evaluate_solution(best_switch_solution, starting_from)
                    > best_possible_start_val
                ):
                    best_possible_start_val = evaluate_solution(
                        best_switch_solution, starting_from
                    )
                    best_possible_start = best_switch_solution[:]
            all_data[starting_from]["best_possible"] = best_possible_start

        pickle.dump(all_data, open(f"data/all_data_{file_identifier}.pkl", "wb"))

        pickle.dump(all_data, open(f"data/all_data_{file_identifier}.pkl", "wb"))
    all_data = pickle.load(open(f"data/all_data_{file_identifier}.pkl", "rb"))

    # in case the change to prevent down-left going to same causes a start that doesn't exist
    for starting_from in all_data:
        for theta_i in range(4):
            sft_i = all_data[starting_from]["best_possible"][theta_i]
            sft_i_start = all_data[starting_from]["best_possible"][
                4 + 4 * theta_i : 4 + 4 * theta_i + 4
            ]
            sft_i_start_str = (
                str(sft_i_start).replace("[", "").replace("]", "").replace(" ", "")
            )
            if sft_i_start_str not in all_data[sft_i]["by_starts"]:
                new_potential = sft_i_start[:]
                for theta_j in range(4):
                    new_potential = (
                        new_potential
                        + all_data[sft_i_start[theta_j]]["best_possible"][:4]
                    )
                all_data[sft_i]["by_starts"][sft_i_start_str] = new_potential[:]
                # possible_values = [[] for i in range(20)]
                # possible_value_dicts = [{} for i in range(20)]
                # for theta_i in range(4):
                #     for j in range(connection_grid.shape[1]):
                #         if connection_grid[starting_from, j, theta_i] > 0 or j == sft_i_start[theta_i]:
                #             possible_values[theta_i].append(j)
                #             possible_value_dicts[theta_i][j] = [j]
                #             for theta_j in range(4):
                #                 possible_value_dicts[4 + 4 * theta_i + theta_j][j] = []
                # # print(possible_values)
                # # print(possible_value_dicts)
                # for theta_i in range(4):
                #     for theta_j in range(4):
                #         for possible_i in possible_values[theta_i]:
                #             for j in range(connection_grid.shape[1]):
                #                 if connection_grid[possible_i, j, theta_j] > 0:
                #                     if j not in possible_values[4 + 4 * theta_i + theta_j]:
                #                         possible_values[4 + 4 * theta_i + theta_j].append(j)
                #                         possible_value_dicts[4 + 4 * theta_i + theta_j][
                #                             possible_i
                #                         ].append(j)
                #
                # all_data[sft_i]["by_starts"][sft_i_start_str] = find_best_for_possible_start(starting_from, sft_i_start,
                #                                                             possible_value_dicts)
            for start in all_data[starting_from]["by_starts"]:

                sft_i = all_data[starting_from]["by_starts"][start][theta_i]
                sft_i_start = all_data[starting_from]["by_starts"][start][
                    4 + 4 * theta_i : 4 + 4 * theta_i + 4
                ]
                sft_i_start_str = (
                    str(sft_i_start).replace("[", "").replace("]", "").replace(" ", "")
                )
                if sft_i_start_str not in all_data[sft_i]["by_starts"]:
                    new_potential = sft_i_start[:]
                    for theta_j in range(4):
                        new_potential = (
                            new_potential
                            + all_data[sft_i_start[theta_j]]["best_possible"][:4]
                        )
                    all_data[sft_i]["by_starts"][sft_i_start_str] = new_potential[:]
    for starting_from in all_data:
        directly_adjacent = []
        for district_i in range(connection_grid.shape[0]):
            if np.max(connection_grid[starting_from, district_i, :]) > 0:
                directly_adjacent.append(district_i)
        best_option = None
        best_option_score = 0
        for option in all_data[starting_from]["by_starts"]:
            option_score = np.sum(
                [
                    connection_grid[
                        starting_from,
                        all_data[starting_from]["by_starts"][option][i],
                        i,
                    ]
                    for i in range(4)
                ]
            ) + 100000 * int(
                evaluate_solution(
                    all_data[starting_from]["by_starts"][option],
                    starting_from,
                    ignore_doubles=True,
                )
                / 100000
            )
            if option_score > best_option_score:
                best_option_score = option_score
                best_option = all_data[starting_from]["by_starts"][option][:]
        #     if starting_from==9:
        #         print(",",option,option_score)
        # if starting_from==9:
        #     print(best_option)
        if best_option is not None:
            all_data[starting_from]["best_possible"] = best_option

    for starting_from in all_data:
        ideal = all_data[starting_from]["best_possible"][:4]
        for sfj in all_data:
            for i in range(4):
                if all_data[sfj]["best_possible"][i] == starting_from:
                    # compare all_data[sfj]["best_possible"][4+4*i:4+4*i+4] with ideal
                    to_compare_initial = np.array(
                        all_data[sfj]["best_possible"][4 + 4 * i : 4 + 4 * i + 4]
                    )
                    to_compare = np.array(to_compare_initial)
                    if np.sum(to_compare != np.array(ideal)) > 0:
                        to_compare[to_compare != np.array(ideal)] = -1
                        if (
                            len(np.unique(to_compare_initial))
                            == len(np.unique(to_compare)) - 1
                        ):
                            all_data[sfj]["best_possible"][
                                4 + 4 * i : 4 + 4 * i + 4
                            ] = ideal
            for start in all_data[sfj]["by_starts"]:
                for i in range(4):
                    if all_data[sfj]["by_starts"][start][i] == starting_from:
                        # compare all_data[sfj]["by_starts"][start][4+4*i:4+4*i+4] with ideal
                        to_compare_initial = np.array(
                            all_data[sfj]["by_starts"][start][4 + 4 * i : 4 + 4 * i + 4]
                        )
                        to_compare = np.array(to_compare_initial)
                        if np.sum(to_compare != np.array(ideal)) > 0:
                            to_compare[to_compare != np.array(ideal)] = -1
                            if (
                                len(np.unique(to_compare_initial))
                                == len(np.unique(to_compare)) - 1
                            ):
                                all_data[sfj]["by_starts"][start][
                                    4 + 4 * i : 4 + 4 * i + 4
                                ] = ideal

    if secondary_movement:
        all_in_direction = [
            [
                np.where(connection_grid[district_i, :, theta_i] > 0)[0]
                for theta_i in range(4)
            ]
            for district_i in range(connection_grid.shape[0])
        ]
        return (
            closest_in_direction,
            all_in_direction,
            connection_grid,
            str(all_data),
            centroid_arr,
        )

    return (
        closest_in_direction,
        [[[el_el] for el_el in el] for el in closest_in_direction],
        connection_grid,
        str(all_data),
        centroid_arr,
    )



if __name__ == "__main__":
    directions_in_conn_grid_order = ["r", "u", "l", "d"]
    directions_in_wanted_order = ["u", "d", "l", "r"]
    direction_dict = {}
    for el in directions_in_wanted_order:
        direction_dict[el] = directions_in_conn_grid_order.index(el)

    shapefile_path = "shapefiles/quadria_districts_with_geom_info.shp"
    shapefile_path = "shapefiles/quintia_districts_with_geom_info.shp"

    shapefile_path = "shapefiles/us_states_48.shp"

    # Load shapefile using geopandas
    gdf = gpd.read_file(shapefile_path)
    names = []
    for i, feature in gdf.iterrows():
        if feature.geometry is None:
            continue
        try:
            name_attribute = (
                feature["STATENAME"] + " " + ordinal(int(feature["DISTRICT"]))
                if int(feature["DISTRICT"]) > 0
                else feature["STATENAME"]
            )
        except:
            if "NAME" in feature:
                name_attribute = feature["NAME"]
            else:
                name_attribute = f"{i}"
        names.append(name_attribute[-1])
    pad_amount = max([len(el) for el in names])

    direction_table, connection_grid = ImprovedRayPlurality(
        shapefile_path, plot=False
    )  # ,stereo_center=np.array([-98, 38]) * math.pi / 180)

    min_presses, optimal_path = get_min_presses(direction_table, return_paths=True)

    for i in range(len(np.where(min_presses == np.max(min_presses))[0])):
        this_from, this_to = (
            np.where(min_presses == np.max(min_presses))[0][i],
            np.where(min_presses == np.max(min_presses))[1][i],
        )
        if optimal_path[this_from][this_to] is not None:
            this_path = optimal_path[this_from][this_to]
            current_loc = this_from
            goes_through = [names[current_loc]]
            for instruction in this_path:
                current_loc = direction_table[current_loc, instruction]
                goes_through.append(directions_in_conn_grid_order[instruction])
                goes_through.append(names[current_loc])
            print(" ".join(goes_through))
    plt.imshow(min_presses)
    plt.show()

    # potential modifications
    best_improvement = None
    best_improvement_score = None
    for i in range(direction_table.shape[0]):
        for theta_i in range(4):
            for j in range(direction_table.shape[0]):
                if connection_grid[i, j, theta_i] > 0:
                    new_direction_table = np.array(direction_table)
                    new_direction_table[i, theta_i] = j
                    if np.max(get_min_presses(new_direction_table)) < np.max(
                        get_min_presses(direction_table)
                    ):
                        if best_improvement is None:
                            best_improvement = [i, theta_i, j]
                            best_improvement_score = [
                                np.max(get_min_presses(new_direction_table)),
                                loss(new_direction_table, connection_grid),
                            ]
                        elif (
                            loss(new_direction_table, connection_grid)
                            > best_improvement_score[1]
                        ):
                            best_improvement = [i, theta_i, j]
                            best_improvement_score = [
                                np.max(get_min_presses(new_direction_table)),
                                loss(new_direction_table, connection_grid),
                            ]

                        print(
                            i,
                            theta_i,
                            j,
                            np.max(get_min_presses(new_direction_table)),
                            loss(new_direction_table, connection_grid),
                        )
    print(best_improvement, best_improvement_score)
    # do the improvement
    direction_table[best_improvement[0], best_improvement[1]] = best_improvement[2]

    min_presses, optimal_path = get_min_presses(direction_table, return_paths=True)
    for i in range(len(np.where(min_presses == np.max(min_presses))[0])):
        this_from, this_to = (
            np.where(min_presses == np.max(min_presses))[0][i],
            np.where(min_presses == np.max(min_presses))[1][i],
        )
        if optimal_path[this_from][this_to] is not None:
            this_path = optimal_path[this_from][this_to]
            current_loc = this_from
            goes_through = [names[current_loc]]
            for instruction in this_path:
                current_loc = direction_table[current_loc, instruction]
                goes_through.append(directions_in_conn_grid_order[instruction])
                goes_through.append(names[current_loc])
            print(" ".join(goes_through))

    plt.imshow(get_min_presses(direction_table))
    plt.show()
    # potential modifications
    finished = False
    while not finished:
        best_improvement = None
        best_improvement_score = None
        for i in range(direction_table.shape[0]):
            for theta_i in range(4):
                for j in range(direction_table.shape[0]):
                    if connection_grid[i, j, theta_i] > 0:
                        new_direction_table = np.array(direction_table)
                        new_direction_table[i, theta_i] = j
                        if np.max(get_min_presses(new_direction_table)) < np.max(
                            get_min_presses(direction_table)
                        ):
                            if best_improvement is None:
                                best_improvement = [i, theta_i, j]
                                best_improvement_score = [
                                    np.max(get_min_presses(new_direction_table)),
                                    loss(new_direction_table, connection_grid),
                                ]
                            elif (
                                loss(new_direction_table, connection_grid)
                                > best_improvement_score[1]
                            ):
                                best_improvement = [i, theta_i, j]
                                best_improvement_score = [
                                    np.max(get_min_presses(new_direction_table)),
                                    loss(new_direction_table, connection_grid),
                                ]

                            print(
                                i,
                                theta_i,
                                j,
                                np.max(get_min_presses(new_direction_table)),
                                loss(new_direction_table, connection_grid),
                            )
        if best_improvement is not None:
            direction_table[
                best_improvement[0], best_improvement[1]
            ] = best_improvement[2]
        else:
            finished = True
        print(best_improvement, best_improvement_score)
    plt.imshow(get_min_presses(direction_table))
    plt.show()
    ardsarsd
    RayPlurality("shapefiles/quintia_districts_with_geom_info.shp", plot=True)
    RayPlurality("shapefiles/quintia_districts_with_geom_info.shp", plot=True)

    connection_grid, centroid_arr, _ = connection_grid_first_collision_method(
        gdf, density_used=10
    )

    direction_table = np.zeros((connection_grid.shape[0], 4))
    direction_strings = [
        " " * pad_amount + "  " + (" " * (pad_amount)).join(directions_in_wanted_order)
    ]
    for i in range(connection_grid.shape[0]):
        line_string = f"{' '*(pad_amount-len(names[i])) + names[i]}:"
        for theta_i_i in range(4):
            theta_i = direction_dict[directions_in_wanted_order[theta_i_i]]
            if (
                connection_grid[i, i, theta_i]
                < np.sum(connection_grid[i, :, theta_i]) / 2
            ):
                connection_grid[i, i, theta_i] = 0
                direction_table[i, theta_i] = np.argmax(connection_grid[i, :, theta_i])
                line_string = f"{line_string} {' '*(pad_amount-len(names[int(direction_table[i,theta_i])])) + names[int(direction_table[i,theta_i])]}"
            else:
                direction_table[i, theta_i] = -1
                line_string = f"{line_string} {'-'*pad_amount}"
        direction_strings.append(line_string)
    direction_strings[1:] = sorted(direction_strings[1:])
    print("\n".join(direction_strings))
    print(connection_grid, centroid_arr)
    plt.imshow(direction_table)
    plt.show()
