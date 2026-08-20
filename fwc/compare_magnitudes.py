"""
For every flag in svg_star_ref.json, print each star's name, whether it's
found in hip_data.csv, and its magnitude (or "--" if missing) -- plus, once
looked up, the star's class id from star_geometry.json.

A star's "name" is normally a real star name (looked up by matching
hip_data.csv's name column), but a few stars in svg_star_ref.json that
lack a popular name in hip_data.csv are instead given a plain number --
the row index of the correct star in hip_data.csv, found by RA/Dec
rather than name. Those are looked up by row position instead.
"""
import csv
import json
import os
from matplotlib import pyplot as plt
import numpy as np
import scipy

import itertools

all_orders = np.array(list(itertools.permutations(range(5))))
all_orders4 = np.array(list(itertools.permutations(range(4))))

print(all_orders)
print(len(all_orders))

HERE = os.path.dirname(os.path.abspath(__file__))


def load_hip_rows():
    with open(os.path.join(HERE, "hip_data.csv"), encoding="utf-8") as f:
        return list(csv.reader(f))


def magnitude_of(rows, mags_by_name, name):
    if isinstance(name, int):
        return float(rows[name][4]) if 0 <= name < len(rows) else None
    return mags_by_name.get(name)


def radec_of(rows, radecs_by_name, name):
    if isinstance(name, int):
        return [float(rows[name][2]),float(rows[name][3])] if 0 <= name < len(rows) else None
    return radecs_by_name.get(name)


def star_class(flag_key, star_name):
    geometry = json.load(open(os.path.join(HERE, "star_geometry.json"), encoding="utf-8"))
    for s in geometry[flag_key]["stars"]:
        if s["name"] == star_name:
            return s.get("class")
    return None


def star_pos(flag_key, star_name):
    geometry = json.load(open(os.path.join(HERE, "star_geometry.json"), encoding="utf-8"))
    for s in geometry[flag_key]["stars"]:
        if s["name"] == star_name:
            return s.get("pos")
    return [0,0]


def scalerotmat(A, B):
    """
    Least-squares fit of M = sR minimizing ||AM^T - B||_F
    (equivalently np.matmul(M, A[:,:,None])[:,:,0] - B).

    Parameters
    ----------
    A, B : (n,2) arrays
        Corresponding vectors.

    Returns
    -------
    M : (2,2) array
        Best-fit scalar times rotation matrix.
    """
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)

    if A.shape != B.shape or A.shape[1] != 2:
        raise ValueError("A and B must both have shape (n,2)")

    # Cross-covariance
    C = B.T @ A

    U, S, Vt = np.linalg.svd(C)

    # Proper rotation (det = +1)
    D = np.eye(2)
    D[1, 1] = np.linalg.det(U @ Vt)

    R = U @ D @ Vt

    # Optimal scale
    s = np.trace(R @ C) / np.sum(A * A)

    return s * R


def make_scalerot(theta, scale):
    c = np.cos(theta)
    s = np.sin(theta)
    return scale * np.array([
        [c, -s],
        [s,  c]
    ])


def scalerot_inv(M):
    """
    Recover (theta, scale) from a 2x2 scaled rotation matrix.
    """
    M = np.asarray(M, dtype=float)

    scale = np.linalg.norm(M[:, 0])

    if scale == 0:
        raise ValueError("Matrix has zero scale.")

    R = M / scale

    theta = np.arctan2(R[1, 0], R[0, 0])

    return theta, scale
from tqdm import tqdm
import time

to_redo = []
to_redo_from_before = ['BR', 'FED', 'NIUE', 'TOK', 'NZMAT', 'ACT', 'CX', 'EUREKA', 'COLONIAL', 'SANTACRUZ', 'TDF', 'SONSOROL', 'ENB', 'NIRE', 'SIMBU', 'WNB', 'WPNG', 'SPLOUGHBLUE', 'SPLOUGHINLA', 'SPLOUGHLABOUR', 'EKARELIAP1', 'EKARELIAP2', 'NT', 'MAG', 'PNG']

to_redo_from_before = ['NIUE', 'TOK', 'NZMAT', 'ACT', 'CX', 'EUREKA', 'COLONIAL', 'SANTACRUZ', 'TDF', 'SONSOROL', 'ENB', 'NIRE', 'SIMBU', 'WNB', 'WPNG', 'SPLOUGHBLUE', 'SPLOUGHINLA', 'SPLOUGHLABOUR', 'EKARELIAP1', 'EKARELIAP2', 'NT', 'MAG', 'PNG']

to_redo_from_before = ['FED', 'EUREKA', 'SONSOROL', 'NIRE', 'WNB', 'MAG', 'PNG']

best_by_name = {}

def main():
    ref = json.load(open(os.path.join(HERE, "svg_star_ref.json"), encoding="utf-8"))
    rows = load_hip_rows()
    mags_by_name = {row[1]: float(row[4]) for row in rows}
    radecs_by_name = {row[1]: [float(row[2]),float(row[3])] for row in rows}

    for flag_key, entry in ref.items():
        if flag_key=="SPLOUGH1914":
            continue
        if flag_key not in to_redo_from_before:
            continue
        # filename = "svg_star_ref.json"

        # with open(filename, "r") as f:
        #     data = json.load(f)

        # if len(data[flag_key]["stars"])>5:
        #     continue

        # crux_star_name_list = []
        # for i in range(len(data[flag_key]["stars"])):
        #     crux_star_name_list.append(data["AU"]["stars"][i]["name"])

        # if data[flag_key]["stars"][0]["name"]==data[flag_key]["stars"][1]["name"]:
        #     for i in range(len(data[flag_key]["stars"])):
        #         data[flag_key]["stars"][i]["name"] = crux_star_name_list[i]
        # with open(filename, "w") as f:
        #     json.dump(data, f, indent=4)
        # continue

        print(flag_key)

        pos_list = []
        name_list = []
        latlon_list = []

        for star in entry["stars"]:
            name = star["name"]
            # print(name)
            mag = magnitude_of(rows, mags_by_name, name)
            latlon = radec_of(rows, radecs_by_name, name)
            mag_str = f"{mag}" if mag is not None else "--"
            cls = star_class(flag_key, name)
            pos_list.append(star_pos(flag_key, name))
            name_list.append(name)
            latlon_list.append(latlon)
            # print(f"  {name}: {'found' if mag is not None else 'missing'}, magnitude {mag_str}, class {cls}")

        pos_list = np.array(pos_list)
        latlon_list = np.array(latlon_list)
        xyz_list = np.vstack([np.cos(latlon_list[:,0])*np.cos(latlon_list[:,1]),
                            np.sin(latlon_list[:,0])*np.cos(latlon_list[:,1]),
                            np.sin(latlon_list[:,1])]).T
        # print(xyz_list)
        xyz_center = np.mean(xyz_list,axis=0)
        right_vec = np.cross(xyz_center,np.array([0,0,1]))
        up_vec = np.cross(xyz_center,right_vec)
        xyz_center /= np.linalg.norm(xyz_center)
        right_vec /= np.linalg.norm(right_vec)
        up_vec /= np.linalg.norm(up_vec)
        rotmat = np.linalg.inv(np.vstack([xyz_center,right_vec,up_vec]).T)
        rotted_xyz = np.matmul(rotmat,xyz_list[:,:,None])[:,:,0]
        ortho = rotted_xyz[:,1:]
        scalerot = scalerotmat(ortho, pos_list-np.mean(pos_list,axis=0)[None,:])
        ortho = np.matmul(scalerot,ortho[:,:,None])[:,:,0]
        ortho += np.mean(pos_list,axis=0)[None,:]


        def do_positioning(x,points_xyz,flip=False):
            points_xyz_use = np.array(points_xyz)
            center_lon = x[0]
            center_lat = x[1]
            center_xyz = np.array([np.cos(center_lon)*np.cos(center_lat),
                            np.sin(center_lon)*np.cos(center_lat),
                            np.sin(center_lat)]).T
            c_right_vec = np.cross(center_xyz,np.array([0,0,1]))
            c_up_vec = np.cross(center_xyz,right_vec)
            center_xyz /= np.linalg.norm(center_xyz)
            c_right_vec /= np.linalg.norm(c_right_vec)
            c_up_vec /= np.linalg.norm(c_up_vec)
            c_rotmat = np.linalg.inv(np.vstack([center_xyz,c_right_vec,c_up_vec]).T)

            rotation = x[2]
            scale = x[3]
            c_scalerot = make_scalerot(rotation,scale)

            new_center = np.array([x[4],x[5]])
            new_pos = new_center+np.matmul(c_scalerot,np.matmul(c_rotmat,points_xyz_use[:,:,None])[:,1:]*(np.array([1,-1]) if flip else np.array([1,1]))[None,:])[:,:,0]
            return new_pos

        to_flip = True
        best_order = None
        best_loss = 100000000000000000
        best_x = None
        all_orders_use = None
        if len(pos_list)==5:
            all_orders_use = all_orders
        elif len(pos_list)==4:
            all_orders_use = all_orders4
        else:
            all_orders_use = [np.arange(len(pos_list))]
        for order in tqdm(all_orders_use[:1]):
            # print(order)
            for flip in [True,False]:
                def loss(x):
                    new_pos = do_positioning(x,xyz_list,flip=flip)
                    return np.mean(np.linalg.norm(new_pos-pos_list[order],axis=1)**2)/x[3]

                x_0 = np.zeros((6,))
                theta,scale_ = scalerot_inv(scalerot)
                x_0[0] = np.arctan2(xyz_center[1],xyz_center[0])
                x_0[1] = np.arcsin(xyz_center[2])
                x_0[2] = theta
                x_0[3] = scale_
                x_0[4] = np.mean(pos_list,axis=0)[0]
                x_0[5] = np.mean(pos_list,axis=0)[1]

                res = scipy.optimize.minimize(loss,x_0)
                # print(res)
                # print(res.x)
                if loss(res.x)<best_loss:
                    to_flip = flip
                    best_loss = loss(res.x)
                    best_order = order
                    best_x = res.x

        new_pos = do_positioning(best_x,xyz_list,flip=to_flip)

        # print(f"best for {flag_key} is {best_order}")
        best_by_name[flag_key] = best_order
        print(best_by_name)
        print(best_loss)


        # filename = "svg_star_ref.json"

        # with open(filename, "r") as f:
        #     data = json.load(f)

        # star_name_list = []
        # for i in range(len(best_order)):
        #     star_name_list.append(data[flag_key]["stars"][i]["name"])
        # for i in range(len(best_order)):
        #     # Swap the names
        #     data[flag_key]["stars"][best_order[i]]["name"] = star_name_list[i]

        # # Save back to JSON
        # with open(filename, "w") as f:
        #     json.dump(data, f, indent=4)

        def loss(x):
            new_pos = do_positioning(x,xyz_list,flip=to_flip)
            return np.mean(np.linalg.norm(new_pos-pos_list[best_order],axis=1)**2)/x[3]

        if loss(best_x)*100>20:
            to_redo.append(flag_key)
            print("to_redo =",to_redo)

        # print(scalerot)
        plt.title(f"{flag_key}: {loss(best_x)*100}")
        plt.plot(pos_list[:,0],pos_list[:,1])
        plt.plot(ortho[:,0],ortho[:,1])
        plt.plot(new_pos[:,0],new_pos[:,1])
        plt.scatter(pos_list[:2,0],pos_list[:2,1])
        plt.scatter(ortho[:2,0],ortho[:2,1])
        plt.scatter(new_pos[:2,0],new_pos[:2,1])





        points = np.asarray(pos_list)
        for (x, y), label in zip(points, name_list):
            plt.text(x, y, label)





        plt.show()

if __name__ == "__main__":
    main()
