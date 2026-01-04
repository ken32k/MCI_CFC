# Import packages
import os
from tqdm import tqdm
import numpy as np
import nibabel as nb
import scipy
from scipy import spatial
import pygeodesic.geodesic as geodesic
import scipy.io as sio
import pandas as pd
import nilearn.surface

CSVD_DIR = "/public/home/baishw/WMH_MCI"

euc_output_dir = f"{CSVD_DIR}/data/euc_distance"
mica_output_dir = f"{CSVD_DIR}/csvd_mica_out/micapipe_v0.2.0"
out_sub_dirs = sorted(os.listdir(mica_output_dir))

fs_dir = f"{CSVD_DIR}/csvd_mica_out/freesurfer"


def get_euclidean_distance(sub, atlas="schaefer-200", matrix="euc"):
    atlas = atlas
    matrix = matrix

    rh_annot = fs_dir + "/sub-" + sub + "/label/rh.schaefer-200_mics.annot"
    lh_annot = fs_dir + "/sub-" + sub + "/label/lh.schaefer-200_mics.annot"

    lh_surf = (
        mica_output_dir
        + "/sub-"
        + sub
        + "/surf/sub-"
        + sub
        + "_hemi-L_surf-fsnative_label-midthickness.surf.gii"
    )
    rh_surf = (
        mica_output_dir
        + "/sub-"
        + sub
        + "/surf/sub-"
        + sub
        + "_hemi-R_surf-fsnative_label-midthickness.surf.gii"
    )

    # Load LEFT surfaces
    # lh = nb.load(rh_surf)
    # vertices_lh = lh.agg_data('NIFTI_INTENT_POINTSET')
    vertices_lh = nilearn.surface.load_surf_mesh(lh_surf).coordinates

    # Load RIGHT surfaces
    # rh = nb.load(rh_surf)
    # vertices_rh = rh.agg_data('NIFTI_INTENT_POINTSET')

    vertices_rh = nilearn.surface.load_surf_mesh(rh_surf).coordinates

    vertices = np.append(vertices_lh, vertices_rh, axis=0)
    # faces = np.append(faces_lh, faces_rh+len(vertices_lh), axis = 0)

    [labels_lh, ctab_lh, names_lh] = nb.freesurfer.io.read_annot(
        lh_annot, orig_ids=True
    )
    [labels_rh, ctab_rh, names_rh] = nb.freesurfer.io.read_annot(
        rh_annot, orig_ids=True
    )

    nativeLength = len(labels_lh) + len(labels_rh)

    parc = np.zeros((nativeLength))
    for x, _ in enumerate(labels_lh):
        parc[x] = np.where(ctab_lh[:, 4] == labels_lh[x])[0][0]
    for x, _ in enumerate(labels_rh):
        parc[x + len(labels_lh)] = np.where(ctab_rh[:, 4] == labels_rh[x])[0][0] + len(
            ctab_lh
        )

    # Find centre vertex
    uparcel = np.unique(parc)
    uparcel = np.delete(uparcel, 200)
    uparcel = np.delete(uparcel, 0)
    voi = np.zeros([1, len(uparcel)])
    center_coords = np.zeros([len(uparcel), 3])

    # print("[Info]..... Finding central vertex for each parcel")
    for n, _ in enumerate(uparcel):
        this_parc = np.where(parc == uparcel[n])[0]
        distances = spatial.distance.pdist(
            np.squeeze(vertices[this_parc, :]), "euclidean"
        )  # Returns condensed matrix of distances
        distancesSq = spatial.distance.squareform(distances)  # convert to square form
        sumDist = np.sum(distancesSq, axis=1)  # sum distance across columns
        # minimum sum distance index
        index = np.where(sumDist == np.min(sumDist))
        voi[0, n] = this_parc[index[0][0]]

        center_coords[n, :] = vertices[this_parc[index[0][0]]]

    ED_matrix = scipy.spatial.distance_matrix(center_coords, center_coords)
    ED_matrix = pd.DataFrame(ED_matrix)
    sio.savemat(
        euc_output_dir + "/sub-" + sub + "_" + matrix + "_" + atlas + ".mat",
        {"ed": ED_matrix},
    )


if __name__ == "__main__":
    for subject in tqdm(out_sub_dirs):
        sub = subject[-3:]
        if not os.path.exists(euc_output_dir + "/sub-" + sub + "_euc_schaefer-200.mat"):
            try:
                get_euclidean_distance(sub)
                print("[Info]..... subject %s: done." % sub, flush=True)
            except Exception as e:
                print("[Error].... subject %s: %s" % (sub, e), flush=True)
        else:
            print("[Info]..... subject %s: Already exist." % sub, flush=True)
