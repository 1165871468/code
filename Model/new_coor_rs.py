#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from pathlib import Path

line_file = Path("../../Output/Model/profile2d/line.txt")
source_file = Path("../../Input/source_3d.txt")
receiver_file = Path("../../Input/receiver_3d.txt")

out_combined = Path("../../Output/Model/profile2d/source_receiver_2d_coordinates.txt")
out_source = Path("../../Output/Model/profile2d/source_new_2d.txt")
out_receiver = Path("../../Output/Model/profile2d/receiver_new_2d.txt")


def load_txt(path):
    return np.loadtxt(path, ndmin=2)


def project_to_2d(point, start, end):
    point = np.asarray(point, dtype=float)
    start = np.asarray(start, dtype=float)
    end = np.asarray(end, dtype=float)

    vector = end - start
    length = np.linalg.norm(vector)
    if length <= 0.0:
        raise ValueError("Profile length is zero. Please check line.txt.")

    unit_vector = vector / length
    x2 = float(np.dot(point - start, unit_vector))
    z2 = float(point[2])

    # Horizontal off-profile distance, used to check whether the point is on the profile.
    vector_xy = end[:2] - start[:2]
    length_xy = np.linalg.norm(vector_xy)
    if length_xy > 0.0:
        unit_xy = vector_xy / length_xy
        projected_xy = start[:2] + np.dot(point[:2] - start[:2], unit_xy) * unit_xy
        off_profile = float(np.linalg.norm(point[:2] - projected_xy))
    else:
        off_profile = float(np.linalg.norm(point[:2] - start[:2]))

    return x2, z2, off_profile


def main():
    line = load_txt(line_file)
    source = load_txt(source_file)
    receiver = load_txt(receiver_file)

    combined_rows = []
    source_rows = []
    receiver_rows = []

    for row in line:
        source_id = int(row[0])
        receiver_id = int(row[1])

        start = row[2:5]
        end = row[5:8]
        profile_length = float(row[8]) if len(row) >= 9 else float(np.linalg.norm(end - start))

        source_point = source[source_id - 1]
        receiver_point = receiver[receiver_id - 1]

        sx2, sz2, soff = project_to_2d(source_point, start, end)
        rx2, rz2, roff = project_to_2d(receiver_point, start, end)

        combined_rows.append([
            source_id, receiver_id, sx2, sz2, rx2, rz2, profile_length, soff, roff
        ])
        source_rows.append([source_id, receiver_id, sx2, sz2, soff])
        receiver_rows.append([source_id, receiver_id, rx2, rz2, roff])

    combined = pd.DataFrame(
        combined_rows,
        columns=[
            "source_id", "receiver_id",
            "source_x2_m", "source_z2_m",
            "receiver_x2_m", "receiver_z2_m",
            "profile_length_m",
            "source_off_profile_m", "receiver_off_profile_m",
        ],
    )
    source_new = pd.DataFrame(
        source_rows,
        columns=["source_id", "receiver_id", "source_x2_m", "source_z2_m", "source_off_profile_m"],
    )
    receiver_new = pd.DataFrame(
        receiver_rows,
        columns=["source_id", "receiver_id", "receiver_x2_m", "receiver_z2_m", "receiver_off_profile_m"],
    )

    combined.to_csv(out_combined, sep="\t", index=False, float_format="%.6f")
    source_new.to_csv(out_source, sep="\t", index=False, float_format="%.6f")
    receiver_new.to_csv(out_receiver, sep="\t", index=False, float_format="%.6f")

    print(f"Done: {out_combined}")
    print(f"Done: {out_source}")
    print(f"Done: {out_receiver}")


if __name__ == "__main__":
    main()
