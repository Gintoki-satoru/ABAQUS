from __future__ import annotations

from pathlib import Path
import csv
from typing import Dict, Tuple, Optional

from odbAccess import openOdb
from abaqusConstants import INTEGRATION_POINT


def find_element_set(odb, set_name: str, instance_name: Optional[str] = None):
    """
    Returns (instanceObj, elementSetObj). Requires instance_name if the set is not uniquely found.
    """
    ra = odb.rootAssembly
    # If instance is provided, resolve directly
    if instance_name is not None:
        if instance_name not in ra.instances:
            raise KeyError(f"Instance '{instance_name}' not found. Available: {list(ra.instances.keys())}")
        inst = ra.instances[instance_name]
        if set_name not in inst.elementSets:
            raise KeyError(f"Element set '{set_name}' not found in instance '{instance_name}'.")
        return inst, inst.elementSets[set_name]
    # Try assembly-level set first (less common for element sets)
    if set_name in ra.elementSets:
        # Assembly-level sets can include multiple instances -> ambiguous for centroid coords
        raise RuntimeError(
            f"Element set '{set_name}' found at assembly-level. "
            f"Please pass instance_name explicitly to avoid ambiguity."
        )
    # Search instances
    matches = []
    for iname, inst in ra.instances.items():
        if set_name in inst.elementSets:
            matches.append((inst, inst.elementSets[set_name], iname))
    if not matches:
        raise KeyError(f"Element set '{set_name}' not found in any instance.")
    if len(matches) > 1:
        names = [m[2] for m in matches]
        raise RuntimeError(
            f"Element set '{set_name}' found in multiple instances: {names}. "
            f"Pass instance_name explicitly."
        )
    inst, eset, _iname = matches[0]
    return inst, eset


def element_centroids_from_odb_instance(inst, element_set) -> Dict[int, Tuple[float, float, float]]:
    """
    Compute centroid of each element in element_set by averaging its node coordinates in the ODB instance.
    Returns dict: elemLabel -> (xc, yc, zc)
    """
    # Node label -> coordinates
    node_coord: Dict[int, Tuple[float, float, float]] = {}
    for nd in inst.nodes:
        node_coord[nd.label] = nd.coordinates  # tuple (x,y,z)
    centroids: Dict[int, Tuple[float, float, float]] = {}
    for e in element_set.elements:
        conn = e.connectivity  # node labels
        xs = ys = zs = 0.0
        ncnt = len(conn)
        for nl in conn:
            x, y, z = node_coord[nl]
            xs += x
            ys += y
            zs += z
        centroids[e.label] = (xs / ncnt, ys / ncnt, zs / ncnt)
    return centroids


def extract_set_stress_and_centroid(
    odb_path: str,
    set_name: str,
    out_csv_path: str,
    instance_name: Optional[str] = None,
    step_name: Optional[str] = None,   # None -> last step
    frame_index: int = -1,             # -1 -> last frame
    average_over_ips: bool = False     # True -> one row per element (IP-averaged)
) -> None:
    odb = openOdb(path=odb_path, readOnly=True)
    # Step/frame
    if step_name is None:
        step_name = list(odb.steps.keys())[-1]
    step = odb.steps[step_name]
    frame = step.frames[frame_index]
    # Stress field
    if "S" not in frame.fieldOutputs:
        odb.close()
        raise RuntimeError("Stress output 'S' not found in the selected frame. Request 'S' in Field Output.")
    Sfield = frame.fieldOutputs["S"].getSubset(position=INTEGRATION_POINT)
    # Resolve element set + instance
    inst, elem_set = find_element_set(odb, set_name=set_name, instance_name=instance_name)
    # Centroids
    elem_centroid = element_centroids_from_odb_instance(inst, elem_set)
    elem_labels_in_set = set(elem_centroid.keys())
    # Subset stress values to the region (fast)
    S_sub = Sfield.getSubset(region=elem_set)
    out_path = Path(out_csv_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # Prepare CSV
    if average_over_ips:
        header = ["elementLabel", "x_c", "y_c", "z_c", "S11", "S22", "S33", "S12", "S13", "S23"]
        accum: Dict[int, Tuple[list, int]] = {}  # elem -> ([sum6], count)
    else:
        header = ["elementLabel", "integrationPoint", "x_c", "y_c", "z_c",
                  "S11", "S22", "S33", "S12", "S13", "S23"]
    with out_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for v in S_sub.values:
            el = v.elementLabel
            if el not in elem_labels_in_set:
                continue
            xc, yc, zc = elem_centroid[el]
            # v.data is tuple; for 3D stresses typically (S11,S22,S33,S12,S13,S23)
            data = v.data
            s = [0.0] * 6
            for i in range(min(6, len(data))):
                s[i] = float(data[i])
            if average_over_ips:
                if el not in accum:
                    accum[el] = ([0.0] * 6, 0)
                sums, cnt = accum[el]
                for i in range(6):
                    sums[i] += s[i]
                accum[el] = (sums, cnt + 1)
            else:
                ip = getattr(v, "integrationPoint", None)
                w.writerow([el, ip, xc, yc, zc, s[0], s[1], s[2], s[3], s[4], s[5]])
        if average_over_ips:
            for el, (sums, cnt) in accum.items():
                if cnt == 0:
                    continue
                xc, yc, zc = elem_centroid[el]
                avg = [val / float(cnt) for val in sums]
                w.writerow([el, xc, yc, zc, avg[0], avg[1], avg[2], avg[3], avg[4], avg[5]])
    odb.close()
    print(f"Wrote: {out_path.resolve()}")


if __name__ == "__main__":
    # ----------- USER INPUTS -----------
    odb_path = r"Job-1.odb"
    set_name = "PLY_01"            # your element set (layer)
    instance_name = "SUPERELLIPSOID_2D-1"     # change to your instance
    out_csv = r"PLY_S_and_centroid.csv"
    extract_set_stress_and_centroid(
        odb_path=odb_path,
        set_name=set_name,
        out_csv_path=out_csv,
        instance_name=instance_name,
        step_name=None,          # last step
        frame_index=-1,          # last frame
        average_over_ips=False 
    )


###############################################          Node Stress         ####################################################

from __future__ import annotations

from pathlib import Path
import csv
from typing import Optional, Dict, Tuple, List

from odbAccess import openOdb
from abaqusConstants import ELEMENT_NODAL, NODAL


def _get_last_step_and_frame(odb, step_name: Optional[str], frame_index: int):
    if step_name is None:
        step_name = list(odb.steps.keys())[-1]
    step = odb.steps[step_name]
    frame = step.frames[frame_index]
    return step_name, frame


def _resolve_node_set(odb, set_name: str, instance_name: Optional[str] = None):
    """
    Returns (instanceObj, nodeSetObj). Requires instance_name if set is ambiguous.
    In ODB, sets live on instances (most common) or assembly.
    """
    ra = odb.rootAssembly
    if instance_name is not None:
        if instance_name not in ra.instances:
            raise KeyError(f"Instance '{instance_name}' not found. Available: {list(ra.instances.keys())}")
        inst = ra.instances[instance_name]
        if set_name in inst.nodeSets:
            return inst, inst.nodeSets[set_name]
        if set_name in ra.nodeSets:
            # assembly-level node set; still usable, but coords need instance disambiguation
            return None, ra.nodeSets[set_name]
        raise KeyError(f"Node set '{set_name}' not found in instance '{instance_name}' nor assembly.")
    # Try assembly-level first
    if set_name in ra.nodeSets:
        return None, ra.nodeSets[set_name]
    # Search instances
    matches = []
    for iname, inst in ra.instances.items():
        if set_name in inst.nodeSets:
            matches.append((inst, inst.nodeSets[set_name], iname))
    if not matches:
        raise KeyError(f"Node set '{set_name}' not found in assembly or any instance.")
    if len(matches) > 1:
        names = [m[2] for m in matches]
        raise RuntimeError(
            f"Node set '{set_name}' found in multiple instances: {names}. "
            f"Pass instance_name explicitly."
        )
    inst, nset, _ = matches[0]
    return inst, nset


def extract_nodeset_S_and_coords(
    odb_path: str,
    node_set_name: str,
    out_csv_path: str,
    instance_name: Optional[str] = None,
    step_name: Optional[str] = None,
    frame_index: int = -1,
    average_per_node: bool = True,
) -> None:
    odb = openOdb(path=odb_path, readOnly=True)
    step_name, frame = _get_last_step_and_frame(odb, step_name, frame_index)
    # --- Stress field (prefer ELEMENT_NODAL) ---
    if "S" not in frame.fieldOutputs:
        odb.close()
        raise RuntimeError("Field output 'S' not found. Request stress output (S) in the analysis.")
    # ELEMENT_NODAL gives stress values at nodes extrapolated from IPs (per element-node)
    S_el_nodal = frame.fieldOutputs["S"].getSubset(position=ELEMENT_NODAL)
    # --- Resolve node set ---
    inst, nset = _resolve_node_set(odb, node_set_name, instance_name)
    # Collect node labels in the set
    # Note: nset.nodes is a sequence of node objects, possibly grouped by instance
    node_labels: List[int] = []
    node_coords: Dict[Tuple[str, int], Tuple[float, float, float]] = {}  # (instName, nodeLabel) -> coords
    if inst is not None:
        # Node set on a specific instance
        inst_name = inst.name
        for nd in nset.nodes:
            node_labels.append(nd.label)
            node_coords[(inst_name, nd.label)] = nd.coordinates
        inst_names = [inst_name]
    else:
        # Assembly-level set: may include nodes from multiple instances
        inst_names = []
        for nd in nset.nodes:
            # In assembly-level, Abaqus node object may carry instance internally but not exposed uniformly.
            # Safer: user should provide instance_name. We'll try best-effort by searching instances by label.
            pass
        odb.close()
        raise RuntimeError(
            "Assembly-level node set detected. Please pass instance_name so coordinates are unambiguous."
        )
    node_label_set = set(node_labels)
    # --- Filter stress values by this node set ---
    # Values have attributes: instance, elementLabel, nodeLabel, data
    # We'll either average per node (recommended) or write per element-node.
    out_path = Path(out_csv_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as f:
        w = csv.writer(f)
        if average_per_node:
            w.writerow(["odb", "step", "frameIndex", "frameValue",
                        "instance", "nodeLabel", "x", "y", "z",
                        "S11", "S22", "S33", "S12", "S13", "S23",
                        "nContrib"])
            # accum: (inst, node) -> sums(6), count
            accum: Dict[Tuple[str, int], Tuple[List[float], int]] = {}
            for sv in S_el_nodal.values:
                nlab = getattr(sv, "nodeLabel", None)
                if nlab is None or nlab not in node_label_set:
                    continue
                iname = sv.instance.name if sv.instance else inst_names[0]
                key = (iname, nlab)
                data = sv.data
                s = [0.0] * 6
                for i in range(min(6, len(data))):
                    s[i] = float(data[i])
                if key not in accum:
                    accum[key] = ([0.0] * 6, 0)
                sums, cnt = accum[key]
                for i in range(6):
                    sums[i] += s[i]
                accum[key] = (sums, cnt + 1)
            for (iname, nlab), (sums, cnt) in sorted(accum.items(), key=lambda x: (x[0][0], x[0][1])):
                x, y, z = node_coords[(iname, nlab)]
                avg = [val / float(cnt) for val in sums]
                w.writerow([Path(odb_path).name, step_name, frame_index, frame.frameValue,
                            iname, nlab, x, y, z,
                            avg[0], avg[1], avg[2], avg[3], avg[4], avg[5],
                            cnt])
        else:
            w.writerow(["odb", "step", "frameIndex", "frameValue",
                        "instance", "elementLabel", "nodeLabel", "x", "y", "z",
                        "S11", "S22", "S33", "S12", "S13", "S23"])
            for sv in S_el_nodal.values:
                nlab = getattr(sv, "nodeLabel", None)
                if nlab is None or nlab not in node_label_set:
                    continue
                iname = sv.instance.name if sv.instance else inst_names[0]
                elab = getattr(sv, "elementLabel", None)
                x, y, z = node_coords[(iname, nlab)]
                data = sv.data
                s = [0.0] * 6
                for i in range(min(6, len(data))):
                    s[i] = float(data[i])
                w.writerow([Path(odb_path).name, step_name, frame_index, frame.frameValue,
                            iname, elab, nlab, x, y, z,
                            s[0], s[1], s[2], s[3], s[4], s[5]])
    odb.close()
    print(f"Wrote: {out_path.resolve()}")


if __name__ == "__main__":
    # -------- USER INPUTS --------
    odb_path = r"Job-1.odb"
    instance_name = "SUPERELLIPSOID_2D-1"
    node_set_name = "PLY_00_IFACE"
    out_csv = r"PLY_stress_nodal.csv"
    extract_nodeset_S_and_coords(
        odb_path=odb_path,
        node_set_name=node_set_name,
        out_csv_path=out_csv,
        instance_name=instance_name,
        step_name=None,       # last step
        frame_index=-1,       # last frame
        average_per_node=True # recommended
    )

