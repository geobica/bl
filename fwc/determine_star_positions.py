"""
Programmatically determine every star's position, size class, point-count,
rotation and color purely from the raw SVG geometry named by svg_star_ref.json's
svg_loc ranges -- no manually-curated pos/class/classes/star_color values are
read or trusted; this script derives all of it from scratch.

For each star:
  - the svg_loc snippet(s) are expanded through any <use>/<g> chain (composing
    transforms, including the ancestor chain in the original document -- found
    via an expat-based char-offset element index) down to leaf <path>/
    <polygon>/<polyline> shapes, each with its final absolute vertices and
    inherited fill color
  - if the combined leaf vertices cleanly alternate between two tight radius
    bands around their centroid, the star is a regular idealized N-pointed
    star: record (pos, outer, inner, points, rotation)
  - otherwise (irregular hand-traced art, or a compound multi-layer star built
    from independently-overlapping pieces that doesn't reduce to one clean
    alternation) the star is recorded as "exact_d": a list of absolute path
    "d" strings, one per leaf shape (kept separate, not merged into one
    multi-subpath path, since merging measurably changes anti-aliasing at
    shared edges) -- these are replayed verbatim by reconstruct_and_verify.py

Regular stars are then clustered per-flag by (outer radius, points, rotation)
into size classes, sorted biggest to smallest. Flags where the SAME star also
has a second, larger same-centroid layer (e.g. EUREKA's black outline behind
a colored fill) get an additional "outline_classes"/"outline_color" pair.

Writes new/star_geometry.json: {key: {"stars": [...], "classes": [[outer,
inner,points,rotation], ...], "star_color", ["outline_classes",
"outline_color"]}}. Deliberately does NOT include "center" -- that lives in
svg_star_ref.json already and reconstruct_and_verify.py doesn't need it.
"""
import json
import math
import os
import re
import statistics
import xml.etree.ElementTree as ET
import xml.parsers.expat

from shapely.geometry import Polygon
from shapely.ops import unary_union

HERE = os.path.dirname(os.path.abspath(__file__))
XLINK_HREF = "{http://www.w3.org/1999/xlink}href"
WRAPPER = ('<svg xmlns="http://www.w3.org/2000/svg" '
           'xmlns:xlink="http://www.w3.org/1999/xlink" '
           'xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" '
           'xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd" '
           'xmlns:dc="http://purl.org/dc/elements/1.1/" '
           'xmlns:cc="http://creativecommons.org/ns#" '
           'xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">{}</svg>')


# ---------------------------------------------------------------- transforms

def parse_transform_list(t):
    """Compose every transform function in an attribute (SVG allows several
    space-separated functions, applied left to right) into one 2x3 affine."""
    M = (1, 0, 0, 1, 0, 0)
    if not t:
        return M
    for fn, args in re.findall(r"(\w+)\(([^)]*)\)", t):
        vals = [float(v) for v in re.split(r"[ ,]+", args.strip()) if v]
        if fn == "translate":
            tx = vals[0]; ty = vals[1] if len(vals) > 1 else 0.0
            F = (1, 0, 0, 1, tx, ty)
        elif fn == "scale":
            sx = vals[0]; sy = vals[1] if len(vals) > 1 else sx
            F = (sx, 0, 0, sy, 0, 0)
        elif fn == "matrix":
            F = tuple(vals)
        elif fn == "rotate":
            deg = vals[0]
            rad = math.radians(deg)
            ca, sa = math.cos(rad), math.sin(rad)
            if len(vals) >= 3:
                cx, cy = vals[1], vals[2]
                F = compose((1, 0, 0, 1, cx, cy),
                             compose((ca, sa, -sa, ca, 0, 0), (1, 0, 0, 1, -cx, -cy)))
            else:
                F = (ca, sa, -sa, ca, 0, 0)
        else:
            continue
        M = compose(M, F)
    return M


def compose(M, F):
    """M ∘ F -- apply F first, then M."""
    a1, b1, c1, d1, e1, f1 = M
    a2, b2, c2, d2, e2, f2 = F
    return (a1 * a2 + c1 * b2, b1 * a2 + d1 * b2,
            a1 * c2 + c1 * d2, b1 * c2 + d1 * d2,
            a1 * e2 + c1 * f2 + e1, b1 * e2 + d1 * f2 + f1)


def apply_affine(M, x, y):
    a, b, c, d, e, f = M
    return (a * x + c * y + e, b * x + d * y + f)


# --------------------------------------------------------------- path parse

PARAM_COUNTS = {"M": 2, "L": 2, "H": 1, "V": 1, "C": 6, "S": 4, "Q": 4, "T": 2, "A": 7, "Z": 0}
_TOKEN_RE = re.compile(r"([MmLlHhVvCcSsQqTtAaZz])|(-?\d*\.\d+(?:[eE]-?\d+)?|-?\d+(?:[eE]-?\d+)?)")


def path_vertices(d):
    """Absolute (x,y) vertices of an SVG path 'd' (M/L/H/V exact; C/S/Q/T/A
    endpoints only, control points ignored -- fine for the straight-edged
    star shapes this covers)."""
    toks = [(m.group(1) or m.group(2), m.group(1) is not None) for m in _TOKEN_RE.finditer(d)]
    i, n = 0, len(toks)
    cur = start = (0.0, 0.0)
    cmd = None
    pts = []
    while i < n:
        tokval, is_cmd = toks[i]
        if is_cmd:
            cmd = tokval
            i += 1
        upper = cmd.upper()
        if upper == "Z":
            # closing back to 'start' duplicates the first vertex -- skip it,
            # every caller here re-closes with its own explicit Z anyway
            cur = start
            continue
        pcount = PARAM_COUNTS[upper]
        vals = [float(toks[i + k][0]) for k in range(pcount)]
        i += pcount
        relative = cmd.islower()
        if upper in ("M", "L", "T"):
            x, y = vals
            if relative:
                x += cur[0]; y += cur[1]
            cur = (x, y)
        elif upper == "H":
            x = vals[0] + cur[0] if relative else vals[0]
            cur = (x, cur[1])
        elif upper == "V":
            y = vals[0] + cur[1] if relative else vals[0]
            cur = (cur[0], y)
        elif upper in ("C", "S", "Q", "A"):
            x, y = vals[-2], vals[-1]
            if relative:
                x += cur[0]; y += cur[1]
            cur = (x, y)
        pts.append(cur)
        if cmd in ("M", "m"):
            start = cur
            cmd = "l" if cmd == "m" else "L"
    return pts


def path_d_absolute(d, M):
    """Re-emit 'd' with every coordinate (including curve control points --
    NOT just vertices/endpoints) transformed by M and written out as an
    absolute path. Used for exact_d fallback output, where curve fidelity
    matters (unlike the flattened vertices path_vertices() gives, which are
    only used for the shapely-union classification step, a fine polygon
    approximation for that purpose but not for verbatim geometry replay)."""
    toks = [(mm.group(1) or mm.group(2), mm.group(1) is not None) for mm in _TOKEN_RE.finditer(d)]
    i, n = 0, len(toks)
    cur = start = (0.0, 0.0)
    cmd = None
    out = []
    while i < n:
        tokval, is_cmd = toks[i]
        if is_cmd:
            cmd = tokval
            i += 1
        upper = cmd.upper()
        if upper == "Z":
            cur = start
            out.append("Z")
            continue
        pcount = PARAM_COUNTS[upper]
        vals = [float(toks[i + k][0]) for k in range(pcount)]
        i += pcount
        relative = cmd.islower()
        if upper == "H":
            x = vals[0] + cur[0] if relative else vals[0]
            cur = (x, cur[1])
            px, py = apply_affine(M, *cur)
            out.append(f"L{px:.4f},{py:.4f}")
        elif upper == "V":
            y = vals[0] + cur[1] if relative else vals[0]
            cur = (cur[0], y)
            px, py = apply_affine(M, *cur)
            out.append(f"L{px:.4f},{py:.4f}")
        else:
            local_pts = list(zip(vals[0::2], vals[1::2]))
            abs_pts = []
            for vx, vy in local_pts:
                if relative:
                    vx += cur[0]; vy += cur[1]
                abs_pts.append((vx, vy))
            cur = abs_pts[-1]
            tpts = [apply_affine(M, x, y) for x, y in abs_pts]
            out.append(upper + " ".join(f"{x:.4f},{y:.4f}" for x, y in tpts))
        if cmd in ("M", "m"):
            start = cur
            cmd = "l" if cmd == "m" else "L"
    return "".join(out) + "Z"


def own_vertices(elem):
    tag = elem.tag.split("}")[-1]
    if tag == "path" and elem.get("d"):
        return elem.get("d"), path_vertices(elem.get("d"))
    if tag in ("polygon", "polyline") and elem.get("points"):
        nums = [float(v) for v in re.split(r"[ ,]+", elem.get("points").strip()) if v]
        pts = list(zip(nums[0::2], nums[1::2]))
        d = "M" + " L".join(f"{x},{y}" for x, y in pts) + " Z"
        return d, pts
    return None, None


def parse_css_classes(text):
    """Extract {classname: fill-color} from a <style>.name{fill:#xxx}</style>
    block (several flags set star colors via CSS class instead of a fill/
    style attribute)."""
    css = {}
    m = re.search(r"<style[^>]*>(.*?)</style>", text, re.S)
    if not m:
        return css
    for name, body in re.findall(r"\.([\w-]+)\s*\{([^}]*)\}", m.group(1)):
        fm = re.search(r"fill:\s*([^;]+)", body)
        if fm:
            css[name] = fm.group(1).strip()
    return css


def own_fill(elem, css):
    f = elem.get("fill")
    if f and f != "none":
        return f
    style = elem.get("style")
    if style:
        m = re.search(r"fill:\s*([^;]+)", style)
        if m and m.group(1).strip() != "none":
            return m.group(1).strip()
    for cls in (elem.get("class") or "").split():
        if cls in css:
            return css[cls]
    return None


# ------------------------------------------------------- char-offset index

def build_element_index(text):
    """Parse with expat to get every element's char-offset [start,end) span
    (covering its full open..close tag range) and parent link, keyed by
    start offset, so a raw svg_loc offset can be matched to its element and
    that element's ancestor chain (for transforms/fill inherited from
    outside the extracted snippet)."""
    encoded = text.encode("utf-8")
    nodes = []
    stack = []

    def start(name, attrs):
        node = {"tag": name, "attrs": attrs, "start_b": parser.CurrentByteIndex,
                "end_b": None, "parent": stack[-1] if stack else None}
        nodes.append(node)
        stack.append(node)

    def end(name):
        node = stack.pop()
        idx = encoded.index(b">", parser.CurrentByteIndex) + 1
        node["end_b"] = idx

    parser = xml.parsers.expat.ParserCreate()
    parser.StartElementHandler = start
    parser.EndElementHandler = end
    parser.Parse(encoded, True)

    def b2c(b):
        return len(encoded[:b].decode("utf-8"))

    by_start = {}
    for node in nodes:
        node["start"] = b2c(node["start_b"])
        node["end"] = b2c(node["end_b"])
        by_start.setdefault(node["start"], node)
    return by_start, nodes


def find_node_at(by_start, nodes, offset):
    node = by_start.get(offset)
    if node is not None:
        return node
    # fall back to nearest containing node (small slack for whitespace)
    best = None
    for node in nodes:
        if node["start"] <= offset <= node["end"] and (best is None or node["start"] > best["start"]):
            best = node
    return best


def ancestor_chain(node):
    chain = []
    while node is not None:
        chain.append(node)
        node = node["parent"]
    chain.reverse()
    return chain


def node_transform(node):
    return parse_transform_list(node["attrs"].get("transform"))


def node_fill(node, css):
    f = node["attrs"].get("fill")
    if f and f != "none":
        return f
    style = node["attrs"].get("style")
    if style:
        m = re.search(r"fill:\s*([^;]+)", style)
        if m and m.group(1).strip() != "none":
            return m.group(1).strip()
    for cls in (node["attrs"].get("class") or "").split():
        if cls in css:
            return css[cls]
    return None


# -------------------------------------------------------------- leaf expand

def expand_leaves(root, elem, M, fill_ctx, css, out, depth=0):
    """Recursively expand elem (own transform not yet applied) into leaf
    (abs_d, abs_vertices, fill) tuples appended to out. fill_ctx is the
    inherited fill color from ancestors/use-site so far."""
    if depth > 20:
        return
    tag = elem.tag.split("}")[-1]
    t = parse_transform_list(elem.get("transform"))
    M2 = compose(M, t)
    f = own_fill(elem, css)
    ctx = f if f is not None else fill_ctx

    if tag in ("path", "polygon", "polyline"):
        d, pts = own_vertices(elem)
        if pts:
            abs_pts = [apply_affine(M2, x, y) for x, y in pts]
            out.append((path_d_absolute(d, M2), abs_pts, ctx))
    elif tag == "use":
        href = elem.get(XLINK_HREF) or elem.get("href")
        if href and href.startswith("#"):
            target = root.find(f".//*[@id='{href[1:]}']")
            if target is not None:
                ux = float(elem.get("x", 0))
                uy = float(elem.get("y", 0))
                M3 = compose(M2, (1, 0, 0, 1, ux, uy))
                expand_leaves(root, target, M3, ctx, css, out, depth + 1)
    else:
        for c in elem:
            expand_leaves(root, c, M2, ctx, css, out, depth + 1)


def resolve_part(root, by_start, nodes, text, start, end, css):
    """Resolve one svg_loc [start,end] range into leaf shapes, seeding the
    initial transform/fill from the element's real ancestor chain in the
    document (so ancestor <g transform=.../fill=...> wrapping the snippet is
    correctly accounted for).

    Two svg_loc shapes occur in this corpus: a whole tag (starts with '<',
    e.g. a <use>/<path>/<g> snippet -- resolved via expand_leaves), or a bare
    subpath fragment of a bigger <path>'s "d" attribute (several flags draw
    background art and multiple stars as subpaths of ONE combined path, and
    svg_loc points at just that star's own "M...Z" text) -- resolved by
    finding the enclosing <path> element and applying ITS OWN + ancestor
    transform/fill directly to the fragment's own vertices."""
    snippet = text[start:end + 1]
    out = []

    if not snippet.lstrip().startswith("<"):
        node = find_node_at(by_start, nodes, start)
        if node is None or node["tag"].split("}")[-1] != "path":
            raise RuntimeError(f"bare path fragment at [{start},{end}] has no enclosing <path>")
        M = (1, 0, 0, 1, 0, 0)
        fill_ctx = None
        for anc in ancestor_chain(node):  # includes the enclosing <path> itself
            M = compose(M, node_transform(anc))
            f = node_fill(anc, css)
            if f is not None:
                fill_ctx = f
        pts = path_vertices(snippet)
        abs_pts = [apply_affine(M, x, y) for x, y in pts]
        out.append((path_d_absolute(snippet, M), abs_pts, fill_ctx))
        return out

    node = find_node_at(by_start, nodes, start)
    M = (1, 0, 0, 1, 0, 0)
    fill_ctx = None
    if node is not None:
        for anc in ancestor_chain(node)[:-1]:  # ancestors only, not the node itself
            M = compose(M, node_transform(anc))
            f = node_fill(anc, css)
            if f is not None:
                fill_ctx = f

    frag = ET.fromstring(WRAPPER.format(snippet))
    for e in frag:
        expand_leaves(root, e, M, fill_ctx, css, out)
    return out


# --------------------------------------------------------- star classify

def union_leaves(leaves):
    """Boolean-union every leaf's polygon (many star constructions draw a
    star as several overlapping triangles/kites -- e.g. one arm mirrored and
    rotated N times -- so the individual leaf vertices include seam points
    that never appear on the true outer silhouette; shapely's union removes
    them, leaving only the real outline). Returns a list of shapely Polygons
    (normally exactly one; more if the leaves aren't actually one connected
    shape) plus the dominant fill color."""
    polys = []
    for _, pts, _ in leaves:
        if len(pts) < 3:
            continue
        p = Polygon(pts)
        if not p.is_valid:
            p = p.buffer(0)
        if not p.is_empty and p.area > 1e-9:
            polys.append(p)
    colors = [f for _, _, f in leaves if f]
    color = statistics.mode(colors) if colors else None
    if not polys:
        return [], color
    union = unary_union(polys)
    geoms = list(union.geoms) if union.geom_type == "MultiPolygon" else [union]
    geoms = [g for g in geoms if g.geom_type == "Polygon" and not g.is_empty]
    return geoms, color


def classify_star(leaves):
    """leaves: list of (abs_d, abs_vertices, fill). Returns either
    ("regular", pos, outer, inner, points, rotation, color) or
    ("exact_d", pos, [abs_d, ...], color). The exact_d fallback always
    replays the leaves' OWN transformed "d" strings (curve control points
    intact) rather than anything shapely-derived -- shapely polygons are
    straight-edged only, so a fallback built from them would silently
    flatten any star drawn with curved/rounded tips."""
    all_pts = [p for _, pts, _ in leaves for p in pts]
    colors = [f for _, _, f in leaves if f]
    color = statistics.mode(colors) if colors else None
    fallback_pos = ((statistics.mean(p[0] for p in all_pts), statistics.mean(p[1] for p in all_pts))
                     if all_pts else (0.0, 0.0))
    fallback = ("exact_d", fallback_pos, [d for d, _, _ in leaves], color)

    geoms, _ = union_leaves(leaves)
    if not geoms:
        return None if not all_pts else fallback

    if len(geoms) != 1 or geoms[0].interiors:
        return fallback

    poly = geoms[0]
    # merge near-collinear vertices left over from unioning touching seams
    tol = 1e-5 * max(poly.bounds[2] - poly.bounds[0], poly.bounds[3] - poly.bounds[1], 1.0)
    poly = poly.simplify(tol, preserve_topology=True)
    coords = list(poly.exterior.coords)[:-1]

    cx, cy = poly.centroid.x, poly.centroid.y
    radii = [math.hypot(x - cx, y - cy) for x, y in coords]
    maxr = max(radii) if radii else 0.0
    if maxr == 0 or len(coords) < 6:
        return fallback

    order = sorted(radii)
    gaps = [(order[i + 1] - order[i], i) for i in range(len(order) - 1)]
    gap, gi = max(gaps, key=lambda g: g[0])
    inner_band, outer_band = order[:gi + 1], order[gi + 1:]
    if not inner_band or not outer_band or len(inner_band) != len(outer_band):
        return fallback

    def cv(band):
        m = statistics.mean(band)
        return (statistics.pstdev(band) / m) if m else 1.0

    if cv(inner_band) > 0.02 or cv(outer_band) > 0.02:
        return fallback

    inner = statistics.mean(inner_band)
    outer = statistics.mean(outer_band)
    points = len(outer_band)
    if points < 3 or points > 16:
        return fallback

    # rotation: angle of a real outer vertex vs the idealized i=0 angle (-90deg)
    mid = (inner + outer) / 2
    outer_vertices = [(x, y) for x, y in coords if math.hypot(x - cx, y - cy) > mid]
    if not outer_vertices:
        return fallback
    ox, oy = outer_vertices[0]
    ang = math.degrees(math.atan2(oy - cy, ox - cx))
    step = 360.0 / points
    rotation = (ang + 90.0) % step
    if rotation > step / 2:
        rotation -= step

    return ("regular", (cx, cy), outer, inner, points, round(rotation, 4), color)


def cluster_classes(shapes):
    """shapes: list of (outer, inner, points, rotation). Group by relative
    outer-radius tolerance (+ matching points/rotation), sorted descending
    by outer radius. Returns (classes, index_for_each_shape)."""
    groups = []  # each: {"outers":[], "inners":[], "points":, "rotation":}
    idx_for = []
    for outer, inner, points, rotation in shapes:
        match = None
        for g in groups:
            gm_outer = statistics.mean(g["outers"])
            if (abs(outer - gm_outer) / gm_outer < 0.05 and g["points"] == points
                    and abs(g["rotation"] - rotation) < 3.0):
                match = g
                break
        if match is None:
            match = {"outers": [], "inners": [], "points": points, "rotation": rotation}
            groups.append(match)
        match["outers"].append(outer)
        match["inners"].append(inner)
        idx_for.append(match)

    groups.sort(key=lambda g: statistics.mean(g["outers"]), reverse=True)
    classes = [[round(statistics.mean(g["outers"]), 4), round(statistics.mean(g["inners"]), 4),
                g["points"], g["rotation"]] for g in groups]
    idx_map = {id(g): i for i, g in enumerate(groups)}
    indices = [idx_map[id(g)] for g in idx_for]
    return classes, indices


# --------------------------------------------------------------------- main

def flatten_ranges(loc):
    return loc if isinstance(loc[0], list) else [loc]


def main():
    ref = json.load(open(os.path.join(HERE, "svg_star_ref.json"), encoding="utf-8"))
    out = {}
    n_regular, n_exact = 0, 0

    for key, entry in ref.items():
        path = os.path.join(HERE, "svg", "original", entry["file"])
        text = open(path, encoding="utf-8").read()
        root = ET.fromstring(text)
        by_start, nodes = build_element_index(text)
        css = parse_css_classes(text)

        # resolve every star into either a "regular" (outer,inner,points,rot,color)
        # tuple or an "exact_d" (fragments, color) tuple, keeping BOTH layers
        # (outline + fill) when a star's parts split cleanly into two same-
        # centroid concentric shapes of different sizes/colors (EUREKA-style)
        resolved = []  # (name, kind, pos, payload, color) where payload is
                        # (outer,inner,points,rot) or [exact_d,...]
        for star in entry["stars"]:
            parts = flatten_ranges(star["svg_loc"])
            leaves = []
            for s, e in parts:
                leaves.extend(resolve_part(root, by_start, nodes, text, s, e, css))
            if not leaves:
                raise RuntimeError(f"{key}/{star['name']}: no resolvable geometry")

            # try treating the whole thing as one star first
            whole = classify_star(leaves)

            # also try splitting into concentric same-centroid layers by
            # leaf-shape grouping when there are >=2 leaves with clearly
            # different sizes AND colors (outline+fill duo) -- otherwise a
            # single leaf-shape star built from an arm+mirror+rotate pattern
            # would wrongly get split leaf-by-leaf
            layers = None
            fills = {f for _, _, f in leaves if f}
            if len(fills) >= 2:
                by_fill = {}
                for d, pts, f in leaves:
                    by_fill.setdefault(f, []).append((d, pts, f))
                sub = [classify_star(v) for v in by_fill.values()]
                if all(s and s[0] == "regular" for s in sub) and len(sub) == 2:
                    sub.sort(key=lambda s: s[2], reverse=True)  # by outer radius
                    layers = sub

            if layers is not None:
                resolved.append((star["name"], "layered", layers[0][1], layers, None))
            elif whole[0] == "regular":
                _, pos, outer, inner, points, rot, color = whole
                resolved.append((star["name"], "regular", pos, (outer, inner, points, rot), color))
                n_regular += 1
            else:
                _, pos, frags, color = whole
                resolved.append((star["name"], "exact_d", pos, frags, color))
                n_exact += 1

        # gather primary-layer shapes for classing + color vote
        primary_shapes, primary_colors = [], []
        outline_shapes, outline_colors = [], []
        for name, kind, pos, payload, color in resolved:
            if kind == "regular":
                primary_shapes.append(payload)
                if color:
                    primary_colors.append(color)
            elif kind == "layered":
                outer_layer, inner_layer = payload
                outline_shapes.append(outer_layer[2:6])
                if outer_layer[6]:
                    outline_colors.append(outer_layer[6])
                primary_shapes.append(inner_layer[2:6])
                if inner_layer[6]:
                    primary_colors.append(inner_layer[6])

        classes, class_idx = cluster_classes(primary_shapes) if primary_shapes else ([], [])
        star_color = statistics.mode(primary_colors) if primary_colors else "#ffffff"

        out_entry_stars = []
        pi = 0
        outline_classes, outline_color = None, None
        if outline_shapes:
            outline_classes, oidx = cluster_classes(outline_shapes)
            outline_color = statistics.mode(outline_colors) if outline_colors else star_color
            oi = 0

        for name, kind, pos, payload, color in resolved:
            rec = {"name": name, "pos": [round(pos[0], 4), round(pos[1], 4)]}
            if kind == "regular":
                rec["class"] = class_idx[pi]
                pi += 1
            elif kind == "layered":
                rec["class"] = class_idx[pi]
                pi += 1
                rec["outline_class"] = oidx[oi]
                oi += 1
            else:
                rec["exact_d"] = payload if len(payload) > 1 else payload[0]
            out_entry_stars.append(rec)

        out_entry = {"stars": out_entry_stars, "classes": classes, "star_color": star_color}
        if outline_classes:
            out_entry["outline_classes"] = outline_classes
            out_entry["outline_color"] = outline_color
        out[key] = out_entry

    with open(os.path.join(HERE, "star_geometry.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)

    print(f"{n_regular} regular stars, {n_exact} exact_d stars, across {len(out)} flags")
    print("wrote new/star_geometry.json")


if __name__ == "__main__":
    main()
