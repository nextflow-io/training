#!/usr/bin/env python3
"""Convert an Excalidraw SVG export into an editable .drawio.svg in the house style.

Usage:
    python3 _scripts/diagrams/excalidraw_to_drawio.py INPUT.svg OUTPUT.drawio.svg

This is a starting point, not a finished diagram. It recovers boxes, text and arrows
from the exported SVG (the Excalidraw scene is usually not embedded), maps the old fills
to roles (process, file, channel, item), snaps to the grid, attaches arrows to boxes and
applies the house colours, fonts and line weights. Then:

1. open the result in VS Code (draw.io extension) and fix the layout by hand;
2. run check.py and fix what it reports;
3. save once in the editor, so draw.io redraws the image.
"""
import math
import os
import re
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import house_style as S  # noqa: E402
from drawio import Diagram  # noqa: E402

NS = "{http://www.w3.org/2000/svg}"
ATTACH = 30  # px: an arrow end this close to a box attaches to it


def snap(v):
    return round(v / S.GRID) * S.GRID


# ---------------------------------------------------------------- read the Excalidraw export
def parse(path):
    root = ET.parse(path).getroot()
    items = []
    for g in [g for g in root if g.tag == NS + "g"]:
        t = g.get("transform", "")
        if "translate" not in t:
            for c in g.iter(NS + "g"):
                if "translate" in c.get("transform", ""):
                    t = c.get("transform")
                    break
        m = re.search(r"translate\(([-\d.e]+) ([-\d.e]+)\)", t)
        r = re.search(r"rotate\(([-\d.e]+) ([-\d.e]+) ([-\d.e]+)\)", t)
        if not m or not r:
            continue
        X, Y, CX, CY = float(m.group(1)), float(m.group(2)), float(r.group(2)), float(r.group(3))
        texts = g.findall(".//" + NS + "text")
        if texts:
            items.append(dict(kind="text", x=X, y=Y, w=2 * CX, h=2 * CY, lines=[x.text or "" for x in texts],
                              size=float(texts[0].get("font-size", "20").rstrip("px")),
                              family=texts[0].get("font-family", "")))
            continue
        paths = g.findall(".//" + NS + "path")
        fills = [p.get("fill") for p in paths if p.get("fill") not in (None, "none")]
        dash = [p.get("stroke-dasharray") for p in paths if p.get("stroke-dasharray")]
        items.append(dict(kind="shape", x=X, y=Y, w=2 * CX, h=2 * CY, fill=fills[0].lower() if fills else None,
                          dash=dash[0] if dash else None, d=paths[0].get("d") if paths else "", tr=t,
                          flipped="scale(1 -1)" in ET.tostring(g, encoding="unicode")))
    return items


def extract(items):
    boxes, texts, arrows, prev = [], [], [], None
    for it in items:
        if it["kind"] == "text":
            texts.append(it)
        elif it["fill"]:
            boxes.append(dict(kind=S.LEGACY_FILLS.get(it["fill"], "file"), x=it["x"], y=it["y"], w=it["w"], h=it["h"]))
        elif it["d"].startswith("M32 0") or it["dash"]:
            boxes.append(dict(kind="group" if it["dash"] else "frame", x=it["x"], y=it["y"], w=it["w"], h=it["h"]))
        else:
            if it["tr"] == prev:  # arrowhead strokes come in their own groups
                continue
            prev = it["tr"]
            parts = re.split(r"(?=M)", it["d"])
            if len(parts) < 2:
                continue
            n = [float(x) for x in re.findall(r"-?\d+(?:\.\d+)?", parts[1])]
            sx, sy, ex, ey = n[0], n[1], n[-2], n[-1]
            if it["flipped"]:
                sy, ey = -sy, -ey
            arrows.append(dict(s=(it["x"] + sx, it["y"] + sy), t=(it["x"] + ex, it["y"] + ey),
                               bent=parts[1].count("C") > 1))
    return boxes, texts, arrows


# ---------------------------------------------------------------- build the model
def convert(inp, outp):
    boxes, texts, arrows = extract(parse(inp))
    for b in boxes:
        b.update(X=snap(b["x"]), Y=snap(b["y"]), W=snap(b["w"]), H=snap(b["h"]), label="", fs=None, parent=None)
        if b["kind"] in ("process", "file") and abs(b["H"] - S.BOX_H) <= 30:
            b["H"] = S.BOX_H
    # boxes inside a channel are items
    chans = [b for b in boxes if b["kind"] == "channel"]
    for b in boxes:
        for c in chans:
            if b is not c and _contains(c, b):
                b["parent"] = c
                if b["kind"] == "process":
                    b["kind"] = "item"
                b["H"] = S.ITEM_H
    # text inside a box is its label; the rest are notes, captions, parameters or code
    free = []
    for t in texts:
        cx, cy = t["x"] + t["w"] / 2, t["y"] + t["h"] / 2
        host = [b for b in boxes if b["kind"] in S.ROLES and b["x"] <= cx <= b["x"] + b["w"] and b["y"] <= cy <= b["y"] + b["h"]
                and not any(o is not b and o.get("parent") is b and _point_in(o, cx, cy) for o in boxes)]
        if host:
            host[-1]["label"] = "\n".join(t["lines"])
        else:
            free.append(t)
    tcells = []
    for t in free:
        text = "\n".join(t["lines"]).strip()
        if not text:
            continue
        if text.startswith("--"):
            role, fs = "param", S.SIZE_PARAM
        elif re.search(r"\w\(.*\)|\{|->", text) and "Cascadia" in t["family"] or text.startswith("channel."):
            role, fs = "code", S.SIZE_CODE
        elif len(t["lines"]) == 1 and t["size"] >= 24 and not _beside_box(t, boxes):
            role, fs = "caption", S.SIZE_CAPTION
        else:
            role, fs = "note", S.SIZE_NOTE
        n = len(text.split("\n"))
        W = S.param_label_width(text, fs) if role in ("param", "code") else \
            math.ceil((max(len(x) for x in text.split("\n")) * fs * 0.55 + 8) / S.GRID) * S.GRID
        H = round(n * fs * 1.2 + 10)
        X = snap(t["x"])
        Y = t["y"] + t["h"] / 2 - H / 2
        if role == "caption":
            X, W = t["x"], t["w"]
        if role == "note":
            near = _beside_box(t, boxes)
            if near:
                X, Y = near["X"] + near["W"] + S.NOTE_GAP, near["Y"] + near["H"] / 2 - H / 2
        tcells.append(dict(kind=role, X=X, Y=Y, W=W, H=H, label=text, parent=None))
    # containment in frames and groups (innermost wins)
    outers = sorted([b for b in boxes if b["kind"] in ("frame", "group")], key=lambda b: b["W"] * b["H"], reverse=True)
    for c in boxes + tcells:
        if c.get("parent"):
            continue
        for f in outers:
            if f is not c and _contains(f, c):
                c["parent"] = f

    d = Diagram()
    order = sorted(boxes, key=lambda b: _depth(b)) + tcells
    for i, c in enumerate(order):
        c["id"] = f"{c['kind']}{i}"
    for c in order:
        d.add_vertex(c["id"], c["label"], S.vertex_style(c["kind"]), c["X"], c["Y"], c["W"], c["H"],
                     parent=c["parent"]["id"] if c.get("parent") else "1")

    # arrows: attach each end to the nearest box (or a frame border), straighten near-axis arrows
    def rect_dist(p, c):
        dx = max(c["X"] - p[0], 0, p[0] - (c["X"] + c["W"]))
        dy = max(c["Y"] - p[1], 0, p[1] - (c["Y"] + c["H"]))
        return math.hypot(dx, dy)

    def pick(p, other):
        best = None
        for c in order:
            if c["kind"] in ("frame", "group"):
                inside = c["X"] < other[0] < c["X"] + c["W"] and c["Y"] < other[1] < c["Y"] + c["H"]
                edge = min(abs(p[0] - c["X"]), abs(p[0] - c["X"] - c["W"]), abs(p[1] - c["Y"]), abs(p[1] - c["Y"] - c["H"]))
                if inside or edge > ATTACH:
                    continue
                dist = edge + 5
            else:
                dist = rect_dist(p, c)
                if dist > ATTACH + (20 if c["kind"] in ("param", "note") else 0):
                    continue
            if best is None or dist < best[0]:
                best = (dist, c)
        return best[1] if best else None

    for i, a in enumerate(arrows):
        src, tgt = pick(a["s"], a["t"]), pick(a["t"], a["s"])
        if not src or not tgt:
            continue
        kind = "param" if src["kind"] == "param" else "data"
        sp, tp = _project(a["s"], src), _project(a["t"], tgt)
        if not a["bent"]:
            if abs(sp[1] - tp[1]) < 15 and sp[2] in "lr" and tp[2] in "lr":
                level = abs((src["Y"] + src["H"] / 2) - (tgt["Y"] + tgt["H"] / 2)) < 15
                y = (tgt["Y"] + tgt["H"] / 2) if level else (sp[1] + tp[1]) / 2
                sp, tp = (sp[0], y, sp[2]), (tp[0], y, tp[2])
            if abs(sp[0] - tp[0]) < 15 and sp[2] in "tb" and tp[2] in "tb":
                x = tgt["X"] + tgt["W"] / 2 if src["kind"] in ("frame", "group") else src["X"] + src["W"] / 2
                sp, tp = (x, sp[1], sp[2]), (x, tp[1], tp[2])
        if kind == "param":
            sp = (src["X"] + src["W"], src["Y"] + src["H"] / 2, "r")
        rel = lambda p, c: ((p[0] - c["X"]) / c["W"], (p[1] - c["Y"]) / c["H"])
        d.add_edge(f"edge{i}", src["id"], tgt["id"], S.edge_style(kind, orthogonal=a["bent"]),
                   exit=rel(sp, src), entry=rel(tp, tgt))
    d.save(outp)
    print(f"wrote {outp}: {len(order)} shapes, {len(d.edges())} arrows ({len(arrows) - len(d.edges())} not attached)")


def _contains(o, c):
    return o["X"] <= c["X"] and c["X"] + c["W"] <= o["X"] + o["W"] and o["Y"] <= c["Y"] and c["Y"] + c["H"] <= o["Y"] + o["H"]


def _point_in(b, x, y):
    return b["x"] <= x <= b["x"] + b["w"] and b["y"] <= y <= b["y"] + b["h"]


def _depth(b):
    return 0 if not b.get("parent") else 1 + _depth(b["parent"])


def _beside_box(t, boxes):
    cy = t["y"] + t["h"] / 2
    for b in boxes:
        if b["kind"] in S.ROLES and b["y"] <= cy <= b["y"] + b["h"] and 0 < t["x"] - (b["x"] + b["w"]) < 80:
            return b
    return None


def _project(p, c):
    """Nearest point on the border of c, with the side as l/r/t/b."""
    x, y, w, h = c["X"], c["Y"], c["W"], c["H"]
    px, py = min(max(p[0], x), x + w), min(max(p[1], y), y + h)
    dist = {"l": px - x, "r": x + w - px, "t": py - y, "b": y + h - py}
    side = min(dist, key=dist.get)
    px = x if side == "l" else x + w if side == "r" else px
    py = y if side == "t" else y + h if side == "b" else py
    return px, py, side


if __name__ == "__main__":
    convert(sys.argv[1], sys.argv[2])
