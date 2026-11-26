#!/usr/bin/env python3
"""Check course diagrams (.drawio.svg) against the house style.

Usage:
    python3 _scripts/diagrams/check.py [--strict] FILE_OR_DIR...

With a directory, checks every *.drawio.svg below it. Errors (missing editable model,
fonts, HTML labels, colours, line weights, grid, crooked or off-centre arrows) make the
exit status 1. Spacing findings are warnings; --strict makes them fail too.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import house_style as S  # noqa: E402
from drawio import Diagram  # noqa: E402

TOL = 0.6  # px


class Report:
    def __init__(self, path):
        self.path, self.errors, self.warnings = path, [], []

    def error(self, msg):
        self.errors.append(msg)

    def warn(self, msg):
        self.warnings.append(msg)


def role_of(d, v):
    st = d.style(v)
    if st.get("shape") in ("partialRectangle", "singleArrow"):
        return "bracket"  # decorations: brackets and "turns into" arrows
    if "text" in st:
        if st.get("align", "center") == "center":
            # captions are bold; centred plain text is a note (e.g. above an item in a channel)
            return "caption" if int(st.get("fontStyle", 0)) & 1 else "note"
        if st.get("fontFamily") == S.CODE_FONT:
            return "param" if any(e.get("source") == v.get("id") for e in d.edges()) else "code"
        return "note"
    if st.get("container") == "1" and st.get("fillColor", "none") == "none":
        return "group" if st.get("dashed") == "1" else "frame"
    fill = st.get("fillColor", "").lower()
    for role, (f, _, _) in S.ROLES.items():
        if fill == f:
            return role
    return "unknown"


def name(d, v):
    return f"'{d.label(v)}'" if d.label(v) else v.get("id")


def check(path):
    r = Report(path)
    raw = open(path, encoding="utf-8").read()
    try:
        d = Diagram(path)
    except Exception as exc:  # noqa: BLE001
        r.error(f"not an editable draw.io diagram: {exc}")
        return r
    for f in S.FORBIDDEN_FONTS:
        if f in raw:
            r.error(f"uses forbidden font {f}")
    cells = d.cells
    roles = {v.get("id"): role_of(d, v) for v in d.vertices()}
    boxes = {cid: d.box(cells[cid]) for cid in roles}
    palette_fills = {f for f, _, _ in S.ROLES.values()} | {"none"}
    palette_strokes = {s for _, s, _ in S.ROLES.values()} | {S.INK, S.GROUP_STROKE}

    # ---------------------------------------------------------------- cell styles
    for v in d.vertices():
        st, role, n = d.style(v), roles[v.get("id")], name(d, v)
        if v.get("value"):
            if st.get("html", "0") != "0":
                r.error(f"{n}: HTML label (set html=0)")
            if st.get("whiteSpace") == "wrap":
                r.error(f"{n}: wrapped label (whiteSpace=wrap makes an HTML label; use explicit line breaks)")
            if st.get("fontFamily", S.LABEL_FONT) not in S.ALLOWED_FONTS:
                r.error(f"{n}: font {st.get('fontFamily')} (use {S.LABEL_FONT} or {S.CODE_FONT})")
        if role == "bracket":
            continue
        if role == "unknown":
            r.error(f"{n}: fill {st.get('fillColor')} is not a palette colour")
            continue
        if "text" in st:
            continue
        if st.get("strokeColor", S.INK).lower() not in palette_strokes:
            r.error(f"{n}: stroke {st.get('strokeColor')} is not a palette colour")
        want = S.W_FRAME if role in ("frame", "group") else S.W_BOX
        if float(st.get("strokeWidth", 1)) != want:
            r.error(f"{n}: outline width {st.get('strokeWidth', 1)} (want {want:g})")
        x, y, _, _ = boxes[v.get("id")]
        w, h = boxes[v.get("id")][2:]
        on = lambda a: abs(a / S.GRID - round(a / S.GRID)) < 0.01
        # left/top edge on the grid, or the centre (boxes centred under/beside another box)
        inside_channel = any(roles.get(a) in ("channel", "item") for a in _ancestors(cells, v.get("id")))
        if role != "item" and not inside_channel and not ((on(x) or on(x + w / 2)) and (on(y) or on(y + h / 2))):
            r.error(f"{n}: position ({x:g}, {y:g}) is off the {S.GRID} px grid (edge or centre)")
        in_channel = False
        p = v.get("parent")
        while p in roles:
            in_channel |= roles[p] in ("channel", "item")
            p = cells[p].get("parent")
        has_kids = any(c.get("parent") == v.get("id") for c in d.vertices())
        if role in ("process", "file") and not in_channel and not has_kids and boxes[v.get("id")][3] != S.BOX_H:
            r.warn(f"{n}: height {boxes[v.get('id')][3]:g} (boxes are {S.BOX_H} high)")
        rows = (boxes[v.get("id")][3] + S.ITEM_GAP) / (S.ITEM_H + S.ITEM_GAP)  # an item may span whole rows
        if role in ("item", "file") and in_channel and not has_kids and abs(rows - round(rows)) > 0.01:
            r.warn(f"{n}: height {boxes[v.get('id')][3]:g} (items in a channel are {S.ITEM_H} high)")
    for v in d.vertices():
        if roles[v.get("id")] == "param":
            want = S.param_label_width(v.get("value"), float(d.style(v).get("fontSize", S.SIZE_PARAM)))
            if abs(boxes[v.get("id")][2] - want) > 1:
                r.warn(f"{name(d, v)}: label width {boxes[v.get('id')][2]:g}, want {want:g} "
                       f"(arrow starts one monospace space after the text)")

    # ---------------------------------------------------------------- edges
    def is_box(cid):
        return roles.get(cid) in ("process", "file", "channel", "item")

    fanout = {}
    for e in d.edges():
        st = d.style(e)
        s, t = e.get("source"), e.get("target")
        label = f"{name(d, cells[s]) if s in cells else '(free)'} -> {name(d, cells[t]) if t in cells else '(free)'}"
        param = st.get("dashed") == "1"
        want = S.W_PARAM_ARROW if param else S.W_DATA_ARROW
        if float(st.get("strokeWidth", 1)) != want:
            r.error(f"{label}: line width {st.get('strokeWidth', 1)} (want {want:g})")
        if st.get("strokeColor", S.INK).lower() != S.INK:
            r.error(f"{label}: line colour {st.get('strokeColor')} (want {S.INK})")
        if param and st.get("dashPattern") != S.PARAM_DASH:
            r.error(f"{label}: dash pattern {st.get('dashPattern')} (want {S.PARAM_DASH})")
        if float(st.get("endSize", 6)) != S.ARROW_SIZE:
            r.error(f"{label}: arrowhead size {st.get('endSize', 6)} (want {S.ARROW_SIZE})")
        pts = d.endpoints(e)
        if None in pts:
            r.error(f"{label}: dangling end")
            continue
        straight = len(pts) == 2
        (x1, y1), (x2, y2) = pts[0], pts[-1]
        if straight and abs(y2 - y1) > 2 * abs(x2 - x1) and abs(y2 - y1) > 20:
            if abs(x2 - x1) > TOL:
                r.error(f"{label}: vertical arrow is off by {x2 - x1:+.1f} px (move the downstream box, not the anchor)")
            if is_box(s) and is_box(t) and (st.get("exitX") != "0.5" or st.get("entryX") != "0.5"):
                r.error(f"{label}: vertical arrow must leave and land at the centre (exitX/entryX 0.5)")
            # vertical gaps
            sb, tb = boxes.get(s), boxes.get(t)
            if sb and tb:
                frames = {cid for cid, ro in roles.items() if ro == "frame"}
                fs = {f for f in frames if _inside(boxes, s, f)}
                ft = {f for f in frames if _inside(boxes, t, f)}
                if s in frames or t in frames or fs != ft:
                    if s in frames:
                        gap = tb[1] - (sb[1] + sb[3])
                    elif t in frames:
                        gap = tb[1] - (sb[1] + sb[3])
                    else:
                        crossed = (ft - fs) or (fs - ft)
                        f = boxes[next(iter(crossed))]
                        gap = f[1] - (sb[1] + sb[3]) if ft - fs else tb[1] - (f[1] + f[3])
                    if abs(gap - S.GAP_FRAME) > TOL and abs(gap - S.GAP_FRAME_PARAM) > TOL:
                        r.warn(f"{label}: gap {gap:g} into/out of a frame (want {S.GAP_FRAME}, or "
                               f"{S.GAP_FRAME_PARAM} when a parameter label sits there)")
                else:
                    gap = tb[1] - (sb[1] + sb[3])
                    if abs(gap - S.GAP_V) > TOL:
                        r.warn(f"{label}: vertical gap {gap:g} (want {S.GAP_V})")
        if straight and st.get("exitX") == "1" and st.get("entryX") == "0" and is_box(s) and is_box(t):
            sb, tb = boxes[s], boxes[t]
            gap = tb[0] - (sb[0] + sb[2])
            want_h = S.GAP_MAP if "course-map" in os.path.basename(path) else S.GAP_H
            if abs(gap - want_h) > TOL:
                r.warn(f"{label}: horizontal gap {gap:g} (want {want_h})")
            if abs((sb[1] + sb[3] / 2) - (tb[1] + tb[3] / 2)) < TOL and (st.get("exitY") != "0.5" or st.get("entryY") != "0.5"):
                r.error(f"{label}: level boxes, so the arrow must leave and land at the centre (exitY/entryY 0.5)")
            fanout.setdefault(s, []).append(t)
        if param and straight and abs(y2 - y1) > TOL:
            r.error(f"{label}: parameter arrow is not level with its label")

    # fan-outs: stacked STACK_GAP apart, stack centred on the source
    for s, ts in fanout.items():
        if len(ts) < 2:
            continue
        col = sorted((boxes[t] for t in ts), key=lambda b: b[1])
        if len({b[0] for b in col}) != 1:
            continue
        gaps = [b[1] - (a[1] + a[3]) for a, b in zip(col, col[1:])]
        if any(abs(g - S.STACK_GAP) > TOL for g in gaps):
            r.warn(f"fan-out from {name(d, cells[s])}: stacked gaps {gaps} (want {S.STACK_GAP})")
        mid = (col[0][1] + col[-1][1] + col[-1][3]) / 2
        sb = boxes[s]
        if abs(mid - (sb[1] + sb[3] / 2)) > TOL:
            r.warn(f"fan-out from {name(d, cells[s])}: stack centre {mid:g} vs source centre {sb[1] + sb[3] / 2:g}")

    # notes sit NOTE_GAP right of their box, centred on it
    for v in d.vertices():
        if roles[v.get("id")] != "note":
            continue
        x, y, w, h = boxes[v.get("id")]
        near = [cid for cid, b in boxes.items() if is_box(cid) and 0 <= x - (b[0] + b[2]) <= 60
                and not any(c.get("parent") == cid for c in d.vertices())
                and not any(roles.get(a) in ("channel", "item", "file", "process") for a in _ancestors(cells, cid))
                and b[1] < y + h / 2 < b[1] + b[3]]
        for cid in near:
            b = boxes[cid]
            if abs(x - (b[0] + b[2] + S.NOTE_GAP)) > TOL or abs((y + h / 2) - (b[1] + b[3] / 2)) > TOL:
                r.warn(f"note {name(d, v)}: should sit {S.NOTE_GAP} px right of {name(d, cells[cid])}, centred on it")

    # captions: CAPTION_GAP under the box/channel above them (CAPTION_GAP_GROUP under a dashed group);
    # course maps stack part titles over subtitles and brackets, so they're exempt
    for v in ([] if "course-map" in os.path.basename(path) else d.vertices()):
        if roles[v.get("id")] != "caption":
            continue
        x, y, w, h = boxes[v.get("id")]
        above = [cid for cid, b in boxes.items() if roles[cid] in ("process", "file", "channel", "group", "frame", "item")
                 and b[1] + b[3] <= y + TOL and b[0] < x + w / 2 < b[0] + b[2]]
        if not above:
            # caption above its shape (top of a captioned stack): same gap, mirrored
            below = [cid for cid, b in boxes.items() if roles[cid] in ("process", "file", "channel", "group", "frame", "item")
                     and b[1] >= y + h - TOL and b[0] < x + w / 2 < b[0] + b[2]]
            if below:
                cid = min(below, key=lambda c: boxes[c][1])
                gap = boxes[cid][1] - (y + h)
                if abs(gap - S.CAPTION_GAP) > TOL:
                    r.warn(f"caption {name(d, v)}: {gap:g} px above {name(d, cells[cid])} (want {S.CAPTION_GAP})")
            continue
        cid = max(above, key=lambda c: boxes[c][1] + boxes[c][3])
        gap = y - (boxes[cid][1] + boxes[cid][3])
        want = S.CAPTION_GAP_GROUP if roles[cid] == "group" else S.CAPTION_GAP
        if roles[cid] != "group" and abs(gap - S.CAPTION_GAP_PROCESS) <= TOL:
            continue  # single-process view
        if abs(gap - want) > TOL:
            r.warn(f"caption {name(d, v)}: {gap:g} px under {name(d, cells[cid])} (want {want})")

    # text must not run into other shapes (text width measured from the font, not the cell)
    def ancestors(cid):
        out = set()
        while cid in cells and cells[cid].get("vertex") == "1":
            cid = cells[cid].get("parent")
            out.add(cid)
        return out

    for v in d.vertices():
        st = d.style(v)
        if not v.get("value") or "text" not in st:
            continue
        x, y, w, h = boxes[v.get("id")]
        fs = float(st.get("fontSize", 12))
        lines = v.get("value").split("\n")
        per = S.MENLO_ADVANCE if st.get("fontFamily") == S.CODE_FONT else 0.55
        tw = max(len(line) for line in lines) * fs * per
        th = len(lines) * fs * 1.2
        if st.get("align") == "left":
            x0 = x + float(st.get("spacingLeft", S.TEXT_SPACING_LEFT))
        else:
            x0 = x + (w - tw) / 2
        y0 = y + (h - th) / 2
        for cid, b in boxes.items():
            if cid == v.get("id") or cid in ancestors(v.get("id")) or "text" in d.style(cells[cid]):
                continue
            if x0 < b[0] + b[2] and b[0] < x0 + tw and y0 < b[1] + b[3] and b[1] < y0 + th:
                r.error(f"text {name(d, v)} runs into {name(d, cells[cid])}")

    # frames: FRAME_PAD on every side
    for fid, ro in roles.items():
        if ro != "frame":
            continue
        fx, fy, fw, fh = boxes[fid]
        kids = [boxes[c] for c in roles if cells[c].get("parent") == fid]
        if not kids:
            continue
        pads = (min(k[0] for k in kids) - fx, min(k[1] for k in kids) - fy,
                fx + fw - max(k[0] + k[2] for k in kids), fy + fh - max(k[1] + k[3] for k in kids))
        if any(abs(p - S.FRAME_PAD) > TOL for p in pads):
            r.warn(f"frame padding left/top/right/bottom {[round(p) for p in pads]} (want {S.FRAME_PAD})")
    return r


def _ancestors(cells, cid):
    out = []
    while cid in cells and cells[cid].get("vertex") == "1":
        cid = cells[cid].get("parent")
        out.append(cid)
    return out


def _inside(boxes, cid, fid):
    a, f = boxes[cid], boxes[fid]
    return cid != fid and f[0] <= a[0] and a[0] + a[2] <= f[0] + f[2] and f[1] <= a[1] and a[1] + a[3] <= f[1] + f[3]


def main(argv):
    strict = "--strict" in argv
    targets = [a for a in argv if a != "--strict"]
    files = []
    for t in targets:
        if os.path.isdir(t):
            for root, _, names in os.walk(t):
                files += [os.path.join(root, n) for n in sorted(names) if n.endswith(".drawio.svg")]
        else:
            files.append(t)
    failed = False
    for f in files:
        r = check(f)
        status = "ok" if not r.errors and not r.warnings else ("FAIL" if r.errors else "warn")
        print(f"{status:4}  {f}")
        for m in r.errors:
            print(f"      error: {m}")
        for m in r.warnings:
            print(f"      warn:  {m}")
        failed |= bool(r.errors) or (strict and bool(r.warnings))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
