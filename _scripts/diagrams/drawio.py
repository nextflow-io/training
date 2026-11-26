"""Read, edit, render and write .drawio.svg files.

A .drawio.svg is a normal SVG image whose root `content` attribute holds the draw.io model.
The draw.io editor (VS Code extension "hediet.vscode-drawio", or diagrams.net) redraws the
image when you save. `Diagram.save()` also redraws it, with a simple built-in renderer, so a
file edited by a script can be viewed and committed without opening the editor.
"""
import base64
import html
import math
import re
import urllib.parse
import xml.etree.ElementTree as ET
import zlib

from house_style import INK, LABEL_FONT, TEXT_SPACING_LEFT


def parse_style(s):
    out = {}
    for part in (s or "").split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
        elif part:
            out[part] = True  # bare keywords such as "text" or "rounded"
    return out


def format_style(d):
    return "".join(f"{k};" if v is True else f"{k}={v};" for k, v in d.items())


def num(v):
    v = round(float(v), 2)
    return f"{v:g}"


class Diagram:
    def __init__(self, path=None):
        self.path = path
        if path:
            raw = open(path, encoding="utf-8").read()
            content = ET.fromstring(raw).get("content")
            if content is None:
                raise ValueError(f"{path}: no embedded draw.io model (content attribute)")
            self.mxfile = ET.fromstring(content)
            d = self.mxfile.find("diagram")
            model = d.find("mxGraphModel")
            if model is None:
                model = ET.fromstring(urllib.parse.unquote(zlib.decompress(base64.b64decode(d.text), -15).decode()))
                d.text = None
                d.append(model)
            self.mxfile.attrib.pop("compressed", None)
            self.model = model
        else:
            self.mxfile = ET.Element("mxfile", host="drawio")
            d = ET.SubElement(self.mxfile, "diagram", id="page-1", name="Page-1")
            self.model = ET.SubElement(d, "mxGraphModel", grid="1", gridSize="10", guides="1", tooltips="1",
                                       connect="1", arrows="1", fold="1", page="0", pageScale="1",
                                       pageWidth="1600", pageHeight="1000", math="0", shadow="0")
            root = ET.SubElement(self.model, "root")
            ET.SubElement(root, "mxCell", id="0")
            ET.SubElement(root, "mxCell", id="1", parent="0")

    # ------------------------------------------------------------ access
    @property
    def root(self):
        return self.model.find("root")

    @property
    def cells(self):
        return {c.get("id"): c for c in self.model.iter("mxCell")}

    def vertices(self):
        return [c for c in self.model.iter("mxCell") if c.get("vertex") == "1"]

    def edges(self):
        return [c for c in self.model.iter("mxCell") if c.get("edge") == "1"]

    def style(self, cell):
        return parse_style(cell.get("style"))

    def set_style(self, cell, **kv):
        st = self.style(cell)
        st.update({k: v for k, v in kv.items() if v is not None})
        for k, v in kv.items():
            if v is None:
                st.pop(k, None)
        cell.set("style", format_style(st))

    def geo(self, cell):
        return cell.find("mxGeometry")

    def box(self, cell):
        """Absolute (x, y, w, h) of a vertex."""
        cells = self.cells
        g = self.geo(cell)
        w, h = float(g.get("width", 0)), float(g.get("height", 0))
        x = y = 0.0
        c = cell
        while c is not None and c.get("vertex") == "1":
            gc = self.geo(c)
            x += float(gc.get("x", 0))
            y += float(gc.get("y", 0))
            c = cells.get(c.get("parent"))
        return x, y, w, h

    def origin(self, cell_id):
        c = self.cells.get(cell_id)
        if c is None or c.get("vertex") != "1":
            return 0.0, 0.0
        x, y, _, _ = self.box(c)
        return x, y

    def is_text(self, cell):
        return "text" in self.style(cell)

    def label(self, cell):
        return (cell.get("value") or "").replace("\n", " ")

    def endpoints(self, edge):
        """Absolute points of an edge: start, waypoints..., end."""
        st = self.style(edge)
        cells = self.cells
        g = self.geo(edge)
        ox, oy = self.origin(edge.get("parent"))
        pts = []
        arr = g.find("Array") if g is not None else None
        way = [(ox + float(p.get("x", 0)), oy + float(p.get("y", 0))) for p in arr] if arr is not None else []

        def end(which):
            cid = edge.get(which)
            if cid in cells:
                x, y, w, h = self.box(cells[cid])
                k = "exit" if which == "source" else "entry"
                if f"{k}X" in st:
                    return x + float(st[f"{k}X"]) * w, y + float(st[f"{k}Y"]) * h
                return None  # floating: resolved below
            p = g.find(f"mxPoint[@as='{which}Point']") if g is not None else None
            return (ox + float(p.get("x", 0)), oy + float(p.get("y", 0))) if p is not None else None

        s, t = end("source"), end("target")

        def floating(cid, toward):
            x, y, w, h = self.box(cells[cid])
            cx, cy = x + w / 2, y + h / 2
            dx, dy = toward[0] - cx, toward[1] - cy
            if abs(dx) * h > abs(dy) * w:
                return (x + w if dx > 0 else x), cy
            return cx, (y + h if dy > 0 else y)

        if s is None and edge.get("source") in cells:
            s = floating(edge.get("source"), way[0] if way else (t or self._centre(edge.get("target"))))
        if t is None and edge.get("target") in cells:
            t = floating(edge.get("target"), way[-1] if way else s)
        return [s] + way + [t]

    def _centre(self, cid):
        x, y, w, h = self.box(self.cells[cid])
        return x + w / 2, y + h / 2

    # ------------------------------------------------------------ building
    def add_vertex(self, cid, value, style, x, y, w, h, parent="1"):
        """Add a vertex at ABSOLUTE coordinates; converted to parent-relative."""
        px, py = self.origin(parent)
        c = ET.SubElement(self.root, "mxCell", id=cid, value=value, style=style, vertex="1", parent=parent)
        ET.SubElement(c, "mxGeometry", x=num(x - px), y=num(y - py), width=num(w), height=num(h), **{"as": "geometry"})
        return c

    def add_edge(self, cid, source, target, style, exit=None, entry=None, points=(), parent="1"):
        st = parse_style(style)
        for k, v in (("exit", exit), ("entry", entry)):
            if v is not None:
                st.update({f"{k}X": f"{v[0]:.5g}", f"{k}Y": f"{v[1]:.5g}", f"{k}Dx": 0, f"{k}Dy": 0})
        c = ET.SubElement(self.root, "mxCell", id=cid, style=format_style(st), edge="1", parent=parent,
                          source=source, target=target)
        g = ET.SubElement(c, "mxGeometry", relative="1", **{"as": "geometry"})
        if points:
            px, py = self.origin(parent)
            arr = ET.SubElement(g, "Array", **{"as": "points"})
            for x, y in points:
                ET.SubElement(arr, "mxPoint", x=num(x - px), y=num(y - py))
        return c

    # ------------------------------------------------------------ rendering
    def render(self, pad=10):
        vs = self.vertices()
        boxes = [self.box(v) for v in vs]
        pts = [p for e in self.edges() for p in self.endpoints(e) if p]
        minx = min([b[0] for b in boxes] + [p[0] for p in pts]) - pad
        miny = min([b[1] for b in boxes] + [p[1] for p in pts]) - pad
        maxx = max([b[0] + b[2] for b in boxes] + [p[0] for p in pts]) + pad
        maxy = max([b[1] + b[3] for b in boxes] + [p[1] for p in pts]) + pad
        W, H = maxx - minx, maxy - miny
        o = [f'<rect x="{num(minx)}" y="{num(miny)}" width="{num(W)}" height="{num(H)}" fill="#ffffff"/>']
        # containers first so their contents draw on top
        depth = lambda c: 0 if c.get("parent") == "1" else 1 + depth(self.cells[c.get("parent")])
        for v in sorted(vs, key=depth):
            st = self.style(v)
            x, y, w, h = self.box(v)
            op = f' opacity="{float(st["opacity"]) / 100:g}"' if "opacity" in st else ""
            if op:
                o.append(f"<g{op}>")
            if st.get("shape") == "partialRectangle":
                sides = {"top": ((x, y), (x + w, y)), "bottom": ((x, y + h), (x + w, y + h)),
                         "left": ((x, y), (x, y + h)), "right": ((x + w, y), (x + w, y + h))}
                for side, ((ax, ay), (bx, by)) in sides.items():
                    if st.get(side, "1") != "0":
                        o.append(f'<line x1="{num(ax)}" y1="{num(ay)}" x2="{num(bx)}" y2="{num(by)}" '
                                 f'stroke="{st.get("strokeColor", INK)}" stroke-width="{st.get("strokeWidth", 1)}"/>')
            elif st.get("shape") == "singleArrow":
                aw, asz = float(st.get("arrowWidth", 0.3)), float(st.get("arrowSize", 0.2))
                cy, half, hx = y + h / 2, aw * h / 2, x + w - asz * w
                pts = [(x, cy - half), (hx, cy - half), (hx, y), (x + w, cy), (hx, y + h), (hx, cy + half), (x, cy + half)]
                o.append('<path d="M' + " L".join(f"{num(a)} {num(b)}" for a, b in pts) + f' Z" fill="{st.get("fillColor", "none")}" '
                         f'stroke="{st.get("strokeColor", INK)}" stroke-width="{st.get("strokeWidth", 1)}" stroke-linejoin="round"/>')
            elif "text" not in st:
                fill = st.get("fillColor", "none")
                stroke = st.get("strokeColor", INK)
                r = float(st.get("arcSize", 16)) / 2 if st.get("absoluteArcSize") == "1" else min(w, h) * 0.1
                r = min(r, h / 2, w / 2)
                dash = f' stroke-dasharray="{st["dashPattern"]}"' if st.get("dashed") == "1" else ""
                o.append(f'<rect x="{num(x)}" y="{num(y)}" width="{num(w)}" height="{num(h)}" rx="{num(r)}" '
                         f'fill="{fill}" stroke="{stroke}" stroke-width="{st.get("strokeWidth", 1)}"{dash}/>')
            lines = (v.get("value") or "").split("\n")
            if not lines[0]:
                if op:
                    o.append("</g>")
                continue
            fs = float(st.get("fontSize", 12))
            fam = st.get("fontFamily", LABEL_FONT)
            fam_css = "Menlo, Consolas, monospace" if fam == "Menlo" else f"{fam}, Arial, sans-serif"
            lh = fs * 1.2
            y0 = y + h / 2 - lh * (len(lines) - 1) / 2
            if st.get("verticalAlign") == "top":  # e.g. a file name as the header of a container
                y0 = y + float(st.get("spacingTop", 10)) + lh / 2
            if st.get("align") == "left":
                anchor, tx = "start", x + float(st.get("spacingLeft", TEXT_SPACING_LEFT))
            else:
                anchor, tx = "middle", x + w / 2
            for i, line in enumerate(lines):
                o.append(f'<text x="{num(tx)}" y="{num(y0 + i * lh)}" font-family="{fam_css}" font-size="{num(fs)}px" '
                         f'fill="{st.get("fontColor", INK)}" text-anchor="{anchor}" dominant-baseline="central" '
                         f'{"font-weight=\"bold\" " if int(st.get("fontStyle", 0)) & 1 else ""}'
                         f'xml:space="preserve">{html.escape(line)}</text>')
            if op:
                o.append("</g>")
        for e in self.edges():
            st = self.style(e)
            p = self.endpoints(e)
            if None in p:
                continue
            sw = float(st.get("strokeWidth", 1))
            L = float(st.get("endSize", 6)) * 1.1 + sw
            (ax, ay), (bx, by) = p[-2], p[-1]
            dl = math.hypot(bx - ax, by - ay) or 1
            ux, uy = (bx - ax) / dl, (by - ay) / dl
            filled = st.get("endFill", "1") == "1"
            back = L * 0.75 if filled else sw / 2
            q = p[:-1] + [(bx - ux * back, by - uy * back)]
            d = f"M{num(q[0][0])} {num(q[0][1])}"
            r = 10 if st.get("rounded") == "1" else 0
            for i in range(1, len(q) - 1):
                (x0, y0), (x1, y1), (x2, y2) = q[i - 1], q[i], q[i + 1]
                l1, l2 = math.hypot(x1 - x0, y1 - y0) or 1, math.hypot(x2 - x1, y2 - y1) or 1
                rr = min(r, l1 / 2, l2 / 2)
                d += f" L{num(x1 - (x1 - x0) / l1 * rr)} {num(y1 - (y1 - y0) / l1 * rr)}"
                d += f" Q{num(x1)} {num(y1)} {num(x1 + (x2 - x1) / l2 * rr)} {num(y1 + (y2 - y1) / l2 * rr)}"
            d += f" L{num(q[-1][0])} {num(q[-1][1])}"
            dash = f' stroke-dasharray="{" ".join(num(float(n) * sw) for n in st.get("dashPattern", "3 3").split())}"' \
                if st.get("dashed") == "1" else ""
            stroke = st.get("strokeColor", INK)
            op = f' opacity="{float(st["opacity"]) / 100:g}"' if "opacity" in st else ""
            o.append(f"<g{op}>")
            o.append(f'<path d="{d}" fill="none" stroke="{stroke}" stroke-width="{num(sw)}"{dash}/>')
            nx, ny = -uy, ux
            hw = L * 0.42
            a = (bx - ux * L + nx * hw, by - uy * L + ny * hw)
            b = (bx - ux * L - nx * hw, by - uy * L - ny * hw)
            if filled:
                c = (bx - ux * L * 0.75, by - uy * L * 0.75)
                o.append(f'<path d="M{num(bx)} {num(by)} L{num(a[0])} {num(a[1])} L{num(c[0])} {num(c[1])} '
                         f'L{num(b[0])} {num(b[1])} Z" fill="{stroke}" stroke="{stroke}" stroke-width="{num(sw)}" '
                         f'stroke-linejoin="round"/>')
            else:
                o.append(f'<path d="M{num(a[0])} {num(a[1])} L{num(bx)} {num(by)} L{num(b[0])} {num(b[1])}" '
                         f'fill="none" stroke="{stroke}" stroke-width="{num(sw)}" stroke-linejoin="round"/>')
            o.append("</g>")
        return minx, miny, W, H, "\n".join(o)

    def save(self, path=None):
        path = path or self.path
        minx, miny, W, H, body = self.render()
        content = html.escape(ET.tostring(self.mxfile, encoding="unicode"), quote=True)
        svg = (f'<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="{num(W)}px" height="{num(H)}px" '
               f'viewBox="{num(minx)} {num(miny)} {num(W)} {num(H)}" content="{content}">\n'
               f'<g>\n{body}\n</g>\n</svg>\n')
        open(path, "w", encoding="utf-8").write(svg)
        self.path = path
