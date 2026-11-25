"""Standard diagram layouts built from the house style.

process_view(): the diagram used wherever the text defines a process. It is deliberately
abstract: one process call, with plain file boxes for its inputs and outputs and no channels
(the workflow diagrams show channels and their contents).

- File input: a file box. Several files received as one input (e.g. `path input_files`)
  sit together in one white group box.
- Without value inputs, everything is on one row: input -> process -> output.
- With value inputs, the file input moves to the row above and reaches the process by a bent
  arrow into its top (the staircase), and the value inputs enter from the left, labelled with
  their input names in the order the process declares them.
- Outputs on the right, level with the process; two outputs are stacked 40 px apart.
- Example values in file names. File contents sit under each file box (above the top box of a
  stack of two outputs), in the note style, centred.

Example:
    d = process_view("translateGreeting", inputs=[("Mateo.txt", "Hello, Mateo!")],
                     outputs=[("es-Mateo.txt", "Hola, Mateo!")])
    d.save("docs/en/docs/hello_nextflow/img/hello-translate-process.drawio.svg")
"""
import house_style as S
from drawio import Diagram

FONT_W = 0.55  # Helvetica advance per px of font size (approximate, for sizing boxes)


def _box_w(text, size, pad=S.BOX_SIDE_PAD, step=20):
    return int(-(-(len(text) * size * FONT_W + 2 * pad) // step) * step)


CONTENTS_GAP = 20  # file contents under (or above) their box


def _contents(d, cid, text, x, w, y, above=False):
    """File contents centred under the box spanning x..x+w whose bottom (or top, if above) is y."""
    if not text:
        return
    lines = text.split("\n")
    h = round(len(lines) * S.SIZE_NOTE * 1.2 + 6)
    tw = max(max(len(t) for t in lines) * S.SIZE_NOTE * FONT_W + 10, w)
    ty = y - CONTENTS_GAP - h if above else y + CONTENTS_GAP
    d.add_vertex(cid, text, S.vertex_style("note").replace("align=left", "align=center"), x + w / 2 - tw / 2, ty, tw, h)


def _input(d, files, x, y):
    """One file box, or a white group holding several files; files are (name, contents). Returns (w, h)."""
    if len(files) == 1:
        name, text = files[0]
        w = _box_w(name, S.SIZE_BOX)
        d.add_vertex("in", name, S.vertex_style("file"), x, y, w, S.BOX_H)
        _contents(d, "in-contents", text, x, w, y + S.BOX_H)
        return w, S.BOX_H
    iw = max(_box_w(f, S.SIZE_ITEM, pad=15, step=10) for f, _ in files)
    w = 2 * S.CHANNEL_PAD + len(files) * iw + (len(files) - 1) * S.ITEM_GAP
    d.add_vertex("in", "", S.vertex_style("item"), x, y, w, S.BOX_H)
    d.set_style(d.cells["in"], container=1, collapsible=0, recursiveResize=0)
    for i, (f, text) in enumerate(files):
        ix = x + S.CHANNEL_PAD + i * (iw + S.ITEM_GAP)
        d.add_vertex(f"in-{i}", f, S.vertex_style("file", S.SIZE_ITEM), ix, y + S.CHANNEL_PAD, iw, S.ITEM_H, parent="in")
        _contents(d, f"in-{i}-contents", text, ix, iw, y + S.BOX_H)
    return w, S.BOX_H


def process_view(process, inputs, outputs, values=(), process_w=None):
    """inputs: [(file name, contents)] received as one input; outputs: [(file name, contents)];
    values: input names, or (input name, example value) to show a value under a single label."""
    d = Diagram()
    pw = process_w or _box_w(process, S.SIZE_BOX)
    x0, y0 = 40, 40
    if not inputs:
        # no file input: the process on the first row, value inputs (if any) from the left
        vals = [v if isinstance(v, tuple) else (v, "") for v in values]
        lw = max([S.param_label_width(v) for v, _ in vals] or [0])
        px = -(-(x0 + lw + S.GAP_H) // 10) * 10 if vals else x0
        py = y0
        d.add_vertex("proc", process, S.vertex_style("process"), px, py, pw, S.BOX_H)
        fracs = {0: [], 1: [0.5], 2: [0.25, 0.75], 3: [0.2, 0.5, 0.8]}[len(vals)]
        for i, ((v, example), f) in enumerate(zip(vals, fracs)):
            y = py + f * S.BOX_H
            w = S.param_label_width(v)
            d.add_vertex(f"val{i}", v, S.vertex_style("param"), px - S.GAP_H - w, y - 19, w, 38)
            d.add_edge(f"e-val{i}", f"val{i}", "proc", S.edge_style("param"), exit=(1, 0.5), entry=(0, f))
            if example:
                _contents(d, f"val{i}-example", example, px - S.GAP_H - w, w - S.MENLO_ADVANCE * S.SIZE_PARAM, y + 19 - 14)
    elif values:
        iw, ih = _input(d, inputs, x0, y0)
        px = -(-(x0 + iw + 30) // 10) * 10
        # the process sits just low enough for its value labels to clear the input's contents
        lines = max([t.count("\n") + 1 for _, t in inputs if t] or [0])
        clear = y0 + ih + (20 + lines * S.SIZE_NOTE * 1.2 + 6 if lines else 0) + 20
        py = clear + 19 - min({1: [0.5], 2: [0.25, 0.75], 3: [0.2, 0.5, 0.8]}[len(values)]) * S.BOX_H
        py = max(py, y0 + ih + S.GAP_V)
        py = -(-py // 10) * 10
        d.add_vertex("proc", process, S.vertex_style("process"), px, py, pw, S.BOX_H)
        d.add_edge("e-in", "in", "proc", S.edge_style(orthogonal=True), exit=(1, 0.5), entry=(0.5, 0),
                   points=[(px + pw / 2, y0 + ih / 2)])
        fracs = {1: [0.5], 2: [0.25, 0.75], 3: [0.2, 0.5, 0.8]}[len(values)]
        values = [v if isinstance(v, tuple) else (v, "") for v in values]
        lx = px - S.GAP_H - max(S.param_label_width(v) for v, _ in values)
        for i, ((v, example), f) in enumerate(zip(values, fracs)):
            y = py + f * S.BOX_H
            w = S.param_label_width(v)
            d.add_vertex(f"val{i}", v, S.vertex_style("param"), lx, y - 19, w, 38)
            d.add_edge(f"e-val{i}", f"val{i}", "proc", S.edge_style("param"), exit=(1, 0.5), entry=(0, f))
            if example:  # an example value, under its label like file contents (single value only)
                _contents(d, f"val{i}-example", example, lx, w - S.MENLO_ADVANCE * S.SIZE_PARAM, y + 19 - 14)
    else:
        iw, ih = _input(d, inputs, x0, y0)
        px, py = x0 + iw + S.GAP_H, y0
        d.add_vertex("proc", process, S.vertex_style("process"), px, py, pw, S.BOX_H)
        d.add_edge("e-in", "in", "proc", S.edge_style(), exit=(1, 0.5), entry=(0, 0.5))
    ox = px + pw + S.GAP_H
    tops = [py] if len(outputs) == 1 else [py + S.BOX_H / 2 - S.BOX_H - S.STACK_GAP / 2, py + S.BOX_H / 2 + S.STACK_GAP / 2]
    exits = [0.5] if len(outputs) == 1 else [0.25, 0.75]
    ow = max(_box_w(name, S.SIZE_BOX) for name, _ in outputs)  # stacked outputs share a width
    for i, ((name, text), y, ex) in enumerate(zip(outputs, tops, exits)):
        w = ow
        d.add_vertex(f"out{i}", name, S.vertex_style("file"), ox, y, w, S.BOX_H)
        d.add_edge(f"e-out{i}", "proc", f"out{i}", S.edge_style(), exit=(1, ex), entry=(0, 0.5))
        top_of_stack = len(outputs) > 1 and i == 0
        _contents(d, f"out{i}-contents", text, ox, w, y if top_of_stack else y + S.BOX_H, above=top_of_stack)
    return d
