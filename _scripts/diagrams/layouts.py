"""Standard diagram layouts built from the house style.

process_view(): the diagram used wherever the text defines or changes a process. It shows
the process on its own, from the process's point of view:

- its channel input on the top row, reaching the process by a bent arrow into its top
  (the staircase), captioned in plain words;
- its value inputs entering from the left, labelled with their input names, in the order
  the process declares them;
- its output channel(s) on the right; two outputs are stacked with the top caption above;
- captions 25 px from their channel (CAPTION_GAP_PROCESS).

Example:
    d = process_view("translateGreeting",
                     channel=Channel(["Raphael.txt", "Celine.txt", "Mateo.txt"], "greetings"),
                     values=["language"],
                     outputs=[Channel(["es-Raphael.txt", "es-Celine.txt", "es-Mateo.txt"], "translated greetings")])
    d.save("docs/en/docs/hello_nextflow/img/hello-translate-process.drawio.svg")
"""
from dataclasses import dataclass, field

import house_style as S
from drawio import Diagram

FONT_W = 0.55  # Helvetica advance per px of font size, for sizing items (approximate)


@dataclass
class Channel:
    items: list
    caption: str
    nested: bool = False  # one element that is a list of the items (after collect())
    note: str = ""  # sample contents, shown to the right of the channel
    item_w: int = field(default=0)

    def width(self):
        if not self.item_w:
            longest = max(len(i) for i in self.items) * S.SIZE_ITEM * FONT_W
            self.item_w = int(-(-(longest + 30) // 10) * 10)
        n = len(self.items)
        w = 2 * S.CHANNEL_PAD + n * self.item_w + (n - 1) * S.ITEM_GAP
        return w + (2 * S.CHANNEL_PAD if self.nested else 0)

    def height(self):
        return 110 if self.nested else 80


def _channel(d, cid, ch, x, y, caption_y, note_lines=0):
    w, h = ch.width(), ch.height()
    d.add_vertex(cid, "", S.vertex_style("channel"), x, y, w, h)
    parent, ix, iy = cid, x + S.CHANNEL_PAD, y + S.CHANNEL_PAD
    if ch.nested:
        inner_w = w - 2 * S.CHANNEL_PAD
        d.add_vertex(cid + "-list", "", S.vertex_style("item"), ix, iy, inner_w, 80, parent=cid)
        d.set_style(d.cells[cid + "-list"], container=1, collapsible=0, recursiveResize=0)
        parent, ix, iy = cid + "-list", ix + S.CHANNEL_PAD, iy + S.CHANNEL_PAD
    for i, item in enumerate(ch.items):
        d.add_vertex(f"{cid}-{i}", item, S.vertex_style("file", S.SIZE_ITEM), ix + i * (ch.item_w + S.ITEM_GAP), iy,
                     ch.item_w, S.ITEM_H, parent=parent)
    d.add_vertex(cid + "-cap", ch.caption, S.vertex_style("caption"), x, caption_y, w, S.CAPTION_H)
    if ch.note:
        lines = ch.note.count("\n") + 1
        nh = round(lines * S.SIZE_NOTE * 1.2 + 10)
        d.add_vertex(cid + "-note", ch.note, S.vertex_style("note"), x + w + S.NOTE_GAP, y + h / 2 - nh / 2,
                     max(len(x) for x in ch.note.split("\n")) * S.SIZE_NOTE * FONT_W + 10, nh)
    return w, h


def process_view(process, channel, outputs, values=(), process_w=None):
    d = Diagram()
    x0, y0 = 40, 40
    cw, ch = _channel(d, "in", channel, x0, y0, y0 + channel.height() + S.CAPTION_GAP_PROCESS)
    pw = process_w or int(-(-(len(process) * S.SIZE_BOX * FONT_W + 2 * S.BOX_SIDE_PAD) // 20) * 20)
    px = x0 + cw + 30
    px = -(-px // 10) * 10
    # leave room under the input caption for the value labels when there are any
    py = y0 + ch + (130 if values else S.GAP_V)
    py += (-(py + S.BOX_H / 2)) % 10 if (py + S.BOX_H / 2) % 10 else 0
    d.add_vertex("proc", process, S.vertex_style("process"), px, py, pw, S.BOX_H)
    d.add_edge("e-in", "in", "proc", S.edge_style(orthogonal=True), exit=(1, 0.5), entry=(0.5, 0),
               points=[(px + pw / 2, y0 + ch / 2)])
    # value inputs: evenly spaced on the left side, labels level with their entry
    n = len(values)
    fracs = {1: [0.5], 2: [0.25, 0.75], 3: [0.2, 0.5, 0.8]}[n] if n else []
    lx = px - S.GAP_H - max([S.param_label_width(v) for v in values] or [0]) + 0
    for i, (v, f) in enumerate(zip(values, fracs)):
        y = py + f * S.BOX_H
        d.add_vertex(f"val{i}", v, S.vertex_style("param"), lx, y - 19, S.param_label_width(v), 38)
        d.add_edge(f"e-val{i}", f"val{i}", "proc", S.edge_style("param"), exit=(1, 0.5), entry=(0, f))
    # outputs: level with the process, or a stack of two centred on it (top caption above)
    ox = px + pw + S.GAP_H
    if len(outputs) == 1:
        o = outputs[0]
        y = py + S.BOX_H / 2 - o.height() / 2
        _channel(d, "out0", o, ox, y, y + o.height() + S.CAPTION_GAP_PROCESS)
        d.add_edge("e-out0", "proc", "out0", S.edge_style(), exit=(1, 0.5), entry=(0, 0.5))
    else:
        top = py + S.BOX_H / 2 - (outputs[0].height() + S.STACK_GAP + outputs[1].height()) / 2
        y1 = top + outputs[0].height() + S.STACK_GAP
        _channel(d, "out0", outputs[0], ox, top, top - S.CAPTION_GAP_PROCESS - S.CAPTION_H)
        _channel(d, "out1", outputs[1], ox, y1, y1 + outputs[1].height() + S.CAPTION_GAP_PROCESS)
        d.add_edge("e-out0", "proc", "out0", S.edge_style(), exit=(1, 0.25), entry=(0, 0.5))
        d.add_edge("e-out1", "proc", "out1", S.edge_style(), exit=(1, 0.75), entry=(0, 0.5))
    return d
