# Course diagrams

Course diagrams are `.drawio.svg` files: ordinary SVG images that also carry the editable draw.io model.
They are inlined into lesson pages with `<figure class="drawio">` and a `--8<--` snippet.
Edit them in VS Code with the [Draw.io Integration](https://marketplace.visualstudio.com/items?itemName=hediet.vscode-drawio) extension, or at [diagrams.net](https://app.diagrams.net).

The scripts here need only Python 3 (no extra packages).

| Script                    | Purpose                                                                                        |
| ------------------------- | ---------------------------------------------------------------------------------------------- |
| `check.py`                | Check diagrams against the house style. Errors fail; spacing issues are warnings (`--strict`). |
| `excalidraw_to_drawio.py` | Turn an old Excalidraw SVG export into a styled `.drawio.svg` starting point.                  |
| `house_style.py`          | Every colour, font, line weight and spacing value. Change the style here.                      |
| `drawio.py`               | Read, edit, render and write `.drawio.svg` files from Python.                                  |

```bash
python3 _scripts/diagrams/check.py docs/en/docs/hello_nextflow/img/
python3 _scripts/diagrams/excalidraw_to_drawio.py old.svg new.drawio.svg
```

## Workflow

1. Generate a starting file (convert an old export, or build one with `drawio.py`).
2. Open it in VS Code and adjust the layout by hand.
3. Run `check.py` and fix what it reports.
4. Save once in the editor so draw.io redraws the image. Files written by `drawio.py` are already redrawn by its own renderer.

After a script edits a file, close and reopen its editor tab before saving, or the editor overwrites the change with the copy it has open.

## House style

**Colours** (fill / outline / text):

| Role                   | Fill      | Outline   | Text      |
| ---------------------- | --------- | --------- | --------- |
| Process                | `#2d799b` | `#055672` | white     |
| File                   | `#f0f9fb` | `#02899a` | `#1e1e1e` |
| Channel                | `#76d2c8` | `#018178` | `#1e1e1e` |
| Value inside a channel | white     | `#018178` | `#1e1e1e` |
| Lines, notes, captions |           | `#1e1e1e` |           |
| Parameter labels, code |           |           | `#1e1e1e` |

The palette stays distinct under protanopia, deuteranopia and tritanopia, and keeps its hues under the site's dark-mode filter.

**Text**: plain SVG text only (`html=0`, no `whiteSpace=wrap`; break lines explicitly).
Helvetica for labels (boxes 28 px, items in channels 22, notes 22, captions 28 bold in sentence case; a caption naming a code object, such as `sayHello.out`, is Menlo bold); Menlo for parameters (26) and code (20).
No handwritten fonts.

**Lines** (diagrams are shown at under half size, so they are drawn heavy):
data arrows 3 px with filled arrowheads (size 12);
parameter arrows 2.5 px, black, dashed `4 2`, open arrowheads;
box outlines 2.5; frames 2.

**Layout** (10 px grid; a box's edge or centre sits on it):

- Boxes are 80 px high; items inside a channel 50.
- Gaps between boxes joined by an arrow: 130 px horizontally, 70 vertically.
- Arrow into or out of a frame: 80 px gap, or 100 when a parameter label sits in that gap.
- Frames have 30 px padding on every side.
- Vertical arrows are exactly vertical, leaving and landing at the centre of both boxes. Fix misalignment by moving the downstream box, not the arrow's anchor.
- Arrows between level boxes leave and land at the middle of the facing sides.
- Fan-outs: targets stacked 40 px apart, the stack centred on the source; arrows leave at evenly spaced points (0.25 / 0.5 / 0.75, or 0.3 / 0.7) and land at each target's centre.
- Several arrows into one side of a box land at evenly spaced points.
- Notes sit 20 px right of their box, centred on it. Captions (36 px high) sit 40 px under the box or channel they describe, or 15 px under a dashed group's outline (the group's own padding supplies the rest).
- Parameter labels sit level with the point where their arrow lands, and end one monospace space before the arrow starts. A parameter that feeds several steps uses one trunk that branches.
