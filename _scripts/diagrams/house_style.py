"""House style for course diagrams (.drawio.svg).

Every rule the checker enforces and every value the converter writes comes from here,
so this file is the single place to change the style.
"""

GRID = 10

# ---------------------------------------------------------------- colours
INK = "#1e1e1e"  # all lines and default text

# role: (fill, stroke, text)
ROLES = {
    "process": ("#2d799b", "#055672", "#ffffff"),
    "file": ("#f0f9fb", "#02899a", INK),
    "channel": ("#76d2c8", "#018178", INK),
    "item": ("#ffffff", "#018178", INK),  # a value (not a file) inside a channel
}
# fills used by the old Excalidraw diagrams, mapped to roles by the converter
LEGACY_FILLS = {
    "#ffec99": "process",
    "#fff9db": "process",
    "#a5d8ff": "file",
    "#e7f5ff": "file",
    "#b2f2bb": "channel",
    "#ffffff": "item",
}

# ---------------------------------------------------------------- fonts
LABEL_FONT = "Helvetica"
CODE_FONT = "Menlo"
ALLOWED_FONTS = {LABEL_FONT, CODE_FONT}
FORBIDDEN_FONTS = ("Excalifont", "Virgil", "Xiaolai", "Comic Shanns", "Cascadia")
SIZE_BOX = 28  # labels in boxes
SIZE_ITEM = 22  # labels of items inside a channel
SIZE_NOTE = 22  # sample file contents next to a box
SIZE_CAPTION = 28  # captions under groups ("input channel"), bold, sentence case
CAPTION_H = 36
CAPTION_GAP = 40  # caption below the box or channel it describes
CAPTION_GAP_PROCESS = 25  # same, in a diagram of a single process definition (tighter, compact view)
CAPTION_GAP_GROUP = 15  # caption below a dashed group outline (the group padding already adds space)
SIZE_PARAM = 26  # parameter labels (--batch trio), Menlo
SIZE_CODE = 20  # code snippets, Menlo
MENLO_ADVANCE = 0.602  # Menlo advance width per em (every glyph is the same width)
TEXT_SPACING_LEFT = 4

# ---------------------------------------------------------------- line weights
# diagrams are drawn ~1700 px wide and shown at <50% on the page, so draw heavy
W_DATA_ARROW = 3
W_PARAM_ARROW = 2.5
W_BOX = 2.5
W_FRAME = 2
ARROW_SIZE = 12
PARAM_DASH = "4 2"
GROUP_DASH = "6 4"  # dashed group outline (e.g. "process calls")
GROUP_STROKE = "#868e96"

# ---------------------------------------------------------------- geometry
BOX_H = 80  # process and file boxes
ITEM_H = 50  # items inside a channel
GAP_H = 130  # horizontal gap for an arrow between two boxes
GAP_V = 70  # vertical gap for an arrow between two boxes
GAP_FRAME = 80  # gap for an arrow into or out of a frame...
GAP_FRAME_PARAM = 100  # ...when a parameter label or trunk sits in that gap
STACK_GAP = 40  # vertical gap between stacked boxes (fan-out outputs, parallel calls)
FRAME_PAD = 30
NOTE_GAP = 20  # note to the right of its box
CHANNEL_PAD = 15  # padding inside a channel around its items
ITEM_GAP = 15  # gap between items in a channel
BOX_SIDE_PAD = 38  # min padding each side of a box label

ARC_BOX = 16
ARC_CHANNEL = 80
ARC_FRAME = 24


def menlo_width(text, size=SIZE_PARAM):
    return len(text) * size * MENLO_ADVANCE


def param_label_width(text, size=SIZE_PARAM):
    """Label cell width so the arrow starts one monospace space after the text."""
    return TEXT_SPACING_LEFT + menlo_width(text, size) + size * MENLO_ADVANCE


def _s(**kv):
    return "".join(f"{k}={v};" for k, v in kv.items())


BASE_TEXT = _s(html=0, whiteSpace="nowrap")


def vertex_style(role, size=None, code=False):
    """code=True: a caption that names a code object (sayHello.out) is set in the code font."""
    """draw.io style string for a vertex of the given role."""
    if role in ROLES:
        fill, stroke, text = ROLES[role]
        arc = ARC_CHANNEL if role == "channel" else ARC_BOX
        st = "rounded=1;" + _s(absoluteArcSize=1, arcSize=arc, fillColor=fill, strokeColor=stroke, strokeWidth=W_BOX,
                               fontFamily=LABEL_FONT, fontSize=size or (SIZE_ITEM if role == "item" else SIZE_BOX),
                               fontColor=text) + BASE_TEXT
        if role == "channel":
            st += _s(container=1, collapsible=0, recursiveResize=0)
        return st
    if role == "frame":
        return "rounded=1;" + _s(absoluteArcSize=1, arcSize=ARC_FRAME, fillColor="none", strokeColor=INK,
                                 strokeWidth=W_FRAME, container=1, collapsible=0, recursiveResize=0) + BASE_TEXT
    if role == "group":
        return "rounded=1;" + _s(absoluteArcSize=1, arcSize=ARC_FRAME, fillColor="none", strokeColor=GROUP_STROKE,
                                 strokeWidth=W_FRAME, dashed=1, dashPattern=GROUP_DASH, container=1, collapsible=0,
                                 recursiveResize=0) + BASE_TEXT
    if role == "note":
        return "text;" + _s(align="left", verticalAlign="middle", spacing=0, spacingLeft=TEXT_SPACING_LEFT,
                            fontFamily=LABEL_FONT, fontSize=size or SIZE_NOTE, fontColor=INK) + BASE_TEXT
    if role == "caption":
        return "text;" + _s(align="center", verticalAlign="middle", fontFamily=LABEL_FONT,
                            fontSize=size or SIZE_CAPTION, fontStyle=1, fontColor=INK).replace(
                                f"fontFamily={LABEL_FONT}", f"fontFamily={CODE_FONT if code else LABEL_FONT}") + BASE_TEXT
    if role == "param":
        return "text;" + _s(align="left", verticalAlign="middle", spacing=0, spacingLeft=TEXT_SPACING_LEFT,
                            fontFamily=CODE_FONT, fontSize=size or SIZE_PARAM, fontColor=INK) + BASE_TEXT
    if role == "code":
        return "text;" + _s(align="left", verticalAlign="middle", spacing=0, spacingLeft=TEXT_SPACING_LEFT,
                            fontFamily=CODE_FONT, fontSize=size or SIZE_CODE, fontColor=INK) + BASE_TEXT
    raise ValueError(role)


def edge_style(kind="data", orthogonal=False):
    st = _s(edgeStyle="orthogonalEdgeStyle" if orthogonal else "none", html=0, strokeColor=INK, endSize=ARROW_SIZE)
    if orthogonal:
        st += _s(rounded=1, arcSize=40)
    if kind == "param":
        return st + _s(strokeWidth=W_PARAM_ARROW, dashed=1, dashPattern=PARAM_DASH, endArrow="open", endFill=0)
    return st + _s(strokeWidth=W_DATA_ARROW, endArrow="classic", endFill=1)
