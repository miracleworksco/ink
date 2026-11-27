# HarfBuzz's shaping of std's vendored fonts (uharfbuzz 0.57, fontTools for
# glyph names), written to shaped.json: the reference font/tests.rs checks
# ink's shaping against. Each text's glyphs: [id, x advance, x offset,
# y offset, cluster (a byte of the text's UTF-8), name], in font units.
# test: not rerun (needs uharfbuzz and fontTools)
#   python3 shaped.py ../../../../std/fonts > shaped.json
import json, sys
import uharfbuzz as hb
from fontTools.ttLib import TTFont
fonts = sys.argv[1]
FILES = {"sans": "Lexend-wght.ttf", "rounded": "Fredoka-wdth-wght.ttf", "serif": "Bitter-wght.ttf", "script": "Pacifico-Regular.ttf"}
CASES = [
    ("sans", {}, ["naïve", "café", "Zoë", "office", "fi fl", "AVATAR", "e\u0301te\u0301", "ni\u0308"]),
    ("sans", {"wght": 900}, ["office", "Wave", "naïve", "x\u0301q\u0308"]),
    ("rounded", {"wght": 700, "wdth": 75}, ["naïve", "nai\u0308ve", "café", "Zoë", "office", "x\u0301q\u0308"]),
    ("serif", {"wght": 300}, ["office", "affluent", "Zoë", "cafe\u0301"]),
    ("script", {}, ["Grace", "Ada Lovelace", "office", "naïve", "Zoë"]),
]
out = []
for key, axes, texts in CASES:
    path = f"{fonts}/{FILES[key]}"
    blob = hb.Blob.from_file_path(path)
    face = hb.Face(blob)
    font = hb.Font(face)
    if axes:
        font.set_variations(axes)
    order = TTFont(path).getGlyphOrder()
    for t in texts:
        buf = hb.Buffer()
        buf.add_str(t)
        buf.guess_segment_properties()
        hb.shape(font, buf, {})
        gs = [[i.codepoint, p.x_advance, p.x_offset, p.y_offset, i.cluster, order[i.codepoint]] for i, p in zip(buf.glyph_infos, buf.glyph_positions)]
        out.append({"font": key, "axes": axes, "text": t, "glyphs": gs})
print(json.dumps(out, ensure_ascii=False, indent=None).replace("}, {", "},\n{"))
