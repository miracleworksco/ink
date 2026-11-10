# Crafted fonts lettering-max must bound or refuse, first set: each
# glyph a rectangle, "wide" and "syllable" 3000 units tall. Written to
# testdata/; font/tests.rs's reproducer table says how each is read.
# test: not rerun (needs fontTools)
#   python3 make-fonts.py
from fontTools.fontBuilder import FontBuilder
from fontTools.pens.ttGlyphPen import TTGlyphPen
from fontTools.feaLib.builder import addOpenTypeFeaturesFromString
from fontTools.varLib.featureVars import addFeatureVariations
from pathlib import Path
out=Path(__file__).parent / 'testdata'

def build(name,fea='',cmap=None,var=False):
    fb=FontBuilder(1000,isTTF=True)
    names=['.notdef','a','b','x','wide','H','jamoL','jamoV','syllable']
    fb.setupGlyphOrder(names)
    fb.setupCharacterMap(cmap or {97:'a',98:'b',72:'H'})
    glyphs={}
    for n in names:
        p=TTGlyphPen(None)
        h=3000 if n in ['wide','syllable'] else 500
        p.moveTo((0,0));p.lineTo((500,0));p.lineTo((500,h));p.lineTo((0,h));p.closePath()
        glyphs[n]=p.glyph()
    fb.setupGlyf(glyphs)
    fb.setupHorizontalMetrics({n:(600,0) for n in names})
    fb.setupHorizontalHeader(ascent=3000,descent=-500)
    fb.setupOS2(sTypoAscender=3000,sTypoDescender=-500,sCapHeight=500)
    fb.setupNameTable({'familyName':'Review','styleName':'Regular','uniqueFontIdentifier':name,'fullName':name,'psName':name})
    fb.setupPost();fb.setupMaxp()
    if var: fb.setupFvar([('wght',100,400,900,'Weight')],[])
    if fea: addOpenTypeFeaturesFromString(fb.font,fea)
    return fb.font

def save(name,fea='',**kw):
    f=build(name,fea,**kw); f.save(out/(name+'.ttf')); return f
save('required','languagesystem latn dflt; feature xxxx { script latn; language dflt required; sub a by wide; } xxxx;')
save('stages','languagesystem latn dflt; lookup EARLY { sub x by wide; } EARLY; lookup LATE { sub a by x; } LATE; feature ccmp { lookup LATE; } ccmp; feature liga { lookup EARLY; } liga;')
save('repeated','languagesystem latn dflt; lookup SHIFT { pos a <0 1000 0 0>; } SHIFT; feature kern { pos a\' lookup SHIFT lookup SHIFT; } kern;')
save('hangul',cmap={0x1100:'jamoL',0x1161:'jamoV',0xac00:'syllable',72:'H'})
f=build('variation','languagesystem latn dflt;',var=True)
addFeatureVariations(f,[([{'wght':(0.5,1.0)}],{'a':'wide'})],featureTag='rvrn')
f.save(out/'variation.ttf')
save('arabic-stages','languagesystem arab dflt; lookup EARLY { sub x by wide; } EARLY; lookup LATE { sub a by x; } LATE; feature ccmp { lookup LATE; } ccmp; feature liga { lookup EARLY; } liga;',cmap={0x628:'a',72:'H'})
save('arabic-context-stages','languagesystem arab dflt; lookup EARLY { sub a x\' by wide; } EARLY; lookup LATE { sub b by x; } LATE; feature ccmp { lookup LATE; } ccmp; feature liga { lookup EARLY; } liga;',cmap={0x628:'a',0x62a:'b',72:'H'})
save('required-context','languagesystem latn dflt; feature xxxx { script latn; language dflt required; sub a b\' by wide; } xxxx;')
save('arabic-context-stages2','languagesystem arab dflt; lookup EARLY { sub x\' a by wide; } EARLY; lookup LATE { sub b by x; } LATE; feature ccmp { lookup LATE; } ccmp; feature liga { lookup EARLY; } liga;',cmap={0x628:'a',0x62a:'b',72:'H'})
save('arabic-ligature-stages','languagesystem arab dflt; lookup EARLY { sub a x by wide; } EARLY; lookup LATE { sub b by x; } LATE; feature ccmp { lookup LATE; } ccmp; feature liga { lookup EARLY; } liga;',cmap={0x628:'a',0x62a:'b',72:'H'})
