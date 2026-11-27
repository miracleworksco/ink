# Crafted fonts lettering-max must bound or refuse, second set (each
# glyph a 500-unit square unless said). Written to testdata/ (the
# required-feature stages font as stages-required.ttf).
# test: not rerun (needs fontTools)
#   python3 make-fonts-2.py
from pathlib import Path
from fontTools.fontBuilder import FontBuilder
from fontTools.pens.ttGlyphPen import TTGlyphPen
from fontTools.feaLib.builder import addOpenTypeFeaturesFromString
from fontTools.varLib.builder import buildVarRegionList, buildVarData, buildVarStore, buildVarDevTable
from fontTools.ttLib import newTable
from fontTools.ttLib.tables import otTables
OUT=Path(__file__).parent / 'testdata'

def build(name,fea='',var=False):
    fb=FontBuilder(1000,isTTF=True)
    names=['.notdef','a','b','x','H']
    fb.setupGlyphOrder(names)
    fb.setupCharacterMap({97:'a',98:'b',120:'x',72:'H'})
    glyphs={}
    for n in names:
        p=TTGlyphPen(None)
        p.moveTo((0,0));p.lineTo((500,0));p.lineTo((500,500));p.lineTo((0,500));p.closePath()
        glyphs[n]=p.glyph()
    fb.setupGlyf(glyphs)
    fb.setupHorizontalMetrics({n:(600,0) for n in names})
    fb.setupHorizontalHeader(ascent=1000,descent=-500)
    fb.setupOS2(sTypoAscender=1000,sTypoDescender=-500,sCapHeight=500)
    fb.setupNameTable({'familyName':'Review2','styleName':'Regular','uniqueFontIdentifier':name,'fullName':name,'psName':name})
    fb.setupPost();fb.setupMaxp()
    if var: fb.setupFvar([('wght',100,400,900,'Weight')],[])
    if fea: addOpenTypeFeaturesFromString(fb.font,fea)
    return fb.font

f=build('fractional','languagesystem latn dflt; feature kern { pos a <0 1 0 0>; } kern;',True)
gdef=f.get('GDEF')
if gdef is None:
    gdef=newTable('GDEF');gdef.table=otTables.GDEF();f['GDEF']=gdef
    for k in ['GlyphClassDef','AttachList','LigCaretList','MarkAttachClassDef']: setattr(gdef.table,k,None)
gdef.table.Version=0x00010003
gdef.table.VarStore=buildVarStore(buildVarRegionList([{'wght':(0,1,1)}],['wght']),[buildVarData([0],[[1]],False)])
v=f['GPOS'].table.LookupList.Lookup[0].SubTable[0].Value
v.YPlacement=0
v.YPlaDevice=buildVarDevTable(0)
f['GPOS'].table.LookupList.Lookup[0].SubTable[0].ValueFormat |= 0x20
f.save(OUT/'fractional.ttf')

for name,fea in [
 ('nested','lookup SHIFT { pos a <0 1000 0 0>; } SHIFT; lookup INNER { pos a a\' lookup SHIFT; } INNER; feature kern { pos a\' lookup INNER; } kern;'),
 ('split','lookup SPLIT { sub x by a a a; } SPLIT; lookup LIG { sub a a by x; } LIG; feature rvrn { lookup LIG; } rvrn; feature liga { lookup SPLIT; } liga;'),
 ('stages','lookup SHIFT { pos a <0 1000 0 0>; } SHIFT; feature kern { lookup SHIFT; } kern; feature xxxx { script latn; language dflt required; lookup SHIFT; } xxxx;'),
]:
 f=build(name,'languagesystem latn dflt; '+fea)
 f.save(OUT/(('stages-required' if name == 'stages' else name)+'.ttf'))
# Nested context ignores bases after its first position. Its caller's
# single-glyph input lets several starting bases reach one ligature.
f=build('nested-skip', '''languagesystem latn dflt;
table GDEF { GlyphClassDef [a H], [b], , ; } GDEF;
lookup SHIFT { pos b <0 1000 0 0>; } SHIFT;
lookup INNER { lookupflag IgnoreBaseGlyphs; pos a' b' lookup SHIFT; } INNER;
feature kern { pos a' lookup INNER; } kern;
''')
f.save(OUT/'nested-skip.ttf')
f=build('avar2','languagesystem latn dflt; feature kern { pos a <0 1 0 0>; } kern;',True)
gdef=newTable('GDEF');gdef.table=otTables.GDEF();f['GDEF']=gdef
for k in ['GlyphClassDef','AttachList','LigCaretList','MarkAttachClassDef']: setattr(gdef.table,k,None)
gdef.table.Version=0x00010003
gdef.table.VarStore=buildVarStore(buildVarRegionList([{'wght':(0,1,1)}],['wght']),[buildVarData([0],[[1000]],False)])
v=f['GPOS'].table.LookupList.Lookup[0].SubTable[0].Value
v.YPlacement=0;v.YPlaDevice=buildVarDevTable(0)
f['GPOS'].table.LookupList.Lookup[0].SubTable[0].ValueFormat|=0x20
avar=newTable('avar');avar.majorVersion=2;avar.minorVersion=0
avar.segments={'wght':{-1:-1,0:0,.5:1,1:1}}
avar.table=otTables.avar();avar.table.VarIdxMap=None;avar.table.VarStore=None
f['avar']=avar
f.save(OUT/'avar2.ttf')
# Legacy AAT insertion: after a, b inserts 31 x's. Each alone is unchanged.
import struct
from fontTools.ttLib.tables.DefaultTable import DefaultTable
p16=lambda *v:struct.pack('>'+'H'*len(v),*v)
p32=lambda *v:struct.pack('>'+'I'*len(v),*v)
classes=p16(1,2)+bytes([4,5])
states=bytes([0,0,0,0,1,0]*2+[0,0,0,0,1,2])
entries=p16(16,0,0xffff,0xffff)+p16(28,0,0xffff,0xffff)+p16(16,31<<5,0,0xffff)
body=p16(6,10,16,34,58)+classes+states+entries+p16(*([3]*31))
sub=p16(8+len(body),5)+p32(1)+body
chain=p32(1,12+len(sub))+p16(0,1)+sub
f=build('mort')
t=DefaultTable('mort');t.data=p32(0x10000,1)+chain;f['mort']=t
f.save(OUT/'mort.ttf')
f=build('large-glyph')
p=TTGlyphPen(None)
p.moveTo((0,0));p.lineTo((20000,0));p.lineTo((20000,20000));p.lineTo((0,20000));p.closePath()
f['glyf']['a']=p.glyph()
f['hmtx'].metrics['a']=(20000,0)
f.save(OUT/'large-glyph.ttf')
for name,fea in [
 ('deep-context','''lookup SHIFT { pos b <0 1000 0 0>; } SHIFT;
 lookup C1 { pos b' lookup SHIFT lookup SHIFT; } C1;
 lookup C2 { pos b' lookup C1 lookup C1; } C2;
 lookup C3 { pos a' b' lookup C2 lookup C2; } C3;
 feature kern { pos a' lookup C3; } kern;'''),
 ('shared-stages','''lookup SPLIT { sub b by b b; } SPLIT;
 lookup CTX { sub a b' lookup SPLIT; } CTX;
 feature rvrn { lookup CTX; } rvrn;
 feature liga { lookup CTX; } liga;
 feature calt { lookup CTX; } calt;'''),
]:
 f=build(name,'languagesystem latn dflt; '+fea)
 f.save(OUT/(('stages-required' if name == 'stages' else name)+'.ttf'))
from fontTools.varLib.featureVars import addFeatureVariationsRaw
f=build('gpos-variation','''languagesystem latn dflt;
lookup SHIFT { pos a <0 1000 0 0>; } SHIFT;
feature kern { pos b <0 0 0 0>; } kern;
''',True)
addFeatureVariationsRaw(f,f['GPOS'].table,[({'wght':(.5,1)},[0])],featureTag='kern')
f.save(OUT/'gpos-variation.ttf')
