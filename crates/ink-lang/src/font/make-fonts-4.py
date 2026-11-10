# Crafted fonts lettering-max must bound or refuse, fourth set: variation
# deltas past f32's whole numbers (2^24), which the shaper turns into f32
# before it rounds them: GPOS y placements of ±16777217, ±16777219 and
# ±33554435 units at wght 900, and HVAR advances of 67108865 to
# 134217739 units at wght 650. Written to testdata/.
# test: not rerun (needs fontTools)
#   python3 make-fonts-4.py
from pathlib import Path
from fontTools.fontBuilder import FontBuilder
from fontTools.pens.ttGlyphPen import TTGlyphPen
from fontTools.feaLib.builder import addOpenTypeFeaturesFromString
from fontTools.ttLib import newTable
from fontTools.ttLib.tables import otTables
from fontTools.varLib.builder import buildVarRegionList,buildVarData,buildVarStore,buildVarDevTable
out=Path(__file__).parent / 'testdata'
cases=[]
def build(name,fea='',var=False,cmap=None):
 fb=FontBuilder(1000,isTTF=True);names=['.notdef','a','b','x','H','mark','space','acute','eacute']
 fb.setupGlyphOrder(names);fb.setupCharacterMap(cmap or {32:'space',72:'H',97:'a',98:'b',120:'x'})
 gs={}
 for n in names:
  p=TTGlyphPen(None)
  p.moveTo((0,0));p.lineTo((500,0));p.lineTo((500,500));p.lineTo((0,500));p.closePath();gs[n]=p.glyph()
 fb.setupGlyf(gs);fb.setupHorizontalMetrics({n:(600,0) for n in names});fb.setupHorizontalHeader(ascent=1000,descent=-500)
 fb.setupOS2(sTypoAscender=1000,sTypoDescender=-500,sCapHeight=500)
 fb.setupNameTable({'familyName':'Review4','styleName':'Regular','uniqueFontIdentifier':name,'fullName':name,'psName':name});fb.setupPost();fb.setupMaxp()
 if var:fb.setupFvar([('wght',100,400,900,'Weight')],[])
 if fea:addOpenTypeFeaturesFromString(fb.font,fea)
 return fb.font

def save(f,name,chars='ab',text='ab',axes='',tracking=0,n=None):
 f.save(out/(name+'.ttf'))
 cases.append(dict(name=name,chars=chars,text=text,axes=axes,tracking=tracking,n=n or len(text)))
def add_store(f,delta,region=(0,1,1)):
 g=newTable('GDEF');g.table=otTables.GDEF();f['GDEF']=g
 for k in ['GlyphClassDef','AttachList','LigCaretList','MarkAttachClassDef']:setattr(g.table,k,None)
 g.table.Version=0x10003;g.table.VarStore=buildVarStore(buildVarRegionList([{'wght':region}],['wght']),[buildVarData([0],[[delta]],False)])
 return g.table.VarStore

for delta in [16777217,16777219,-16777217,-16777219,33554435,-33554435]:
 name='gpos-float-'+str(delta)
 f=build(name,'languagesystem latn dflt; feature kern { pos a <0 1 0 0>; } kern;',True)
 add_store(f,delta);s=f['GPOS'].table.LookupList.Lookup[0].SubTable[0];s.Value.YPlacement=0;s.Value.YPlaDevice=buildVarDevTable(0);s.ValueFormat|=0x20
 save(f,name,'a','a','(wght 900)')


for delta in [67108865,67108869,67108873,134217739]:
 name='hvar-float-'+str(delta)
 f=build(name,var=True);vs=add_store(f,delta)
 t=newTable('HVAR');t.table=otTables.HVAR();t.table.Version=0x10000;t.table.VarStore=vs
 t.table.AdvWidthMap=otTables.VarIdxMap();t.table.AdvWidthMap.mapping={g:0 for g in f.getGlyphOrder()}
 t.table.LsbMap=None;t.table.RsbMap=None;f['HVAR']=t
 save(f,name,'ab','abab','(wght 650)')
