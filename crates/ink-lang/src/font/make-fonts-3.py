# Crafted fonts lettering-max must bound or refuse, third set: the two
# font readers disagreeing (avar 1.1, HVAR overflow and versions, axis
# normalization, GDEF 1.4, variation arithmetic) and layout formats with
# wider offsets. Reuses make-fonts-2.py's builder. Written to testdata/.
# test: not rerun (needs fontTools)
#   python3 make-fonts-3.py
from pathlib import Path
from fontTools.ttLib import newTable
from fontTools.ttLib.tables import otTables
from fontTools.varLib.builder import buildVarRegionList,buildVarData,buildVarStore,buildVarDevTable
r=Path(__file__).parent / 'testdata'
ns={'__file__':__file__}; exec((Path(__file__).parent/'make-fonts-2.py').read_text().split("f=build('fractional'")[0],ns)
build=ns['build']
def store(delta):return buildVarStore(buildVarRegionList([{'wght':(0,1,1)}],['wght']),[buildVarData([0],[[delta]],False)])
def gdef(f,delta):
 g=newTable('GDEF');g.table=otTables.GDEF();f['GDEF']=g
 for k in ['GlyphClassDef','AttachList','LigCaretList','MarkAttachClassDef']:setattr(g.table,k,None)
 g.table.Version=0x10003;g.table.VarStore=store(delta)
for name,delta in [('hvar-negative',-10000),('hvar-overflow',70000)]:
 f=build(name,var=True)
 t=newTable('HVAR');t.table=otTables.HVAR();t.table.Version=0x10000;t.table.VarStore=store(delta)
 # All glyphs use item zero.
 t.table.AdvWidthMap=otTables.VarIdxMap();t.table.AdvWidthMap.mapping={g:0 for g in f.getGlyphOrder()}
 t.table.LsbMap=None;t.table.RsbMap=None;f['HVAR']=t;f.save(r/(name+'.ttf'))
for name,delta in [('gpos-long',70000),('gpos-small',1000)]:
 f=build(name,'languagesystem latn dflt; feature kern { pos a <0 1 0 0>; } kern;',True)
 gdef(f,delta);s=f['GPOS'].table.LookupList.Lookup[0].SubTable[0];s.Value.YPlacement=0;s.Value.YPlaDevice=buildVarDevTable(0);s.ValueFormat|=0x20;f.save(r/(name+'.ttf'))
 # ttf-parser requires avar 1.0; read-fonts accepts major==1?
 a=newTable('avar');a.majorVersion=1;a.minorVersion=1;a.segments={'wght':{-1:-1,0:0,.5:1,1:1}};f['avar']=a;f.save(r/(name+'-avar11.ttf'))
for name,kind in [('gpos-normalization','GPOS'),('hvar-normalization','HVAR')]:
 f=build(name,'languagesystem latn dflt; feature kern { pos a <0 1 0 0>; } kern;' if kind=='GPOS' else '',True)
 vs=buildVarStore(buildVarRegionList([{'wght':(.5,8193/16384,8194/16384)}],['wght']),[buildVarData([0],[[1000]],False)])
 if kind=='GPOS':
  gdef(f,1000);f['GDEF'].table.VarStore=vs;s=f['GPOS'].table.LookupList.Lookup[0].SubTable[0];s.Value.YPlacement=0;s.Value.YPlaDevice=buildVarDevTable(0);s.ValueFormat|=0x20
 else:
  t=newTable('HVAR');t.table=otTables.HVAR();t.table.Version=0x10000;t.table.VarStore=vs;t.table.AdvWidthMap=otTables.VarIdxMap();t.table.AdvWidthMap.mapping={g:0 for g in f.getGlyphOrder()};t.table.LsbMap=None;t.table.RsbMap=None;f['HVAR']=t
 f.save(r/(name+'.ttf'))
from fontTools.ttLib import TTFont
for name,tag,ver in [('gdef14','GDEF',0x10004),('gdef11','GDEF',0x10001),('hvar11','HVAR',0x10001)]:
 f=TTFont(r/('gpos-small.ttf' if tag=='GDEF' else 'hvar-overflow.ttf'))
 f[tag].table.Version=ver
 if ver==0x10004:
  for k in ['GlyphClassDef2','AttachList2','LigCaretList2','MarkAttachClassDef2','MarkGlyphSetsDef2']:setattr(f[tag].table,k,None)
 f.save(r/(name+'.ttf'))
f=build('gpos-wide','languagesystem latn dflt; feature kern { pos a <0 1000 0 0>; } kern;')
f['GPOS'].table.LookupList.Lookup[0].SubTable[0].Format=3
f.save(r/'gpos-wide.ttf')
f=build('gsub-wide','languagesystem latn dflt; feature liga { sub a b\' by x x x x; } liga;')
for look in f['GSUB'].table.LookupList.Lookup:
 if look.LookupType==2:
  for s in look.SubTable:s.Format=2
f.save(r/'gsub-wide.ttf')
import types
f=build('gsub-wide','languagesystem latn dflt; feature liga { sub a b\' by x x x x; } liga;')
for look in f['GSUB'].table.LookupList.Lookup:
 if look.LookupType==2:
  for s in look.SubTable:
   def prewrite(s,font):
    cov=otTables.Coverage();cov.glyphs=sorted(s.mapping,key=font.getGlyphID);s.Format=2
    return {'Coverage':cov,'Sequence':[s.makeSequence_(s.mapping[g],True) for g in cov.glyphs]}
   s.preWrite=types.MethodType(prewrite,s)
f.save(r/'gsub-wide.ttf')
f=build('gpos-wide-context',"languagesystem latn dflt; lookup SHIFT { pos b <0 1000 0 0>; } SHIFT; feature kern { pos a b' lookup SHIFT; } kern;")
for look in f['GPOS'].table.LookupList.Lookup:
 if look.LookupType in [7,8]:
  for s in look.SubTable:s.Format+=3
f.save(r/'gpos-wide-context.ttf')
f=build('gpos-scalar','languagesystem latn dflt; feature kern { pos a <0 1 0 0>; } kern;',True)
gdef(f,700000);f['GDEF'].table.VarStore=buildVarStore(buildVarRegionList([{'wght':(0,.75,.75)}],['wght']),[buildVarData([0],[[700000]],False)])
s=f['GPOS'].table.LookupList.Lookup[0].SubTable[0];s.Value.YPlacement=0;s.Value.YPlaDevice=buildVarDevTable(0);s.ValueFormat|=0x20
f.save(r/'gpos-scalar.ttf')
f=build('gpos-wide-context',"languagesystem latn dflt; lookup SHIFT { pos b <0 1000 0 0>; } SHIFT; feature kern { pos a' b' lookup SHIFT; } kern;")
for look in f['GPOS'].table.LookupList.Lookup:
 if look.LookupType==7:
  for s in look.SubTable:s.Format+=3
f.save(r/'gpos-wide-context.ttf')
