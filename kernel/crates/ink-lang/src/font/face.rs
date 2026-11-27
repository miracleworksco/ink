//! A font as ink reads it: read-fonts, the reader the shaper (harfrust) is
//! built on, at the one instance the shaper sets text with. Outlines
//! (skrifa, on the same read-fonts, its points' fractions kept as HarfBuzz
//! keeps them), advances, the character map, line metrics and variation
//! deltas all come from that instance, at its normalized coordinates, so
//! what ink draws, measures and bounds is what the shaper places.

use read_fonts::model::glyph::outline::PathSink;
use read_fonts::tables::gdef::Gdef;
use read_fonts::tables::variations::{DeltaSetIndex, ItemVariationStore};
use read_fonts::types::{F2Dot14, GlyphId, Tag};
use read_fonts::{FontData, FontRead, FontRef};
use skrifa::MetadataProvider;
use skrifa::instance::{LocationRef, Size};
use skrifa::outline::DrawSettings;
use skrifa::outline::pen::PathStyle;

/// A font's bytes and its instance at the axes set (the shaper's).
pub(crate) struct Face<'a> {
    pub font: &'a harfrust::Font,
    raw: FontRef<'a>,
}

/// Receives an outline's path: move, line, quadratic and cubic segments,
/// close; in font units.
pub(crate) trait Pen {
    fn move_to(&mut self, x: f32, y: f32);
    fn line_to(&mut self, x: f32, y: f32);
    fn quad_to(&mut self, x1: f32, y1: f32, x: f32, y: f32);
    fn curve_to(&mut self, x1: f32, y1: f32, x2: f32, y2: f32, x: f32, y: f32);
    fn close(&mut self);
}

struct Sink<'p>(&'p mut dyn Pen);

impl PathSink for Sink<'_> {
    fn move_to(&mut self, x: f32, y: f32) {
        self.0.move_to(x, y);
    }
    fn line_to(&mut self, x: f32, y: f32) {
        self.0.line_to(x, y);
    }
    fn quad_to(&mut self, x1: f32, y1: f32, x: f32, y: f32) {
        self.0.quad_to(x1, y1, x, y);
    }
    fn curve_to(&mut self, x1: f32, y1: f32, x2: f32, y2: f32, x: f32, y: f32) {
        self.0.curve_to(x1, y1, x2, y2, x, y);
    }
    fn close(&mut self) {
        self.0.close();
    }
}

impl<'a> Face<'a> {
    /// The face of `bytes` at instance `font` (made from the same bytes).
    pub fn new(bytes: &'a [u8], font: &'a harfrust::Font) -> Result<Face<'a>, String> {
        let raw = FontRef::new(bytes).map_err(|e| format!("not a font ink reads: {e}"))?;
        Ok(Face { font, raw })
    }

    pub fn units_per_em(&self) -> u16 {
        self.font.units_per_em()
    }

    pub fn number_of_glyphs(&self) -> u16 {
        self.font.num_glyphs().min(u16::MAX as u32) as u16
    }

    /// Character `c`'s glyph, as the font's character map maps it.
    pub fn glyph_index(&self, c: char) -> Option<u16> {
        self.font.charmap().map_unicode(c as u32).map(|g| g.to_u32() as u16)
    }

    /// Glyph `g`'s outline at the instance, drawn into `pen` in font units,
    /// unrounded.
    pub fn outline_glyph(&self, g: u16, pen: &mut dyn Pen) {
        if let Some(o) = self.raw.outline_glyphs().get(GlyphId::new(g as u32)) {
            let settings = DrawSettings::unhinted(Size::unscaled(), LocationRef::new(self.coords())).with_path_style(PathStyle::HarfBuzz);
            let _ = o.draw(settings, &mut Sink(pen));
        }
    }

    /// Glyph `g`'s advance at the instance, exactly (the shaper rounds it).
    pub fn advance(&self, g: u16) -> f64 {
        self.font.glyph_metrics().h_advance_exact(GlyphId::new(g as u32)).to_f64()
    }

    /// The instance's normalized coordinates.
    pub fn coords(&self) -> &[F2Dot14] {
        self.font.normalized_coords()
    }

    /// Table `tag`'s bytes.
    pub fn table(&self, tag: &[u8; 4]) -> Option<&'a [u8]> {
        self.raw.table_data(Tag::new(tag)).map(|d| d.as_bytes())
    }

    /// GDEF's item variation store, as the shaper reads it: None without a
    /// GDEF or a store, an error where it has one that cannot be read.
    pub fn var_store(&self) -> Result<Option<ItemVariationStore<'a>>, String> {
        let Some(t) = self.table(b"GDEF") else { return Ok(None) };
        let gdef = Gdef::read(FontData::new(t)).map_err(|e| format!("its GDEF table cannot be read: {e}"))?;
        gdef.item_var_store().transpose().map_err(|e| format!("its GDEF variation store cannot be read: {e}"))
    }

    /// The delta of variation index (outer, inner) at the instance, as the
    /// shaper computes it (before it rounds): 0 at the default instance.
    pub fn delta(store: &ItemVariationStore, coords: &[F2Dot14], outer: u16, inner: u16) -> f64 {
        store.compute_delta(DeltaSetIndex { outer, inner }, coords).map_or(0.0, |d| d.to_f64())
    }
}
