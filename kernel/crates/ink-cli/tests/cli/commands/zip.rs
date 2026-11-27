//! A zip archive read from its central directory, each entry inflated
//! (RFC 1951) and its CRC-32 checked: the 3MF container, read without
//! ink-io's writer.

/// Each entry's name and contents; panics on a bad CRC or a method other
/// than stored or deflate.
pub fn entries(data: &[u8]) -> Vec<(String, Vec<u8>)> {
    let u16_at = |o: usize| u16::from_le_bytes(data[o..o + 2].try_into().unwrap()) as usize;
    let u32_at = |o: usize| u32::from_le_bytes(data[o..o + 4].try_into().unwrap()) as usize;
    let end = (0..=data.len() - 22).rev().find(|&o| u32_at(o) == 0x0605_4b50).expect("no end of central directory");
    let (n, mut cd) = (u16_at(end + 10), u32_at(end + 16));
    let mut out = vec![];
    for _ in 0..n {
        assert_eq!(u32_at(cd), 0x0201_4b50, "central directory entry");
        let (method, crc, csize, usize_) = (u16_at(cd + 10), u32_at(cd + 16), u32_at(cd + 20), u32_at(cd + 24));
        let (nlen, xlen, clen, off) = (u16_at(cd + 28), u16_at(cd + 30), u16_at(cd + 32), u32_at(cd + 42));
        let name = String::from_utf8(data[cd + 46..cd + 46 + nlen].to_vec()).unwrap();
        assert_eq!(u32_at(off), 0x0403_4b50, "{name}: local header");
        let start = off + 30 + u16_at(off + 26) + u16_at(off + 28);
        let raw = &data[start..start + csize];
        let body = match method {
            0 => raw.to_vec(),
            8 => inflate(raw),
            m => panic!("{name}: compression method {m}"),
        };
        assert_eq!(body.len(), usize_, "{name}: size");
        assert_eq!(crc32(&body) as usize, crc, "{name}: CRC");
        out.push((name, body));
        cd += 46 + nlen + xlen + clen;
    }
    out
}

pub fn crc32(data: &[u8]) -> u32 {
    let mut c = !0u32;
    for &b in data {
        c ^= b as u32;
        for _ in 0..8 {
            c = if c & 1 != 0 { 0xEDB8_8320 ^ (c >> 1) } else { c >> 1 };
        }
    }
    !c
}

struct Bits<'a> {
    data: &'a [u8],
    pos: usize,
}

impl Bits<'_> {
    fn bit(&mut self) -> u32 {
        let b = (self.data[self.pos >> 3] >> (self.pos & 7)) & 1;
        self.pos += 1;
        b as u32
    }
    fn bits(&mut self, n: u32) -> u32 {
        (0..n).fold(0, |v, k| v | (self.bit() << k))
    }
}

/// A canonical Huffman code: counts per length and symbols by code.
struct Huff {
    counts: [u16; 16],
    symbols: Vec<u16>,
}

impl Huff {
    fn new(lengths: &[u8]) -> Huff {
        let mut counts = [0u16; 16];
        lengths.iter().for_each(|&l| counts[l as usize] += 1);
        counts[0] = 0;
        let mut offs = [0u16; 16];
        for l in 1..16 {
            offs[l] = offs[l - 1] + counts[l - 1];
        }
        let mut symbols = vec![0; lengths.len()];
        for (s, &l) in lengths.iter().enumerate() {
            if l != 0 {
                symbols[offs[l as usize] as usize] = s as u16;
                offs[l as usize] += 1;
            }
        }
        Huff { counts, symbols }
    }
    fn decode(&self, b: &mut Bits) -> usize {
        let (mut code, mut first, mut index) = (0i32, 0i32, 0i32);
        for l in 1..16 {
            code |= b.bit() as i32;
            let count = self.counts[l] as i32;
            if code - count < first {
                return self.symbols[(index + code - first) as usize] as usize;
            }
            index += count;
            first = (first + count) << 1;
            code <<= 1;
        }
        panic!("bad Huffman code");
    }
}

const LEN_BASE: [usize; 29] = [3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258];
const LEN_EXTRA: [u32; 29] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0];
const DIST_BASE: [usize; 30] = [1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577];
const DIST_EXTRA: [u32; 30] = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13];

/// Raw deflate, every block kind.
pub fn inflate(data: &[u8]) -> Vec<u8> {
    let mut b = Bits { data, pos: 0 };
    let mut out = Vec::new();
    loop {
        let last = b.bit();
        match b.bits(2) {
            0 => {
                b.pos = (b.pos + 7) & !7;
                let p = b.pos >> 3;
                let len = u16::from_le_bytes([data[p], data[p + 1]]) as usize;
                assert_eq!(len as u16, !u16::from_le_bytes([data[p + 2], data[p + 3]]), "stored block length");
                out.extend_from_slice(&data[p + 4..p + 4 + len]);
                b.pos = (p + 4 + len) << 3;
            }
            kind @ (1 | 2) => {
                let (lit, dist) = if kind == 1 {
                    let mut l = [8u8; 288];
                    l[144..256].fill(9);
                    l[256..280].fill(7);
                    (Huff::new(&l), Huff::new(&[5; 30]))
                } else {
                    let (hlit, hdist, hclen) = (b.bits(5) as usize + 257, b.bits(5) as usize + 1, b.bits(4) as usize + 4);
                    const ORDER: [usize; 19] = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];
                    let mut cl = [0u8; 19];
                    for &o in &ORDER[..hclen] {
                        cl[o] = b.bits(3) as u8;
                    }
                    let clh = Huff::new(&cl);
                    let mut lens = Vec::with_capacity(hlit + hdist);
                    while lens.len() < hlit + hdist {
                        match clh.decode(&mut b) {
                            s @ 0..=15 => lens.push(s as u8),
                            16 => {
                                let prev = *lens.last().expect("repeat with no length before");
                                let n = 3 + b.bits(2);
                                lens.extend((0..n).map(|_| prev));
                            }
                            17 => {
                                let n = 3 + b.bits(3);
                                lens.extend((0..n).map(|_| 0));
                            }
                            _ => {
                                let n = 11 + b.bits(7);
                                lens.extend((0..n).map(|_| 0));
                            }
                        }
                    }
                    (Huff::new(&lens[..hlit]), Huff::new(&lens[hlit..]))
                };
                loop {
                    let s = lit.decode(&mut b);
                    if s < 256 {
                        out.push(s as u8);
                    } else if s == 256 {
                        break;
                    } else {
                        let len = LEN_BASE[s - 257] + b.bits(LEN_EXTRA[s - 257]) as usize;
                        let d = dist.decode(&mut b);
                        let back = DIST_BASE[d] + b.bits(DIST_EXTRA[d]) as usize;
                        let from = out.len() - back;
                        for k in 0..len {
                            out.push(out[from + k]);
                        }
                    }
                }
            }
            k => panic!("deflate block type {k}"),
        }
        if last == 1 {
            return out;
        }
    }
}
