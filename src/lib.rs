//! Encoding and decoding of the base (17, 9, 5) cyclic code used by
//! [DMR](https://en.wikipedia.org/wiki/Digital_mobile_radio) and
//! [P25](https://en.wikipedia.org/wiki/Project_25).
//!
//! The generator polynomial for the base code is given by
//!
//! > g(x) = x<sup>8</sup> + x<sup>5</sup> + x<sup>4</sup> + x<sup>3</sup> + 1
//!
//! It can detect up to 4 errors or correct up to 2 errors.
//!
//! ## DMR "quadrature residue" code
//!
//! The DMR air interface extends this code to (18, 9, 6) with an extra parity check bit
//! in the LSB, then shortens it to (16, 7, 6) by deleting two MSB data bits. The extra
//! parity bit is computed over the data bits using the mask `1010111`.
//!
//! ## P25 "shortened cyclic" code
//!
//! The P25 air interface shortens this code to (16, 8, 5) by deleting the MSB data bit.
//!
//! ## References
//!
//! The decoding algorithm is based on the algorithm described in Lin and Costello's
//! *Error Control Coding* (1983) and Roman's *Coding and Information Theory* (1992),
//! p345.

extern crate binfield_matrix;

use binfield_matrix::{matrix_mul, matrix_mul_systematic};

/// Encode the given 9 data bits into a 17-bit codeword.
pub fn encode(data: u16) -> u32 {
    assert_eq!(data >> 9, 0);
    matrix_mul_systematic(data, GEN)
}

/// Try to decode the given 17-bit word to the nearest codeword, correcting up to 2
/// errors.
///
/// If decoding was successful, return `Some((data, err))`, where `data` is the 9 data
/// bits and `err` is the number of corrected bits. Otherwise, return `None` to indicate
/// an unrecoverable error.
pub fn decode(word: u32) -> Option<(u16, usize)> {
    assert_eq!(word >> 17, 0);

    // Go through a full cycle of the codeword, so the data bits end up in their original
    // position.
    let (fixed, word) = (0..17).fold((Some(0), word), |(fixed, word), _| {
        let syndrome: u8 = matrix_mul(word, PAR);

        if syndrome == 0 {
            return (fixed, rotate_17(word));
        }

        match PATTERNS[syndrome as usize] {
            0 => (None, rotate_17(word)),
            pat => (Some(pat.count_ones() as usize), rotate_17(word ^ pat)),
        }
    });

    fixed.map(|err| ((word >> 8) as u16, err))
}

/// Transpose of the generator matrix, without the identity part.
const GEN: &[u16] = &[
    0b100111100,
    0b010011110,
    0b001001111,
    0b100011011,
    0b110110001,
    0b111100100,
    0b011110010,
    0b001111001,
];

/// Transpose of parity-check matrix.
///
/// This was derived in the standard way from the generator matrix.
const PAR: &[u32] = &[
    0b10011110010000000,
    0b01001111001000000,
    0b00100111100100000,
    0b10001101100010000,
    0b11011000100001000,
    0b11110010000000100,
    0b01111001000000010,
    0b00111100100000001,
];

/// Maps each 8-bit syndrome to an error pattern.
///
/// If a syndrome is invalid, the pattern is zero. Because the code is cyclic, we only
/// need to store patterns for syndromes where the LSB is set.
const PATTERNS: [u32; 256] = [
    0,
    0b00000000000000001,
    0,
    0b00000000000000011,
    0,
    0b00000000000000101,
    0, 0, 0,
    0b00000000000001001,
    0, 0, 0, 0, 0, 0, 0,
    0b00000000000010001,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0b00000000000100001,
    0, 0, 0, 0,
    0b00100000000000001,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0b00000000100000001,
    0, 0, 0, 0, 0, 0, 0, 0,
    0b00000000001000001,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0b01000000000000001,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0b00000001000000001,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0b00000000010000001,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0b00010000000000001,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0b10000000000000001,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,
    0b00001000000000001,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0b00000010000000001,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0b00000100000000001,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

/// Cyclically rotate the word right as if it was 17 bits long.
fn rotate_17(word: u32) -> u32 {
    let lsb = word & 1;
    word >> 1 | lsb << 16
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_encode() {
        assert_eq!(encode(0b000000000), 0b000000000_00000000);
        assert_eq!(encode(0b111111111), 0b111111111_11111111);
        assert_eq!(encode(0b100000001), 0b100000001_10100101);
        assert_eq!(encode(0b000001001), 0b000001001_11001000);
        assert_eq!(encode(0b000001011), 0b000001011_10111010);
        assert_eq!(encode(0b000001010), 0b000001010_10000011);
        assert_eq!(encode(0b000001000), 0b000001000_11110001);
    }

    #[test]
    fn test_decode() {
        // Exhaustively test loopback of all possible input words.
        for w in 0..1<<9 {
            assert_eq!(decode(encode(w)), Some((w, 0)));
        }

        let w = encode(0b1010101);
        assert_eq!(w, 0b1010101_00100001);
        assert_eq!(decode(w), Some((0b1010101, 0)));

        // Exhaustively test one-bit errors.
        for i in 0..17 {
            assert_eq!(decode(w ^ 1 << i), Some((0b1010101, 1)));
        }

        // Exhaustively test two-bit errors.
        for (i, j) in (0..17).zip(0..17) {
            if i != j {
                assert_eq!(decode(w ^ (1 << i) ^ (1 << j)), Some((0b1010101, 2)));
            }
        }
    }

    #[test]
    fn test_rotate_17() {
        assert_eq!(rotate_17(0b00000000000000000), 0b00000000000000000);
        assert_eq!(rotate_17(0b10000000000000000), 0b01000000000000000);
        assert_eq!(rotate_17(0b01000000000000000), 0b00100000000000000);
        assert_eq!(rotate_17(0b00100000000000000), 0b00010000000000000);
        assert_eq!(rotate_17(0b00010000000000000), 0b00001000000000000);
        assert_eq!(rotate_17(0b00001000000000000), 0b00000100000000000);
        assert_eq!(rotate_17(0b00000100000000000), 0b00000010000000000);
        assert_eq!(rotate_17(0b00000010000000000), 0b00000001000000000);
        assert_eq!(rotate_17(0b00000001000000000), 0b00000000100000000);
        assert_eq!(rotate_17(0b00000000100000000), 0b00000000010000000);
        assert_eq!(rotate_17(0b00000000010000000), 0b00000000001000000);
        assert_eq!(rotate_17(0b00000000001000000), 0b00000000000100000);
        assert_eq!(rotate_17(0b00000000000100000), 0b00000000000010000);
        assert_eq!(rotate_17(0b00000000000010000), 0b00000000000001000);
        assert_eq!(rotate_17(0b00000000000001000), 0b00000000000000100);
        assert_eq!(rotate_17(0b00000000000000100), 0b00000000000000010);
        assert_eq!(rotate_17(0b00000000000000010), 0b00000000000000001);
        assert_eq!(rotate_17(0b00000000000000001), 0b10000000000000000);
        assert_eq!(rotate_17(0b01111111111111111), 0b10111111111111111);

        let mut word = 0b11100011001010101;

        for _ in 0..17 {
            word = rotate_17(word);
        }

        assert_eq!(word, 0b11100011001010101);
    }
}
