use sp1_lib::{secp256k1::Secp256k1Point, utils::{AffinePoint as SP1AffinePointTrait, WeierstrassPoint}};

use crate::{
    arithmetic::{scalar::Scalar, FieldElement},
    Secp256k1,
    FieldBytes,
};
// use crate::{CompressedPoint, EncodedPoint, FieldBytes, PublicKey, Scalar};
use core::{ops::{Mul, Neg, Add, Sub, MulAssign, AddAssign, SubAssign}, iter::Sum};
use core::convert::From;
use elliptic_curve::{
    ops::LinearCombination,
    group::{prime::PrimeCurveAffine, Curve, Group, GroupEncoding},
    point::{AffineCoordinates, DecompactPoint, DecompressPoint},
    sec1::{self, FromEncodedPoint, ToEncodedPoint},
    subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption},
    zeroize::DefaultIsZeroes,
    Error, Result,
    rand_core::RngCore,
    ff::Field,
};

use alloc::vec::Vec;

pub use affine::Sp1AffinePoint;
pub use projective::Sp1ProjectivePoint;

mod affine {
    use crate::AffinePoint;

    use super::*;

    /// Elliptic curve point in affine coordinates.
    // type AffinePoint: 'static
    // + AffineCoordinates<FieldRepr = FieldBytes<Self>>
    // + Copy
    // + ConditionallySelectable
    // + ConstantTimeEq
    // + Debug
    // + Default
    // + DefaultIsZeroes
    // + Eq
    // + PartialEq
    // + Sized
    // + Send
    // + Sync;

    #[derive(Clone, Copy, Debug)]
    pub struct Sp1AffinePoint {
        pub(crate) x: FieldElement,
        pub(crate) y: FieldElement,
        pub(crate) is_infinity: u8,
    }

    impl Sp1AffinePoint {
        pub(super) fn as_zkvm_point(&self) -> Secp256k1Point {
            <Secp256k1Point as From<Sp1AffinePoint>>::from(self.clone())
        }

        pub(super) fn generator() -> Self {
            Sp1AffinePoint::from(Secp256k1Point(WeierstrassPoint::Affine(Secp256k1Point::GENERATOR)))
        }

        pub(super) const fn identity() -> Self {
            Sp1AffinePoint {
                x: FieldElement::ZERO,
                y: FieldElement::ZERO,
                is_infinity: 1,
            }
        }
    }

    impl From<Sp1AffinePoint> for Secp256k1Point {
        fn from(p: Sp1AffinePoint) -> Self {
            let mut bytes = [0u8; 64];

            // Returns the bytes in BE format.
            let mut x_bytes = p.x.to_bytes().as_slice().to_vec();
            x_bytes.reverse();

            let mut y_bytes = p.y.to_bytes().as_slice().to_vec();
            y_bytes.reverse();

            bytes[..32].copy_from_slice(&x_bytes);
            bytes[32..].copy_from_slice(&y_bytes);

            Secp256k1Point::from_le_bytes(&bytes)
        }
    }

    impl From<Secp256k1Point> for Sp1AffinePoint {
        fn from(p: Secp256k1Point) -> Self {
            let bytes = p.to_le_bytes();

            let mut x_bytes: [u8; 32] = bytes[..32].try_into().unwrap();
            x_bytes.reverse();

            let mut y_bytes: [u8; 32] = bytes[32..].try_into().unwrap();
            y_bytes.reverse();

            // Needs to be in BE format.
            Sp1AffinePoint {
                x: FieldElement::from_bytes(&x_bytes.try_into().unwrap()).unwrap(),
                y: FieldElement::from_bytes(&y_bytes.try_into().unwrap()).unwrap(),
                is_infinity: 0,
            }
        }
    }

    impl DecompressPoint<Secp256k1> for Sp1AffinePoint {
        fn decompress(x: &FieldBytes, y_is_odd: Choice) -> CtOption<Self> {
            let point: Option<AffinePoint> = AffinePoint::decompress(x, y_is_odd).into();

            // In the zkvm, were not concerned with constant time operations.
            if let Some(point) = point {
                return CtOption::new(Self {
                    x: point.x,
                    y: point.y,
                    is_infinity: point.is_identity().unwrap_u8(),
                }, Choice::from(1));
            }

            CtOption::new(Sp1AffinePoint::identity(), Choice::from(0))
        }
    }

    impl AffineCoordinates for Sp1AffinePoint {
        type FieldRepr = FieldBytes;

        fn x(&self) -> FieldBytes {
            self.x.to_bytes()
        }

        fn y_is_odd(&self) -> Choice {
            Choice::from(self.y.is_odd())
        }
    }

    impl ConditionallySelectable for Sp1AffinePoint {
        fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
            Sp1AffinePoint {
                x: FieldElement::conditional_select(&a.x, &b.x, choice),
                y: FieldElement::conditional_select(&a.y, &b.y, choice),
                is_infinity: u8::conditional_select(&a.is_infinity, &b.is_infinity, choice),
            }
        }
    }

    impl ConstantTimeEq for Sp1AffinePoint {
        fn ct_eq(&self, other: &Self) -> Choice {
            self.x.ct_eq(&other.x)
                & self.y.ct_eq(&other.y)
                & self.is_infinity.ct_eq(&other.is_infinity)
        }
    }

    impl PartialEq for Sp1AffinePoint {
        fn eq(&self, other: &Self) -> bool {
            self.ct_eq(other).into()
        }
    }

    impl Eq for Sp1AffinePoint {}

    impl Default for Sp1AffinePoint {
        fn default() -> Self {
            Sp1AffinePoint::identity()
        }
    }

    impl DefaultIsZeroes for Sp1AffinePoint {}
}

/// In our case, we actually only care about affine points.
///
/// So this type is purely to satisfy trait bounds.
mod projective {
    use elliptic_curve::ops::MulByGenerator;

    use super::*;

    // type ProjectivePoint: ConditionallySelectable
    // + ConstantTimeEq
    // + Default
    // + DefaultIsZeroes
    // + From<Self::AffinePoint>
    // + Into<Self::AffinePoint>
    // + LinearCombination
    // + MulByGenerator
    // + group::Curve<AffineRepr = Self::AffinePoint>
    // + group::Group<Scalar = Self::Scalar>;

    /// While the underlying `RustCrypto` algorithm uses projective points, we need
    /// affine points for our syscalls.
    ///
    /// So this type is merely a wrapper around the `Secp256k1Point` type.
    #[derive(Clone, Copy, Debug)]
    pub struct Sp1ProjectivePoint {
        inner: Sp1AffinePoint,
    }

    impl Sp1ProjectivePoint {
        pub(super) const fn identity() -> Self {
            Sp1ProjectivePoint {
                inner: Sp1AffinePoint::identity(),
            }
        }
    }

    impl From<Sp1AffinePoint> for Sp1ProjectivePoint {
        fn from(p: Sp1AffinePoint) -> Self {
            Sp1ProjectivePoint { inner: p }
        }
    }

    impl From<Sp1ProjectivePoint> for Sp1AffinePoint {
        fn from(p: Sp1ProjectivePoint) -> Self {
            p.inner
        }
    }

    impl Group for Sp1ProjectivePoint {
        type Scalar = Scalar;

        fn identity() -> Self {
            Sp1ProjectivePoint::identity()
        }

        fn random(rng: impl RngCore) -> Self {
            Self::generator() * Scalar::random(rng)
        }

        fn double(&self) -> Self {
            *self + self
        }

        fn generator() -> Self {
            Self {
                inner: Sp1AffinePoint::generator(),
            }
        }

        fn is_identity(&self) -> Choice {
            self.inner.is_infinity.into()
        }
    }

    impl Curve for Sp1ProjectivePoint {
        type AffineRepr = Sp1AffinePoint;

        fn to_affine(&self) -> Self::AffineRepr {
            self.inner
        }
    }

    impl MulByGenerator for Sp1ProjectivePoint {}

    impl LinearCombination for Sp1ProjectivePoint {
        fn lincomb(x: &Self, k: &Self::Scalar, y: &Self, l: &Self::Scalar) -> Self {
            let x = x.inner.as_zkvm_point();
            let y = y.inner.as_zkvm_point();

            let a_bits_le = be_bytes_to_le_bits(&k.to_bytes().as_slice().try_into().unwrap());
            let b_bits_le = be_bytes_to_le_bits(&l.to_bytes().as_slice().try_into().unwrap());

            let sp1_point = Secp256k1Point::multi_scalar_multiplication(&a_bits_le, x, &b_bits_le, y);

            Sp1ProjectivePoint {
                inner: Sp1AffinePoint::from(sp1_point),
            }
        }
    }

    impl Neg for Sp1ProjectivePoint {
        type Output = Sp1ProjectivePoint;

        fn neg(self) -> Self::Output {
            Self {
                inner: Sp1AffinePoint {
                    x: self.inner.x,
                    y: -self.inner.y,
                    is_infinity: self.inner.is_infinity,
                },
            }
        }
    }

    impl Add<Sp1ProjectivePoint> for Sp1ProjectivePoint {
        type Output = Sp1ProjectivePoint;

        fn add(self, rhs: Sp1ProjectivePoint) -> Self::Output {
            let mut sp1_point = self.inner.as_zkvm_point();
            
            sp1_point.add_assign(&rhs.inner.as_zkvm_point());

            Sp1ProjectivePoint {
                inner: Sp1AffinePoint::from(sp1_point),
            }
        }
    }

    impl Sub<Sp1ProjectivePoint> for Sp1ProjectivePoint {
        type Output = Sp1ProjectivePoint;

        fn sub(self, rhs: Sp1ProjectivePoint) -> Self::Output {
            self + rhs.neg()
        }
    }

    impl Add<&Sp1ProjectivePoint> for Sp1ProjectivePoint {
        type Output = Sp1ProjectivePoint;

        fn add(self, rhs: &Sp1ProjectivePoint) -> Self::Output {
            let mut sp1_point = self.inner.as_zkvm_point();
            
            sp1_point.add_assign(&rhs.inner.as_zkvm_point());

            Sp1ProjectivePoint {
                inner: Sp1AffinePoint::from(sp1_point),
            }
        }
    }

    impl Sub<&Sp1ProjectivePoint> for Sp1ProjectivePoint {
        type Output = Sp1ProjectivePoint;

        fn sub(self, rhs: &Sp1ProjectivePoint) -> Self::Output {
            self + (*rhs).neg()
        }
    }

    impl Mul<Scalar> for Sp1ProjectivePoint {
        type Output = Sp1ProjectivePoint;

        fn mul(self, rhs: Scalar) -> Self::Output {
            let mut sp1_point = self.inner.as_zkvm_point();
            let scalar_bytes_be = rhs.to_bytes().as_slice().to_vec();

            sp1_point.mul_assign(&be_bytes_to_le_words(scalar_bytes_be));

            Sp1ProjectivePoint {
                inner: Sp1AffinePoint::from(sp1_point),
            }
        }
    }

    impl Mul<&Scalar> for Sp1ProjectivePoint {
        type Output = Sp1ProjectivePoint;

        fn mul(self, rhs: &Scalar) -> Self::Output {
            let mut sp1_point = self.inner.as_zkvm_point();
            let scalar_bytes_be = rhs.to_bytes().as_slice().to_vec();

            sp1_point.mul_assign(&be_bytes_to_le_words(scalar_bytes_be));

            Sp1ProjectivePoint {
                inner: Sp1AffinePoint::from(sp1_point),
            }
        }
    }

    impl Sum for Sp1ProjectivePoint {
        fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
            iter.fold(Sp1ProjectivePoint::identity(), |a, b| a + b)
        }
    }
    
    impl<'a> Sum<&'a Sp1ProjectivePoint> for Sp1ProjectivePoint {
        fn sum<I: Iterator<Item = &'a Sp1ProjectivePoint>>(iter: I) -> Self {
            iter.cloned().sum()
        }
    }

    impl MulAssign<Scalar> for Sp1ProjectivePoint {
        fn mul_assign(&mut self, rhs: Scalar) {
            *self = self.mul(rhs);
        }
    }

    impl MulAssign<&Scalar> for Sp1ProjectivePoint {
        fn mul_assign(&mut self, rhs: &Scalar) {
            *self = self.mul(rhs);
        }
    }

    impl AddAssign<Sp1ProjectivePoint> for Sp1ProjectivePoint {
        fn add_assign(&mut self, rhs: Sp1ProjectivePoint) {
            *self = self.add(rhs);
        }
    }

    impl AddAssign<&Sp1ProjectivePoint> for Sp1ProjectivePoint {
        fn add_assign(&mut self, rhs: &Sp1ProjectivePoint) {
            *self = self.add(rhs);
        }
    }

    impl SubAssign<Sp1ProjectivePoint> for Sp1ProjectivePoint {
        fn sub_assign(&mut self, rhs: Sp1ProjectivePoint) {
            *self = self.sub(rhs);
        }
    }

    impl SubAssign<&Sp1ProjectivePoint> for Sp1ProjectivePoint {
        fn sub_assign(&mut self, rhs: &Sp1ProjectivePoint) {
            *self = self.sub(rhs);
        }
    }

    impl Default for Sp1ProjectivePoint {
        fn default() -> Self {
            Sp1ProjectivePoint {
                inner: Sp1AffinePoint::identity(),
            }
        }
    }

    impl Add<Sp1AffinePoint> for Sp1ProjectivePoint {
        type Output = Sp1ProjectivePoint;

        fn add(self, rhs: Sp1AffinePoint) -> Self::Output {
            self + Sp1ProjectivePoint { inner: rhs }
        }
    }

    impl Add<&Sp1AffinePoint> for Sp1ProjectivePoint {
        type Output = Sp1ProjectivePoint;

        fn add(self, rhs: &Sp1AffinePoint) -> Self::Output {
            self + Sp1ProjectivePoint { inner: *rhs }
        }
    }

    impl Sub<Sp1AffinePoint> for Sp1ProjectivePoint {
        type Output = Sp1ProjectivePoint;

        fn sub(self, rhs: Sp1AffinePoint) -> Self::Output {
            self - Sp1ProjectivePoint { inner: rhs }
        }
    }

    impl Sub<&Sp1AffinePoint> for Sp1ProjectivePoint {
        type Output = Sp1ProjectivePoint;

        fn sub(self, rhs: &Sp1AffinePoint) -> Self::Output {
            self - Sp1ProjectivePoint { inner: *rhs }
        }
    }

    impl AddAssign<Sp1AffinePoint> for Sp1ProjectivePoint {
        fn add_assign(&mut self, rhs: Sp1AffinePoint) {
            *self = self.add(rhs);
        }
    }

    impl AddAssign<&Sp1AffinePoint> for Sp1ProjectivePoint {
        fn add_assign(&mut self, rhs: &Sp1AffinePoint) {
            *self = self.add(rhs);
        }
    }

    impl SubAssign<Sp1AffinePoint> for Sp1ProjectivePoint {
        fn sub_assign(&mut self, rhs: Sp1AffinePoint) {
            *self = self.sub(rhs);
        }
    }

    impl SubAssign<&Sp1AffinePoint> for Sp1ProjectivePoint {
        fn sub_assign(&mut self, rhs: &Sp1AffinePoint) {
            *self = self.sub(rhs);
        }
    }

    impl DefaultIsZeroes for Sp1ProjectivePoint {}

    impl ConditionallySelectable for Sp1ProjectivePoint {
        fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
            Sp1ProjectivePoint {
                inner: Sp1AffinePoint::conditional_select(&a.inner, &b.inner, choice),
            }
        }
    }

    impl ConstantTimeEq for Sp1ProjectivePoint {
        fn ct_eq(&self, other: &Self) -> Choice {
            self.inner.ct_eq(&other.inner)
        }
    }

    impl PartialEq for Sp1ProjectivePoint {
        fn eq(&self, other: &Self) -> bool {
            self.ct_eq(other).into()
        }
    }

    impl Eq for Sp1ProjectivePoint {}
}

/// WARNING: The values in this type are UNTRUSTED.
///
/// The values must be constrained by the caller, by either checking the sqrt is correct and canon,
/// or checking the NQR property.
pub(crate) enum SqrtReturn {
    /// The square root was found, this is the bytes in BE.
    Found(Vec<u8>),

    /// The square root was not found.
    /// This is instead the square root of the product of
    /// a non-quadratic residue and the original value.
    NotFound(Vec<u8>),
}

/// Call the sp1 sqrt hook.
///
/// This hook takes in a field element and returns the square root of the element (with respect to the modulus).
///
/// If the element is not a quadratic residue, it returns the square root of the product of
/// the element and the nqr.
///
/// - `x`: The field element to square root.
/// - `modulus`: The modulus to square root with respect to.
/// - `nqr`: The non-quadratic residue wrt the modulus.
pub(crate) fn call_sqrt_hook(x: &[u8], modulus: &[u8], nqr: &[u8]) -> SqrtReturn {
    todo!()
}

#[inline]
fn be_bytes_to_le_words(mut bytes: Vec<u8>) -> [u32; 16] {
    bytes.reverse();

    bytes
        .chunks(4)
        .map(|b| u32::from_le_bytes(b.try_into().unwrap()))
        .collect::<Vec<_>>()
        .try_into()
        .unwrap()
}

/// Convert big-endian bytes with the most significant bit first to little-endian bytes with the least significant bit first.
#[inline]
fn be_bytes_to_le_bits(be_bytes: &[u8; 32]) -> [bool; 256] {
    let mut bits = [false; 256];
    // Reverse the byte order to little-endian.
    for (i, &byte) in be_bytes.iter().rev().enumerate() {
        for j in 0..8 {
            // Flip the bit order so the least significant bit is now the first bit of the chunk.
            bits[i * 8 + j] = ((byte >> j) & 1) == 1;
        }
    }
    bits
}