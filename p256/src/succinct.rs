use sp1_lib::{secp256r1::Secp256r1Point, utils::{AffinePoint as SP1AffinePointTrait, WeierstrassPoint}};

use crate::{
    arithmetic::{scalar::Scalar, field::FieldElement},
    NistP256,
    FieldBytes,
    CompressedPoint,
};
// use crate::{CompressedPoint, EncodedPoint, FieldBytes, PublicKey, Scalar};
use core::{
    ops::{Mul, Add, Sub, MulAssign, AddAssign, SubAssign, Neg}, 
    iter::Sum
};
use core::convert::From;
use elliptic_curve::{
    ops::LinearCombination,
    group::{Curve, Group, GroupEncoding},
    point::{AffineCoordinates, DecompactPoint, DecompressPoint},
    sec1::{self, FromEncodedPoint, ToEncodedPoint},
    subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption},
    zeroize::DefaultIsZeroes,
    rand_core::RngCore,
    ff::Field,
};

use alloc::vec::Vec;

pub use affine::Sp1AffinePoint;
pub use projective::Sp1ProjectivePoint;

// a = -3
const EQUATION_A: FieldElement = FieldElement::neg(&FieldElement::from_u64(3));

const EQUATION_B: FieldElement =
    FieldElement::from_hex("5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b");

#[allow(missing_docs)]
mod affine {
    use sp1_lib::utils::WeierstrassAffinePoint;

    use crate::EncodedPoint;

    use super::*;

    #[derive(Clone, Copy)]
    pub struct Sp1AffinePoint {
        pub(crate) point: Secp256r1Point
    }

    impl core::fmt::Debug for Sp1AffinePoint {
        fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
            write!(f, "Sp1AffinePoint {{ point }}")
        }
    }

    impl Sp1AffinePoint {
        pub(crate) fn from_field_elements_unchecked(x: FieldElement, y: FieldElement) -> Self {
            let mut x_slice = x.to_bytes();
            let x_slice = x_slice.as_mut_slice();
            x_slice.reverse();

            let mut y_slice = y.to_bytes();
            let y_slice = y_slice.as_mut_slice();
            y_slice.reverse();

            Sp1AffinePoint {
                point: <Secp256r1Point as SP1AffinePointTrait<16>>::from(&x_slice, &y_slice)
            }
        }

        /// # Panics 
        /// - if we have non canon represntations of the field elements.
        pub(crate) fn field_elements(&self) -> (FieldElement, FieldElement) {
            if self.is_identity().into() {
                return (FieldElement::ZERO, FieldElement::ZERO);
            }

            let bytes = self.point.to_le_bytes();
            
            let mut x_bytes: [u8; 32] = bytes[..32].try_into().unwrap();
            x_bytes.reverse();

            let mut y_bytes: [u8; 32] = bytes[32..].try_into().unwrap();
            y_bytes.reverse();

            let x = FieldElement::from_bytes(&x_bytes.try_into().unwrap()).unwrap();
            let y = FieldElement::from_bytes(&y_bytes.try_into().unwrap()).unwrap();
            (x, y)
        }

        pub(super) const fn generator() -> Self {
            Sp1AffinePoint {
                point: Secp256r1Point(WeierstrassPoint::Affine(Secp256r1Point::GENERATOR)),
            }
        }

        pub(super) const fn identity() -> Self {
            Sp1AffinePoint {
                point: Secp256r1Point(WeierstrassPoint::Infinity),
            }
        }

        pub(crate) fn is_identity(&self) -> Choice {
            Choice::from(self.point.is_infinity() as u8)
        }
    }

    impl From<Sp1AffinePoint> for Secp256r1Point {
        fn from(p: Sp1AffinePoint) -> Self {
            p.point
        }
    }

    impl From<Secp256r1Point> for Sp1AffinePoint {
        fn from(p: Secp256r1Point) -> Self {
            Sp1AffinePoint { point: p }
        }
    }

    impl FromEncodedPoint<NistP256> for Sp1AffinePoint {
        fn from_encoded_point(point: &EncodedPoint) -> CtOption<Self> {
            match point.coordinates() {
                sec1::Coordinates::Identity => CtOption::new(Self::identity(), 1.into()),
                sec1::Coordinates::Compact { x } => Self::decompact(x),
                sec1::Coordinates::Compressed { x, y_is_odd } => {
                    Sp1AffinePoint::decompress(x, Choice::from(y_is_odd as u8))
                }
                sec1::Coordinates::Uncompressed { x, y } => {
                    let x = FieldElement::from_bytes(x);
                    let y = FieldElement::from_bytes(y);

                    x.and_then(|x| {
                        y.and_then(|y| {
                            // Check that the point is on the curve
                            let lhs = (y * &y).neg();
                            let rhs = (&x * &x * &x) + (EQUATION_A * &x) + &EQUATION_B;

                            let point = Self::from_field_elements_unchecked(x, y);

                            CtOption::new(point, (lhs + &rhs).is_zero())
                        })
                    })
                }
            }
        }
    }

    impl ToEncodedPoint<NistP256> for Sp1AffinePoint {
        fn to_encoded_point(&self, compress: bool) -> EncodedPoint {
            // If the point is the identity point, we can just return the identity point.
            if self.is_identity().into() {
                return EncodedPoint::identity();
            }

            let (x, y) = self.field_elements();

            EncodedPoint::from_affine_coordinates(&x.to_bytes(), &y.to_bytes(), compress)
        }
    }

    impl DecompressPoint<NistP256> for Sp1AffinePoint {
        fn decompress(x_bytes: &FieldBytes, y_is_odd: Choice) -> CtOption<Self> {
            FieldElement::from_bytes(x_bytes).and_then(|x| {
                let alpha =  (x * &x * &x) + (EQUATION_A * &x) + &EQUATION_B;
                let beta = alpha.sqrt();
    
                beta.map(|beta| {
                    let y = FieldElement::conditional_select(
                        &beta.neg(),
                        &beta,
                        beta.is_odd().ct_eq(&y_is_odd),
                    );
    
                    Sp1AffinePoint::from_field_elements_unchecked(x, y)
                })
            })
        }
    }

    impl DecompactPoint<NistP256> for Sp1AffinePoint {
        fn decompact(x_bytes: &FieldBytes) -> CtOption<Self> {
            Self::decompress(x_bytes, Choice::from(0))
        }
    }
    

    impl AffineCoordinates for Sp1AffinePoint {
        type FieldRepr = FieldBytes;

        fn x(&self) -> FieldBytes {
            let (x, _) = self.field_elements();

            x.to_bytes()
        }

        fn y_is_odd(&self) -> Choice {
            let (_, y) = self.field_elements();

            Choice::from(y.is_odd())
        }
    }

    impl ConditionallySelectable for Sp1AffinePoint {
        fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
            // In the vm, we dont care about constant time selection.
            if choice.into() {
                *b
            } else {
                *a
            }
        }
    }

    impl ConstantTimeEq for Sp1AffinePoint {
        fn ct_eq(&self, other: &Self) -> Choice {
            let (x1, y1) = self.field_elements();
            let (x1, y1) = (x1, y1);

            let (x2, y2) = other.field_elements();
            let (x2, y2) = (x2, y2);

            x1.ct_eq(&x2) & y1.ct_eq(&y2)   
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

    impl GroupEncoding for Sp1AffinePoint {
        type Repr = CompressedPoint;
    
        fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
            EncodedPoint::from_bytes(bytes)
                .map(|point| CtOption::new(point, Choice::from(1)))
                .unwrap_or_else(|_| {
                    // SEC1 identity encoding is technically 1-byte 0x00, but the
                    // `GroupEncoding` API requires a fixed-width `Repr`
                    let is_identity = bytes.ct_eq(&Self::Repr::default());
                    CtOption::new(EncodedPoint::identity(), is_identity)
                })
                .and_then(|point| Self::from_encoded_point(&point))
        }
    
        fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
            // No unchecked conversion possible for compressed points
            Self::from_bytes(bytes)
        }
    
        fn to_bytes(&self) -> Self::Repr {
            let encoded = self.to_encoded_point(true);
            let mut result = CompressedPoint::default();
            result[..encoded.len()].copy_from_slice(encoded.as_bytes());
            result
        }
    }
}

/// In our case, we actually only care about affine points.
///
/// So this type is purely to satisfy trait bounds.
mod projective {
    use elliptic_curve::{group::{cofactor::CofactorGroup, prime::PrimeGroup}, ops::MulByGenerator};

    use super::*;

    /// While the underlying `RustCrypto` algorithm uses projective points, we need
    /// affine points for our syscalls.
    ///
    /// So this type is merely a wrapper around the `Secp256r1Point` type.
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

        pub(crate) fn to_affine(self) -> Sp1AffinePoint {
            self.inner
        }

        fn zkvm_point(&self) -> Secp256r1Point {
            self.inner.point
        }
    }

    impl From<Sp1AffinePoint> for Sp1ProjectivePoint {
        fn from(p: Sp1AffinePoint) -> Self {
            Sp1ProjectivePoint { inner: p }
        }
    }

    impl From<Secp256r1Point> for Sp1ProjectivePoint {
        fn from(p: Secp256r1Point) -> Self {
            Sp1ProjectivePoint { inner: Sp1AffinePoint::from(p) }
        }
    }

    impl From<&Sp1AffinePoint> for Sp1ProjectivePoint {
        fn from(p: &Sp1AffinePoint) -> Self {
            Sp1ProjectivePoint { inner: *p }
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
            self.inner.is_identity()
        }
    }

    /// Technically this can also implement [`elliptic_curve::PrimeCurve`], but we dont need it.
    impl Curve for Sp1ProjectivePoint {
        type AffineRepr = Sp1AffinePoint;

        fn to_affine(&self) -> Self::AffineRepr {
            self.inner
        }
    }


    impl MulByGenerator for Sp1ProjectivePoint {}

    impl LinearCombination for Sp1ProjectivePoint {
        fn lincomb(x: &Self, k: &Self::Scalar, y: &Self, l: &Self::Scalar) -> Self {
            let x = x.zkvm_point();
            let y = y.zkvm_point();

            let a_bits_le = be_bytes_to_le_bits(k.to_bytes().as_slice());
            let b_bits_le = be_bytes_to_le_bits(&l.to_bytes().as_slice());

            let sp1_point = Secp256r1Point::multi_scalar_multiplication(&a_bits_le, x, &b_bits_le, y);

            Sp1ProjectivePoint {
                inner: Sp1AffinePoint::from(sp1_point),
            }
        }
    }

    impl Neg for Sp1ProjectivePoint {
        type Output = Sp1ProjectivePoint;

        fn neg(self) -> Self::Output {
            if self.is_identity().into() {
                return self;
            }

            let point = self.to_affine();
            let (x, y) = point.field_elements();

            Sp1AffinePoint::from_field_elements_unchecked(x, y.neg()).into()
        }
    }

    impl Add<Sp1ProjectivePoint> for Sp1ProjectivePoint {
        type Output = Sp1ProjectivePoint;

        fn add(self, rhs: Sp1ProjectivePoint) -> Self::Output {
            let mut sp1_point = self.zkvm_point();
            
            sp1_point.add_assign(&rhs.zkvm_point());

            sp1_point.into()
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
            let mut sp1_point = self.zkvm_point();
            
            sp1_point.add_assign(&rhs.zkvm_point());

            sp1_point.into()
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
            let mut sp1_point = self.zkvm_point();
            let mut scalar_bytes_be = rhs.to_bytes();

            sp1_point.mul_assign(&be_bytes_to_le_words(scalar_bytes_be.as_mut_slice()));

            sp1_point.into()
        }
    }

    impl Mul<&Scalar> for Sp1ProjectivePoint {
        type Output = Sp1ProjectivePoint;

        fn mul(self, rhs: &Scalar) -> Self::Output {
            let mut sp1_point = self.zkvm_point();
            let mut scalar_bytes_be = rhs.to_bytes();
            
            sp1_point.mul_assign(&be_bytes_to_le_words(scalar_bytes_be.as_mut_slice()));

            sp1_point.into()
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
            self.inner.point.mul_assign(&be_bytes_to_le_words(rhs.to_bytes().as_mut_slice()));
        }
    }

    impl MulAssign<&Scalar> for Sp1ProjectivePoint {
        fn mul_assign(&mut self, rhs: &Scalar) {
            self.inner.point.mul_assign(&be_bytes_to_le_words(rhs.to_bytes().as_mut_slice()));
        }
    }

    impl AddAssign<Sp1ProjectivePoint> for Sp1ProjectivePoint {
        fn add_assign(&mut self, rhs: Sp1ProjectivePoint) {
            self.inner.point.add_assign(&rhs.inner.point);
        }
    }
        
    impl AddAssign<&Sp1ProjectivePoint> for Sp1ProjectivePoint {
        fn add_assign(&mut self, rhs: &Sp1ProjectivePoint) {
            self.inner.point.add_assign(&rhs.inner.point);
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
            self.inner.point.add_assign(&rhs.point);
        }
    }

    impl AddAssign<&Sp1AffinePoint> for Sp1ProjectivePoint {
        fn add_assign(&mut self, rhs: &Sp1AffinePoint) {
            self.inner.point.add_assign(&rhs.point);
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

    // Traits for hash2curve
    impl GroupEncoding for Sp1ProjectivePoint {
        type Repr = CompressedPoint;
    
        fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
            <Sp1AffinePoint as GroupEncoding>::from_bytes(bytes).map(Into::into)
        }
    
        fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
            // No unchecked conversion possible for compressed points
            Self::from_bytes(bytes)
        }
    
        fn to_bytes(&self) -> Self::Repr {
            self.inner.to_bytes()
        }
    }

    impl PrimeGroup for Sp1ProjectivePoint {}

    /// The curve has prime order, so the cofactor is 1.
    impl CofactorGroup for Sp1ProjectivePoint {
        type Subgroup = Self;

        fn clear_cofactor(&self) -> Self {
            *self
        }

        fn into_subgroup(self) -> CtOption<Self> {
            CtOption::new(self, Choice::from(1))
        }

        fn is_torsion_free(&self) -> Choice {
            Choice::from(1)
        }
    }
}

/// Call the sp1 sqrt hook.
///
/// This hook takes in a field element and returns the square root of the element (with respect to the modulus).
///
/// If the element is not a quadratic residue, it returns the square root of the product of
/// the element and the nqr.
///
/// - `x`: The field element to square root.
/// - `modulus`: The HEX encoded modulus to square root with respect to.
/// - `nqr`: The non-quadratic residue wrt the modulus.
pub(crate) fn call_sqrt_hook(x: &[u8], modulus: &'static str, nqr: &[u8]) -> (u8, Vec<u8>) {
    let mut buf = Vec::new();
    buf.extend_from_slice(&32_u32.to_be_bytes());
    buf.extend_from_slice(x);
    buf.extend_from_slice(&hex::decode(modulus).unwrap());
    buf.extend_from_slice(nqr);

    sp1_lib::unconstrained! {
        sp1_lib::io::write(
            sp1_lib::io::FD_FP_SQRT,
            buf.as_slice()
        );
    }

    let status: u8 = sp1_lib::io::read_vec().first().copied().expect("sqrt hook should have a status");
    let result = sp1_lib::io::read_vec();

    (status, result)
}

/// Call the sp1 inverse hook.
///
/// This hook takes in a field element and returns the inverse of the element (with respect to the modulus).
///
/// - `x`: The field element to inverse.
/// - `modulus`: The HEX encoded modulus to inverse with respect to.
pub(crate) fn call_inv_hook(x: &[u8], modulus: &'static str) -> Vec<u8> {
    let mut buf = Vec::new();
    buf.extend_from_slice(&32_u32.to_be_bytes());
    buf.extend_from_slice(x);
    buf.extend_from_slice(&hex::decode(modulus).unwrap());

    sp1_lib::unconstrained! {
        sp1_lib::io::write(sp1_lib::io::FD_FP_INV, buf.as_slice());
    }

    sp1_lib::io::read_vec()
}

#[inline]
fn be_bytes_to_le_words(bytes: &mut [u8]) -> [u32; 16] {
    bytes.reverse();

    bytes
        .chunks(4)
        .map(|b| u32::from_le_bytes(b.try_into().unwrap()))
        .collect::<Vec<_>>()
        .try_into()
        .unwrap()
}

/// Convert big-endian bytes with the most significant bit first to little-endian bytes with the least significant bit first.
/// Panics: If the bytes have len > 32.
#[inline]
fn be_bytes_to_le_bits(be_bytes: &[u8]) -> [bool; 256] {
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