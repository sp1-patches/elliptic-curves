//! Pure Rust implementation of group operations on secp256r1.
//!
//! Curve parameters can be found in [NIST SP 800-186] § G.1.2: Curve P-256.
//!
//! [NIST SP 800-186]: https://csrc.nist.gov/publications/detail/sp/800-186/final

pub(crate) mod field;

#[cfg(all(not(target_os = "zkvm"), feature = "hash2curve"))]
mod hash2curve;

pub(crate) mod scalar;
pub(crate) mod util;

use self::field::FieldElement;
use crate::NistP256;

#[cfg(not(target_os = "zkvm"))]
use {
    elliptic_curve::PrimeCurveArithmetic,
    primeorder::{point_arithmetic, PrimeCurveParams},
};

use elliptic_curve::CurveArithmetic;

#[cfg(not(target_os = "zkvm"))]
mod native_types {
    use super::*;

    /// Elliptic curve point in affine coordinates.
    pub type AffinePoint = primeorder::AffinePoint<NistP256>;

    /// Elliptic curve point in projective coordinates.
    pub type ProjectivePoint = primeorder::ProjectivePoint<NistP256>;

    pub type Scalar = crate::arithmetic::scalar::Scalar;
}

#[cfg(not(target_os = "zkvm"))]
pub use native_types::*;

#[cfg(target_os = "zkvm")]
mod succinct_types {
    use super::{NistP256, FieldElement};
    use elliptic_curve::{FieldBytes, subtle::CtOption};

    impl sp1_lib::ecdsa::ECDSACurve for NistP256 {
        // a = -3
        const EQUATION_A: FieldElement = FieldElement::neg(&FieldElement::from_u64(3));

        const EQUATION_B: FieldElement =
            FieldElement::from_hex("5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b");

        type FieldElement = FieldElement;

        type ScalarImpl = crate::arithmetic::scalar::Scalar;

        type SP1AffinePoint = sp1_lib::secp256r1::Secp256r1Point;
    }

    impl sp1_lib::ecdsa::Field<NistP256> for FieldElement {
        fn from_bytes(bytes: &FieldBytes<NistP256>) -> CtOption<Self> {
            FieldElement::from_bytes(bytes)
        }

        fn to_bytes(self) -> FieldBytes<NistP256> {
            FieldElement::to_bytes(self)
        }

        #[inline]
        fn normalize(self) -> Self {
            self
        }
    }

    /// Elliptic curve point in affine coordinates.
    ///
    /// For use inside the SP1 zkvm.
    pub type AffinePoint = sp1_lib::ecdsa::AffinePoint<NistP256>;

    /// Elliptic curve point in projective coordinates.
    /// 
    /// For use inside the SP1 zkvm.
    pub type ProjectivePoint = sp1_lib::ecdsa::ProjectivePoint<NistP256>;

    /// The actual scalar type used in the SP1 zkvm.
    pub type Scalar = sp1_lib::ecdsa::Scalar<NistP256>;
}

#[cfg(target_os = "zkvm")]
pub use succinct_types::*;

impl CurveArithmetic for NistP256 {
    type AffinePoint = AffinePoint;
    type ProjectivePoint = ProjectivePoint;
    type Scalar = Scalar;
}

#[cfg(not(target_os = "zkvm"))]
impl PrimeCurveArithmetic for NistP256 {
    type CurveGroup = ProjectivePoint;
}

/// Adapted from [NIST SP 800-186] § G.1.2: Curve P-256.
///
/// [NIST SP 800-186]: https://csrc.nist.gov/publications/detail/sp/800-186/final
#[cfg(not(target_os = "zkvm"))]
impl PrimeCurveParams for NistP256 {
    type FieldElement = FieldElement;
    type PointArithmetic = point_arithmetic::EquationAIsMinusThree;

    /// a = -3
    const EQUATION_A: FieldElement = FieldElement::from_u64(3).neg();

    const EQUATION_B: FieldElement =
        FieldElement::from_hex("5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b");

    /// Base point of P-256.
    ///
    /// Defined in NIST SP 800-186 § G.1.2:
    ///
    /// ```text
    /// Gₓ = 6b17d1f2 e12c4247 f8bce6e5 63a440f2 77037d81 2deb33a0 f4a13945 d898c296
    /// Gᵧ = 4fe342e2 fe1a7f9b 8ee7eb4a 7c0f9e16 2bce3357 6b315ece cbb64068 37bf51f5
    /// ```
    const GENERATOR: (FieldElement, FieldElement) = (
        FieldElement::from_hex("6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296"),
        FieldElement::from_hex("4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5"),
    );
}
