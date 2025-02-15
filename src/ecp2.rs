/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements.  See the NOTICE file
distributed with this work for additional information
regarding copyright ownership.  The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied.  See the License for the
specific language governing permissions and limitations
under the License.
*/

use super::big;
use super::big::Big;
use super::ecp;
use super::fp2::FP2;
use super::rom;
use crate::types::{CurvePairingType, SexticTwist, SignOfX};
use crate::std::{string::String, fmt, str::SplitWhitespace, format};

/// Elliptic Curve Point over Fp2
///
/// A projective elliptic curve point defined over Fp2.
/// (X, Y, Z)
#[derive(Clone)]
pub struct ECP2 {
    x: FP2,
    y: FP2,
    z: FP2,
}

impl PartialEq for ECP2 {
    fn eq(&self, other: &ECP2) -> bool {
        self.equals(other)
    }
}

impl Eq for ECP2 {}

impl fmt::Display for ECP2 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "ECP2: [ {}, {}, {} ]", self.x, self.y, self.z)
    }
}

impl fmt::Debug for ECP2 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "ECP2: [ {}, {}, {} ]", self.x, self.y, self.z)
    }
}

#[allow(non_snake_case)]
impl ECP2 {
    /// New
    ///
    /// Creates a new elliptic curve point at infinity: (0, 1, 0).
    #[inline(always)]
    pub fn new() -> ECP2 {
        ECP2 {
            x: FP2::new(),
            y: FP2::new_int(1),
            z: FP2::new(),
        }
    }

    /// New Fp2's
    ///
    /// Constructs from (x,y).
    /// Set to infinity if not on curve.
    #[allow(non_snake_case)]
    #[inline(always)]
    pub fn new_fp2s(x: FP2, y: FP2) -> ECP2 {
        let mut E = ECP2 {
            x,
            y,
            z: FP2::new_int(1),
        };
        E.x.norm();

        let rhs = ECP2::rhs(&E.x);
        let mut y2 = E.y.clone();
        y2.sqr();
        if !y2.equals(&rhs) {
            E.inf();
        }
        E
    }

    /// New Fp2
    ///
    /// Constructs point from x, by calculating y.
    /// Set to infinity if not on curve.
    #[inline(always)]
    pub fn new_fp2(ix: &FP2) -> ECP2 {
        let mut E = ECP2::new();
        E.x = ix.clone();
        E.y.one();
        E.z.one();
        E.x.norm();
        let mut rhs = ECP2::rhs(&E.x);
        if rhs.sqrt() {
            E.y = rhs;
        } else {
            E.inf();
        }
        return E;
    }

    /// New Porjective
    ///
    /// Constructs point from (X, Y, Z) with no gaurentee of correctness.
    /// Do not use this on untrusted input!
    #[inline(always)]
    pub fn new_projective(x: FP2, y: FP2, z: FP2) -> ECP2 {
        ECP2 { x, y, z }
    }

    /* Test this=O? */
    pub fn is_infinity(&self) -> bool {
        self.x.is_zilch() && self.z.is_zilch()
    }

    /* set self=O */
    pub fn inf(&mut self) {
        self.x.zero();
        self.y.one();
        self.z.zero();
    }

    /* set self=-self */
    pub fn neg(&mut self) {
        self.y.norm();
        self.y.neg();
        self.y.norm();
    }

    /* Conditional move of Q to self dependant on d */
    pub fn cmove(&mut self, Q: &ECP2, d: isize) {
        self.x.cmove(&Q.x, d);
        self.y.cmove(&Q.y, d);
        self.z.cmove(&Q.z, d);
    }

    /* return 1 if b==c, no branching */
    fn teq(b: i32, c: i32) -> isize {
        let mut x = b ^ c;
        x -= 1; // if x=0, x now -1
        return ((x >> 31) & 1) as isize;
    }

    /* Constant time select from pre-computed table */
    pub fn selector(&mut self, W: &[ECP2], b: i32) {
        let m = b >> 31;
        let mut babs = (b ^ m) - m;

        babs = (babs - 1) / 2;

        self.cmove(&W[0], ECP2::teq(babs, 0)); // conditional move
        self.cmove(&W[1], ECP2::teq(babs, 1));
        self.cmove(&W[2], ECP2::teq(babs, 2));
        self.cmove(&W[3], ECP2::teq(babs, 3));
        self.cmove(&W[4], ECP2::teq(babs, 4));
        self.cmove(&W[5], ECP2::teq(babs, 5));
        self.cmove(&W[6], ECP2::teq(babs, 6));
        self.cmove(&W[7], ECP2::teq(babs, 7));

        let mut MP = self.clone();
        MP.neg();
        self.cmove(&MP, (m & 1) as isize);
    }

    /* Test if P == Q */
    pub fn equals(&self, Q: &ECP2) -> bool {
        let mut a = self.x.clone();
        let mut b = Q.x.clone();

        a.mul(&Q.z);
        b.mul(&self.z);
        if !a.equals(&b) {
            return false;
        }
        a = self.getpy();
        a.mul(&Q.z);
        b = Q.getpy();
        b.mul(&self.z);
        if !a.equals(&b) {
            return false;
        }

        return true;
    }

    /* set to Affine - (x,y,z) to (x,y) */
    pub fn affine(&mut self) {
        if self.is_infinity() {
            return;
        }
        let one = FP2::new_int(1);
        if self.z.equals(&one) {
            return;
        }
        self.z.inverse();

        self.x.mul(&self.z);
        self.x.reduce();
        self.y.mul(&self.z);
        self.y.reduce();
        self.z = one.clone();
    }

    /* extract affine x as FP2 */
    pub fn getx(&self) -> FP2 {
        let mut W = self.clone();
        W.affine();
        W.x.clone()
    }

    /* extract affine y as FP2 */
    pub fn gety(&self) -> FP2 {
        let mut W = self.clone();
        W.affine();
        W.y.clone()
    }

    /* extract projective x */
    pub fn getpx(&self) -> FP2 {
        self.x.clone()
    }
    /* extract projective y */
    pub fn getpy(&self) -> FP2 {
        self.y.clone()
    }
    /* extract projective z */
    pub fn getpz(&self) -> FP2 {
        self.z.clone()
    }

    /* convert to byte array */
    pub fn to_bytes(&self, b: &mut [u8]) {
        let mut t: [u8; big::MODBYTES as usize] = [0; big::MODBYTES as usize];
        let mb = big::MODBYTES as usize;
        let mut W = self.clone();

        W.affine();
        W.x.geta().to_bytes(&mut t);
        for i in 0..mb {
            b[i] = t[i]
        }
        W.x.getb().to_bytes(&mut t);
        for i in 0..mb {
            b[i + mb] = t[i]
        }

        W.y.geta().to_bytes(&mut t);
        for i in 0..mb {
            b[i + 2 * mb] = t[i]
        }
        W.y.getb().to_bytes(&mut t);
        for i in 0..mb {
            b[i + 3 * mb] = t[i]
        }
    }

    /// From Bytes
    ///
    /// Converts byte array to point.
    /// Pancis if insufficient bytes are given.
    #[inline(always)]
    pub fn from_bytes(b: &[u8]) -> ECP2 {
        let mut t: [u8; big::MODBYTES as usize] = [0; big::MODBYTES as usize];
        let mb = big::MODBYTES as usize;

        for i in 0..mb {
            t[i] = b[i]
        }
        let ra = Big::from_bytes(&t);
        for i in 0..mb {
            t[i] = b[i + mb]
        }
        let rb = Big::from_bytes(&t);
        let rx = FP2::new_bigs(ra, rb);

        for i in 0..mb {
            t[i] = b[i + 2 * mb]
        }
        let ra = Big::from_bytes(&t);
        for i in 0..mb {
            t[i] = b[i + 3 * mb]
        }
        let rb = Big::from_bytes(&t);
        let ry = FP2::new_bigs(ra, rb);

        ECP2::new_fp2s(rx, ry)
    }

    /// To String
    ///
    /// Converts `ECP2` to a hex string.
    pub fn to_string(&self) -> String {
        let mut W = self.clone();
        W.affine();
        if W.is_infinity() {
            return String::from("infinity");
        }
        return format!("({},{})", W.x.to_string(), W.y.to_string());
    }

    /// To Hex
    ///
    /// Converts each projective (X, Y, Z) to a hex string, separated by spaces.
    pub fn to_hex(&self) -> String {
        format!(
            "{} {} {}",
            self.x.to_hex(),
            self.y.to_hex(),
            self.z.to_hex()
        )
    }

    /// From Hex Iterator
    #[inline(always)]
    pub fn from_hex_iter(iter: &mut SplitWhitespace) -> ECP2 {
        ECP2 {
            x: FP2::from_hex_iter(iter),
            y: FP2::from_hex_iter(iter),
            z: FP2::from_hex_iter(iter),
        }
    }

    /// From Hex
    #[inline(always)]
    pub fn from_hex(val: String) -> ECP2 {
        let mut iter = val.split_whitespace();
        return ECP2::from_hex_iter(&mut iter);
    }

    /* Calculate RHS of twisted curve equation x^3+B/i */
    pub fn rhs(x: &FP2) -> FP2 {
        let mut r = x.clone();
        r.sqr();
        let mut b = FP2::new_big(Big::new_ints(&rom::CURVE_B));
        if ecp::SEXTIC_TWIST == SexticTwist::DType {
            b.div_ip();
        }
        if ecp::SEXTIC_TWIST == SexticTwist::MType {
            b.norm();
            b.mul_ip();
            b.norm();
        }

        r.mul(x);
        r.add(&b);

        r.reduce();
        return r;
    }

    /* self+=self */
    pub fn dbl(&mut self) -> isize {
        let mut iy = self.y.clone();
        if ecp::SEXTIC_TWIST == SexticTwist::DType {
            iy.mul_ip();
            iy.norm();
        }

        let mut t0 = self.y.clone(); //***** Change
        t0.sqr();
        if ecp::SEXTIC_TWIST == SexticTwist::DType {
            t0.mul_ip();
        }
        let mut t1 = iy.clone();
        t1.mul(&self.z);
        let mut t2 = self.z.clone();
        t2.sqr();

        self.z = t0.clone();
        self.z.add(&t0);
        self.z.norm();
        self.z.dbl();
        self.z.dbl();
        self.z.norm();

        t2.imul(3 * rom::CURVE_B_I);
        if ecp::SEXTIC_TWIST == SexticTwist::MType {
            t2.mul_ip();
            t2.norm();
        }
        let mut x3 = t2.clone();
        x3.mul(&self.z);

        let mut y3 = t0.clone();

        y3.add(&t2);
        y3.norm();
        self.z.mul(&t1);
        t1 = t2.clone();
        t1.add(&t2);
        t2.add(&t1);
        t2.norm();
        t0.sub(&t2);
        t0.norm(); //y^2-9bz^2
        y3.mul(&t0);
        y3.add(&x3); //(y^2+3z*2)(y^2-9z^2)+3b.z^2.8y^2
        t1 = self.x.clone();
        t1.mul(&iy); //
        self.x = t0.clone();
        self.x.norm();
        self.x.mul(&t1);
        self.x.dbl(); //(y^2-9bz^2)xy2

        self.x.norm();
        self.y = y3;
        self.y.norm();

        return 1;
    }

    /* self+=Q - return 0 for add, 1 for double, -1 for O */
    pub fn add(&mut self, Q: &ECP2) -> isize {
        let b = 3 * rom::CURVE_B_I;
        let mut t0 = self.x.clone();
        t0.mul(&Q.x); // x.Q.x
        let mut t1 = self.y.clone();
        t1.mul(&Q.y); // y.Q.y

        let mut t2 = self.z.clone();
        t2.mul(&Q.z);
        let mut t3 = self.x.clone();
        t3.add(&self.y);
        t3.norm(); //t3=X1+Y1
        let mut t4 = Q.x.clone();
        t4.add(&Q.y);
        t4.norm(); //t4=X2+Y2
        t3.mul(&t4); //t3=(X1+Y1)(X2+Y2)
        t4 = t0.clone();
        t4.add(&t1); //t4=X1.X2+Y1.Y2

        t3.sub(&t4);
        t3.norm();
        if ecp::SEXTIC_TWIST == SexticTwist::DType {
            t3.mul_ip();
            t3.norm(); //t3=(X1+Y1)(X2+Y2)-(X1.X2+Y1.Y2) = X1.Y2+X2.Y1
        }
        t4 = self.getpy();
        t4.add(&self.z);
        t4.norm(); //t4=Y1+Z1
        let mut x3 = Q.y.clone();
        x3.add(&Q.z);
        x3.norm(); //x3=Y2+Z2

        t4.mul(&x3); //t4=(Y1+Z1)(Y2+Z2)
        x3 = t1.clone(); //
        x3.add(&t2); //X3=Y1.Y2+Z1.Z2

        t4.sub(&x3);
        t4.norm();
        if ecp::SEXTIC_TWIST == SexticTwist::DType {
            t4.mul_ip();
            t4.norm(); //t4=(Y1+Z1)(Y2+Z2) - (Y1.Y2+Z1.Z2) = Y1.Z2+Y2.Z1
        }
        x3 = self.getpx();
        x3.add(&self.z);
        x3.norm(); // x3=X1+Z1
        let mut y3 = Q.x.clone();
        y3.add(&Q.z);
        y3.norm(); // y3=X2+Z2
        x3.mul(&y3); // x3=(X1+Z1)(X2+Z2)
        y3 = t0.clone();
        y3.add(&t2); // y3=X1.X2+Z1+Z2
        y3.rsub(&x3);
        y3.norm(); // y3=(X1+Z1)(X2+Z2) - (X1.X2+Z1.Z2) = X1.Z2+X2.Z1

        if ecp::SEXTIC_TWIST == SexticTwist::DType {
            t0.mul_ip();
            t0.norm(); // x.Q.x
            t1.mul_ip();
            t1.norm(); // y.Q.y
        }
        x3 = t0.clone();
        x3.add(&t0);
        t0.add(&x3);
        t0.norm();
        t2.imul(b);
        if ecp::SEXTIC_TWIST == SexticTwist::MType {
            t2.mul_ip();
            t2.norm();
        }
        let mut z3 = t1.clone();
        z3.add(&t2);
        z3.norm();
        t1.sub(&t2);
        t1.norm();
        y3.imul(b);
        if ecp::SEXTIC_TWIST == SexticTwist::MType {
            y3.mul_ip();
            y3.norm();
        }
        x3 = y3.clone();
        x3.mul(&t4);
        t2 = t3.clone();
        t2.mul(&t1);
        x3.rsub(&t2);
        y3.mul(&t0);
        t1.mul(&z3);
        y3.add(&t1);
        t0.mul(&t3);
        z3.mul(&t4);
        z3.add(&t0);

        self.x = x3;
        self.x.norm();
        self.y = y3;
        self.y.norm();
        self.z = z3;
        self.z.norm();

        return 0;
    }

    /* set this-=Q */
    pub fn sub(&mut self, Q: &ECP2) -> isize {
        let mut NQ = Q.clone();
        NQ.neg();
        let d = self.add(&NQ);
        return d;
    }

    /* set this*=q, where q is Modulus, using Frobenius */
    pub fn frob(&mut self, x: &FP2) {
        let mut x2 = x.clone();
        x2.sqr();
        self.x.conj();
        self.y.conj();
        self.z.conj();
        self.z.reduce();
        self.x.mul(&x2);
        self.y.mul(&x2);
        self.y.mul(x);
    }

    /// Multiplication
    ///
    /// Return e * self
    #[inline(always)]
    pub fn mul(&self, e: &Big) -> ECP2 {
        if self.is_infinity() {
            return ECP2::new();
        }

        let mut W: [ECP2; 8] = [
            ECP2::new(),
            ECP2::new(),
            ECP2::new(),
            ECP2::new(),
            ECP2::new(),
            ECP2::new(),
            ECP2::new(),
            ECP2::new(),
        ];

        const CT: usize = 1 + (big::NLEN * (big::BASEBITS as usize) + 3) / 4;
        let mut w: [i8; CT] = [0; CT];

        /* precompute table */
        let mut Q = self.clone();
        Q.dbl();

        W[0] = self.clone();

        for i in 1..8 {
            W[i] = W[i - 1].clone();
            W[i].add(&Q);
        }

        /* make exponent odd - add 2P if even, P if odd */
        let mut t = e.clone();
        let s = t.parity();
        t.inc(1);
        t.norm();
        let ns = t.parity();
        let mut mt = t.clone();
        mt.inc(1);
        mt.norm();
        t.cmove(&mt, s);
        Q.cmove(&self, ns);
        let C = Q.clone();

        let nb = 1 + (t.nbits() + 3) / 4;

        /* convert exponent to signed 4-bit window */
        for i in 0..nb {
            w[i] = (t.lastbits(5) - 16) as i8;
            t.dec(w[i] as isize);
            t.norm();
            t.fshr(4);
        }
        w[nb] = (t.lastbits(5)) as i8;

        let mut P = W[((w[nb] as usize) - 1) / 2].clone();
        for i in (0..nb).rev() {
            Q.selector(&W, w[i] as i32);
            P.dbl();
            P.dbl();
            P.dbl();
            P.dbl();
            P.add(&Q);
        }
        P.sub(&C);
        P.affine();
        P
    }

    /// Multiply 4 Points
    ///
    /// P = u0 * Q0 + u1 * Q1 + u2 * Q2 + u3 * Q3
    /// Bos & Costello https://eprint.iacr.org/2013/458.pdf
    /// Faz-Hernandez & Longa & Sanchez  https://eprint.iacr.org/2013/158.pdf
    /// Side channel attack secure
    ///
    /// Panics if 4 points and 4 scalars are not given.
    #[inline(always)]
    pub fn mul4(Q: &mut [ECP2], u: &[Big]) -> ECP2 {
        let mut P = ECP2::new();

        let mut T: [ECP2; 8] = [
            ECP2::new(),
            ECP2::new(),
            ECP2::new(),
            ECP2::new(),
            ECP2::new(),
            ECP2::new(),
            ECP2::new(),
            ECP2::new(),
        ];

        let mut mt = Big::new();

        let mut t: [Big; 4] = [u[0].clone(), u[1].clone(), u[2].clone(), u[3].clone()];

        const CT: usize = 1 + big::NLEN * (big::BASEBITS as usize);
        let mut w: [i8; CT] = [0; CT];
        let mut s: [i8; CT] = [0; CT];

        for i in 0..4 {
            t[i].norm();
        }

        T[0] = Q[0].clone();
        let mut W = T[0].clone();
        T[1] = W.clone();
        T[1].add(&Q[1]); // Q[0]+Q[1]
        T[2] = W.clone();
        T[2].add(&Q[2]);
        W = T[1].clone(); // Q[0]+Q[2]
        T[3] = W.clone();
        T[3].add(&Q[2]);
        W = T[0].clone(); // Q[0]+Q[1]+Q[2]
        T[4] = W.clone();
        T[4].add(&Q[3]);
        W = T[1].clone(); // Q[0]+Q[3]
        T[5] = W.clone();
        T[5].add(&Q[3]);
        W = T[2].clone(); // Q[0]+Q[1]+Q[3]
        T[6] = W.clone();
        T[6].add(&Q[3]);
        W = T[3].clone(); // Q[0]+Q[2]+Q[3]
        T[7] = W.clone();
        T[7].add(&Q[3]); // Q[0]+Q[1]+Q[2]+Q[3]

        // Make it odd
        let pb = 1 - t[0].parity();
        t[0].inc(pb);
        t[0].norm();

        // Number of bits
        mt.zero();
        for i in 0..4 {
            mt.or(&t[i]);
        }

        let nb = 1 + mt.nbits();

        // Sign pivot

        s[nb - 1] = 1;
        for i in 0..nb - 1 {
            t[0].fshr(1);
            s[i] = (2 * t[0].parity() - 1) as i8;
        }

        // Recoded exponent
        for i in 0..nb {
            w[i] = 0;
            let mut k = 1;
            for j in 1..4 {
                let bt = s[i] * (t[j].parity() as i8);
                t[j].fshr(1);
                t[j].dec((bt >> 1) as isize);
                t[j].norm();
                w[i] += bt * (k as i8);
                k = 2 * k;
            }
        }

        // Main loop
        P.selector(&T, (2 * w[nb - 1] + 1) as i32);
        for i in (0..nb - 1).rev() {
            P.dbl();
            W.selector(&T, (2 * w[i] + s[i]) as i32);
            P.add(&W);
        }

        // apply correction
        W = P.clone();
        W.sub(&Q[0]);
        P.cmove(&W, pb);
        P.affine();

        return P;
    }

    /// Map It
    ///
    /// Maps bytes to a curve point using hash and test.
    /// Not conformant to hash-to-curve standards.
    #[allow(non_snake_case)]
    #[inline(always)]
    pub fn mapit(h: &[u8]) -> ECP2 {
        let q = Big::new_ints(&rom::MODULUS);
        let mut x = Big::from_bytes(h);
        x.rmod(&q);
        let mut Q: ECP2;
        let one = Big::new_int(1);

        loop {
            let X = FP2::new_bigs(one.clone(), x.clone());
            Q = ECP2::new_fp2(&X);
            if !Q.is_infinity() {
                break;
            }
            x.inc(1);
            x.norm();
        }
        Q.clear_cofactor();
        Q
    }

    pub fn clear_cofactor(&mut self) {
        let mut X = FP2::new_bigs(Big::new_ints(&rom::FRA), Big::new_ints(&rom::FRB));
        if ecp::SEXTIC_TWIST == SexticTwist::MType {
            X.inverse();
            X.norm();
        }
        let x = Big::new_ints(&rom::CURVE_BNX);

        if ecp::CURVE_PAIRING_TYPE == CurvePairingType::Bn {
            let mut T = self.mul(&x);
            if ecp::SIGN_OF_X == SignOfX::NegativeX {
                T.neg();
            }
            let mut K = T.clone();
            K.dbl();
            K.add(&T);

            K.frob(&X);
            self.frob(&X);
            self.frob(&X);
            self.frob(&X);
            self.add(&T);
            self.add(&K);
            T.frob(&X);
            T.frob(&X);
            self.add(&T);
        }
        if ecp::CURVE_PAIRING_TYPE == CurvePairingType::Bls {
            let mut xQ = self.mul(&x);
            let mut x2Q = xQ.mul(&x);

            if ecp::SIGN_OF_X == SignOfX::NegativeX {
                xQ.neg();
            }
            x2Q.sub(&xQ);
            x2Q.sub(&self);

            xQ.sub(&self);
            xQ.frob(&X);

            self.dbl();
            self.frob(&X);
            self.frob(&X);

            self.add(&x2Q);
            self.add(&xQ);
        }

        self.affine();
    }

    /// Generator
    ///
    /// Returns the generator of the group.
    #[inline(always)]
    pub fn generator() -> ECP2 {
        return ECP2::new_fp2s(
            FP2::new_bigs(
                Big::new_ints(&rom::CURVE_PXA),
                Big::new_ints(&rom::CURVE_PXB),
            ),
            FP2::new_bigs(
                Big::new_ints(&rom::CURVE_PYA),
                Big::new_ints(&rom::CURVE_PYB),
            ),
        );
    }
}
