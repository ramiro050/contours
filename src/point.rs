use std::ops::{Add, Sub};

#[derive(PartialEq, Copy, Clone, Debug)]
pub struct Point(pub usize, pub usize);
impl Point {
    pub fn origin() -> Point {
        Point(0, 0)
    }

    pub fn right(&self) -> Point {
        Point(self.0 + 1, self.1)
    }

    pub fn left(&self) -> Point {
        Point(self.0 - 1, self.1)
    }

    pub fn up(&self) -> Point {
        Point(self.0, self.1 + 1)
    }

    pub fn down(&self) -> Point {
        Point(self.0, self.1 - 1)
    }

    pub fn tuple(&self) -> (usize, usize) {
        (self.0, self.1)
    }

    pub fn get_adjacent_pts(&self) -> Vec<Point> {
        vec![self.right(), self.left(), self.up(), self.down()]
    }

    pub fn is_in_rect(&self, p_min: Point, p_max: Point) -> bool {
        p_min.0 <= self.0 && self.0 < p_max.0 &&
            p_min.1 <= self.1 && self.1 < p_max.1
    }
}

impl Add for Point {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0, self.1 + other.1)
    }
}

impl Sub for Point {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0, self.1 - other.1)
    }
}

#[cfg(test)]
mod tests {
    use quickcheck as qc;

    impl qc::Arbitrary for super::Point {
        fn arbitrary<G: qc::Gen>(g: &mut G) -> super::Point {
            super::Point(usize::arbitrary(g), usize::arbitrary(g))
        }
    }

    qc::quickcheck! {
        fn prop_add_pts_assoc(p: super::Point,
                              q: super::Point,
                              r: super::Point) -> bool {
            (p + q) + r == p + (q + r)
        }

        fn prop_add_pts_commut(p: super::Point, q: super::Point) -> bool {
            p + q == q + p
        }

        fn prop_sub_pts(p: super::Point, q: super::Point) -> bool {
            p + q - q == p
        }

        fn prop_adjacent_pts(p: super::Point) -> bool {
            let sum: super::Point = p
                .get_adjacent_pts()
                .iter()
                .fold(super::Point(0, 0), |acc, &x| acc + x);

            super::Point(p.0 * 4, p.1 * 4) == sum
        }
    }
}