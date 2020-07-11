//! Point
//!
//! The following defines a structure and a set of methods and functions
//! to operate on 2D coordinates.

use std::ops::Add;
use std::convert::identity;

/// Structure representing 2D array indexes
#[derive(PartialEq, Eq, Hash, Copy, Clone, Debug)]
pub struct Point(pub usize, pub usize);
impl Point {
    /// Initializes a `Point` at the origin
    pub fn origin() -> Point {
        Point(0, 0)
    }

    /// Create a new point to the right of `self`
    pub fn right(&self) -> Option<Point> {
        Some(Point(self.0 + 1, self.1))
    }

    /// Create a new point to the left of `self`. If the coordinates become
    /// negative, then `None` is returned.
    pub fn left(&self) -> Option<Point> {
        if self.0 == 0 { None }
        else { Some(Point(self.0 - 1, self.1)) }
    }

    /// Create a new point above `self`
    pub fn up(&self) -> Option<Point> {
        Some(Point(self.0, self.1 + 1))
    }

    /// Create a new point below `self`. If the coordinates become
    /// negative, then `None` is returned.
    pub fn down(&self) -> Option<Point> {
        if self.1 == 0 { None }
        else { Some(Point(self.0, self.1 - 1)) }
    }

    /// Turns `Point` into a tuple
    pub fn tuple(&self) -> (usize, usize) {
        (self.0, self.1)
    }

    /// Gets a vector of all the points adjacent to `self` that do
    /// not have negative coordinates
    pub fn get_adjacent_pts(&self) -> Vec<Point> {
        vec![self.right(), self.left(), self.up(), self.down()]
            .into_iter()
            .filter_map(identity)
            .collect()
    }

    /// Checks if `self` is in the rectangle specified by the opposite
    /// corners `p_min` and `p_max`.
    ///
    /// *Note:* The bound made by `p_min` is inclusive  and the bound
    /// made by `p_max` is exclusive
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

        fn prop_full_circle(p: super::Point) -> bool {
            let clockwise = Some(p) == p
                .up()
                .map(|x| x.right()).flatten()
                .map(|x| x.down()).flatten()
                .map(|x| x.left()).flatten();
            let counterclockwise = Some(p) == p
                .right()
                .map(|x| x.up()).flatten()
                .map(|x| x.left()).flatten()
                .map(|x| x.down()).flatten();

            clockwise && counterclockwise
        }

        fn prop_neighbor_of_neighbor(p: super::Point) -> bool {
            p
                .get_adjacent_pts()
                .iter()
                .flat_map(|&x| x.get_adjacent_pts())
                .any(|x| x == p)
        }
    }
}
