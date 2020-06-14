//! Contours
//!
//! `countours` is a collection of utilities for calculating a nesting
//! of contours in a 2D map.

use peroxide::prelude::*;
use std::collections::VecDeque;
mod point;
use point::{Point};

fn main() {
    let a = py_matrix(vec![
        vec![1, 2],
        vec![3, 4]]);

    a.print();

    let mut ft = zeros(2, 2);
    let sp = vec![Point::origin()];
}

struct Region {
    region: Matrix,
    boundary: Vec<Point>
}

impl Region {
    fn is_flooded_at(&self, p: Point) -> bool {
        self.region[p.tuple()] > 0.0
    }
}

fn flood<F>(terrain: &Matrix, flood_reg: &mut Region, can_flood: F) -> Vec<Point>
where F: Fn(&Matrix, Point) -> bool {
    let flood_reg_shape = (flood_reg.region.row, flood_reg.region.col);
    assert_eq!((terrain.row, terrain.col), flood_reg_shape);

    let max_pt = Point(terrain.row, terrain.col);

    let mut queue: VecDeque<Point> = flood_reg.boundary
        .iter()
        .flat_map(|&x| x.get_adjacent_pts())
        .filter(|&x| x.is_in_rect(Point::origin(), max_pt)
                && flood_reg.is_flooded_at(x))
        .collect();

    let mut flooded_pts: Vec<Point> =  Vec::new();
    let mut new_boundary: Vec<Point> = Vec::new();

    while !queue.is_empty() {
        let curr_pt = queue.pop_front().unwrap_or(Point::origin());
        let adj_pts: Vec<Point> = curr_pt
            .get_adjacent_pts()
            .into_iter()
            .filter(|&x| x.is_in_rect(Point::origin(), max_pt)
                    && !flood_reg.is_flooded_at(x))
            .collect();

        for adj_pt in adj_pts {
            if can_flood(terrain, adj_pt) {
                flood_reg.region[adj_pt.tuple()] = 1.0;

                flooded_pts.push(adj_pt);
                queue.push_back(adj_pt);
            } else {
                new_boundary.push(adj_pt);
            }
        }
    }

    flood_reg.boundary = new_boundary;
    return flooded_pts;
}

#[cfg(test)]
mod tests {
    #[test]
    fn eye() {
        let mut a = super::zeros(10, 10);
        a[(0, 0)] = 1.0;
        let boundary = vec![super::Point(0, 0)];

        let mut region = super::Region { region: a, boundary: boundary};

        let terrain = super::py_matrix(vec![vec![1; 10]; 10]);

        super::flood(&terrain, &mut region, |_, _| true);
        assert_eq!(1, 1);
    }
}
