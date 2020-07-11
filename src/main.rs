//! Contours
//!
//! `countours` is a collection of utilities for calculating a nesting
//! of contours in a 2D map.

use peroxide::prelude::*;
use std::collections::{VecDeque, HashSet};
use std::cmp::Ordering::{self, Less, Greater};
mod point;
use point::{Point};
use gnuplot::{Figure};

fn main() {
    let a = match Matrix::read("../data/sampleMap.csv", false, ',') {
        Ok(m) => Terrain::new(m),
        _ => panic!("Unable to read csv")
    };
    let mut flood_reg = Region::zeros_like(&a.matrix);
    let mut boundary = Boundary::new(vec![Point::origin()]);

    let h = -0.22553230822086334;
    let (_, new_boundary_pts) = flood(&a, &mut flood_reg, &boundary, h);

    boundary.pts = new_boundary_pts;

    println!("Boundary len: {}", boundary.pts.len());

    let (lower_h, upper_h) = find_saddle(&a, &mut flood_reg, &boundary);

    println!("Found heights: {:?}, {:?}", lower_h, upper_h);

    let mut fig = Figure::new();
    fig
        .axes2d()
        .image(flood_reg.matrix.data, flood_reg.matrix.row,
               flood_reg.matrix.col, None, &[]);

    //fig.show();
}

#[derive(Copy, Clone)]
pub enum Curvature {Pos, Neg}

pub struct Terrain {
    matrix: Matrix,
    maxs: Vec<Point>,
    mins: Vec<Point>
}

impl Terrain {
    pub fn new(m: Matrix) -> Terrain {
        Terrain {
            matrix: m,
            maxs: Vec::new(),
            mins: Vec::new()
        }
    }

    pub fn zeros_like(m: &Matrix) -> Terrain {
        Terrain{
            matrix: zeros(m.row, m.col),
            maxs: Vec::new(),
            mins: Vec::new()
        }
    }

    pub fn shape(&self) -> (usize, usize) {
        (self.matrix.row, self.matrix.col)
    }

    pub fn can_flood_at(&self, p: Point, height_bound: f64,
                        curvature: Curvature) -> bool {
        let terrain_val = self.matrix[p.tuple()];

        match curvature {
            Curvature::Pos => terrain_val >= height_bound,
            Curvature::Neg => terrain_val < height_bound
        }
    }
}

pub struct Boundary {
    pts: Vec<Point>,
    curvature: Curvature
}

impl Boundary {
    pub fn new(pts: Vec<Point>) -> Boundary {
        Boundary {
            pts: pts,
            curvature: Curvature::Pos
        }
    }
}

pub struct Region {
    matrix: Matrix
}

impl Region {
    pub fn new(m: Matrix) -> Region {
        Region { matrix: m }
    }

    pub fn zeros_like(m: &Matrix) -> Region {
        Region{ matrix: zeros(m.row, m.col) }
    }

    pub fn is_flooded_at(&self, p: Point) -> bool {
        self.matrix[p.tuple()] > 0.0
    }

    pub fn flood_at(&mut self, p: Point) {
        self.matrix[p.tuple()] = 1.0;
    }

    pub fn unflood_at(&mut self, p: Point) {
        self.matrix[p.tuple()] = 0.0;
    }

    pub fn shape(&self) -> (usize, usize) {
        (self.matrix.row, self.matrix.col)
    }

    pub fn is_pt_in_region(&self, p: Point) -> bool {
        p.is_in_rect(Point::origin(), Point(self.matrix.row, self.matrix.col))
    }
}

pub fn flood(terrain: &Terrain, flood_reg: &mut Region,
             boundary: &Boundary, h: f64) -> (Vec<Point>, Vec<Point>) {
    assert_eq!(terrain.shape(), flood_reg.shape());

    let mut queue: VecDeque<Point> = VecDeque::new();
    let mut flooded_pts: Vec<Point> =  Vec::new();
    let mut new_boundary_pts: Vec<Point> = Vec::new();

    for &boundary_pt in &boundary.pts {
        if terrain.can_flood_at(boundary_pt, h, boundary.curvature) {
            flood_reg.flood_at(boundary_pt);
            flooded_pts.push(boundary_pt);
            queue.push_back(boundary_pt);
        } else {
            new_boundary_pts.push(boundary_pt);
        }
    }

    while !queue.is_empty() {
        let curr_pt = queue.pop_front().unwrap();
        let adj_pts: Vec<Point> = curr_pt
            .get_adjacent_pts()
            .into_iter()
            .filter(|&x| flood_reg.is_pt_in_region(x)
                    && !flood_reg.is_flooded_at(x))
            .collect();

        for adj_pt in adj_pts {
            if terrain.can_flood_at(adj_pt, h, boundary.curvature) {
                flood_reg.flood_at(adj_pt);
                flooded_pts.push(adj_pt);
                queue.push_back(adj_pt);
            } else {
                new_boundary_pts.push(adj_pt);
            }
        }
    }

    return (flooded_pts, new_boundary_pts);
}

pub fn trace_boundary(flood_reg: &Region, start: Point) -> HashSet<Point> {
    let mut visited: HashSet<Point> = HashSet::new();
    if flood_reg.is_flooded_at(start) { return visited }

    let mut queue: VecDeque<Point> = VecDeque::new();

    if start
        .get_adjacent_pts()
        .into_iter()
        .any(|x| flood_reg.is_pt_in_region(x)
             && flood_reg.is_flooded_at(x)) {
            visited.insert(start);
            queue.push_back(start);
    } else {
        return visited;
    }

    while !queue.is_empty() {
        let curr_pt = queue.pop_front().unwrap();
        let adj_pts: Vec<Point> = curr_pt
            .get_adjacent_pts()
            .into_iter()
            .filter(|&x| !visited.contains(&x)
                    && flood_reg.is_pt_in_region(x)
                    && !flood_reg.is_flooded_at(x)
                    && x
                    .get_adjacent_pts()
                    .into_iter()
                    .any(|x| flood_reg.is_pt_in_region(x)
                         && flood_reg.is_flooded_at(x)))
            .collect();

        for adj_pt in adj_pts {
            visited.insert(adj_pt);
            queue.push_back(adj_pt);
        }
    }

    visited
}

pub fn cmp(terrain: &Terrain, flood_reg: &mut Region,
           boundary: &Boundary, h: f64) -> Ordering {
    let (flooded_pts, new_boundary_pts) = flood(terrain, flood_reg, boundary, h);

    let extrema = match boundary.curvature {
        Curvature::Pos => &terrain.maxs,
        Curvature::Neg => &terrain.mins
    };

    if extrema
        .into_iter()
        .any(|&x| flood_reg.is_flooded_at(x)) {
            for pt in flooded_pts {
                flood_reg.unflood_at(pt);
            }
            return Greater;
        }

    let first_pt = match new_boundary_pts.first() {
        Some(&first) => first,
        None => panic!() // Deal with this without crashing
    };

    let circle = trace_boundary(flood_reg, first_pt);
    if circle.len() > new_boundary_pts.len() {
        for pt in flooded_pts {
            flood_reg.unflood_at(pt);
        }
        return Greater;
    }

    return Less;
}

fn increasing<T>(x: &T, y: &T) -> Ordering where T: PartialOrd {
    x.partial_cmp(y).unwrap()
}

fn decreasing<T>(x: &T, y: &T) -> Ordering where T: PartialOrd {
    y.partial_cmp(x).unwrap()
}

pub fn find_saddle(terrain: &Terrain, flood_reg: &mut Region,
                   boundary: &Boundary) -> (Option<f64>, Option<f64>) {
    let sort_fn: fn(&f64, &f64) -> Ordering =
        match boundary.curvature {
            Curvature::Pos => decreasing,
            Curvature::Neg => increasing
        };

    let mut heights: Vec<f64> = (&boundary.pts)
        .into_iter()
        .map(|p| terrain.matrix[p.tuple()])
        .collect();

    heights.sort_by(sort_fn);

    let cmp_fn = |&h| cmp(terrain, flood_reg, boundary, h);
    let h_index = match heights.binary_search_by(cmp_fn) {
        Err(i) => i,
        _      => panic!("Comparison function returned equality!")
    };

    let upper_h = if h_index == heights.len() {None} else {Some(heights[h_index])};
    let lower_h = if h_index == 0 {None} else {Some(heights[h_index - 1])};

    (lower_h, upper_h)
}

#[cfg(test)]
mod tests {
    use quickcheck as qc;

    qc::quickcheck! {
        fn prop_flood_all(i: usize, j: usize) -> bool {
            let start = super::Point(i % 50, j % 50);
            let region = super::zeros(50, 50);
            let mut flood_reg = super::Region::new(region);

            let terrain = super::Terrain::zeros_like(&flood_reg.matrix);

            let mut boundary = super::Boundary::new(vec![start]);
            boundary.curvature = super::Curvature::Neg;

            super::flood(&terrain, &mut flood_reg, &boundary, 100.0);

            flood_reg.matrix.data.into_iter().all(|x| x == 1.0)
        }

        fn prop_trace_line(i: usize) -> bool {
            let mut region = super::zeros(50, 50);
            for j in 0..50 {
                region[(i % 50, j)] = 1.0;
            }

            let flood_reg = super::Region::new(region);

            let start: super::Point;
            if i % 50 == 0 { start = super::Point(i % 50 + 1, 0) }
            else { start = super::Point(i % 50 - 1, 0) }

            super::trace_boundary(&flood_reg, start).len() == 50
        }

        fn prop_trace_no_boundary(i: usize, j: usize) -> bool {
            let start = super::Point(i % 50, j % 50);
            let region = super::zeros(50, 50);

            let flood_reg = super::Region::new(region);

            super::trace_boundary(&flood_reg, start).len() == 0
        }
    }
}
