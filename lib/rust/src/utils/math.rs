use std::ops::{Add, Mul, Sub};

use num::{Float};

fn chmin<T>(a: &mut T, b: T) -> bool
    where T: PartialOrd
{
    if *a > b {
        *a = b;
        true
    } else {
        false
    }
}

fn chmax<T>(a: &mut T, b: T) -> bool
    where T: PartialOrd
{
    if *a < b {
        *a = b;
        true
    } else {
        false
    }
}

struct Vector2<T> {
    x: T,
    y: T,
}

impl<T> Vector2<T> 
    where T: Float
{
    fn new(x: T, y: T) -> Self {
        Self { x, y }
    }

    fn add(&self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }

    fn sub(&self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }

    fn dot(&self, other: Self) -> T {
        self.x * other.x + self.y * other.y
    }

    fn cross(&self, other: Self) -> T {
        self.x * other.y - self.y * other.x
    }

    fn norm(&self) -> T {
        self.x * self.x + self.y * self.y
    }

    fn abs(&self) -> T {
        self.norm().sqrt()
    }

    fn unit(&self) -> Self {
        let abs = self.abs();
        Self {
            x: self.x / abs,
            y: self.y / abs,
        }
    }

    fn rotate(&self, theta: T) -> Self {
        let (s, c) = (theta.sin(), theta.cos());
        Self {
            x: self.x * c - self.y * s,
            y: self.x * s + self.y * c,
        }
    }

    fn rotate90(&self) -> Self {
        Self {
            x: -self.y,
            y: self.x,
        }
    }

    fn rotate270(&self) -> Self {
        Self {
            x: self.y,
            y: -self.x,
        }
    }

    fn rotate180(&self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
        }
    }
}

impl<T> Add for Vector2<T>
    where T: Float
{
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Self::Output {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl<T> Sub for Vector2<T> 
    where T: Float
{
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl<T> Mul<T> for Vector2<T> 
    where T: Float
{
    type Output = Self;

    fn mul(self, scalar: T) -> Self {
        Self {
            x: self.x * scalar,
            y: self.y * scalar,
        }
    }
}

impl<T> Mul for Vector2<T> 
    where T: Float
{
    type Output = T;

    fn mul(self, other: Self) -> T {
        self.x * other.x + self.y * other.y
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vector2() {
        let v1 = Vector2::new(1.0, 2.0);
        let v2 = Vector2::new(3.0, 4.0);
        let v3 = v1.add(v2);
        assert_eq!(v3.x, 4.0);
        assert_eq!(v3.y, 6.0);
    }

    #[test]
    fn test_vector2_add() {
        let v1 = Vector2::new(1.0, 2.0);
        let v2 = Vector2::new(3.0, 4.0);
        let v3 = v1 + v2;
        assert_eq!(v3.x, 4.0);
        assert_eq!(v3.y, 6.0);
    }

    #[test]
    fn test_vector2_sub() {
        let v1 = Vector2::new(1.0, 2.0);
        let v2 = Vector2::new(3.0, 4.0);
        let v3 = v1 - v2;
        assert_eq!(v3.x, -2.0);
        assert_eq!(v3.y, -2.0);
    }

    #[test]
    fn test_vector2_mul() {
        let v1 = Vector2::new(1.0, 2.0);
        let v2 = v1 * 2.0;
        assert_eq!(v2.x, 2.0);
        assert_eq!(v2.y, 4.0);
    }

    #[test]
    fn test_vector2_dot() {
        let v1 = Vector2::new(1.0, 2.0);
        let v2 = Vector2::new(3.0, 4.0);
        let dot = v1 * v2;
        assert_eq!(dot, 11.0);
    }

    #[test]
    fn test_vector2_cross() {
        let v1 = Vector2::new(1.0, 2.0);
        let v2 = Vector2::new(3.0, 4.0);
        let cross = v1.cross(v2);
        assert_eq!(cross, -2.0);
    }

}