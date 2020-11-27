use kiss3d::nalgebra as na;

use kiss3d::light::Light;
use kiss3d::scene::SceneNode;
use kiss3d::window::{State, Window};
use na::Translation3;
use na::Vector3;

const DT: f32 = 0.005;
const GRAVITY: f32 = -0.1;

struct AppState {
    particles: Vec<Particle>,
    grid_nodes: Vec<Vec<Vec<GridNode>>>,
}

impl State for AppState {
    fn step(&mut self, _: &mut Window) {
        for point in self.particles.iter_mut() {
            // We will use Forward Euler Method
            point.position += DT * point.velocity;
            point.velocity += DT * Vector3::new(0.0, GRAVITY, 0.0);
            let translation = Translation3::from(point.position);
            if let Some(scene_node) = point.scene_node.as_mut() {
                scene_node.set_local_translation(translation);
            } else {
                panic!("Scene node has not been set for particle");
            }
        }
    }
}

impl AppState {
    fn interpolate_to_grid(&mut self) {
        for point in &self.particles {
            panic!("Not implemented!");
        }
    }
}

fn main() {
    let mut window = Window::new("Balls");

    // TODO: Grid should be constant
    let grid_nodes = setup_grid(15, 15, 15);
    let mut particles = setup_particles(15.0, 15.0, 15.0, 10);

    for p in &mut particles {
        // TODO: sphere_size should be constant
        let c = window.add_sphere(0.05);
        p.scene_node = Some(c);
    }

    let mut c = window.add_cube(15.0, 0.03, 15.0);
    c.set_color(0.0, 0.0, 1.0);
    c.prepend_to_local_translation(&Translation3::new(7.5, 0.0, 7.5));
    window.set_light(Light::StickToCamera);

    let state = AppState {
        grid_nodes,
        particles,
    };
    window.render_loop(state);
}

/// Interpolated nodes to a single point in the grid
struct GridNode {
    mass: f32,
    velocity: Vector3<f32>,
}

struct Particle {
    position: Vector3<f32>,
    velocity: Vector3<f32>,
    mass: f32,
    scene_node: Option<SceneNode>,
    // volume:
}

fn setup_grid(x: i32, y: i32, z: i32) -> Vec<Vec<Vec<GridNode>>> {
    let mut res = Vec::new();
    for _ in 0..x {
        let mut ys = Vec::new();
        for _ in 0..y {
            let mut zs = Vec::new();
            for _ in 0..z {
                zs.push(GridNode {
                    mass: 0.0,
                    velocity: Vector3::from_element(0.0),
                })
            }
            ys.push(zs);
        }
        res.push(ys);
    }

    res
}

fn setup_particles(x: f32, y: f32, z: f32, num_particles_1d: i32) -> Vec<Particle> {
    let size = 0.1;
    let x_offset = x / 2.0 - (num_particles_1d / 2) as f32 * size;
    let y_offset = y / 2.0 - (num_particles_1d / 2) as f32 * size;
    let z_offset = z / 2.0 - (num_particles_1d / 2) as f32 * size;
    let mut particles = Vec::new();
    // Only 2D for now
    for i in 0..num_particles_1d {
        for j in 0..num_particles_1d {
            // TODO: mass, size should be a constant
            let particle = Particle {
                position: Vector3::new(
                    x_offset + size * i as f32,
                    y_offset + size * j as f32,
                    z_offset,
                ),
                velocity: Vector3::new(0.0, 0.0, 0.0),
                mass: 5.0,
                scene_node: None,
            };
            particles.push(particle);
        }
    }

    particles
}
