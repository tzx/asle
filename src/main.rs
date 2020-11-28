use kiss3d::nalgebra as na;

use kiss3d::light::Light;
use kiss3d::scene::SceneNode;
use kiss3d::window::{State, Window};
use na::Matrix3;
use na::Translation3;
use na::Vector3;

const DT: f32 = 0.005;
const GRAVITY: f32 = -5.0;
const grid_spacing: f32 = 0.2;

// Parameters
const E_0: f32 = 1.4e5;
const NU: f32 = 0.2;
const XI: f32 = 10.0;
const THETA_C: f32 = 2.5e-2;
const THETA_S: f32 = 7.5e-3;

// From Wikipedia: Lame parameters
const MU_0: f32 = E_0 / (2.0 * (1.0 + NU));
const LAMBDA_0: f32 = E_0 * NU / ((1.0 + NU) * (1.0 - 2.0 * NU));

// For particle velocities
const ALPHA: f32 = 0.95;

struct AppState {
    particles: Vec<Particle>,
    grid_nodes: Vec<Vec<Vec<GridNode>>>,
}

impl State for AppState {
    fn step(&mut self, _: &mut Window) {
        // for point in self.particles.iter_mut() {
        //     // We will use Forward Euler Method
        //     point.position += DT * point.velocity;
        //     point.velocity += DT * Vector3::new(0.0, GRAVITY, 0.0);
        //     let translation = Translation3::from(point.position);
        //     if let Some(scene_node) = point.scene_node.as_mut() {
        //         scene_node.set_local_translation(translation);
        //     } else {
        //         panic!("Scene node has not been set for particle");
        //     }
        // }
        self.reset_grid();
        self.interpolate_to_grid();
        self.compute_particle_volumes();
        self.compute_force_hat_elastic();
        self.compute_grid_forces();
        self.compute_grid_velocities();
        self.update_deformation_gradient();
        self.update_particle_velocities();
        self.update_particle_positions();
        for point in self.particles.iter_mut() {
            let translation = Translation3::from(point.position);
            if let Some(scene_node) = point.scene_node.as_mut() {
                scene_node.set_local_translation(translation);
            } else {
                panic!("Scene node has not been set for particle");
            }
        }
    }
}

fn get_bounds(
    position: Vector3<f32>,
    x_max: usize,
    y_max: usize,
    z_max: usize,
) -> (usize, usize, usize, usize, usize, usize) {
    let i_low = std::cmp::max((position.x / grid_spacing).floor() as i32 - 2, 0) as usize;
    let j_low = std::cmp::max((position.y / grid_spacing).floor() as i32 - 2, 0) as usize;
    let k_low = std::cmp::max((position.z / grid_spacing).floor() as i32 - 2, 0) as usize;
    let i_high =
        std::cmp::min((position.x / grid_spacing).floor() as i32 + 2, x_max as i32) as usize;
    let j_high =
        std::cmp::min((position.y / grid_spacing).floor() as i32 + 2, y_max as i32) as usize;
    let k_high =
        std::cmp::min((position.z / grid_spacing).floor() as i32 + 2, z_max as i32) as usize;

    (i_low, i_high, j_low, j_high, k_low, k_high)
}

impl AppState {
    fn reset_grid(&mut self) {
        self.grid_nodes = setup_grid(100, 100, 100);
    }

    fn interpolate_to_grid(&mut self) {
        for point in &self.particles {
            let (i_low, i_high, j_low, j_high, k_low, k_high) = get_bounds(
                point.position,
                self.grid_nodes.len(),
                self.grid_nodes[0].len(),
                self.grid_nodes[0][0].len(),
            );

            // Scatter to the grid nodes
            for i in i_low..i_high {
                for j in j_low..j_high {
                    for k in k_low..k_high {
                        let weight =
                            b_spline(Vector3::new(i as f32, j as f32, k as f32), point.position);
                        self.grid_nodes[i][j][k].mass += weight * point.mass;
                        self.grid_nodes[i][j][k].velocity += weight * point.mass * point.velocity;
                    }
                }
            }
        }
        // Normalized grid velocity
        for r in &mut self.grid_nodes {
            for c in r {
                for grid_node in c {
                    grid_node.velocity /= grid_node.mass;
                }
            }
        }
    }

    fn compute_particle_volumes(&mut self) {
        for point in &mut self.particles {
            let (i_low, i_high, j_low, j_high, k_low, k_high) = get_bounds(
                point.position,
                self.grid_nodes.len(),
                self.grid_nodes[0].len(),
                self.grid_nodes[0][0].len(),
            );

            let mut density = 0.0;
            for i in i_low..i_high {
                for j in j_low..j_high {
                    for k in k_low..k_high {
                        let weight =
                            b_spline(Vector3::new(i as f32, j as f32, k as f32), point.position);
                        density += weight * self.grid_nodes[i][j][k].mass;
                    }
                }
            }
            let h3 = grid_spacing.powi(3);
            density /= h3;
            if density == 0.0 {
                panic!("Density is 0!");
            }
            point.volume = point.mass / density;
        }
    }

    fn compute_force_hat_elastic(&mut self) {
        for point in &mut self.particles {
            let (i_low, i_high, j_low, j_high, k_low, k_high) = get_bounds(
                point.position,
                self.grid_nodes.len(),
                self.grid_nodes[0].len(),
                self.grid_nodes[0][0].len(),
            );
            let mut sum = Matrix3::<f32>::zeros();
            for i in i_low..i_high {
                for j in j_low..j_high {
                    for k in k_low..k_high {
                        let weight_gradient = b_spline_gradient(
                            Vector3::new(i as f32, j as f32, k as f32),
                            point.position,
                        );
                        let velocity = self.grid_nodes[i][j][k].velocity;
                        sum += DT * velocity * weight_gradient.transpose();
                    }
                }
            }
            point.f_hat_ep =
                (Matrix3::<f32>::identity() + sum) * point.deformation_gradient_elastic;
        }
    }

    fn compute_grid_forces(&mut self) {
        for particle in &mut self.particles {
            let j_p = particle.deformation_gradient_plastic.determinant();
            let j_e = particle.deformation_gradient_elastic.determinant();

            let svd = particle.deformation_gradient_elastic.svd(true, true);
            let u = svd.u.unwrap();
            let v = svd.v_t.unwrap();
            let r_e = u * v;

            let lame_mu = MU_0 * (XI * (1.0 - j_p)).exp();
            let lame_lambda = LAMBDA_0 * (XI * (1.0 - j_p)).exp();
            // Taking derivative leads to:
            let sigma = 2.0 * lame_mu * (particle.f_hat_ep - r_e) * particle.f_hat_ep.transpose()
                + lame_lambda * (j_e - 1.0) * j_e * Matrix3::<f32>::identity();

            let force_unweighted = particle.volume * sigma;

            let (i_low, i_high, j_low, j_high, k_low, k_high) = get_bounds(
                particle.position,
                self.grid_nodes.len(),
                self.grid_nodes[0].len(),
                self.grid_nodes[0][0].len(),
            );
            for i in i_low..i_high {
                for j in j_low..j_high {
                    for k in k_low..k_high {
                        let weight_gradient = b_spline_gradient(
                            Vector3::new(i as f32, j as f32, k as f32),
                            particle.position,
                        );
                        self.grid_nodes[i][j][k].force -= force_unweighted * weight_gradient;
                    }
                }
            }
        }
    }

    fn compute_grid_velocities(&mut self) {
        for r in &mut self.grid_nodes {
            for c in r {
                for grid_node in c {
                    if grid_node.mass <= 0.0 {
                        panic!("mass should not be 0!")
                    }
                    grid_node.next_velocity =
                        grid_node.velocity + DT * grid_node.force / grid_node.mass;

                    // TODO: TODO: DO COLLISIONS
                }
            }
        }
    }

    fn update_deformation_gradient(&mut self) {
        for point in &mut self.particles {
            let (i_low, i_high, j_low, j_high, k_low, k_high) = get_bounds(
                point.position,
                self.grid_nodes.len(),
                self.grid_nodes[0].len(),
                self.grid_nodes[0][0].len(),
            );

            let mut velocity_gradient = Matrix3::<f32>::zeros();
            for i in i_low..i_high {
                for j in j_low..j_high {
                    for k in k_low..k_high {
                        let weight_gradient = b_spline_gradient(
                            Vector3::new(i as f32, j as f32, k as f32),
                            point.position,
                        );
                        let velocity = self.grid_nodes[i][j][k].next_velocity;
                        velocity_gradient += velocity * weight_gradient.transpose();
                    }
                }
            }

            let f_hat_elastic_next = (Matrix3::<f32>::identity() + DT * velocity_gradient)
                * point.deformation_gradient_elastic;
            let f_next = f_hat_elastic_next * point.deformation_gradient_plastic;

            let svd = f_hat_elastic_next.svd(true, true);
            let u = svd.u.unwrap();
            let v_t = svd.v_t.unwrap();
            let s_hat_vec: Vector3<f32> = svd.singular_values;

            let s_vec = Vector3::new(
                na::clamp(s_hat_vec.x, 1.0 - THETA_C, 1.0 + THETA_S),
                na::clamp(s_hat_vec.y, 1.0 - THETA_C, 1.0 + THETA_S),
                na::clamp(s_hat_vec.z, 1.0 - THETA_C, 1.0 + THETA_S),
            );
            let s_vec_inv = Vector3::new(s_vec.x.recip(), s_vec.y.recip(), s_vec.z.recip());
            let s = Matrix3::from_diagonal(&s_vec);
            let s_inv = Matrix3::from_diagonal(&s_vec_inv);

            point.deformation_gradient_elastic = u * s * v_t;
            point.deformation_gradient_plastic = v_t.transpose() * s_inv * u.transpose() * f_next;
        }
    }

    fn update_particle_velocities(&mut self) {
        for particle in &mut self.particles {
            let mut velocity_pic = Vector3::<f32>::zeros();
            let mut velocity_flip = particle.velocity;
            let (i_low, i_high, j_low, j_high, k_low, k_high) = get_bounds(
                particle.position,
                self.grid_nodes.len(),
                self.grid_nodes[0].len(),
                self.grid_nodes[0][0].len(),
            );

            for i in i_low..i_high {
                for j in j_low..j_high {
                    for k in k_low..k_high {
                        let weight = b_spline(
                            Vector3::new(i as f32, j as f32, k as f32),
                            particle.position,
                        );
                        let node = &self.grid_nodes[i][j][k];
                        velocity_pic += node.next_velocity * weight;
                        velocity_flip += (node.next_velocity - node.velocity) * weight;
                    }
                }
            }
            particle.velocity = (1.0 - ALPHA) * velocity_pic + ALPHA * velocity_flip;
            // Update due to acceleration from gravity
            particle.velocity += Vector3::new(0.0, GRAVITY, 0.0) * DT;
        }
    }

    fn update_particle_positions(&mut self) {
        for particle in &mut self.particles {
            particle.position += particle.velocity * DT;
        }
    }
}

fn b_spline(grid_index: Vector3<f32>, evaluation_position: Vector3<f32>) -> f32 {
    let scaled = (evaluation_position - grid_index * grid_spacing) / grid_spacing;
    n_fn(scaled.x) * n_fn(scaled.y) * n_fn(scaled.z)
}

fn b_spline_gradient(grid_index: Vector3<f32>, evaluation_position: Vector3<f32>) -> Vector3<f32> {
    let scaled = (evaluation_position - grid_index * grid_spacing) / grid_spacing;
    let n_fn_x = n_fn(scaled.x);
    let n_fn_y = n_fn(scaled.y);
    let n_fn_z = n_fn(scaled.z);

    let dx = n_fn_y * n_fn_z * n_fn_derivative(scaled.x) / grid_spacing;
    let dy = n_fn_x * n_fn_z * n_fn_derivative(scaled.y) / grid_spacing;
    let dz = n_fn_x * n_fn_y * n_fn_derivative(scaled.z) / grid_spacing;

    Vector3::new(dx, dy, dz)
}

fn n_fn(x: f32) -> f32 {
    let x = x.abs();
    if x < 1.0 {
        return 0.5 * x.powi(3) - x.powi(2) + 2.0 / 3.0;
    } else if x < 2.0 {
        return -1.0 / 6.0 * x.powi(3) + x.powi(2) - 2.0 * x + 4.0 / 3.0;
    } else {
        return 0.0;
    }
}

fn n_fn_derivative(x: f32) -> f32 {
    let x_abs = x.abs();
    let signed_x = if x > 0.0 { 1.0 } else { -1.0 };
    if x_abs < 1.0 {
        return 1.5 * x * x * signed_x - 2.0 * x;
    } else if x_abs < 2.0 {
        return -0.5 * x * x * signed_x + 2.0 * x - 2.0 * signed_x;
    } else {
        return 0.0;
    }
}

fn main() {
    let mut window = Window::new("Balls");

    // TODO: Grid should be constant
    let grid_nodes = setup_grid(100, 100, 100);
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
    next_velocity: Vector3<f32>,
    force: Vector3<f32>,
}

struct Particle {
    position: Vector3<f32>,
    velocity: Vector3<f32>,
    mass: f32,
    scene_node: Option<SceneNode>,
    volume: f32,
    deformation_gradient_plastic: Matrix3<f32>,
    deformation_gradient_elastic: Matrix3<f32>,
    f_hat_ep: Matrix3<f32>,
}

fn setup_grid(x: i32, y: i32, z: i32) -> Vec<Vec<Vec<GridNode>>> {
    let mut res = Vec::new();
    for _ in 0..x {
        let mut ys = Vec::new();
        for _ in 0..y {
            let mut zs = Vec::new();
            for _ in 0..z {
                zs.push(GridNode {
                    // TODO: Mass should be a constant
                    mass: 5.0,
                    velocity: Vector3::from_element(0.0),
                    next_velocity: Vector3::from_element(0.0),
                    force: Vector3::from_element(0.0),
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
                volume: 1.0,
                deformation_gradient_elastic: Matrix3::<f32>::identity(), // TODO: these really should be optional
                deformation_gradient_plastic: Matrix3::<f32>::identity(),
                f_hat_ep: Matrix3::<f32>::identity(),
            };
            particles.push(particle);
        }
    }

    particles
}
