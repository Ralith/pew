use std::f32;

use {Cartridge, Barrel, rifling_travel};

#[derive(Debug, Copy, Clone)]
pub struct Firing {
    /// m/s
    bullet_linear_velocity: f32,
    /// rad/s
    bullet_angular_velocity: f32,
    /// J
    muzzle_blast_energy: f32,
}

impl Firing {
    pub fn new(cartridge: &Cartridge, barrel: &Barrel) -> Self {
        const R: f32 = 8.3144598;
        const STANDARD_TEMPERATURE: f32 = 293.15;

        // U = cNT; T = U/(cN)
        let initial_temperature = STANDARD_TEMPERATURE + cartridge.powder.heat_of_combustion / cartridge.powder.products.specific_heat_capacity_cv;
        let quantity_per_mass = 1. / cartridge.powder.products.molar_mass;
        // P = nRT/V, factor out volume
        let initial_pressure = cartridge.charge * quantity_per_mass * R * initial_temperature / barrel.chamber_volume();
        let gamma = cartridge.powder.products.heat_capacity_ratio;
        let k = initial_pressure * cartridge.body_volume().powf(gamma);
        let final_volume = cartridge.body_volume() + f32::consts::PI * barrel.bore_radius * barrel.bore_radius * barrel.bore_length;
        let final_pressure = k / final_volume.powf(gamma);
        let n = quantity_per_mass * cartridge.charge;
        // T = PV/(nR)
        let final_temperature = final_pressure * final_volume / (n * R);

        let work = k * (final_volume.powf(1.-gamma) - cartridge.body_volume().powf(1.-gamma)) / (1.-gamma);
        let bullet_mass = cartridge.bullet.mass();

        // KE = KEl + KEr
        // KEr = 0.5 * I * w^2
        // KEl = 0.5 * m * v^2
        // 2 pi v / travel = w

        // KE = 0.5 * I * (2 pi v / travel)^2 + 0.5 * m * v^2
        // Solve for v:
        let travel = rifling_travel(barrel.rifling_angle, barrel.bore_radius);
        let v = f32::consts::SQRT_2 * travel *
            (work
             / (bullet_mass * travel*travel
                + 4. * f32::consts::PI * f32::consts::PI * cartridge.bullet.inertia_around_axis())).sqrt();

        let ke_r = work - 0.5 * bullet_mass * v * v;
        Firing {
            bullet_linear_velocity: v,
            // w = sqrt(2KEr/I)
            bullet_angular_velocity: f32::consts::SQRT_2 * (ke_r/cartridge.bullet.inertia_around_axis()).sqrt(),
            // U = cNT
            muzzle_blast_energy: cartridge.powder.products.specific_heat_capacity_cv * cartridge.charge * final_temperature,
        }
    }
}
