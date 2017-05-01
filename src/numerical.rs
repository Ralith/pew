use std::f32;

use {Cartridge, Barrel, R, STANDARD_TEMPERATURE};

#[derive(Debug, Copy, Clone)]
pub struct Firing {
    /// m/s
    pub bullet_linear_velocity: f32,
    /// rad/s
    pub bullet_angular_velocity: f32,
    /// J
    pub muzzle_blast_energy: f32,
    /// J
    pub cartridge_heat: f32,
    /// J
    pub barrel_heat: f32,
}

pub struct Sim<'a> {
    cartridge: &'a Cartridge,
    barrel: &'a Barrel,
    /// Position of base of bullet along the bore
    x: f32,
    /// Linear velocity of the bullet
    v: f32,
    /// Angular velocity of the bullet
    w: f32,
    /// K
    gas_temperature: f32,
    /// kg
    gas_mass: f32,
    /// J accumulators
    cartridge_heat: f32,
    barrel_heat: f32,
}

impl<'a> Sim<'a> {
    pub fn new(cartridge: &'a Cartridge, barrel: &'a Barrel) -> Self {
        Sim {
            barrel: barrel,
            cartridge: cartridge,
            x: 0.0,
            v: 0.0,
            w: 0.0,
            gas_temperature: STANDARD_TEMPERATURE + cartridge.powder.heat_of_combustion / cartridge.powder.products.specific_heat_capacity_cv,
            gas_mass: cartridge.charge,
            cartridge_heat: 0.0,
            barrel_heat: 0.0,
        }
    }

    fn gas_energy(&self) -> f32 {
        self.cartridge.powder.products.specific_heat_capacity_cv * self.gas_mass * self.gas_temperature
    }

    pub fn run(&mut self, step_size: f32) -> Firing {
        // If atmospheric pressure handling is implemented, we'll need a termination condition for "bullet unable to
        // escape barrel" here.
        while self.x < self.barrel.bore_length {
            self.step(step_size);
        }

        Firing {
            bullet_linear_velocity: self.v,
            bullet_angular_velocity: 0.0,
            muzzle_blast_energy: self.cartridge.powder.products.specific_heat_capacity_cv * self.gas_mass * self.gas_temperature,
            cartridge_heat: self.cartridge_heat,
            barrel_heat: self.barrel_heat,
        }
    }

    fn step(&mut self, dt: f32) {
        // Constants
        let gas_volume = self.barrel.chamber_volume() + f32::consts::PI * self.barrel.bore_radius * self.barrel.bore_radius * self.x;
        let exposed_barrel_area = 2. * f32::consts::PI * self.barrel.bore_radius * self.x;
        let bullet_base_area = f32::consts::PI * self.barrel.bore_radius * self.barrel.bore_radius;
        let gamma = self.cartridge.powder.products.heat_capacity_ratio;

        // P = nRT/V
        let pressure = self.gas_mass * R * self.gas_temperature / (gas_volume * self.cartridge.powder.products.molar_mass);
        let k = pressure * gas_volume.powf(gamma);
        let total_force = pressure * bullet_base_area;
        let acceleration = total_force / self.cartridge.bullet.mass(); // TODO: angular

        const CONVECTIVE_HEAT_TRANSFER: f32 = 250_000.; // W/(m^2 K), guessed

        let temp_diff = self.gas_temperature - STANDARD_TEMPERATURE;

        let barrel_heat_flux_density = CONVECTIVE_HEAT_TRANSFER * temp_diff;
        let barrel_heat_rate = exposed_barrel_area * barrel_heat_flux_density;

        let cartridge_heat_flux_density = CONVECTIVE_HEAT_TRANSFER * temp_diff;
        let cartridge_heat_rate = self.cartridge.body_internal_surface_area() * cartridge_heat_flux_density;

        self.x += self.v * dt;
        self.v += acceleration * dt;

        self.barrel_heat += barrel_heat_rate * dt;
        self.cartridge_heat += cartridge_heat_rate * dt;

        // Adiabatic expansion of pre-existing gas
        let adiabatic_temp = {
            let next_volume = self.barrel.chamber_volume() + f32::consts::PI * self.barrel.bore_radius * self.barrel.bore_radius * self.x;
            // k/v^gamma = p
            let next_pressure = k / next_volume.powf(gamma);
            // PV/(nR) = T
            next_volume * next_pressure * self.cartridge.powder.products.molar_mass / (self.gas_mass * R)
        };
        // Temperature change from convection to barrel/cartridge walls
        let convection_temp_change = (cartridge_heat_rate + barrel_heat_rate) * dt / (self.cartridge.powder.products.specific_heat_capacity_cv * self.gas_mass);

        self.gas_temperature = adiabatic_temp - convection_temp_change;
    }
}

