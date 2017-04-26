//! Unless otherwise documented, units are: m, kg, s

use std::f32;

#[derive(Debug, Copy, Clone)]
pub struct Barrel {
    /// caliber/2
    pub bore_radius: f32,
    /// distance traveled by the bullet from firing to exiting the barrel
    pub bore_length: f32,
    /// meters per turn
    pub rifling_twist: f32,
    pub chamber_body_radius: f32,
    pub chamber_body_length: f32,
}

#[derive(Debug, Copy, Clone)]
pub struct Firing {
    /// m/s
    pub bullet_linear_velocity: f32,
    /// rad/s
    pub bullet_angular_velocity: f32,
    /// J
    pub muzzle_blast_energy: f32,
}

impl Firing {
    pub fn adiabatic(cartridge: &Cartridge, barrel: &Barrel) -> Self {
        const R: f32 = 8.3144598;
        const STANDARD_TEMPERATURE: f32 = 293.15;

        // U = cNT; T = U/(cN)
        let powder_mass = cartridge.internal_volume() * cartridge.powder.density;
        let initial_temperature = STANDARD_TEMPERATURE + (cartridge.powder.heat_of_combustion * powder_mass) / (cartridge.powder.products.specific_heat_capacity_cv * powder_mass);
        let quantity_per_mass = 1. / cartridge.powder.products.molar_mass;
        // P = nRT/V, factor out volume
        let initial_pressure = cartridge.powder.density * quantity_per_mass * R * initial_temperature;
        let gamma = cartridge.powder.products.heat_capacity_ratio;
        let k = initial_pressure * cartridge.internal_volume().powf(gamma);
        let final_volume = cartridge.internal_volume() + f32::consts::PI * barrel.bore_radius * barrel.bore_radius * barrel.bore_length;
        let final_pressure = k / final_volume.powf(gamma);
        let n = quantity_per_mass * powder_mass;
        // T = PV/(nR)
        let final_temperature = final_pressure * final_volume / (n * R);

        let work = k * (final_volume.powf(1.-gamma) - cartridge.internal_volume().powf(1.-gamma)) / (1.-gamma);
        let bullet_mass = cartridge.bullet.mass();

        // KE = KEl + KEr
        // KEr = 0.5 * I * w^2
        // KEl = 0.5 * m * v^2
        // 2 pi v / twist = w

        // KE = 0.5 * I * (2 pi v / twist)^2 + 0.5 * m * v^2
        // Solve for v:
        let v = f32::consts::SQRT_2 * barrel.rifling_twist *
            (work
             / (bullet_mass * barrel.rifling_twist*barrel.rifling_twist
                + 4. * f32::consts::PI * f32::consts::PI * cartridge.bullet.inertia_around_axis())).sqrt();

        let ke_r = work - 0.5 * bullet_mass * v * v;
        Firing {
            bullet_linear_velocity: v,
            // w = sqrt(2KEr/I)
            bullet_angular_velocity: f32::consts::SQRT_2 * (ke_r/cartridge.bullet.inertia_around_axis()).sqrt(),
            // U = cNT
            muzzle_blast_energy: 1000. * cartridge.powder.products.specific_heat_capacity_cv * powder_mass * final_temperature,
        }

        // // We model gunpowder as yielding mostly triatomic gas when burned, so the heat capacity ratio is:
        // const GAMMA: f32 = 1.+2./7.;

        // // g/mol
        // const NITROCELLULOSE_MOLAR_MASS: f32 = 297.1333;
        // const O2_MOLAR_MASS: f32 = 31.99880;
        // const CO2_MOLAR_MASS: f32 = 44.0095;
        // const N2_MOLAR_MASS: f32 = 28.0134;
        // const H2O_MOLAR_MASS: f32 = 18.01528;

        // // nitrocellulose combustion: 4 C6H7N3O11 + 9 O2 ⟶ 24 CO2 + 6 N2 + 14 H2O
        // const NITROCELLULOSE_QUANTITY: f32 = 4;
        // const O2_QUANTITY: f32 = 9.;
        // const CO2_QUANTITY: f32 = 24.;
        // const N2_QUANTITY: f32 = 6.;
        // const H2O_QUANTITY: f32 = 14.;
        // const MIXTURE_QUANTITY: f32 = CO2_QUANTITY + N2_QUANTITY + H2O_QUANTITY;

        // // kg
        // const TOTAL_MASS: f32 = NITROCELLULOSE_QUANTITY * NITROCELLULOSE_MOLAR_MASS + O2_QUANTITY * O2_MOLAR_MASS * .001;

        // const MIXTURE_PER_MASS: f32 = MIXTURE_QUANTITY / TOTAL_MASS;

        // let powder_mass = cartridge.internal_volume() * POWDER_DENSITY;
        // let initial_temperature = NITROCELLULOSE_HEAT_OF_COMBUSTION * powder_mass
        // let initial_pressure = POWDER_DENSITY * MIXTURE_PER_MASS * R * initial_temperature; // P = nRT/V, factor out volume
    }
}

#[derive(Debug, Copy, Clone)]
pub struct Suppressor {
    pub length: f32,
    pub volume: f32,
}

#[derive(Debug, Copy, Clone)]
pub enum Feed {
    Single,
    InternalMagazine,
    Magazine,
    Belt,
}

#[derive(Debug, Copy, Clone)]
pub enum Action {
    Manual,
    Semi,
    Full,
}

#[derive(Debug, Copy, Clone)]
pub enum MagazineType {
    Box,
    Drum,
    DualDrum,
}

#[derive(Debug, Copy, Clone)]
pub struct Magazine {
    pub ty: MagazineType,
    pub capacity: u16,
    pub cartridge_radius: f32,
    pub cartridge_length: f32,
}

#[derive(Debug, Copy, Clone)]
pub struct Cartridge {
    pub bullet: Bullet,
    pub powder: Powder,
    pub body_internal_radius: f32,
    pub body_length: f32,
}

impl Cartridge {
    pub fn internal_volume(&self) -> f32 {
        f32::consts::PI * self.body_internal_radius * self.body_internal_radius * self.body_length
    }
}

#[derive(Debug, Copy, Clone)]
pub struct Powder {
    /// kg/m^3
    pub density: f32,
    /// kJ/kg
    pub heat_of_combustion: f32,

    pub products: Gas,
}

#[derive(Debug, Copy, Clone)]
pub struct Gas {
    pub heat_capacity_ratio: f32,
    /// kJ/(kg K)
    pub specific_heat_capacity_cv: f32,
    /// kg/mol
    pub molar_mass: f32,
}

#[derive(Debug, Copy, Clone)]
pub struct Bullet {
    pub radius: f32,
    pub ogive_length: f32,
    pub cylinder_length: f32,
    pub tail_length: f32,
    pub density: f32,
}

const TAIL_ANGLE: f32 = 7. * f32::consts::PI / 180.;

impl Bullet {
    pub fn volume(&self) -> f32 {
        let tail_volume = ConicFrustum::from_angle(self.radius, self.tail_length, TAIL_ANGLE).volume();
        let cylinder_volume = f32::consts::PI * self.radius * self.radius * self.cylinder_length;
        let nose_volume = TangentOgive::new(self.ogive_length, self.radius).volume();
        tail_volume + cylinder_volume + nose_volume
    }

    pub fn mass(&self) -> f32 {
        self.density * self.volume()
    }

    pub fn inertia_around_axis(&self) -> f32 {
        let tail_moment = self.density * ConicFrustum::from_angle(self.radius, self.tail_length, TAIL_ANGLE).inertia_around_axis();
        // 1/2 m r^2 = 1/2 pi r^2 h r^2
        let cylinder_mass = f32::consts::PI * self.radius * self.radius * self.cylinder_length;
        let cylinder_moment = self.density * 0.5 * cylinder_mass * self.radius * self.radius;
        let nose_moment = self.density * TangentOgive::new(self.ogive_length, self.radius).inertia_around_axis();
        tail_moment + cylinder_moment + nose_moment
    }
}

#[derive(Debug, Copy, Clone)]
pub struct ConicFrustum {
    pub r1: f32,
    pub r2: f32,
    pub h: f32,
}

impl ConicFrustum {
    pub fn from_angle(r: f32, h: f32, angle: f32) -> Self {
        ConicFrustum {
            r1: r,
            r2: r - (h * angle.tan()),
            h: h,
        }
    }

    pub fn volume(&self) -> f32 {
        (1./3.) * f32::consts::PI * self.h * (self.r1 * self.r1 + self.r1 * self.r2 + self.r2 * self.r2)
    }

    pub fn inertia_around_axis(&self) -> f32 {
        // r(t) := r1*(1-t) + r2*t;
        // integrate(1/2 * (d * %pi * r(t)^2) * r(t)^2, t, 0, 1);
        let r2_2 = self.r2 * self.r2;
        let r1_2 = self.r1 * self.r1;
        self.h * (f32::consts::PI*(r2_2*r2_2 + self.r1*r2_2*self.r2 + r1_2*r2_2 + r1_2*self.r1*self.r2 + r1_2*r1_2))/10.
    }
}

#[derive(Debug, Copy, Clone)]
pub struct TangentOgive {
    pub length: f32,
    pub radius: f32,
}

impl TangentOgive {
    pub fn new(length: f32, radius: f32) -> Self {
        TangentOgive {
            length: length,
            radius: radius,
        }
    }

    pub fn rho(&self) -> f32 {
        (self.radius * self.radius + self.length * self.length) / (2. * self.radius)
    }

    pub fn volume(&self) -> f32 {
        let rho = self.rho();
        let rho2 = rho * rho;
        // integrate(pi*(sqrt(rho^2-(L-x)^2)+R-rho)^2, x, 0, L);
        f32::consts::PI*((-(self.length*rho-self.length*self.radius)*(rho2-self.length*self.length).sqrt())-(self.length/rho).asin()*rho2*rho
            +(6.*self.length*rho2-6.*self.length*self.radius*rho+3.*self.length*self.radius*self.radius-self.length*self.length*self.length)/3.
            +self.radius*(self.length/rho).asin()*rho2)
    }

    pub fn inertia_around_axis(&self) -> f32 {
        // r(x) := sqrt(rho^2-(L-x)^2)+R-rho;
        // integrate(1/2 * (d * %pi * r(x)^2) * r(x)^2, x, 0, L);
        let rho = self.rho();
        let a = (self.length/rho).asin();
        (f32::consts::PI
         *((120.*self.length*rho.powi(4)-240.*self.length*self.radius*rho.powi(3)+(180.*self.length*self.radius.powi(2)-40.*self.length.powi(3))*rho.powi(2)
            +(60.*self.length.powi(3)*self.radius-60.*self.length*self.radius.powi(3))*rho+15.*self.length*self.radius.powi(4)-30.*self.length.powi(3)*self.radius.powi(2)+3.*self.length.powi(5))
           /15.
           -(7.*a*rho.powi(5)-15.*self.radius*a*rho.powi(4)
             +(rho.powi(2)-self.length.powi(2)).sqrt()
             *(9.*self.length*rho.powi(3)-17.*self.length*self.radius*rho.powi(2)+(12.*self.length*self.radius.powi(2)-2.*self.length.powi(3))*rho
               -4.*self.length*self.radius.powi(3)+2.*self.length.powi(3)*self.radius)
             +12.*self.radius.powi(2)*a*rho.powi(3)-4.*self.radius.powi(3)*a*rho.powi(2))
           /2.))
            /2.
    }
}

#[test]
fn ogive_volume_sanity() {
    let length = 1.0;
    let radius = 0.5;
    assert!(TangentOgive { length: length, radius: radius }.volume() < f32::consts::PI * radius * radius * length);
}
