extern crate pew;

use pew::{adiabatic, numerical, Barrel, Cartridge, Bullet, Powder, Gas, rifling_angle};

fn main() {
    let standard_powder = Powder {
        density: 777., // Plausible: http://www.tacticoolproducts.com/powder.pdf
        heat_of_combustion: 10_000_000., // Plausible: http://nvlpubs.nist.gov/nistpubs/jres/44/jresv44n4p387_A1b.pdf
        products: Gas {
            heat_capacity_ratio: 1.0 + 2./5.,
            specific_heat_capacity_cv: 1600.,
            molar_mass: 28.0134e-3,
        }
    };

    const LEAD_DENSITY: f32 = 11340.;

    let m2_barrel = Barrel {
        bore_radius: 12.7e-3 / 2.,
        bore_length: 1.56,
        rifling_angle: rifling_angle(0.380, 12.7e-3 / 2.),
        chamber_body_radius: 20.42e-3,
        chamber_body_length: 76.35e-3,
    };
    let fiftybmg = Cartridge {
        body_internal_radius: m2_barrel.chamber_body_radius - 0.3e-3,
        body_length: m2_barrel.chamber_body_length,
        bullet: Bullet {
            radius: 12.7e-3 / 2.,
            ogive_length: 30.48e-3,
            cylinder_length: 22.86e-3,
            tail_length: 7.62e-3,
            density: 8500.,
        },
        powder: standard_powder,
        charge: 0.015098146,
    };
    println!("M2\nadiabatic: {:?}\nnumerical: {:?}\n",
             adiabatic::Firing::new(&fiftybmg, &m2_barrel),
             numerical::Sim::new(&fiftybmg, &m2_barrel).run(1e-6));

    // AR10
    let heavy_rifle_barrel = Barrel {
        bore_radius: 7.62e-3 / 2.,
        bore_length: 0.528,
        rifling_angle: rifling_angle(0.254, 7.62e-3 / 2.),
        chamber_body_radius: 0.0119,
        chamber_body_length: 0.04165,
    };
    let m80 = Cartridge {
        body_internal_radius: heavy_rifle_barrel.chamber_body_radius - 0.0003,
        body_length: heavy_rifle_barrel.chamber_body_length,
        bullet: Bullet {
            radius:          7.62e-3 / 2.,
            ogive_length:    0.016256,
            cylinder_length: 0.0095631,
            tail_length:     0.004318,
            density:         LEAD_DENSITY,
        },
        powder: standard_powder,
        charge: 0.0029807499,
    };
    println!("AR-10\nadiabatic: {:?}\nnumerical: {:?}\n",
             adiabatic::Firing::new(&m80, &heavy_rifle_barrel),
             numerical::Sim::new(&m80, &heavy_rifle_barrel).run(1e-6));

    // M16A4
    let m16_barrel = Barrel {
        bore_radius: 5.56e-3 / 2.,
        bore_length: 508e-3,
        rifling_angle: rifling_angle(177.8e-3, 5.56e-3 / 2.),
        chamber_body_radius: 9e-3,
        chamber_body_length: 39.55e-3,
    };
    let m855 = Cartridge {
        body_internal_radius: m16_barrel.chamber_body_radius - 0.3e-3,
        body_length: m16_barrel.chamber_body_length,
        bullet: Bullet {
            radius: 5.56e-3 / 2.,
            ogive_length: 15.5e-3,
            cylinder_length: 3.59e-3,
            tail_length: 3.91e-3,
            density: LEAD_DENSITY,
        },
        powder: standard_powder,
        charge: 0.0015551738,
    };
    println!("M16A4\nadiabatic: {:?}\nnumerical: {:?}\n",
             adiabatic::Firing::new(&m855, &m16_barrel),
             numerical::Sim::new(&m855, &m16_barrel).run(1e-6));

    // 9mm sig p320
    let pistol_barrel = Barrel {
        bore_radius: 9e-3 / 2.,
        bore_length: 0.120,
        rifling_angle: rifling_angle(0.254, 0.009/2.),
        chamber_body_radius: 0.00993/2.,
        chamber_body_length: 0.01915,
    };
    let pistol_cartridge = Cartridge {
        body_internal_radius: 9e-3 / 2.,
        body_length: pistol_barrel.chamber_body_length,
        bullet: Bullet {
            radius: 9e-3 / 2.,
            ogive_length: 0.01054,
            cylinder_length: 0.005,
            tail_length: 0.001,
            density: LEAD_DENSITY,
        },
        powder: standard_powder,
        charge: 0.00032399455,
    };
    println!("SIG P320\nadiabatic: {:?}\nnumerical: {:?}",
             adiabatic::Firing::new(&pistol_cartridge, &pistol_barrel),
             numerical::Sim::new(&pistol_cartridge, &pistol_barrel).run(1e-6));
}
