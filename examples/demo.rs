extern crate pew;

use pew::{adiabatic, numerical, Barrel, Cartridge, Bullet, Powder, Gas, rifling_angle};

fn main() {
    let standard_powder = Powder {
        density: 777., // Plausible: http://www.tacticoolproducts.com/powder.pdf
        heat_of_combustion: 10_000_000., // Plausible: http://nvlpubs.nist.gov/nistpubs/jres/44/jresv44n4p387_A1b.pdf
        // Products are modeled as pure N2
        products: Gas {
            heat_capacity_ratio: 1. + 2./5.,
            specific_heat_capacity_cv: 743.,
            molar_mass: 28.0134e-3,
        }
    };
    // AR10
    let heavy_rifle_barrel = Barrel {
        bore_radius: 7.62e-3 / 2.,
        bore_length: 0.528,
        rifling_angle: rifling_angle(0.254, 0.00762/2.),
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
            density:         11340.,
        },
        powder: standard_powder,
        charge: 0.0029807499,
    };
    println!("AR-10\nadiabatic: {:?}\nnumerical: {:?}\n",
             adiabatic::Firing::new(&m80, &heavy_rifle_barrel),
             numerical::Sim::new(&m80, &heavy_rifle_barrel).run(1e-6));

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
            density: 11340.,
        },
        powder: standard_powder,
        charge: 0.00032399455,
    };
    println!("SIG P320\nadiabatic: {:?}\nnumerical: {:?}",
             adiabatic::Firing::new(&pistol_cartridge, &pistol_barrel),
             numerical::Sim::new(&pistol_cartridge, &pistol_barrel).run(1e-6));
}
