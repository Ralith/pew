extern crate pew;

use pew::{Barrel, Cartridge, Bullet, Powder, Gas, Firing};

fn main() {
    let standard_powder = Powder {
        density: 777., // Plausible: http://www.tacticoolproducts.com/powder.pdf
        heat_of_combustion: 10000., // Plausible: http://nvlpubs.nist.gov/nistpubs/jres/44/jresv44n4p387_A1b.pdf

        // Products are modeled as pure N2
        products: Gas {
            heat_capacity_ratio: 1. + 2./5.,
            specific_heat_capacity_cv: 0.743,
            molar_mass: 28.0134 / 1000.,
        }
    };
    // AR10
    let heavy_rifle_barrel = Barrel {
        bore_radius: 0.00762/2.,
        bore_length: 0.528,
        rifling_twist: 0.254,
        chamber_body_radius: 0.0119,
        chamber_body_length: 0.04165,
    };
    let m80 = Cartridge {
        body_internal_radius: heavy_rifle_barrel.chamber_body_radius - 0.0003,
        body_length: heavy_rifle_barrel.chamber_body_length,
        bullet: Bullet {
            radius:          0.00762 / 2.,
            ogive_length:    0.016256,
            cylinder_length: 0.0095631,
            tail_length:     0.004318,
            density:         11340.,
        },
        powder: standard_powder,
    };
    println!("AR-10: {:?}", Firing::adiabatic(&m80, &heavy_rifle_barrel));

    // 9mm sig p320
    let pistol_barrel = Barrel {
        bore_radius: 0.009/2.,
        bore_length: 0.120,
        rifling_twist: 0.254,
        chamber_body_radius: 0.00993/2.,
        chamber_body_length: 0.01915,
    };
    let pistol_cartridge = Cartridge {
        body_internal_radius: 0.009/2.,
        body_length: pistol_barrel.chamber_body_length/2.,
        bullet: Bullet {
            radius: 0.009/2.,
            ogive_length: 0.01054,
            cylinder_length: 0.005,
            tail_length: 0.001,
            density: 11340.,
        },
        powder: standard_powder,
    };
    println!("SIG P320: {:?}", Firing::adiabatic(&pistol_cartridge, &pistol_barrel));
}
