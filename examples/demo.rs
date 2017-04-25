extern crate pew;

use pew::{Barrel, Cartridge, Bullet, Powder, Gas};

fn main() {
    let barrel = Barrel {
        bore_radius: 0.00762/2.,
        bore_length: 0.4064,
        rifling_twist: 0.3048,
        chamber_body_radius: 0.0119,
        chamber_body_length: 0.04165,
    };
    let m80a1 = Cartridge {
        body_internal_radius: barrel.chamber_body_radius - 0.0003,
        body_length: barrel.chamber_body_length,
        bullet: Bullet {
            radius:          0.00762 / 2.,
            ogive_length:    0.016256,
            cylinder_length: 0.0095631,
            tail_length:     0.004318,
            density:         11340.,
        },
        powder: Powder {
            density: 777., // Plausible: http://www.tacticoolproducts.com/powder.pdf
            heat_of_combustion: 10000., // Plausible: http://nvlpubs.nist.gov/nistpubs/jres/44/jresv44n4p387_A1b.pdf

            // Products are modeled as pure N2
            products: Gas {
                heat_capacity_ratio: 1. + 2./5.,
                specific_heat_capacity_cv: 0.743,
                molar_mass: 28.0134 / 1000.,
            }
        }
    };
    println!("{:?}", barrel.fire(&m80a1));
}
