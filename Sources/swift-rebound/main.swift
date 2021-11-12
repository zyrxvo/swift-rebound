print("Welcome to swift-rebound!")

let G = 1.0
let m = 0.0
let a = 1.0
let e = 0.0
let inc = 0.0
let Omega = 0.0
let omega = 0.0
let f = PI

let sim = Simulation()
sim.add(m: 1)
sim.add(m: 1e-6, a: 1)
print(sim)

//do {
//    let p = try Particle(m: 1)
//    print(p)
//
//
//    let particle = try Particle(m: m, a: a, primary: p, G: G)
//    print(particle)
//    let orbit = try tools_particle_to_orbit(G: G, p: particle, primary: p)
//    print(orbit)
//    let particle2 = try tools_orbit_to_particle(G: G, primary: p, m: m, orbit: orbit)
//    print(particle2)
//
//}
//catch {
//    print(error)
//}
