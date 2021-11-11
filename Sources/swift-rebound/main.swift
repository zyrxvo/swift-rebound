print("Hello, world!")
var p = Particle(m: 1)
print(p)

let G = 1.0
let m = 0.0
let a = 1.0
let e = 0.0
let inc = 0.0
let Omega = 0.0
let omega = 0.0
let f = 0.0

do {
    let res = try tools_orbit_to_particle_err(G: G, primary: p, m: m, a: a, e: e, inc: inc, Omega: Omega, omega: omega, f: f)
    print(res)
    let res2 = tools_particle_to_orbit_err(G: G, p: res, primary: p)
    print(res2)
}
catch {
    print(error)
}