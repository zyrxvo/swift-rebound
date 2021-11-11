print("Hello, world!")
var p = Particle()
print(p)
p = ParticleNAN()
print(p)
p = Particle(m: 1)

let G = 1.0
let m = 0.0
let a = 1.0
let e = 0.0
let inc = 1.0
let Omega = 0.0
let omega = 0.0
let f = 0.0

var res = tools_orbit_to_particle_err(G: G, primary: p, m: m, a: a, e: e, inc: inc, Omega: Omega, omega: omega, f: f)
print(res.0)
var res2 = tools_particle_to_orbit_err(G: G, p: res.0, primary: p)
print(res2)