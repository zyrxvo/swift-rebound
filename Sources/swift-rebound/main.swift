print("Hello, world!")
var p = Particle()
print(p)
p = ParticleNAN()
print(p)
p = Particle(m: 1)

var res = tools_orbit_to_particle_err(G: 0, primary: p, m: 0, a: 0, e: 0, inc: 0, Omega: 0, omega: 0, f: 0)
print(res)