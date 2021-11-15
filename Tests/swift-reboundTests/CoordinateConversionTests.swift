@testable import swift_rebound
import XCTest

class CoordinateConversionTests: XCTestCase {
    var primary : Particle!
    var orbit : Orbit!
    let G = 1.0

    override func setUp() {
        super.setUp()
        do {
            primary = try Particle(m: 1)
        }
        catch {
            print(error)
        }
       
        orbit = Orbit()
    }

    override func tearDown() {
        primary = nil
        orbit = nil
        super.tearDown()
    }

    func test_radial_orbit_initialization() throws {
        let expectedError = OrbitError.radialOrbit
        var err: OrbitError?
        XCTAssertThrowsError(try tools_orbital_elements_to_particle(G: G, primary: primary, m: 0.0, a: 1.0, e: 1.0, inc: 0.0, Omega: 0.0, omega: 0.0, f: 0.0)) {
            thrownError in err = thrownError as? OrbitError
        }
        XCTAssertEqual(expectedError, err)
        XCTAssertEqual(expectedError.description, err?.description)
    }

    func test_negative_eccentricity_initialization() throws {
        let expectedError = OrbitError.negativeEccentricity
        var err: OrbitError?
        XCTAssertThrowsError(try tools_orbital_elements_to_particle(G: G, primary: primary, m: 0.0, a: 1.0, e: -0.1, inc: 0.0, Omega: 0.0, omega: 0.0, f: 0.0)) {
            thrownError in err = thrownError as? OrbitError
        }
        XCTAssertEqual(expectedError, err)
        XCTAssertEqual(expectedError.description, err?.description)
    }

    func test_bound_orbit_initialization() throws {
        let expectedError = OrbitError.boundOrbitError
        var err: OrbitError?
        XCTAssertThrowsError(try tools_orbital_elements_to_particle(G: G, primary: primary, m: 0.0, a: 1.0, e: 2.5, inc: 0.0, Omega: 0.0, omega: 0.0, f: 0.0)) {
            thrownError in err = thrownError as? OrbitError
        }
        XCTAssertEqual(expectedError, err)
        XCTAssertEqual(expectedError.description, err?.description)
    }

    func test_unbound_orbit_initialization() throws {
        let expectedError = OrbitError.unboundOrbitError
        var err: OrbitError?
        XCTAssertThrowsError(try tools_orbital_elements_to_particle(G: G, primary: primary, m: 0.0, a: -1.0, e: 0.25, inc: 0.0, Omega: 0.0, omega: 0.0, f: 0.0)) {
            thrownError in err = thrownError as? OrbitError
        }
        XCTAssertEqual(expectedError, err)
        XCTAssertEqual(expectedError.description, err?.description)
    }

    func test_unbound_orbit_incorrect_true_anomaly_initialization() throws {
        let expectedError = OrbitError.invalidTrueAnomaly
        var err: OrbitError?
        XCTAssertThrowsError(try tools_orbital_elements_to_particle(G: G, primary: primary, m: 0.0, a: -1.0, e: 2.5, inc: 0.0, Omega: 0.0, omega: 0.0, f: 3.14159265)) {
            thrownError in err = thrownError as? OrbitError
        }
        XCTAssertEqual(expectedError, err)
        XCTAssertEqual(expectedError.description, err?.description)
    }

    func test_massless_primary_initialization() throws {
        let expectedError = OrbitError.masslessPrimary
        var err: OrbitError?
        primary = try Particle(m: 0)

        XCTAssertThrowsError(try tools_orbital_elements_to_particle(G: G, primary: primary, m: 0.0, a: 1.0, e: 0.0, inc: 0.0, Omega: 0.0, omega: 0.0, f: 0.0)) {
            thrownError in err = thrownError as? OrbitError
        }
        XCTAssertEqual(expectedError, err)
        XCTAssertEqual(expectedError.description, err?.description)
    }
    
    func test_massless_primary_orbit_initialization() throws {
        let expectedError = OrbitError.masslessPrimary
        var err: OrbitError?
        let particle = try Particle(a: 1, primary: primary)
        primary = try Particle(m: 0)

        XCTAssertThrowsError(try tools_particle_to_orbit(G: G, p: particle, primary: primary)) {
            thrownError in err = thrownError as? OrbitError
        }
        XCTAssertEqual(expectedError, err)
        XCTAssertEqual(expectedError.description, err?.description)
    }
    
    func test_orbit_coincides_with_primary_initialization() throws {
        let expectedError = OrbitError.coincidesWithPrimary
        var err: OrbitError?
        let particle = try Particle()

        XCTAssertThrowsError(try tools_particle_to_orbit(G: G, p: particle, primary: primary)) {
            thrownError in err = thrownError as? OrbitError
        }
        XCTAssertEqual(expectedError, err)
        XCTAssertEqual(expectedError.description, err?.description)
    }
    
    func test_convert_from_particle_to_orbit_and_back() throws {
        let mass = 0.0
        let particle = try Particle(m: mass, x: 1, vy: 1)

        XCTAssertNoThrow(try tools_particle_to_orbit(G: G, p: particle, primary: primary))
        let orb = try tools_particle_to_orbit(G: 1.0, p: particle, primary: primary)
        XCTAssertNoThrow(try tools_orbit_to_particle(G: 1.0, primary: primary, m: mass, orbit: orb))
        let part = try tools_orbit_to_particle(G: 1.0, primary: primary, m: mass, orbit: orb)
        XCTAssertEqual(particle, part)
    }

}
