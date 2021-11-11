@testable import swift_rebound
import XCTest

class CoordinateConversionTests: XCTestCase {
    var primary : Particle!
    var orbit : Orbit!

    override func setUp() {
        super.setUp()
        primary = Particle(m: 1)
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
        XCTAssertThrowsError(try tools_orbital_elements_to_particle(G: 1.0, primary: primary, m: 0.0, a: 1.0, e: 1.0, inc: 0.0, Omega: 0.0, omega: 0.0, f: 0.0)) {
            thrownError in err = thrownError as? OrbitError
        }
        XCTAssertEqual(expectedError, err)
        XCTAssertEqual(expectedError.description, err?.description)
    }

    func test_negative_eccentricity_initialization() throws {
        let expectedError = OrbitError.negativeEccentricity
        var err: OrbitError?
        XCTAssertThrowsError(try tools_orbital_elements_to_particle(G: 1.0, primary: primary, m: 0.0, a: 1.0, e: -0.1, inc: 0.0, Omega: 0.0, omega: 0.0, f: 0.0)) {
            thrownError in err = thrownError as? OrbitError
        }
        XCTAssertEqual(expectedError, err)
        XCTAssertEqual(expectedError.description, err?.description)
    }

    func test_bound_orbit_initialization() throws {
        let expectedError = OrbitError.boundOrbitError
        var err: OrbitError?
        XCTAssertThrowsError(try tools_orbital_elements_to_particle(G: 1.0, primary: primary, m: 0.0, a: 1.0, e: 2.5, inc: 0.0, Omega: 0.0, omega: 0.0, f: 0.0)) {
            thrownError in err = thrownError as? OrbitError
        }
        XCTAssertEqual(expectedError, err)
        XCTAssertEqual(expectedError.description, err?.description)
    }

    func test_unbound_orbit_initialization() throws {
        let expectedError = OrbitError.unboundOrbitError
        var err: OrbitError?
        XCTAssertThrowsError(try tools_orbital_elements_to_particle(G: 1.0, primary: primary, m: 0.0, a: -1.0, e: 0.25, inc: 0.0, Omega: 0.0, omega: 0.0, f: 0.0)) {
            thrownError in err = thrownError as? OrbitError
        }
        XCTAssertEqual(expectedError, err)
        XCTAssertEqual(expectedError.description, err?.description)
    }

    func test_unbound_orbit_incorrect_true_anomaly_initialization() throws {
        let expectedError = OrbitError.invalidTrueAnomaly
        var err: OrbitError?
        XCTAssertThrowsError(try tools_orbital_elements_to_particle(G: 1.0, primary: primary, m: 0.0, a: -1.0, e: 2.5, inc: 0.0, Omega: 0.0, omega: 0.0, f: 3.14159265)) {
            thrownError in err = thrownError as? OrbitError
        }
        XCTAssertEqual(expectedError, err)
        XCTAssertEqual(expectedError.description, err?.description)
    }

    func test_massless_primary_initialization() throws {
        let expectedError = OrbitError.masslessPrimary
        var err: OrbitError?
        primary = Particle(m: 0)

        XCTAssertThrowsError(try tools_orbital_elements_to_particle(G: 1.0, primary: primary, m: 0.0, a: 1.0, e: 0.0, inc: 0.0, Omega: 0.0, omega: 0.0, f: 0.0)) {
            thrownError in err = thrownError as? OrbitError
        }
        XCTAssertEqual(expectedError, err)
        XCTAssertEqual(expectedError.description, err?.description)
    }

}
