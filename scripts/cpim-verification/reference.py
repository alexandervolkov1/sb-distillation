"""Generate SI CPIM references from main.typ at 93c4a9e09dae96427635c2e8130a058c982b6931."""

import argparse
import csv
from dataclasses import dataclass
from itertools import product
from math import pi, sqrt
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq


LENGTH = 0.5  # m
MASS = 0.5  # kg
R = 8.314462618
KB = 1.380649e-23
AVOGADRO = 6.02214076e23
M1, M2 = 0.48704, 0.00201588  # kg/mol
SIGMA1, SIGMA2 = 0.600e-9, 0.289e-9  # m
TORR = 101325 / 760
GAUGE = {100: 0.02010, 200: 0.03680, 300: 0.05110}
STATIONS = np.linspace(0, LENGTH, 11)
STANDARD = (2e-11, 2e-12, 2e-14)  # relative, absolute, and Brent flux tolerances
TIGHT = (2e-13, 2e-14, 2e-16)


def liquid_pressure(temperature):
    # Legacy Antoine law, with temperature in K and result in Pa.
    return 1e5 * 10 ** (2.26041 - 4475.449 / (temperature - 152.352))


@dataclass
class Tube:
    temperature: float  # K
    diameter: float  # m
    feed: float  # ml/h at 293.15 K and 101325 Pa
    outlet: float  # Pa, zero trap drop

    def __post_init__(self):
        t = self.temperature
        self.rt = R * t
        self.area = pi * (self.diameter / 2) ** 2
        self.throughput = 101325 * self.feed * 1e-6 / 3600
        self.hydrogen_flow = self.throughput / (R * 293.15)
        self.n2 = self.hydrogen_flow / self.area
        self.saturation = liquid_pressure(t)
        self.dk1, self.dk2 = (
            self.diameter / 3 * sqrt(8 * self.rt / (pi * mass)) for mass in (M1, M2)
        )
        self.ak = 3 / self.diameter * sqrt(pi * self.rt / 8)
        reduced_mass = 1 / (AVOGADRO * (1 / M1 + 1 / M2))
        self.k12 = (
            3 * KB * t / (16 * ((SIGMA1 + SIGMA2) / 2) ** 2)
            * sqrt(2 * KB * t / (pi * reduced_mass))
        )
        self.mu1 = 5 / (16 * SIGMA1**2) * sqrt(M1 / AVOGADRO * KB * t / pi)
        self.mu2 = 8.76e-6 * (t / 293.85) ** 1.5 * (293.85 + 72) / (t + 72)
        self.phi12, self.phi21 = (
            (1 + sqrt(mu_ratio) * mass_ratio**-0.25) ** 2 / sqrt(8 * (1 + mass_ratio))
            for mu_ratio, mass_ratio in ((self.mu1 / self.mu2, M1 / M2),
                                         (self.mu2 / self.mu1, M2 / M1))
        )
        self.molecular_flux = self.dk1 * self.saturation / (self.rt * LENGTH)

    def rhs(self, n1, continuum):
        # Integrate total pressure and Sb4 mole fraction in distance from the outlet.
        weighted_flux = sqrt(M1) * n1 + sqrt(M2) * self.n2
        viscous_factor = self.diameter**2 / (32 * self.rt)
        diffusion_offset1 = 0 if continuum else 1 / self.dk1
        diffusion_offset2 = 0 if continuum else 1 / self.dk2
        molecular_conductance = 0 if continuum else 1 / self.ak

        def gradient(_distance, state):
            pressure, x1 = state
            x2 = 1 - x1
            mu = (x1 * self.mu1 / (x1 + x2 * self.phi12)
                  + x2 * self.mu2 / (x2 + x1 * self.phi21))
            conductance = viscous_factor * pressure * (x1 * sqrt(M1) + x2 * sqrt(M2)) / mu
            mutual = pressure / self.k12
            composition = x2 * n1 * (mutual + diffusion_offset1)
            composition -= x1 * self.n2 * (mutual + diffusion_offset2)
            return weighted_flux / (conductance + molecular_conductance), self.rt * composition / pressure

        return gradient

    def integrate(self, n1, accuracy, stop_at_saturation=False, continuum=False):
        def saturation(_s, state):
            return state[0] * state[1] - self.saturation

        def hydrogen_zero(_s, state):
            return 1 - state[1]

        saturation.terminal = True
        saturation.direction = 1
        hydrogen_zero.terminal = True
        hydrogen_zero.direction = -1
        events = (hydrogen_zero, saturation) if stop_at_saturation else (hydrogen_zero,)
        rtol, atol, _ = accuracy
        result = solve_ivp(
            self.rhs(n1, continuum), (0.0, LENGTH), (self.outlet, 0.0),
            method="DOP853", rtol=rtol, atol=(atol, atol / self.saturation),
            dense_output=True, events=events,
        )
        if not result.success or len(result.t_events[0]):
            raise ValueError(f"Inadmissible reference trial at N1={n1}: {result.message}")
        return result

    def solve(self, accuracy, continuum=False):
        def objective(n1):
            result = self.integrate(n1, accuracy, stop_at_saturation=True, continuum=continuum)
            distance = result.t[-1]
            if len(result.t_events[1]):
                return self.saturation * (1 - distance / LENGTH)
            pressure, x1 = result.y[:, -1]
            return pressure * x1 - self.saturation

        # An early saturation crossing localizes an overshoot without classifying
        # a trial whose integration failed. The objective vanishes at the inlet.
        bounded = False
        for scale in (1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64):
            upper = scale * self.molecular_flux
            if objective(upper) >= 0:
                bounded = True
                break
        if not bounded:
            raise ValueError("Reference search failed to bracket saturation")
        n1 = brentq(objective, 0.0, upper, xtol=accuracy[2], rtol=9e-16)
        return n1, self.integrate(n1, accuracy, continuum=continuum)


def profile(result, z):
    pressure, x1 = result.sol(LENGTH - np.asarray(z))
    return np.array([pressure * x1, pressure * (1 - x1), x1, 1 - x1]).T


def operating_tube(temp_c, diameter_mm, feed):
    return Tube(temp_c + 273.15, diameter_mm * 1e-3, feed, GAUGE[feed] * TORR)


def lim_tubes():
    for diameter in (1e-5, 1e-6, 1e-7):
        yield Tube(973.15, diameter, 200 * (diameter / 0.01) ** 3, 0.2)
    for outlet in (100, 1000, 10000):
        yield Tube(973.15, 0.010, 200, outlet)


def edn(value):
    if isinstance(value, dict):
        return "{" + "\n ".join(f":{key} {edn(entry)}" for key, entry in value.items()) + "}"
    if isinstance(value, np.ndarray) and value.ndim == 2:
        return "[" + "\n  ".join(edn(row) for row in value) + "]"
    if isinstance(value, (np.ndarray, list, tuple)):
        return "[" + " ".join(edn(item) for item in value) + "]"
    if isinstance(value, (float, np.floating)):
        return repr(float(value))
    return str(value)


def fixture():
    tube = operating_tube(700, 10, 200)
    n1, result = tube.solve(TIGHT)
    mass_flow = M1 * tube.area * n1
    trials = []
    for flux in (0.0, 0.5 * tube.molecular_flux, tube.molecular_flux, 8 * tube.molecular_flux):
        dp_ds, dx_ds = tube.rhs(flux, False)(0, (tube.outlet, 0))
        trials.append({
            "n1": flux,
            "outlet-derivatives": [-tube.outlet * dx_ds, tube.outlet * dx_ds - dp_ds],
            "profile": np.column_stack((STATIONS, profile(tube.integrate(flux, TIGHT), STATIONS))),
        })
    return {
        "area": tube.area, "throughput": tube.throughput, "hydrogen-flow": tube.hydrogen_flow,
        "n2": tube.n2, "saturation": tube.saturation, "outlet": tube.outlet,
        "molecular-flux": tube.molecular_flux, "n1": n1, "mass-flow": mass_flow,
        "time": MASS / mass_flow, "profile": np.column_stack((STATIONS, profile(result, STATIONS))),
        "trials": trials,
        "solid-pressure-at-773.15": TORR * 10 ** (10.571 - 10300 / 773.15),
    }


def grid():
    largest = np.zeros(3)
    rows = []
    for point in product(range(680, 721, 5), range(6, 13), GAUGE):
        tube = operating_tube(*point)
        n1, result = tube.solve(TIGHT)
        coarse_n1, coarse_result = tube.solve(STANDARD)
        sample = profile(result, STATIONS)
        coarse_sample = profile(coarse_result, STATIONS)
        largest = np.maximum(largest, [
            abs(n1 - coarse_n1) / n1, *np.max(abs(sample[:, :2] - coarse_sample[:, :2]), axis=0)
        ])
        if abs(sample[0, 0] - tube.saturation) > 1e-8 or np.any(sample[:, 1] <= 0):
            raise ValueError(f"Reference failed boundary or physical check at {point}")
        rows.append([*point, n1, sample[0, 1], M1 * tube.area * n1, sample[:, :2]])
    print("Reference refinement maxima [relative N1, p1 Pa, p2 Pa]:", largest)
    return rows


def pressure_errors(path):
    """Compare Clojure nodes and common stations against dense SciPy profiles."""
    comparisons = {}
    references = {}
    series = {}
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            point = tuple(int(row[key]) for key in ("temp", "diameter", "feed"))
            key = (*point, row["variant"], row["mode"])
            tube, n1, result = references.get(point, (None, None, None))
            if tube is None:
                tube = operating_tube(*point)
                n1, result = tube.solve(TIGHT)
                references[point] = tube, n1, result
            trial_flux = float(row["n1"])
            if row["mode"] == "fixed":
                n1 = trial_flux
                fixed_key = (*point, n1)
                if fixed_key not in references:
                    references[fixed_key] = tube.integrate(n1, TIGHT)
                result = references[fixed_key]
            z = float(row["z"])
            expected = profile(result, z)[:2]
            actual = np.array([float(row["p1"]), float(row["p2"])])
            errors = comparisons.setdefault(key, np.zeros(7))
            errors[0] = abs(trial_flux / n1 - 1)  # Also the relative mass-flow error.
            if z == 0:
                errors[1] = abs(actual[1] - expected[1])
            offset = 2 if row["kind"] == "node" else 4
            errors[offset:offset + 2] = np.maximum(errors[offset:offset + 2], abs(actual - expected))
            errors[6] = abs(float(row["residual"]))
            nodes, stations, _, _, _ = series.setdefault(key, ([], [], tube, n1, result))
            (nodes if row["kind"] == "node" else stations).append(z)
    print("Celsius, mm, ml/h, variant, mode: rel N1, inlet p2, node p1/p2,"
          " station p1/p2, residual, interpolation-only p1/p2 (Pa)")
    for key, errors in comparisons.items():
        nodes, stations, tube, n1, result = series[key]
        z = np.array(sorted(nodes))
        p, x1 = result.sol(LENGTH - z)
        gradient = tube.rhs(n1, False)
        dp, dx1 = -np.array([gradient(LENGTH - zi, (pressure, x))
                             for zi, pressure, x in zip(z, p, x1)]).T
        dp1 = x1 * dp + p * dx1
        partials = np.column_stack((p * x1, p * (1 - x1)))
        cubic = CubicHermiteSpline(z, partials, np.column_stack((dp1, dp - dp1)))
        interpolation = np.max(abs(cubic(stations) - profile(result, stations)[:, :2]), axis=0)
        print(*key, *(f"{value:.5e}" for value in (*errors, *interpolation)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", type=Path)
    parser.add_argument("--study", type=Path)
    parser.add_argument("--pressures", type=Path)
    args = parser.parse_args()
    if args.fixture:
        args.fixture.write_text(
            ";; Generated by scripts/cpim-verification/reference.py using DOP853 and Brent.\n"
            + edn(fixture()) + "\n", encoding="utf-8",
        )
    if args.study:
        rows = grid()
        limits = []
        for index, tube in enumerate(lim_tubes()):
            n1, result = tube.solve(TIGHT)
            if index < 3:
                comparison_n1 = tube.molecular_flux
                z_from_outlet = LENGTH - STATIONS
                comparison_p = np.column_stack((
                    tube.saturation * z_from_outlet / LENGTH,
                    tube.outlet + tube.rt * tube.n2 / tube.dk2 * z_from_outlet,
                ))
            else:
                comparison_n1, comparison_result = tube.solve(TIGHT, continuum=True)
                comparison_p = profile(comparison_result, STATIONS)[:, :2]
            limits.append([tube.diameter, tube.feed, tube.outlet, n1, comparison_n1,
                           profile(result, STATIONS)[:, :2], comparison_p])
        args.study.write_text(edn({"grid": rows, "limits": limits}) + "\n", encoding="utf-8")
    if args.pressures:
        pressure_errors(args.pressures)


if __name__ == "__main__":
    main()
