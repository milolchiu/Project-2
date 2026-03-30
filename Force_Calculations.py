import math

# Constants
g = 9.81
airDensity = 1.225
dragCoefficient = 0.47
dynamicViscosity = 1.81e-5
sigma = -1e-6
epsilon_0 = 8.854e-12
e = 1.602e-19
H = 7
dt = 0.001      
airVelocity = 1.0
W = 3.5
air_cycle_time = H / airVelocity
tower_volume = H * W * W

def simulate_particle(diam, particleConcentration, particleDensity, y0, x0, particleCharge):
    R = diam / 2
    V = (4/3) * math.pi * R**3
    m = particleDensity * V
    mass_conc = particleConcentration * 1e-9
    number_conc = mass_conc / m
    x, y = x0, y0
    vy = 0.0
    t = 0.0
    xs, ys = [], []
    while t < air_cycle_time and x > 0:
        Fdy = -6 * math.pi * dynamicViscosity * R * vy
        Fg = -m * g
        Fb = airDensity * V * g
        E_plate = sigma / (2 * epsilon_0)
        F_plate = particleCharge * E_plate
        E_cloud = (particleCharge * number_conc / (2 * epsilon_0)) * (2 * x - W)
        F_cloud = particleCharge * E_cloud
        vx = (F_plate - F_cloud) / (6 * math.pi * dynamicViscosity * R)
        ay = (Fg + Fb + Fdy) / m
        vy += ay * dt
        x += vx * dt
        y += vy * dt
        t += dt
        xs.append(x)
        ys.append(y)
    hit_plate = x <= 0
    return xs, ys, hit_plate, t

def get_removal_percentage(diam, particleConcentration, particleDensity, particleCharge, n_samples=50):
    if particleCharge == 0:
        return 0.0, []   # gas molecules cannot be ESP-captured
    x_starts = [W * (i + 0.5) / n_samples for i in range(n_samples)]
    y_start = H * (5/7)
    captured = 0
    all_trajectories = []
    for x0 in x_starts:
        xs, ys, hit_plate, t = simulate_particle(diam, particleConcentration, particleDensity, y_start, x0, particleCharge)
        all_trajectories.append((xs, ys, hit_plate))
        if hit_plate:
            captured += 1
    return (captured / n_samples) * 100, all_trajectories

def number_removed_per_cycle(diam, particleConcentration, particleDensity, particleCharge):
    R = diam / 2
    V_particle = (4/3) * math.pi * R**3
    m = particleDensity * V_particle
    mass_conc = particleConcentration * 1e-9
    number_conc = mass_conc / m
    total_particles = number_conc * tower_volume
    removal_pct, _ = get_removal_percentage(diam, particleConcentration, particleDensity, particleCharge)
    return (removal_pct / 100) * total_particles

def number_removed_per_day(diam, particleConcentration, particleDensity, particleCharge):
    cycles_per_day = 24 * 3600 / air_cycle_time
    return number_removed_per_cycle(diam, particleConcentration, particleDensity, particleCharge) * cycles_per_day

pollutants = [
    {"name": "sulfur_dioxide",      "diam (nm)": 0.36,  "ug": 4.5,   "density": 2.26,  "charge_e": 0},
    {"name": "nitrogen_oxide",      "diam (nm)": 0.32,  "ug": 65,    "density": 1.34,  "charge_e": 0},
    {"name": "nitrogen_dioxide",    "diam (nm)": 0.33,  "ug": 63.2,  "density": 1.88,  "charge_e": 0},
    {"name": "carbon_monoxide",     "diam (nm)": 0.38,  "ug": 557.5, "density": 1.145, "charge_e": 0},
    {"name": "ozone_particles",     "diam (nm)": 0.30,  "ug": 7.1,   "density": 1.96,  "charge_e": 0},
    {"name": "suspended_particles", "diam (nm)": 10000, "ug": 17.1,  "density": 2000,  "charge_e": 71000},
    {"name": "fine_particles",      "diam (nm)": 2500,  "ug": 23.7,  "density": 1600,  "charge_e": 8000},
]

cycles_per_day = 24 * 3600 / air_cycle_time

print("=" * 90)
print(f"  Cycle time: {air_cycle_time:.2f} s  |  Cycles per day: {cycles_per_day:.1f}  |  Tower volume: {tower_volume:.1f} m^3")
print("=" * 90)
print(f"  {'Pollutant':<25} {'Removed/cycle':>15} {'Particles/cycle':>18} {'Particles/day':>18}")
print("-" * 90)
for p in pollutants:
    if p["name"] == "ozone_particles":
        removal_pct = 0.0
        per_day = 0.0
        per_cycle = per_day / cycles_per_day
    else:
        diam_m = p["diam (nm)"] * 1e-9
        particleCharge = p["charge_e"] * e
        removal_pct, _ = get_removal_percentage(diam_m, p["ug"], p["density"], particleCharge)
        per_cycle = number_removed_per_cycle(diam_m, p["ug"], p["density"], particleCharge)
        per_day = number_removed_per_day(diam_m, p["ug"], p["density"], particleCharge)
        print(f"  {p['name']:<25} {removal_pct:>14.1f}% {per_cycle:>18.3e} {per_day:>18.3e}")
print("=" * 90)

