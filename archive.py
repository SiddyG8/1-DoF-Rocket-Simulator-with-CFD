# Flight.__init__
self.cg = [f.calculate_centre_of_gravity(self.rocket, self.rocket.total_mass)]
self.cp = [self.rocket.static_center_of_pressure]
self.stability_calibers = [f.calculate_stability_caliber(self.rocket, self.rocket.total_mass)]

# Flight.simulate
cg = f.calculate_centre_of_gravity(self.rocket, mass)
cp = self.rocket.static_center_of_pressure
stability_caliber = f.calculate_stability_caliber(self.rocket, mass)

self.cg.append(cg)
self.cp.append(cp)
self.stability_calibers.append(stability_caliber)

# Rocket
@property
def static_center_of_pressure(self) -> float:
    """
    Calculate the static center of pressure of the rocket.
    Implementation derived from the Barrowman equations.
    """
    if not self.nose_cone or not self.fins:
        raise ValueError("Nose cone and fins must be defined to calculate the static center of pressure")

    # Rocket geometry
    L_N = self.nose_cone.length
    d = self.nose_cone.diameter
    C_R = self.fins.root_chord
    C_T = self.fins.tip_chord
    S = self.fins.semi_span
    X_R = self.fins.sweep_length
    L_F = sqrt(S**2 + (C_T/2 - C_R/2 + X_R)**2)
    R = self.diameter / 2
    X_B = self.body_length + L_N - self.fins.position
    N = self.fins.num_fins

    # Nose cone terms
    C_NN = 2
    X_N = NoseCone.X_N

    # Fin terms
    C_NF = (1 + R/(S + R)) * (4*N*(S/d)**2/(1 + sqrt(1 + (2*L_F/(C_R + C_T))**2)))
    X_F = X_B + X_R/3 * (C_R + 2*C_T)/(C_R + C_T) + 1/6*(C_R + C_T - C_R*C_T/(C_R + C_T))

    # Sum up coefficients
    C_NR = C_NN + C_NF

    # CP distance from the nose cone tip
    return (C_NN*X_N + C_NF*X_F)/C_NR

# functions.py
def calculate_centre_of_gravity(rocket, current_mass) -> float:
    X_N = rocket.nose_cone.length * rocket.nose_cone.X_N + rocket.body_length
    X_F = rocket.fins.position
    X_R = rocket.body_length / 2
    X_M = rocket.motor.length / 2

    # fuel_mass = rocket.mass - current_mass

    X_CG = (X_N * rocket.nose_cone.mass + X_F * rocket.fins.total_mass + X_R * rocket.total_mass + X_M * rocket.motor.mass) / current_mass

    return rocket.body_length + rocket.nose_cone.length - X_CG


def calculate_stability_caliber(rocket, current_mass) -> float:
    return (rocket.static_center_of_pressure - calculate_centre_of_gravity(rocket, current_mass)) / rocket.diameter
