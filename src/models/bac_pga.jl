# ============================================================
# Funciones mecánicas adimensionales para rod3D_grow
#
# Convención usada en estas rutinas:
#   - Las funciones devuelven VELOCIDADES, no fuerzas.
#   - La fuerza geométrica se divide por L = l + d antes de aplicar
#     la movilidad traslacional.
#   - El torque se construye con esa misma fuerza/L.
#   - Por ello, el prefactor angular 12*mu_rot/L^2 produce en total
#     la dependencia correcta 12*tau_raw/L^3.
# ============================================================

@inline function du_dtheta(theta, phi)
    return (
        -cos(phi) * sin(theta),
         cos(phi) * cos(theta),
         0.0,
    )
end

@inline function du_dphi(theta, phi)
    return (
        -sin(phi) * cos(theta),
        -sin(phi) * sin(theta),
         cos(phi),
    )
end

@inline function rod_dir3d(theta, phi)
    cp = cos(phi)
    return (
        cos(theta) * cp,
        sin(theta) * cp,
        sin(phi),
    )
end

@inline function rod_poles3d_ordered(x, y, z, l, theta, phi)
    ux, uy, uz = rod_dir3d(theta, phi)
    half_l = 0.5 * l

    return (
        (x + half_l * ux, y + half_l * uy, z + half_l * uz),
        (x - half_l * ux, y - half_l * uy, z - half_l * uz),
    )
end

# ------------------------------------------------------------
# Repulsión célula-célula
# ------------------------------------------------------------
function repulsiveForces_rods3d_asym(
    x, y, z, d, l, theta, phi,
    x2, y2, z2, d2, l2, theta2, phi2,
    eta, E, A;
    mu_rot = nothing,
)
    vix = 0.0
    viy = 0.0
    viz = 0.0
    omega_theta = 0.0
    omega_phi   = 0.0

    mu_par  = 1.0 / eta
    mu_perp = 1.0 / (eta * A)
    mu_rot === nothing && (mu_rot = mu_perp)

    ux, uy, uz = rod_dir3d(theta, phi)

    xi, yi, zi, xj, yj, zj =
        CBMMetrics.rodIntersection3d_(
            x, y, z, l, theta, phi,
            x2, y2, z2, l2, theta2, phi2,
        )

    dx = xi - xj
    dy = yi - yj
    dz = zi - zj
    rij2 = dx^2 + dy^2 + dz^2

    contact_dist = 0.5 * (d + d2)

    if rij2 > 0.0 && rij2 < contact_dist^2
        rij = sqrt(rij2)
        overlap = contact_dist - rij

        nx = dx / rij
        ny = dy / rij
        nz = dz / rij

        L = l + d

        # Fuerza geométrica dividida por L:
        # dr/dt = M*F/L
        d_eff = 0.5 * (d + d2)
        Fmag = E * sqrt(d_eff * overlap^3) / L

        Fx = Fmag * nx
        Fy = Fmag * ny
        Fz = Fmag * nz

        u_dot_F = ux * Fx + uy * Fy + uz * Fz

        vix = mu_perp * Fx + (mu_par - mu_perp) * u_dot_F * ux
        viy = mu_perp * Fy + (mu_par - mu_perp) * u_dot_F * uy
        viz = mu_perp * Fz + (mu_par - mu_perp) * u_dot_F * uz

        rx = xi - x
        ry = yi - y
        rz = zi - z

        tau_x = ry * Fz - rz * Fy
        tau_y = rz * Fx - rx * Fz
        tau_z = rx * Fy - ry * Fx

        duθx, duθy, duθz = du_dtheta(theta, phi)
        duφx, duφy, duφz = du_dphi(theta, phi)

        sθx = uy * duθz - uz * duθy
        sθy = uz * duθx - ux * duθz
        sθz = ux * duθy - uy * duθx

        sφx = uy * duφz - uz * duφy
        sφy = uz * duφx - ux * duφz
        sφz = ux * duφy - uy * duφx

        # Correcto aquí: tau ya contiene el factor 1/L.
        scale = 12.0 * mu_rot / L^2

        omega_theta = scale * (
            tau_x * sθx +
            tau_y * sθy +
            tau_z * sθz
        )

        omega_phi = scale * (
            tau_x * sφx +
            tau_y * sφy +
            tau_z * sφz
        )
    end

    return vix, viy, viz, omega_theta, omega_phi
end
# ------------------------------------------------------------
# Adhesión célula-sustrato mediante potencial Yukawa
#
# h = zp - zwall:
# distancia entre el centro del polo y el sustrato.
#
# La adhesión actúa únicamente cuando:
#     d/2 < h < rcut_h
#
# Por tanto, h nunca se evalúa cerca de cero y el potencial
# no presenta singularidad en la región activa.
#
# U(h) = -eps_cs * exp(-gamma*h) / h
# ------------------------------------------------------------
function substrateAttraction_rods_yukawa3d(
    x, y, z, d, l, theta, phi,
    eta, eps_cs, gamma, A;
    zwall = 0.0,
    rcut_h = 5.0,
    mu_rot = nothing,
)
    vx = 0.0
    vy = 0.0
    vz = 0.0

    omega_theta = 0.0
    omega_phi   = 0.0

    mu_par  = 1.0 / eta
    mu_perp = 1.0 / (eta * A)

    if mu_rot === nothing
        mu_rot = mu_perp
    end

    ux, uy, uz = rod_dir3d(theta, phi)

    pole_plus, pole_minus =
        rod_poles3d_ordered(
            x, y, z,
            l, theta, phi,
        )

    Fgx = 0.0
    Fgy = 0.0
    Fgz = 0.0

    tau_x = 0.0
    tau_y = 0.0
    tau_z = 0.0

    L = l + d
    h_contact = d / 2

    for (xp, yp, zp) in (pole_plus, pole_minus)

        # Distancia centro del polo–sustrato
        h = zp - zwall

        # Solo actúa fuera del contacto y dentro del cutoff
        if h_contact < h < rcut_h

            expterm = exp(-gamma * h)

            # U(h) = -eps_cs exp(-gamma*h)/h
            #
            # Fz = -dU/dh
            #    = -eps_cs exp(-gamma*h)
            #      (gamma/h + 1/h²)
            #
            # Fz < 0: atracción hacia el sustrato.
            Fz_p = -eps_cs * expterm *
                   (
                       gamma / h +
                       1.0 / h^2
                   ) / L

            Fgz += Fz_p

            rx = xp - x
            ry = yp - y
            rz = zp - z

            tau_x += ry * Fz_p
            tau_y -= rx * Fz_p
        end
    end

    # Movilidad traslacional anisótropa
    u_dot_F =
        ux * Fgx +
        uy * Fgy +
        uz * Fgz

    vx =
        mu_perp * Fgx +
        (mu_par - mu_perp) * u_dot_F * ux

    vy =
        mu_perp * Fgy +
        (mu_par - mu_perp) * u_dot_F * uy

    vz =
        mu_perp * Fgz +
        (mu_par - mu_perp) * u_dot_F * uz

    # Proyección angular
    duθx, duθy, duθz =
        du_dtheta(theta, phi)

    duφx, duφy, duφz =
        du_dphi(theta, phi)

    sθx = uy * duθz - uz * duθy
    sθy = uz * duθx - ux * duθz
    sθz = ux * duθy - uy * duθx

    sφx = uy * duφz - uz * duφy
    sφy = uz * duφx - ux * duφz
    sφz = ux * duφy - uy * duφx

    # El torque ya contiene el factor 1/L.
    scale = 12.0 * mu_rot / L^2

    omega_theta = scale * (
        tau_x * sθx +
        tau_y * sθy +
        tau_z * sθz
    )

    omega_phi = scale * (
        tau_x * sφx +
        tau_y * sφy +
        tau_z * sφz
    )

    return (
        vx,
        vy,
        vz,
        omega_theta,
        omega_phi,
    )
end

# ------------------------------------------------------------
# Repulsión célula-sustrato
# ------------------------------------------------------------
function substrateRepulsion_rods3d(
    x, y, z, d, l, theta, phi,
    eta, Ebv, A;
    zwall = 0.0,
    mu_rot = nothing,
)
    vx = 0.0
    vy = 0.0
    vz = 0.0
    omega_theta = 0.0
    omega_phi   = 0.0

    mu_par  = 1.0 / eta
    mu_perp = 1.0 / (eta * A)
    mu_rot === nothing && (mu_rot = mu_perp)

    ux, uy, uz = rod_dir3d(theta, phi)
    pole_plus, pole_minus =
        rod_poles3d_ordered(x, y, z, l, theta, phi)

    Fgx = 0.0
    Fgy = 0.0
    Fgz = 0.0

    tau_x = 0.0
    tau_y = 0.0
    tau_z = 0.0

    L = l + d

    for (xp, yp, zp) in (pole_plus, pole_minus)
        overlap = zwall + d / 2 - zp

        if overlap > 0.0
            Fz_p = Ebv * sqrt(d * overlap^3) / L
            Fgz += Fz_p

            rx = xp - x
            ry = yp - y
            rz = zp - z

            tau_x += ry * Fz_p
            tau_y -= rx * Fz_p
        end
    end

    u_dot_F = ux * Fgx + uy * Fgy + uz * Fgz

    vx = mu_perp * Fgx + (mu_par - mu_perp) * u_dot_F * ux
    vy = mu_perp * Fgy + (mu_par - mu_perp) * u_dot_F * uy
    vz = mu_perp * Fgz + (mu_par - mu_perp) * u_dot_F * uz

    duθx, duθy, duθz = du_dtheta(theta, phi)
    duφx, duφy, duφz = du_dphi(theta, phi)

    sθx = uy * duθz - uz * duθy
    sθy = uz * duθx - ux * duθz
    sθz = ux * duθy - uy * duθx

    sφx = uy * duφz - uz * duφy
    sφy = uz * duφx - ux * duφz
    sφz = ux * duφy - uy * duφx

    scale = 12.0 * mu_rot / L^2

    omega_theta = scale * (
        tau_x * sθx +
        tau_y * sθy +
        tau_z * sθz
    )

    omega_phi = scale * (
        tau_x * sφx +
        tau_y * sφy +
        tau_z * sφz
    )

    return vx, vy, vz, omega_theta, omega_phi
end

# ------------------------------------------------------------
# Repulsión célula-pad en un punto del eje
# ------------------------------------------------------------
function padRepulsion_rods_point_3d(
    x, y, z,
    xp, yp, zp,
    d, l, theta, phi,
    eta, Ebp, hPad, A;
    mu_rot = nothing,
)
    vx = 0.0
    vy = 0.0
    vz = 0.0
    omega_theta = 0.0
    omega_phi   = 0.0
    Fz_geom = 0.0

    mu_par  = 1.0 / eta
    mu_perp = 1.0 / (eta * A)
    mu_rot === nothing && (mu_rot = mu_perp)

    ux, uy, uz = rod_dir3d(theta, phi)

    Fgx = 0.0
    Fgy = 0.0
    Fgz = 0.0

    tau_x = 0.0
    tau_y = 0.0
    tau_z = 0.0

    overlap = zp + d / 2 - hPad
    L = l + d

    if overlap > 0.0
        # Densidad lineal geométrica antes de dividir por L.
        Fz_geom = -Ebp * sqrt(d * overlap^3)

        # Contribución a velocidad por unidad de longitud:
        # M*f/L
        Fgz = Fz_geom / L

        rx = xp - x
        ry = yp - y
        rz = zp - z

        tau_x = ry * Fgz
        tau_y = -rx * Fgz
    end

    u_dot_F = ux * Fgx + uy * Fgy + uz * Fgz

    vx = mu_perp * Fgx + (mu_par - mu_perp) * u_dot_F * ux
    vy = mu_perp * Fgy + (mu_par - mu_perp) * u_dot_F * uy
    vz = mu_perp * Fgz + (mu_par - mu_perp) * u_dot_F * uz

    duθx, duθy, duθz = du_dtheta(theta, phi)
    duφx, duφy, duφz = du_dphi(theta, phi)

    sθx = uy * duθz - uz * duθy
    sθy = uz * duθx - ux * duθz
    sθz = ux * duθy - uy * duθx

    sφx = uy * duφz - uz * duφy
    sφy = uz * duφx - ux * duφz
    sφz = ux * duφy - uy * duφx

    scale = 12.0 * mu_rot / L^2

    omega_theta = scale * (
        tau_x * sθx +
        tau_y * sθy +
        tau_z * sθz
    )

    omega_phi = scale * (
        tau_x * sφx +
        tau_y * sφy +
        tau_z * sφz
    )

    return (
        vx,
        vy,
        vz,
        omega_theta,
        omega_phi,
        Fz_geom,
    )
end
