# --- derivadas del eje u respecto a los ángulos ---
# u = (cosφ cosθ, cosφ sinθ, sinφ)
@inline function du_dtheta(theta, phi)
    (-cos(phi)*sin(theta),  cos(phi)*cos(theta),  0.0)
end
@inline function du_dphi(theta, phi)
    (-sin(phi)*cos(theta), -sin(phi)*sin(theta),  cos(phi))
end

@inline dot3(a,b) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]


# --- Utilidades de orientación y polos (orden fijo: + y −) ---
@inline function rod_dir3d(theta, phi)
    cx = cos(theta); sx = sin(theta)
    cp = cos(phi);   sp = sin(phi)
    return (cx*cp, sx*cp, sp)               # u = (cosθ cosφ, sinθ cosφ, sinφ)
end


@inline function rod_poles3d_ordered(x,y,z,l,theta,phi)
    ux,uy,uz = rod_dir3d(theta, phi)
    hx,hy,hz = 0.5l*ux, 0.5l*uy, 0.5l*uz
    return (x+hx, y+hy, z+hz), (x-hx, y-hy, z-hz)  # (+, −)
end



# function repulsiveForces_rods3d_asym(
#     x,y,z,d,l,theta,phi,
#     x2,y2,z2,d2,l2,theta2,phi2,
#     eta, E, A,mu_rot = nothing)

#     Fijx = 0.0; Fijy = 0.0; Fijz = 0.0
#     Wijtheta = 0.0; Wijphi = 0.0

#     # Puntos más cercanos
#     xiAux, yiAux, ziAux,  xjAux, yjAux, zjAux =
#         CBMMetrics.rodIntersection3d_(x, y, z, l, theta, phi,
#                                      x2, y2, z2, l2, theta2, phi2)

#     dx = xiAux - xjAux
#     dy = yiAux - yjAux
#     dz = ziAux - zjAux
#     rij = sqrt(dx*dx + dy*dy + dz*dz)

#     # Dirección del rod i
#     ux =  cos(theta)*cos(phi)
#     uy =  sin(theta)*cos(phi)
#     uz =  sin(phi)
#     # (por si acaso) — debería estar ya normalizado
#     # un = inv(sqrt(ux^2+uy^2+uz^2)); ux*=un; uy*=un; uz*=un

#     if rij > 0.0 && rij < (d + d2)/2
#         hAux = (d + d2)/2 - rij
#         nijx = dx / rij
#         nijy = dy / rij
#         nijz = dz / rij

#         # Fuerza "cruda" de contacto SIN /eta
#         Fraw_mag = E * sqrt(d2 * hAux^3) / (l + d)
#         Fx_raw = Fraw_mag * nijx
#         Fy_raw = Fraw_mag * nijy
#         Fz_raw = Fraw_mag * nijz

#         # === 1) TRASLACIÓN con movilidad anisótropa ===
#         mu_par  = 1.0 / eta  # μ∥
#         mu_perp = 1.0 / (eta*A)         # μ⊥
#         u_dot_Fraw = ux*Fx_raw + uy*Fy_raw + uz*Fz_raw
#         Fasx = mu_perp*Fx_raw + (mu_par - mu_perp)*u_dot_Fraw*ux
#         Fasy = mu_perp*Fy_raw + (mu_par - mu_perp)*u_dot_Fraw*uy
#         Fasz = mu_perp*Fz_raw + (mu_par - mu_perp)*u_dot_Fraw*uz

#         Fijx = Fasx;  Fijy = Fasy;  Fijz = Fasz

#         # === 2) ROTACIÓN con fuerza ISOTRÓPICA (no anisotropa) ===
#         # fuerza isotrópica = Fraw / eta
#         Fx_iso = Fx_raw
#         Fy_iso = Fy_raw
#         Fz_iso = Fz_raw 

#         rx = (xiAux - x); ry = (yiAux - y); rz = (ziAux - z)
#         τx = ry*Fz_iso - rz*Fy_iso
#         τy = rz*Fx_iso - rx*Fz_iso
#         τz = rx*Fy_iso - ry*Fx_iso

#         # Proyección a (θ, φ)
#         duθx, duθy, duθz = du_dtheta(theta, phi)
#         duφx, duφy, duφz = du_dphi(theta,  phi)

#         sθx = uy*duθz - uz*duθy
#         sθy = uz*duθx - ux*duθz
#         sθz = ux*duθy - uy*duθx

#         sφx = uy*duφz - uz*duφy
#         sφy = uz*duφx - ux*duφz
#         sφz = ux*duφy - uy*duφx

#         # Escala rotacional del caso 3D isotrópico
#         scale = 12.0 / ((l + d)^2)*mu_perp

#         Wijtheta = scale * (τx*sθx + τy*sθy + τz*sθz)
#         Wijphi   = scale * (τx*sφx + τy*sφy + τz*sφz)
#     end

#     return Fijx, Fijy, Fijz, Wijtheta, Wijphi
# end

function repulsiveForces_rods3d_asym(
    x, y, z, d, l, theta, phi,
    x2, y2, z2, d2, l2, theta2, phi2,
    eta, E, A; mu_rot = nothing)

    # -----------------------------
    # Salidas: contribuciones a las velocidades
    # -----------------------------
    vix = 0.0
    viy = 0.0
    viz = 0.0
    ωtheta = 0.0
    ωphi   = 0.0

    # -----------------------------
    # Movilidades
    # A > 1  => moverse perpendicular al eje largo cuesta más
    # -----------------------------
    mu_par  = 1.0 / eta
    mu_perp = 1.0 / (eta * A)

    # Si no se especifica, tomamos una movilidad rotacional efectiva
    # igual a la perpendicular (aproximación razonable y simple)
    if mu_rot === nothing
        mu_rot = mu_perp
    end

    # -----------------------------
    # Dirección unitaria del rod i
    # u(theta, phi) = (cosθ cosφ, sinθ cosφ, sinφ)
    # -----------------------------
    ux = cos(theta) * cos(phi)
    uy = sin(theta) * cos(phi)
    uz = sin(phi)

    # -----------------------------
    # Puntos más cercanos entre ambos segmentos
    # -----------------------------
    xiAux, yiAux, ziAux, xjAux, yjAux, zjAux =
        CBMMetrics.rodIntersection3d_(
            x, y, z, l, theta, phi,
            x2, y2, z2, l2, theta2, phi2
        )

    # Vector desde j hacia i en el punto de contacto
    dx = xiAux - xjAux
    dy = yiAux - yjAux
    dz = ziAux - zjAux
    rij = sqrt(dx*dx + dy*dy + dz*dz)

    # -----------------------------
    # Contacto repulsivo
    # -----------------------------
    contact_dist = 0.5 * (d + d2)

    if rij > 0.0 && rij < contact_dist

        # Solapamiento efectivo
        hAux = contact_dist - rij

        # Normal unitaria j -> i
        nijx = dx / rij
        nijy = dy / rij
        nijz = dz / rij

        # -----------------------------
        # 1) Fuerza geométrica de contacto
        # -----------------------------
        # Mantengo tu ley original, pero SIN movilidad aquí
        # porque en overdamped primero va la "fuerza",
        # y luego la movilidad la convierte en velocidad.
        Fmag = E * sqrt(d2 * hAux^3) / (l + d)

        Fx = Fmag * nijx
        Fy = Fmag * nijy
        Fz = Fmag * nijz

        # -----------------------------
        # 2) Velocidad traslacional anisótropa
        # v = M_t * F
        # M_t = mu_perp I + (mu_par - mu_perp) u u^T
        # -----------------------------
        u_dot_F = ux*Fx + uy*Fy + uz*Fz

        vix = mu_perp * Fx + (mu_par - mu_perp) * u_dot_F * ux
        viy = mu_perp * Fy + (mu_par - mu_perp) * u_dot_F * uy
        viz = mu_perp * Fz + (mu_par - mu_perp) * u_dot_F * uz

        # -----------------------------
        # 3) Torque geométrico respecto al centro de i
        # tau = r_contact x F
        # -----------------------------
        rx = xiAux - x
        ry = yiAux - y
        rz = ziAux - z

        τx = ry*Fz - rz*Fy
        τy = rz*Fx - rx*Fz
        τz = rx*Fy - ry*Fx

        # -----------------------------
        # 4) Proyección del torque sobre (theta, phi)
        # -----------------------------
        duθx, duθy, duθz = du_dtheta(theta, phi)
        duφx, duφy, duφz = du_dphi(theta, phi)

        # sθ = u × (du/dθ)
        sθx = uy*duθz - uz*duθy
        sθy = uz*duθx - ux*duθz
        sθz = ux*duθy - uy*duθx

        # sφ = u × (du/dφ)
        sφx = uy*duφz - uz*duφy
        sφy = uz*duφx - ux*duφz
        sφz = ux*duφy - uy*duφx

        # -----------------------------
        # 5) Velocidad angular generalizada
        # ω_gen ~ mu_rot * torque_projected
        # Mantengo tu prefactor geométrico
        # -----------------------------
        scale = 12.0 * mu_rot / ((l + d)^2)

        ωtheta = scale * (τx*sθx + τy*sθy + τz*sθz)
        ωphi   = scale * (τx*sφx + τy*sφy + τz*sφz)
    end

    return vix, viy, viz, ωtheta, ωphi
end

function substrateAttraction_rods_yukawa3d(
    x, y, z, d, l, theta, phi,
    eta, eps_cs, gamma, A;
    rcut_h = 5.0, mu_rot = nothing
)
    Fx = 0.0; Fy = 0.0; Fz = 0.0
    Wθ = 0.0; Wφ = 0.0

    # Movilidades overdamped
    mu_par  = 1.0 / eta
    mu_perp = 1.0 / (eta * A)

    if mu_rot === nothing
        mu_rot = mu_perp
    end

    # Dirección del rod
    ux = cos(theta) * cos(phi)
    uy = sin(theta) * cos(phi)
    uz = sin(phi)

    a = 0.0  # plano en z=0

    # Polos del rod
    (x1, y1, z1), (x2, y2, z2) = rod_poles3d_ordered(x, y, z, l, theta, phi)

    # Acumuladores de fuerza geométrica y torque geométrico
    Fgx = 0.0; Fgy = 0.0; Fgz = 0.0
    τx = 0.0;  τy = 0.0;  τz = 0.0

    for pole in 1:2
        xp, yp, zp = pole == 1 ? (x1, y1, z1) : (x2, y2, z2)

        # MISMA definición que en tu versión que funcionaba
        s = zp
        h = s - a + d/2

        if h > 0.0 && h < rcut_h
            # MISMA ley que antes
            f    = exp(gamma * h) / h
            dUdh = eps_cs * f * (gamma - 1.0 / h)

            # Fuerza geométrica (sin movilidad aún)
            Fz_p = -dUdh / (l + d)
            Fx_p = 0.0
            Fy_p = 0.0

            Fgx += Fx_p
            Fgy += Fy_p
            Fgz += Fz_p

            # Torque geométrico respecto al CM
            rx = xp - x
            ry = yp - y
            rz = zp - z

            τx += ry * Fz_p - rz * Fy_p
            τy += rz * Fx_p - rx * Fz_p
            τz += rx * Fy_p - ry * Fx_p
        end
    end

    # Movilidad traslacional anisótropa
    u_dot_F = ux * Fgx + uy * Fgy + uz * Fgz

    Fx = mu_perp * Fgx + (mu_par - mu_perp) * u_dot_F * ux
    Fy = mu_perp * Fgy + (mu_par - mu_perp) * u_dot_F * uy
    Fz = mu_perp * Fgz + (mu_par - mu_perp) * u_dot_F * uz

    # Proyección angular
    duθx, duθy, duθz = du_dtheta(theta, phi)
    duφx, duφy, duφz = du_dphi(theta, phi)

    sθx = uy * duθz - uz * duθy
    sθy = uz * duθx - ux * duθz
    sθz = ux * duθy - uy * duθx

    sφx = uy * duφz - uz * duφy
    sφy = uz * duφx - ux * duφz
    sφz = ux * duφy - uy * duφx

    scale = 12.0 * mu_rot / ((l + d)^2)
    Wθ = scale * (τx * sθx + τy * sθy + τz * sθz)
    Wφ = scale * (τx * sφx + τy * sφy + τz * sφz)

    return Fx, Fy, Fz, Wθ, Wφ
end


function substrateRepulsion_rods3d(
    x, y, z, d, l, theta, phi,
    eta, Ebv, A;
    zwall = 0.0, mu_rot = nothing
)
    # --------------------------------
    # Salidas: contribuciones a velocidades
    # --------------------------------
    vx = 0.0
    vy = 0.0
    vz = 0.0
    Wθ = 0.0
    Wφ = 0.0

    # --------------------------------
    # Movilidades
    # --------------------------------
    mu_par  = 1.0 / eta
    mu_perp = 1.0 / (eta * A)

    if mu_rot === nothing
        mu_rot = mu_perp
    end

    # --------------------------------
    # Dirección del rod
    # --------------------------------
    ux = cos(theta) * cos(phi)
    uy = sin(theta) * cos(phi)
    uz = sin(phi)

    # --------------------------------
    # Polos del rod
    # --------------------------------
    (x1, y1, z1), (x2, y2, z2) = rod_poles3d_ordered(x, y, z, l, theta, phi)

    # --------------------------------
    # Acumuladores geométricos
    # --------------------------------
    Fxg = 0.0
    Fyg = 0.0
    Fzg = 0.0

    τx = 0.0
    τy = 0.0
    τz = 0.0

    # --------------------------------
    # Repulsión del sustrato en cada polo
    # --------------------------------
    for pole in 1:2
        xp, yp, zp = pole == 1 ? (x1, y1, z1) : (x2, y2, z2)

        # Solapamiento con el plano z = zwall
        δ = zwall + d/2 - zp

        if δ > 0.0
            # Fuerza geométrica repulsiva tipo Hertz
            Fh = Ebv * sqrt(d * δ^3) / (l + d)

            # Fuerza vertical (+z)
            Fx_p = 0.0
            Fy_p = 0.0
            Fz_p = Fh

            Fxg += Fx_p
            Fyg += Fy_p
            Fzg += Fz_p

            # Torque geométrico respecto al CM
            rx = xp - x
            ry = yp - y
            rz = zp - z

            τx += ry * Fz_p - rz * Fy_p
            τy += rz * Fx_p - rx * Fz_p
            τz += rx * Fy_p - ry * Fx_p
        end
    end

    # --------------------------------
    # Velocidad traslacional anisótropa
    # v = M_t * F
    # --------------------------------
    u_dot_F = ux * Fxg + uy * Fyg + uz * Fzg

    vx = mu_perp * Fxg + (mu_par - mu_perp) * u_dot_F * ux
    vy = mu_perp * Fyg + (mu_par - mu_perp) * u_dot_F * uy
    vz = mu_perp * Fzg + (mu_par - mu_perp) * u_dot_F * uz

    # --------------------------------
    # Proyección del torque a (θ, φ)
    # --------------------------------
    duθx, duθy, duθz = du_dtheta(theta, phi)
    duφx, duφy, duφz = du_dphi(theta, phi)

    # sθ = u × du/dθ
    sθx = uy * duθz - uz * duθy
    sθy = uz * duθx - ux * duθz
    sθz = ux * duθy - uy * duθx

    # sφ = u × du/dφ
    sφx = uy * duφz - uz * duφy
    sφy = uz * duφx - ux * duφz
    sφz = ux * duφy - uy * duφx

    scale = 12.0 * mu_rot / ((l + d)^2)

    Wθ = scale * (τx * sθx + τy * sθy + τz * sθz)
    Wφ = scale * (τx * sφx + τy * sφy + τz * sφz)

    return vx, vy, vz, Wθ, Wφ
end



function padRepulsion_rods_point_3d(
    x, y, z, xp, yp, zp, d, l, theta, phi,
    eta, Ebp, hPad, A; mu_rot = nothing
)
    Fx = 0.0; Fy = 0.0; Fz = 0.0
    Wθ = 0.0; Wφ = 0.0
    Fz_geom = 0.0

    mu_par  = 1.0 / eta
    mu_perp = 1.0 / (eta * A)

    if mu_rot === nothing
        mu_rot = mu_perp
    end

    ux = cos(theta) * cos(phi)
    uy = sin(theta) * cos(phi)
    uz = sin(phi)

    Fgx = 0.0; Fgy = 0.0; Fgz = 0.0
    τx = 0.0; τy = 0.0; τz = 0.0

    δ = zp + d/2 - hPad

    if δ > 0.0
        Fz_geom = -Ebp * sqrt(d * δ^3) / (l + d)

        Fgx = 0.0
        Fgy = 0.0
        Fgz = Fz_geom

        rx = xp - x
        ry = yp - y
        rz = zp - z

        τx = ry * Fgz
        τy = -rx * Fgz
        τz = 0.0
    end

    u_dot_F = ux * Fgx + uy * Fgy + uz * Fgz

    Fx = mu_perp * Fgx + (mu_par - mu_perp) * u_dot_F * ux
    Fy = mu_perp * Fgy + (mu_par - mu_perp) * u_dot_F * uy
    Fz = mu_perp * Fgz + (mu_par - mu_perp) * u_dot_F * uz

    duθx, duθy, duθz = du_dtheta(theta, phi)
    duφx, duφy, duφz = du_dphi(theta, phi)

    sθx = uy*duθz - uz*duθy
    sθy = uz*duθx - ux*duθz
    sθz = ux*duθy - uy*duθx

    sφx = uy*duφz - uz*duφy
    sφy = uz*duφx - ux*duφz
    sφz = ux*duφy - uy*duφx

    scale = 12.0 * mu_rot / ((l + d)^2)
    Wθ = scale * (τx*sθx + τy*sθy + τz*sθz)
    Wφ = scale * (τx*sφx + τy*sφy + τz*sφz)

    return Fx, Fy, Fz, Wθ, Wφ, Fz_geom
end

