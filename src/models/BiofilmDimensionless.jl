# Funciones mecánicas adimensionales usadas por rod3D_grow

@inline function rod_dir3d(theta, phi)
    cp = cos(phi)
    return (cos(theta) * cp, sin(theta) * cp, sin(phi))
end

@inline function du_dtheta(theta, phi)
    return (-cos(phi) * sin(theta), cos(phi) * cos(theta), 0.0)
end

@inline function du_dphi(theta, phi)
    return (-sin(phi) * cos(theta), -sin(phi) * sin(theta), cos(phi))
end

@inline function rod_poles3d_ordered(x, y, z, l, theta, phi)
    ux, uy, uz = rod_dir3d(theta, phi)
    h = 0.5 * l
    return ((x + h*ux, y + h*uy, z + h*uz),
            (x - h*ux, y - h*uy, z - h*uz))
end

@inline function anisotropic_velocity(Fx, Fy, Fz, ux, uy, uz, eta, A)
    mu_par  = 1.0 / eta
    mu_perp = 1.0 / (eta * A)
    uF = ux*Fx + uy*Fy + uz*Fz
    vx = mu_perp*Fx + (mu_par - mu_perp)*uF*ux
    vy = mu_perp*Fy + (mu_par - mu_perp)*uF*uy
    vz = mu_perp*Fz + (mu_par - mu_perp)*uF*uz
    return vx, vy, vz
end

@inline function angular_rates_from_torque(taux, tauy, tauz, theta, phi, L, mu_rot)
    ux, uy, uz = rod_dir3d(theta, phi)
    duθx, duθy, duθz = du_dtheta(theta, phi)
    duφx, duφy, duφz = du_dphi(theta, phi)

    sθx = uy*duθz - uz*duθy
    sθy = uz*duθx - ux*duθz
    sθz = ux*duθy - uy*duθx

    sφx = uy*duφz - uz*duφy
    sφy = uz*duφx - ux*duφz
    sφz = ux*duφy - uy*duφx

    scale = 12.0 * mu_rot / L^3
    Wθ = scale * (taux*sθx + tauy*sθy + tauz*sθz)
    Wφ = scale * (taux*sφx + tauy*sφy + tauz*sφz)
    return Wθ, Wφ
end

function repulsiveForces_rods3d_asym(
    x,y,z,d,l,theta,phi,
    x2,y2,z2,d2,l2,theta2,phi2,
    eta,Ebb,A; mu_rot=1.0/(eta*A))

    xi, yi, zi, xj, yj, zj = CBMMetrics.rodIntersection3d_(
        x,y,z,l,theta,phi, x2,y2,z2,l2,theta2,phi2)

    dx, dy, dz = xi-xj, yi-yj, zi-zj
    rij2 = dx^2 + dy^2 + dz^2
    contact = 0.5*(d+d2)
    if rij2 <= 0.0 || rij2 >= contact^2
        return 0.0,0.0,0.0,0.0,0.0
    end

    rij = sqrt(rij2)
    overlap = contact - rij
    nx, ny, nz = dx/rij, dy/rij, dz/rij
    deff = 0.5*(d+d2)
    Fmag = Ebb * sqrt(deff * overlap^3)
    Fx, Fy, Fz = Fmag*nx, Fmag*ny, Fmag*nz

    L = l + d
    ux, uy, uz = rod_dir3d(theta, phi)
    vx, vy, vz = anisotropic_velocity(Fx/L, Fy/L, Fz/L, ux, uy, uz, eta, A)

    rx, ry, rz = xi-x, yi-y, zi-z
    taux = ry*Fz - rz*Fy
    tauy = rz*Fx - rx*Fz
    tauz = rx*Fy - ry*Fx
    Wθ, Wφ = angular_rates_from_torque(taux,tauy,tauz,theta,phi,L,mu_rot)

    return vx,vy,vz,Wθ,Wφ
end

function substrateAttraction_rods_yukawa3d(
    x,y,z,d,l,theta,phi,
    eta,eps_cs,gamma_cs,A;
    mu_rot=1.0/(eta*A), rcut_h=5.0)

    p1, p2 = rod_poles3d_ordered(x,y,z,l,theta,phi)
    Fx=0.0; Fy=0.0; Fz=0.0
    taux=0.0; tauy=0.0; tauz=0.0

    for (xp,yp,zp) in (p1,p2)
        h = zp - d/2
        if 0.0 < h < rcut_h
            Fzp = -eps_cs * exp(-gamma_cs*h) * (gamma_cs/h + 1.0/h^2)
            Fz += Fzp
            rx, ry, rz = xp-x, yp-y, zp-z
            taux += ry*Fzp
            tauy -= rx*Fzp
        end
    end

    L = l + d
    ux, uy, uz = rod_dir3d(theta, phi)
    vx, vy, vz = anisotropic_velocity(Fx/L,Fy/L,Fz/L,ux,uy,uz,eta,A)
    Wθ, Wφ = angular_rates_from_torque(taux,tauy,tauz,theta,phi,L,mu_rot)
    return vx,vy,vz,Wθ,Wφ
end

function substrateRepulsion_rods3d(
    x,y,z,d,l,theta,phi,
    eta,Ebv,A; mu_rot=1.0/(eta*A))

    p1, p2 = rod_poles3d_ordered(x,y,z,l,theta,phi)
    Fx=0.0; Fy=0.0; Fz=0.0
    taux=0.0; tauy=0.0; tauz=0.0

    for (xp,yp,zp) in (p1,p2)
        overlap = d/2 - zp
        if overlap > 0.0
            Fzp = Ebv * sqrt(d * overlap^3)
            Fz += Fzp
            rx, ry, rz = xp-x, yp-y, zp-z
            taux += ry*Fzp
            tauy -= rx*Fzp
        end
    end

    L = l + d
    ux, uy, uz = rod_dir3d(theta, phi)
    vx, vy, vz = anisotropic_velocity(Fx/L,Fy/L,Fz/L,ux,uy,uz,eta,A)
    Wθ, Wφ = angular_rates_from_torque(taux,tauy,tauz,theta,phi,L,mu_rot)
    return vx,vy,vz,Wθ,Wφ
end

function padRepulsion_rods_point_3d(
    x,y,z,xp,yp,zp,d,l,theta,phi,
    eta,Ebp,hPad,A; mu_rot=1.0/(eta*A))

    overlap = zp + d/2 - hPad
    if overlap <= 0.0
        return 0.0,0.0,0.0,0.0,0.0,0.0
    end

    Fz_geom = -Ebp * sqrt(d * overlap^3)
    Fx, Fy, Fz = 0.0, 0.0, Fz_geom

    L = l + d
    ux, uy, uz = rod_dir3d(theta, phi)
    vx, vy, vz = anisotropic_velocity(Fx/L,Fy/L,Fz/L,ux,uy,uz,eta,A)

    rx, ry, rz = xp-x, yp-y, zp-z
    taux = ry*Fz
    tauy = -rx*Fz
    tauz = 0.0
    Wθ, Wφ = angular_rates_from_torque(taux,tauy,tauz,theta,phi,L,mu_rot)

    return vx,vy,vz,Wθ,Wφ,Fz_geom
end
