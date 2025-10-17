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

"""
    repulsiveForces(
        x, y, d, l, theta, etab, type,
        x2, y2, d2, l2, theta2, etag, type2,
        Ebb, Ebg, Egg
    ) -> (Fijx, Fijy, Wij)

Computes the repulsive interaction forces between two agents (either rods or gels), based on their type and geometry. If both are of the same type and are gels, direct contact is assumed; otherwise, virtual contact points are calculated.

### Parameters
- `x, y`: coordinates of the first agent.
- `d, l, theta`: diameter, length, and orientation angle of the first agent.
- `etab, type`: effective viscosity and type of the first agent (`0` for rod, `1` for gel).
- `x2, y2`: coordinates of the second agent.
- `d2, l2, theta2`: diameter, length, and orientation angle of the second agent.
- `etag, type2`: effective viscosity and type of the second agent.
- `Ebb, Ebg, Egg`: elastic moduli for rod-rod, rod-gel, and gel-gel interactions.

### Returns
- `Fijx, Fijy`: components of the repulsive force between the agents.
- `Wij`: radial torque-like contribution (only non-zero when the first agent is a rod).
"""
function repulsiveForces(
    x,y,d,l,theta,etab,type,
    x2,y2,d2,l2,theta2,etag,type2, Ebb, Ebg, Egg )

    Fijx = 0.
    Fijy = 0.
    Wij = 0.

    if type==type2 && type==1
        xiAux,yiAux,xjAux,yjAux=x,y,x2,y2
    else
        #Function that finds the virtual spheres of contact between both rods
        xiAux,yiAux,xjAux,yjAux = CBMMetrics.rodIntersection(x,y,l,theta,x2,y2,l2,theta2)
    end

    #Compute distance between virtual spheres
    rij = sqrt((xiAux-xjAux)^2 +(yiAux-yjAux)^2)
    if rij > 0. && rij < (d+d2)/2. #If it is smaller than a diameter compute forces
        #Compute auxiliar
        hAux = (d+d2)/2. - rij
        #Compute direction
        nijx = (xiAux-xjAux)/rij
        nijy = (yiAux-yjAux)/rij
        if type == 0
            if type2 == 0
                FnAux = Ebb * sqrt(d2* hAux^3)  / (etab * (l + d))
            else
                FnAux = Ebg * sqrt(d2* hAux^3) / (etab * (l + d))
            end
            Fijx = FnAux * nijx
            Fijy = FnAux * nijy
            Wij = ((xiAux-x)*Fijy - (yiAux-y)*Fijx) * 12. / (etab * (l + d)^3)
        else
            if type2 == 0
                FnAux = Ebg * sqrt(d2* hAux^3) / (etag * d)
            else
                FnAux = Egg * sqrt(d2* hAux^3) / (etag * d)
            end
            Fijx = FnAux * nijx
            Fijy = FnAux * nijy
        end
    end

    return Fijx, Fijy, Wij

end





function repulsiveForces_rods(
    x,y,d,l,theta,
    x2,y2,d2,l2,theta2,eta, E)

    Fijx = 0.
    Fijy = 0.
    Wij = 0.


    #Function that finds the virtual spheres of contact between both rods
    xiAux,yiAux,xjAux,yjAux = CBMMetrics.rodIntersection(x,y,l,theta,x2,y2,l2,theta2)
   

    #Compute distance between virtual spheres
    rij = sqrt((xiAux-xjAux)^2 +(yiAux-yjAux)^2)
    if rij > 0. && rij < (d+d2)/2. #If it is smaller than a diameter compute forces
        #Compute auxiliar
        hAux = (d+d2)/2. - rij
        #Compute direction
        nijx = (xiAux-xjAux)/rij
        nijy = (yiAux-yjAux)/rij
    
        FnAux = E * sqrt(d2* hAux^3)  / (eta * (l + d))
           
        Fijx = FnAux * nijx
        Fijy = FnAux * nijy
        Wij = ((xiAux-x)*Fijy - (yiAux-y)*Fijx) * 12. / ((l + d)^2)
      
    end

    return Fijx, Fijy, Wij

end



function repulsiveForces_rods3d(
    x,y,z,d,l,theta,phi,
    x2,y2,z2,d2,l2,theta2,phi2, eta, E)

    Fijx = 0.0; Fijy = 0.0; Fijz = 0.0
    Wijtheta = 0.0; Wijphi = 0.0

    # 1) Puntos de "contacto" (más cercanos) sobre cada segmento
    xiAux, yiAux, ziAux,  xjAux, yjAux, zjAux =
        CBMMetrics.rodIntersection3d(x, y, z, l, theta, phi,
                          x2, y2, z2, l2, theta2, phi2)

    # 2) Distancia y normal j->i
    dx = xiAux - xjAux
    dy = yiAux - yjAux
    dz = ziAux - zjAux
    rij = sqrt(dx*dx + dy*dy + dz*dz)

    if rij > 0.0 && rij < (d + d2)/2
        hAux = (d + d2)/2 - rij
        nijx = dx / rij
        nijy = dy / rij
        nijz = dz / rij

        # Misma ley de contacto que en 2D
        FnAux = E * sqrt(d2 * hAux^3) / (eta * (l + d))

        # Fuerza sobre i debida a j
        Fijx = FnAux * nijx
        Fijy = FnAux * nijy
        Fijz = FnAux * nijz

        # 3) Torque en el centro de i: τ = r_i × F_ij
        rx = (xiAux - x); ry = (yiAux - y); rz = (ziAux - z)
        τx = ry*Fijz - rz*Fijy
        τy = rz*Fijx - rx*Fijz
        τz = rx*Fijy - ry*Fijx

        # 4) Proyección en las direcciones generalizadas (d u / dθ, d u / dφ)
        # u(θ,φ)
        ux =  cos(theta)*cos(phi)
        uy =  sin(theta)*cos(phi)
        uz =  sin(phi)

        # du/dθ y du/dφ (tus funciones)
        duθx, duθy, duθz = du_dtheta(theta, phi)
        duφx, duφy, duφz = du_dphi(theta,  phi)

        # sθ = u × du/dθ
        sθx = uy*duθz - uz*duθy
        sθy = uz*duθx - ux*duθz
        sθz = ux*duθy - uy*duθx

        # sφ = u × du/dφ
        sφx = uy*duφz - uz*duφy
        sφy = uz*duφx - ux*duφz
        sφz = ux*duφy - uy*duφx

        # Proyección correcta del par
        scale = 12.0 / ((l + d)^2)
        Wijtheta = scale * (τx*sθx + τy*sθy + τz*sθz)
        Wijphi   = scale * (τx*sφx + τy*sφy + τz*sφz)

    end

    return Fijx, Fijy, Fijz, Wijtheta, Wijphi
end


function repulsiveForces_as_rods(
    x,y,d,l,theta,
    x2,y2,d2,l2,theta2,eta, E, A)

    Fasx = 0.
    Fasy = 0.
    Wij = 0.


    #Function that finds the virtual spheres of contact between both rods
    xiAux,yiAux,xjAux,yjAux = CBMMetrics.rodIntersection(x,y,l,theta,x2,y2,l2,theta2)
   
    c=cos(theta)
    s=sin(theta)
    #Compute distance between virtual spheres
    rij = sqrt((xiAux-xjAux)^2 +(yiAux-yjAux)^2)

    mu_par  = 1/(eta*A)   # μ∥
    mu_perp = A/eta       # μ⊥
    
    if rij > 0. && rij < (d+d2)/2. #If it is smaller than a diameter compute forces
        #Compute auxiliar
        hAux = (d+d2)/2. - rij
        #Compute direction
        nijx = (xiAux-xjAux)/rij
        nijy = (yiAux-yjAux)/rij
    
        FnAux = E * sqrt(d2* hAux^3)  /  (l + d)
           
        Fijx = FnAux * nijx
        Fijy = FnAux * nijy
        Fasx = (mu_par*c^2 + mu_perp*s^2)*Fijx + (mu_par - mu_perp)*c*s*Fijy
        Fasy = (mu_par - mu_perp)*c*s*Fijx     + (mu_par*s^2 + mu_perp*c^2)*Fijy
        Wij = ((xiAux-x)*Fijy - (yiAux-y)*Fijx) * 12. / (eta/A*(l + d)^2)
      
    end

    return Fasx, Fasy, Wij

end

function attractiveForces_rods(
    x,y,d,l,theta,
    x2,y2,d2,l2,theta2,eta, eps)

    Fijx = 0.
    Fijy = 0.
    Wij = 0.

    #Function that finds the virtual spheres of contact between both rods
    xiAux,yiAux,xjAux,yjAux = CBMMetrics.rodIntersection(x,y,l,theta,x2,y2,l2,theta2)
   

    #Compute distance between virtual spheres
    rij = sqrt((xiAux-xjAux)^2 +(yiAux-yjAux)^2)
    if rij > (d+d2)/2. && rij < 5*(d+d2)/2. #If it is smaller than a diameter compute forces
       
        #Compute direction
        nijx = (xiAux-xjAux)/rij
        nijy = (yiAux-yjAux)/rij
    
        FnAux = -24*eps* ((d+d2)/2)^6/rij^7 / (eta * (l + d))
           
        Fijx = FnAux * nijx
        Fijy = FnAux * nijy
        Wij = ((xiAux-x)*Fijy - (yiAux-y)*Fijx) * 12. / ((l + d)^2)
      
    end

    return Fijx, Fijy, Wij

end

function attractiveForces_rods_yukawa(
    x,y,d,l,theta,
    x2,y2,d2,l2,theta2,
    eta, eps, gamma)

    Fijx = 0.0
    Fijy = 0.0
    Wij  = 0.0

    # Centros de contacto virtual
    xiAux, yiAux, xjAux, yjAux = CBMMetrics.rodIntersection(x,y,l,theta, x2,y2,l2,theta2)

    # Distancia entre "esferas" virtuales
    rij = sqrt((xiAux - xjAux)^2 + (yiAux - yjAux)^2)
    σ = 0.5*(d + d2)
    rcut = 5σ

    if rij > σ && rij < rcut
        # Dirección ij
        nijx = (xiAux - xjAux)/rij
        nijy = (yiAux - yjAux)/rij

        # Fuerza Yukawa atractiva (eps>0 → atracción)
        # u(r) = eps * exp(-γ (r-σ)) / r   (desplazo para que r≈σ no sea enorme)
        # F = -du/dr * n  ⇒ magnitud:
        # du/dr = eps * exp(-γ (r-σ)) * (-γ*r - 1) / r^2
        # ⇒ -du/dr = eps * exp(-γ (r-σ)) * (γ*r + 1) / r^2  (repulsivo si eps>0)
        # Hacemos atractiva invirtiendo el signo:
        Fmag = - eps * exp(-gamma*(rij - σ)) * (gamma*rij + 1.0) / (rij^2)

        # (opcional) reescala como en tu versión:
        Fmag /= (eta * (l + d))

        Fijx = Fmag * nijx
        Fijy = Fmag * nijy

        # Mismo torque por palanca
        Wij = ((xiAux - x)*Fijy - (yiAux - y)*Fijx) * 12.0 / ((l + d)^2)
    end

    return Fijx, Fijy, Wij
end


function attractiveForces_rods_yukawa3d(
    x,y,z,d,l,theta,phi,
    x2,y2,z2,d2,l2,theta2,phi2,
    eta, eps, gamma)

    Fijx = 0.0; Fijy = 0.0; Fijz = 0.0
    Wijtheta = 0.0; Wijphi = 0.0

    # 1) Puntos "más cercanos" sobre cada segmento (ajusta el namespace si toca)
    xiAux, yiAux, ziAux,  xjAux, yjAux, zjAux =
        CBMMetrics.rodIntersection3d(x, y, z, l, theta, phi,
                                     x2, y2, z2, l2, theta2, phi2)

    # 2) Vector y distancia j->i
    dx = xiAux - xjAux
    dy = yiAux - yjAux
    dz = ziAux - zjAux
    rij = sqrt(dx*dx + dy*dy + dz*dz)

    σ    = 0.5*(d + d2)
    rcut = 5*σ

    if rij > σ && rij < rcut
        # Dirección normal n̂ (j -> i)
        nijx = dx/rij
        nijy = dy/rij
        nijz = dz/rij

        # 3) Fuerza Yukawa ATRACTIVA:
        #    u(r) = eps * exp(-γ (r-σ)) / r
        #    F = -∂u/∂r n̂  => magnitud atractiva con signo negativo delante
        Fmag = - eps * exp(-gamma*(rij - σ)) * (gamma*rij + 1.0) / (rij^2)

        # Reescalado efectivo (consistente con tus otras fuerzas)
        Fmag /= (eta * (l + d))

        # Fuerza sobre i
        Fijx = Fmag * nijx
        Fijy = Fmag * nijy
        Fijz = Fmag * nijz

        # 4) Par en el CM de i: τ = r × F
        rx = (xiAux - x); ry = (yiAux - y); rz = (ziAux - z)
        τx = ry*Fijz - rz*Fijy
        τy = rz*Fijx - rx*Fijz
        τz = rx*Fijy - ry*Fijx

        # 5) Proyección CORRECTA a (θ, φ):
        #    Wq = τ · (u × ∂u/∂q)
        ux =  cos(theta)*cos(phi)
        uy =  sin(theta)*cos(phi)
        uz =  sin(phi)

        duθ = du_dtheta(theta, phi)   # (-cosφ sinθ,  cosφ cosθ, 0)
        duφ = du_dphi(theta,  phi)    # (-sinφ cosθ, -sinφ sinθ, cosφ)

        # sθ = u × du/dθ
        sθx = uy*duθ[3] - uz*duθ[2]
        sθy = uz*duθ[1] - ux*duθ[3]
        sθz = ux*duθ[2] - uy*duθ[1]

        # sφ = u × du/dφ
        sφx = uy*duφ[3] - uz*duφ[2]
        sφy = uz*duφ[1] - ux*duφ[3]
        sφz = ux*duφ[2] - uy*duφ[1]

        scale = 12.0 / ((l + d)^2)
        Wijtheta = scale * (τx*sθx + τy*sθy + τz*sθz)
        Wijphi   = scale * (τx*sφx + τy*sφy + τz*sφz)
    end

    return Fijx, Fijy, Fijz, Wijtheta, Wijphi
end


# --- Atracción bacteria–sustrato (Yukawa) para un rod en 3D ---
#   U_cs(h) = eps_cs * exp(gamma*(h - h0)) / (h - h0)
#   Fz = -dU/dh  (dirección -z; aquí el signo sale de la derivada)
#
# Parámetros:
#   eta      : factor de reescalado viscoso (como en tu código)
#   eps_cs   : profundidad (debe ser NEGATIVA si quieres atracción estricta como en el paper)
#   gamma    : parámetro de decaimiento
#   h0       : posición efectiva del "sustrato" en la métrica de h (paper)
#   rcut_h   : (opcional) corte de alcance en h para eficiencia
#
function substrateAttraction_rods_yukawa3d(
    x,y,z,d,l,theta,phi,
    eta, eps_cs, gamma; rcut_h = 5.0
)
    a = 0 # plano en z=0
    (x1,y1,z1), (x2,y2,z2) = rod_poles3d_ordered(x,y,z,l,theta,phi)

    # Acumuladores
    Fx = 0.0; Fy = 0.0; Fz = 0.0
    τx = 0.0; τy = 0.0; τz = 0.0

    # Recorremos polos (aproximación por extremos del esferocilindro)
    for pole in 1:2
        xp, yp, zp = pole == 1 ? (x1,y1,z1) : (x2,y2,z2)

        # gap h = distancia "libre" polo–plano menos el radio a
        #   s = zp  (altura del polo)
        #   h = s - a
        s  = zp
        h  = s - a + d/2

        # Solo aporta si está por encima del plano efectivo (h > h0)
        # y dentro de un corte razonable
        if h > 0
            Δh = h 
            if Δh < rcut_h
                # f(h) = exp(gamma*(h-h0)) / (h-h0)
                # dU/dh = eps_cs * f * (gamma - 1/(h-h0))
                # F = -dU/dh  (hacia -z si dU/dh > 0)
                f     = exp(gamma * Δh) / (Δh)
                dUdh  = eps_cs * f * (gamma - 1.0/Δh)
                Fmag  = - dUdh / (eta * (l + d))   # reescalado consistente

                # Fuerza solo en z (hacia el sustrato si Fmag < 0)
                Fx_p, Fy_p, Fz_p = 0.0, 0.0, Fmag
                Fx += Fx_p; Fy += Fy_p; Fz += Fz_p

                # Par respecto al CM: r × F
                rx = xp - x; ry = yp - y; rz = zp - z
                τx += ry*Fz_p - rz*Fy_p
                τy += rz*Fx_p - rx*Fz_p
                τz += rx*Fy_p - ry*Fx_p
            end
        end
    end

    # --- Proyección del par a (θ, φ) ---
    ux =  cos(theta)*cos(phi)
    uy =  sin(theta)*cos(phi)
    uz =  sin(phi)

    duθ = du_dtheta(theta, phi)   # (-cosφ sinθ,  cosφ cosθ, 0)
    duφ = du_dphi(theta,  phi)    # (-sinφ cosθ, -sinφ sinθ, cosφ)

    # sθ = u × du/dθ
    sθx = uy*duθ[3] - uz*duθ[2]
    sθy = uz*duθ[1] - ux*duθ[3]
    sθz = ux*duθ[2] - uy*duθ[1]

    # sφ = u × du/dφ
    sφx = uy*duφ[3] - uz*duφ[2]
    sφy = uz*duφ[1] - ux*duφ[3]
    sφz = ux*duφ[2] - uy*duφ[1]

    # Mismo factor de escala que usas en el resto del código
    scale = 12.0 / ((l + d)^2)
    Wθ = scale * (τx*sθx + τy*sθy + τz*sθz)
    Wφ = scale * (τx*sφx + τy*sφy + τz*sφz)

    return Fx, Fy, Fz, Wθ, Wφ
end

function substrateRepulsion_rods3d(
    x,y,z,d,l,theta,phi,
    eta, Ebv
)
    a = 0
    (x1,y1,z1), (x2,y2,z2) = rod_poles3d_ordered(x,y,z,l,theta,phi)

    # Acumuladores netos
    Fx = 0.0; Fy = 0.0; Fz = 0.0
    τx = 0.0; τy = 0.0; τz = 0.0

    for pole in 1:2
        xp, yp, zp = pole == 1 ? (x1,y1,z1) : (x2,y2,z2)

        s = zp
        δ = a + d/2 - s

        # --- Repulsión suelo (esfera–plano) ---
        if δ > 0.0
            Fh = Ebv * sqrt(d * δ^3) / (eta * (l + d))
            # Fuerza en el polo: (0,0,Fh)
            Fx += 0.0; Fy += 0.0; Fz += Fh
            # Par respecto al CM: r × F
            rx = xp - x; ry = yp - y; rz = zp - z
            τx += ry*Fh - rz*0.0
            τy += rz*0.0 - rx*Fh
            τz += rx*0.0 - ry*0.0
        end

    end

    # --- Proyección correcta del par a (θ, φ) ---
    ux =  cos(theta)*cos(phi)
    uy =  sin(theta)*cos(phi)
    uz =  sin(phi)

    duθ = du_dtheta(theta, phi)   # (-cosφ sinθ,  cosφ cosθ, 0)
    duφ = du_dphi(theta,  phi)    # (-sinφ cosθ, -sinφ sinθ, cosφ)

    # sθ = u × du/dθ
    sθx = uy*duθ[3] - uz*duθ[2]
    sθy = uz*duθ[1] - ux*duθ[3]
    sθz = ux*duθ[2] - uy*duθ[1]

    # sφ = u × du/dφ
    sφx = uy*duφ[3] - uz*duφ[2]
    sφy = uz*duφ[1] - ux*duφ[3]
    sφz = ux*duφ[2] - uy*duφ[1]

    scale = 12.0 / ((l + d)^2)
    Wθ = scale * (τx*sθx + τy*sθy + τz*sθz)
    Wφ = scale * (τx*sφx + τy*sφy + τz*sφz)

    return Fx, Fy, Fz, Wθ, Wφ
end


function padRepulsion_rods3d(
    x,y,z,d,l,theta,phi,
    eta, Ebp, hPad
)
    a = hPad
    (x1,y1,z1), (x2,y2,z2) = rod_poles3d_ordered(x,y,z,l,theta,phi)

    # Acumuladores netos
    Fx = 0.0; Fy = 0.0; Fz = 0.0
    τx = 0.0; τy = 0.0; τz = 0.0

    for pole in 1:2
        xp, yp, zp = pole == 1 ? (x1,y1,z1) : (x2,y2,z2)

        s = zp
        δ = s + d/2 -a

        # --- Repulsión pad (esfera–plano) ---
        if δ > 0.0
            Fh = -Ebp * sqrt(d * δ^3) / (eta * (l + d)) # negativo: hacia -z
            # Fuerza en el polo: (0,0,Fh)
            Fx += 0.0; Fy += 0.0; Fz +=  Fh
            # Par respecto al CM: r × F
            rx = xp - x; ry = yp - y; rz = zp - z
            τx += ry*Fh - rz*0.0
            τy += rz*0.0 - rx*Fh
            τz += rx*0.0 - ry*0.0
        end

    end

    # --- Proyección correcta del par a (θ, φ) ---
    ux =  cos(theta)*cos(phi)
    uy =  sin(theta)*cos(phi)
    uz =  sin(phi)

    duθ = du_dtheta(theta, phi)   # (-cosφ sinθ,  cosφ cosθ, 0)
    duφ = du_dphi(theta,  phi)    # (-sinφ cosθ, -sinφ sinθ, cosφ)

    # sθ = u × du/dθ
    sθx = uy*duθ[3] - uz*duθ[2]
    sθy = uz*duθ[1] - ux*duθ[3]
    sθz = ux*duθ[2] - uy*duθ[1]

    # sφ = u × du/dφ
    sφx = uy*duφ[3] - uz*duφ[2]
    sφy = uz*duφ[1] - ux*duφ[3]
    sφz = ux*duφ[2] - uy*duφ[1]

    scale = 12.0 / ((l + d)^2)
    Wθ = scale * (τx*sθx + τy*sθy + τz*sθz)
    Wφ = scale * (τx*sφx + τy*sφy + τz*sφz)

    return Fx, Fy, Fz, Wθ, Wφ
end

# --- Fuerzas de pared y pili en 3D para un rod ---

function wallForces_rods_bond3d(
    x,y,z,d,l,theta,phi,
    eta, Ebv,
    kb, Δ0, Δmax, T_eng,hax,
    xpili1::Float64, ypili1::Float64, t_pili1::Int,
    xpili2::Float64, ypili2::Float64, t_pili2::Int
)
    a = 0.5*d
    (x1,y1,z1), (x2,y2,z2) = rod_poles3d_ordered(x,y,z,l,theta,phi)

    # Acumuladores netos
    Fx = 0.0; Fy = 0.0; Fz = 0.0
    τx = 0.0; τy = 0.0; τz = 0.0

    for pole in 1:2
        xp, yp, zp = pole == 1 ? (x1,y1,z1) : (x2,y2,z2)
        xpili = pole == 1 ? xpili1 : xpili2
        ypili = pole == 1 ? ypili1 : ypili2
        t_pili = pole == 1 ? t_pili1 : t_pili2

        s = zp
        δ = a - s

        # --- Repulsión suelo (esfera–plano) ---
        if δ > 0.0
            Fh = Ebv * sqrt(d * δ^3) / (eta * (l + d))
            # Fuerza en el polo: (0,0,Fh)
            Fx += 0.0; Fy += 0.0; Fz += Fh
            # Par respecto al CM: r × F
            rx = xp - x; ry = yp - y; rz = zp - z
            τx += ry*Fh - rz*0.0
            τy += rz*0.0 - rx*Fh
            τz += rx*0.0 - ry*0.0
        end

        # --- PILI ---
        if t_pili < 0 # Enganchado
            xa, ya, za = xpili, ypili, 0.0
            dx = xp - xa; dy = yp - ya; dz = zp - za
            Δ = sqrt(dx*dx + dy*dy + dz*dz)
            if Δ > 1e-12
                Fspr = -kb * (Δ - Δ0)
                if (Δ > Δmax) 
                    t_pili = 0; xpili = NaN; ypili = NaN
                else
                    Fx_p = Fspr * (dx/Δ) / (eta * (l + d))
                    Fy_p = Fspr * (dy/Δ) / (eta * (l + d))
                    Fz_p = Fspr * (dz/Δ) / (eta * (l + d))
                    Fx += Fx_p; Fy += Fy_p; Fz += Fz_p
                    rx = xp - x; ry = yp - y; rz = zp - z
                    τx += ry*Fz_p - rz*Fy_p
                    τy += rz*Fx_p - rx*Fz_p
                    τz += rx*Fy_p - ry*Fx_p
                end
            end

        else  # No enganchado, contando tiempo
            if (s - a) < hax
                if t_pili >= T_eng
                    t_pili = -1; xpili = xp; ypili = yp
                else
                    t_pili += 1
                end
            else
                t_pili = 0
            end
        end

        if pole==1
            xpili1=xpili; ypili1=ypili; t_pili1=t_pili
        else
            xpili2=xpili; ypili2=ypili; t_pili2=t_pili
        end
    end

    # --- Proyección correcta del par a (θ, φ) ---
    ux =  cos(theta)*cos(phi)
    uy =  sin(theta)*cos(phi)
    uz =  sin(phi)

    duθ = du_dtheta(theta, phi)   # (-cosφ sinθ,  cosφ cosθ, 0)
    duφ = du_dphi(theta,  phi)    # (-sinφ cosθ, -sinφ sinθ, cosφ)

    # sθ = u × du/dθ
    sθx = uy*duθ[3] - uz*duθ[2]
    sθy = uz*duθ[1] - ux*duθ[3]
    sθz = ux*duθ[2] - uy*duθ[1]

    # sφ = u × du/dφ
    sφx = uy*duφ[3] - uz*duφ[2]
    sφy = uz*duφ[1] - ux*duφ[3]
    sφz = ux*duφ[2] - uy*duφ[1]

    scale = 12.0 / ((l + d)^2)
    Wθ = scale * (τx*sθx + τy*sθy + τz*sθz)
    Wφ = scale * (τx*sφx + τy*sφy + τz*sφz)

    return Fx, Fy, Fz, Wθ, Wφ,
           xpili1, ypili1, t_pili1,
           xpili2, ypili2, t_pili2
end
