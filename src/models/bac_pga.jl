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

# Helpers mínimos
@inline dot3(a,b) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]

# d u / dtheta y d u / dphi para u = (cosφ cosθ, cosφ sinθ, sinφ)
@inline function du_dtheta(theta, phi)
    (-cos(phi)*sin(theta),  cos(phi)*cos(theta),  0.0)
end
@inline function du_dphi(theta, phi)
    (-sin(phi)*cos(theta), -sin(phi)*sin(theta),  cos(phi))
end

function repulsiveForces_rods3d(
    x,y,z,d,l,theta,phi,
    x2,y2,z2,d2,l2,theta2,phi2, eta, E)

    Fijx = 0.0; Fijy = 0.0; Fijz = 0.0
    Wijtheta = 0.0; Wijphi = 0.0

    # 1) Puntos de "contacto" (más cercanos) sobre cada segmento
    xiAux, yiAux, ziAux,  xjAux, yjAux, zjAux =
        rodIntersection3d(x, y, z, l, theta, phi,
                          x2, y2, z2, l2, theta2, phi2; separation=0.99)

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
        duθ = du_dtheta(theta, phi)
        duφ = du_dphi(theta, phi)

        scale = 12.0 / ((l + d)^2)  # mismo factor que usabas en 2D
        Wijtheta = scale * (τx*duθ[1] + τy*duθ[2] + τz*duθ[3])
        Wijphi   = scale * (τx*duφ[1] + τy*duφ[2] + τz*duφ[3])
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

# --- derivadas del eje u respecto a los ángulos ---
# u = (cosφ cosθ, cosφ sinθ, sinφ)
@inline function du_dtheta(theta, phi)
    (-cos(phi)*sin(theta),  cos(phi)*cos(theta),  0.0)
end
@inline function du_dphi(theta, phi)
    (-sin(phi)*cos(theta), -sin(phi)*sin(theta),  cos(phi))
end

@inline dot3(a,b) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]

function attractiveForces_rods_yukawa3d(
    x,y,z,d,l,theta,phi,
    x2,y2,z2,d2,l2,theta2,phi2,
    eta, eps, gamma)

    Fijx = 0.0; Fijy = 0.0; Fijz = 0.0
    Wijtheta = 0.0; Wijphi = 0.0

    # 1) Puntos de "contacto" (más cercanos) sobre cada segmento
    #    Sustituye por tu espacio de nombres si procede: CBMMetrics.rodIntersection3d(...)
    xiAux, yiAux, ziAux,  xjAux, yjAux, zjAux =
        rodIntersection3d(x, y, z, l, theta, phi,
                          x2, y2, z2, l2, theta2, phi2; separation=0.99)

    # 2) Distancia entre "esferas" virtuales y dirección j->i
    dx = xiAux - xjAux
    dy = yiAux - yjAux
    dz = ziAux - zjAux
    rij = sqrt(dx*dx + dy*dy + dz*dz)

    σ = 0.5*(d + d2)
    rcut = 5σ

    if rij > σ && rij < rcut
        # Dirección normal
        nijx = dx/rij
        nijy = dy/rij
        nijz = dz/rij

        # 3) Fuerza Yukawa ATRACTIVA (eps>0 -> atracción)
        #    u(r) = eps * exp(-γ (r-σ)) / r
        #    F = -∂u/∂r * n   => magnitud:
        #    ∂u/∂r = eps * exp(-γ(r-σ)) * (-γ*r - 1) / r^2
        #    -∂u/∂r = eps * exp(-γ(r-σ)) * (γ*r + 1) / r^2  (esto sería REPULSIVA)
        #    Para que sea atractiva, invertimos el signo:
        Fmag = - eps * exp(-gamma*(rij - σ)) * (gamma*rij + 1.0) / (rij^2)

        # Reescalado como en tu versión
        Fmag /= (eta * (l + d))

        # Fuerza sobre i
        Fijx = Fmag * nijx
        Fijy = Fmag * nijy
        Fijz = Fmag * nijz

        # 4) Torque en el centro de i: τ = r_i × F_ij
        rx = (xiAux - x); ry = (yiAux - y); rz = (ziAux - z)
        τx = ry*Fijz - rz*Fijy
        τy = rz*Fijx - rx*Fijz
        τz = rx*Fijy - ry*Fijx

        # 5) Proyección del torque en las direcciones generalizadas (θ, φ)
        duθ = du_dtheta(theta, phi)
        duφ = du_dphi(theta,  phi)

        scale = 12.0 / ((l + d)^2)
        Wijtheta = scale * (τx*duθ[1] + τy*duθ[2] + τz*duθ[3])
        Wijphi   = scale * (τx*duφ[1] + τy*duφ[2] + τz*duφ[3])
    end

    return Fijx, Fijy, Fijz, Wijtheta, Wijphi
end
