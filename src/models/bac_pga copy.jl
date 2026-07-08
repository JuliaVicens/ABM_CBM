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

# """
#     repulsiveForces(
#         x, y, d, l, theta, etab, type,
#         x2, y2, d2, l2, theta2, etag, type2,
#         Ebb, Ebg, Egg
#     ) -> (Fijx, Fijy, Wij)

# Computes the repulsive interaction forces between two agents (either rods or gels), based on their type and geometry. If both are of the same type and are gels, direct contact is assumed; otherwise, virtual contact points are calculated.

# ### Parameters
# - `x, y`: coordinates of the first agent.
# - `d, l, theta`: diameter, length, and orientation angle of the first agent.
# - `etab, type`: effective viscosity and type of the first agent (`0` for rod, `1` for gel).
# - `x2, y2`: coordinates of the second agent.
# - `d2, l2, theta2`: diameter, length, and orientation angle of the second agent.
# - `etag, type2`: effective viscosity and type of the second agent.
# - `Ebb, Ebg, Egg`: elastic moduli for rod-rod, rod-gel, and gel-gel interactions.

# ### Returns
# - `Fijx, Fijy`: components of the repulsive force between the agents.
# - `Wij`: radial torque-like contribution (only non-zero when the first agent is a rod).
# """
# function repulsiveForces(
#     x,y,d,l,theta,etab,type,
#     x2,y2,d2,l2,theta2,etag,type2, Ebb, Ebg, Egg )

#     Fijx = 0.
#     Fijy = 0.
#     Wij = 0.

#     if type==type2 && type==1
#         xiAux,yiAux,xjAux,yjAux=x,y,x2,y2
#     else
#         #Function that finds the virtual spheres of contact between both rods
#         xiAux,yiAux,xjAux,yjAux = CBMMetrics.rodIntersection(x,y,l,theta,x2,y2,l2,theta2)
#     end

#     #Compute distance between virtual spheres
#     rij = sqrt((xiAux-xjAux)^2 +(yiAux-yjAux)^2)
#     if rij > 0. && rij < (d+d2)/2. #If it is smaller than a diameter compute forces
#         #Compute auxiliar
#         hAux = (d+d2)/2. - rij
#         #Compute direction
#         nijx = (xiAux-xjAux)/rij
#         nijy = (yiAux-yjAux)/rij
#         if type == 0
#             if type2 == 0
#                 FnAux = Ebb * sqrt(d2* hAux^3)  / (etab * (l + d))
#             else
#                 FnAux = Ebg * sqrt(d2* hAux^3) / (etab * (l + d))
#             end
#             Fijx = FnAux * nijx
#             Fijy = FnAux * nijy
#             Wij = ((xiAux-x)*Fijy - (yiAux-y)*Fijx) * 12. / (etab * (l + d)^3)
#         else
#             if type2 == 0
#                 FnAux = Ebg * sqrt(d2* hAux^3) / (etag * d)
#             else
#                 FnAux = Egg * sqrt(d2* hAux^3) / (etag * d)
#             end
#             Fijx = FnAux * nijx
#             Fijy = FnAux * nijy
#         end
#     end

#     return Fijx, Fijy, Wij

# end





# function repulsiveForces_rods(
#     x,y,d,l,theta,
#     x2,y2,d2,l2,theta2,eta, E)

#     Fijx = 0.
#     Fijy = 0.
#     Wij = 0.


#     #Function that finds the virtual spheres of contact between both rods
#     xiAux,yiAux,xjAux,yjAux = CBMMetrics.rodIntersection(x,y,l,theta,x2,y2,l2,theta2)
   

#     #Compute distance between virtual spheres
#     rij = sqrt((xiAux-xjAux)^2 +(yiAux-yjAux)^2)
#     if rij > 0. && rij < (d+d2)/2. #If it is smaller than a diameter compute forces
#         #Compute auxiliar
#         hAux = (d+d2)/2. - rij
#         #Compute direction
#         nijx = (xiAux-xjAux)/rij
#         nijy = (yiAux-yjAux)/rij
    
#         FnAux = E * sqrt(d2* hAux^3)  / (eta * (l + d))
           
#         Fijx = FnAux * nijx
#         Fijy = FnAux * nijy
#         Wij = ((xiAux-x)*Fijy - (yiAux-y)*Fijx) * 12. / ((l + d)^2)

#     # elseif (sqrt((xiAux+l*cos(theta)-(xjAux+l2*cos(theta2)))^2 +(yiAux+l*sin(theta)-(yjAux+l2*sin(theta2)))^2))< (d+d2)/2.
#     #     println("Error in distance calculation between rods")
#     # elseif (sqrt((xiAux-l*cos(theta)-(xjAux-l2*cos(theta2)))^2 +(yiAux-l*sin(theta)-(yjAux-l2*sin(theta2)))^2))< (d+d2)/2.
#     #     println("Error in distance calculation between rods")
#     # elseif (sqrt((xiAux-l*cos(theta)-(xjAux+l2*cos(theta2)))^2 +(yiAux-l*sin(theta)-(yjAux+l2*sin(theta2)))^2))< (d+d2)/2.
#     #     println("Error in distance calculation between rods")
#     # elseif (sqrt((xiAux+l*cos(theta)-(xjAux-l2*cos(theta2)))^2 +(yiAux+l*sin(theta)-(yjAux-l2*sin(theta2)))^2))< (d+d2)/2.
#     #     println("Error in distance calculation between rods")
      
#     end

#     return Fijx, Fijy, Wij

# end


# function repulsiveForces_rods_(
#     x,y,d,l,theta,
#     x2,y2,d2,l2,theta2,eta, E)

#     Fijx = 0.0
#     Fijy = 0.0
#     Wij  = 0.0

#     # --- helpers y pre-cálculos ---
 
#     Rsum = (d + d2)/2.0

#     # medias longitudes (las puntas reales están a ± l/2)
#     hx  = (l/2.0)*cos(theta)
#     hy  = (l/2.0)*sin(theta)
#     hx2 = (l2/2.0)*cos(theta2)
#     hy2 = (l2/2.0)*sin(theta2)

#     # puntas de cada rod
#     # rod 1
#     x_plus  = x + hx;  y_plus  = y + hy
#     x_minus = x - hx;  y_minus = y - hy
#     # rod 2
#     x2_plus  = x2 + hx2;  y2_plus  = y2 + hy2
#     x2_minus = x2 - hx2;  y2_minus = y2 - hy2



#     # =========================================
#     # 1) Caso general: puntos "virtuales" vía rodIntersection
#     # =========================================
#     xiAux, yiAux, xjAux, yjAux = CBMMetrics.rodIntersection_(x,y,l,theta, x2,y2,l2,theta2)

#     rij = sqrt((xiAux - xjAux)^2 + (yiAux - yjAux)^2)
#     if rij > 0.0 && rij < Rsum
#         hAux = Rsum - rij
#         nijx = (xiAux - xjAux)/rij
#         nijy = (yiAux - yjAux)/rij

#         FnAux = E * sqrt(d * hAux^3) / (eta*(l+d))
#         Fijx  = FnAux * nijx
#         Fijy  = FnAux * nijy
#         Wij   = ((xiAux - x)*Fijy - (yiAux - y)*Fijx) * 12.0 / (l^2)

#     # =========================================
#     # 2) Puntas: + con + (punta delantera de ambos)
#     # =========================================
#     elseif sqrt((x_plus - x2_plus)^2 + (y_plus - y2_plus)^2) < Rsum
#         xiAux, yiAux = x_plus,  y_plus
#         xjAux, yjAux = x2_plus, y2_plus
#         rij = sqrt((xiAux - xjAux)^2 + (yiAux - yjAux)^2)
#         hAux = Rsum - rij
#         nijx = (xiAux - xjAux)/rij
#         nijy = (yiAux - yjAux)/rij

#         FnAux = E * sqrt(d * hAux^3) / (eta*(l+d))
#         Fijx  = FnAux * nijx
#         Fijy  = FnAux * nijy
#         Wij   = ((xiAux - x)*Fijy - (yiAux - y)*Fijx) * 12.0 / (l^2)
#     # =========================================
#     # 3) Puntas: - con - (colas de ambos)
#     # =========================================
#     elseif sqrt((x_minus - x2_minus)^2 + (y_minus - y2_minus)^2) < Rsum
#         xiAux, yiAux = x_minus,  y_minus
#         xjAux, yjAux = x2_minus, y2_minus
#         rij = sqrt((xiAux - xjAux)^2 + (yiAux - yjAux)^2)
#         hAux = Rsum - rij
#         nijx = (xiAux - xjAux)/rij
#         nijy = (yiAux - yjAux)/rij

#         FnAux = E * sqrt(d * hAux^3) / (eta*(l+d))
#         Fijx  = FnAux * nijx
#         Fijy  = FnAux * nijy
#         Wij   = ((xiAux - x)*Fijy - (yiAux - y)*Fijx) * 12.0 / (l^2)
#     # =========================================
#     # 4) Puntas cruzadas: - de 1 con + de 2
#     # =========================================
#     elseif sqrt((x_minus - x2_plus)^2 + (y_minus - y2_plus)^2) < Rsum
#         xiAux, yiAux = x_minus,  y_minus
#         xjAux, yjAux = x2_plus,  y2_plus
#         rij = sqrt((xiAux - xjAux)^2 + (yiAux - yjAux)^2)
#         hAux = Rsum - rij
#         nijx = (xiAux - xjAux)/rij
#         nijy = (yiAux - yjAux)/rij

#         FnAux = E * sqrt(d * hAux^3) / (eta*(l+d))
#         Fijx  = FnAux * nijx
#         Fijy  = FnAux * nijy
#         Wij   = ((xiAux - x)*Fijy - (yiAux - y)*Fijx) * 12.0 / (l^2)
#     # =========================================
#     # 5) Puntas cruzadas: + de 1 con - de 2
#     # =========================================
#     elseif sqrt((x_plus - x2_minus)^2 + (y_plus - y2_minus)^2) < Rsum
#         xiAux, yiAux = x_plus,   y_plus
#         xjAux, yjAux = x2_minus, y2_minus
#         rij = sqrt((xiAux - xjAux)^2 + (yiAux - yjAux)^2)
#         hAux = Rsum - rij
#         nijx = (xiAux - xjAux)/rij
#         nijy = (yiAux - yjAux)/rij

#         FnAux = E * sqrt(d * hAux^3) / (eta*(l+d))
#         Fijx  = FnAux * nijx
#         Fijy  = FnAux * nijy
#         Wij   = ((xiAux - x)*Fijy - (yiAux - y)*Fijx) * 12.0 / (l^2)
#     end

#     return Fijx, Fijy, Wij
# end


function repulsiveForces_rods3d(
    x,y,z,d,l,theta,phi,
    x2,y2,z2,d2,l2,theta2,phi2, eta, E)

    Fijx = 0.0; Fijy = 0.0; Fijz = 0.0
    Wijtheta = 0.0; Wijphi = 0.0

    # 1) Puntos de "contacto" (más cercanos) sobre cada segmento
    xiAux, yiAux, ziAux,  xjAux, yjAux, zjAux =
        CBMMetrics.rodIntersection3d_(x, y, z, l, theta, phi,
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
    xiAux,yiAux,xjAux,yjAux = CBMMetrics.rodIntersection_(x,y,l,theta,x2,y2,l2,theta2)
   
    c=cos(theta)
    s=sin(theta)
    #Compute distance between virtual spheres
    rij = sqrt((xiAux-xjAux)^2 +(yiAux-yjAux)^2)

    mu_par  = 1/(eta)   # μ∥
    mu_perp = 1/(eta*A)       # μ⊥
    
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
        Wij = ((xiAux-x)*Fijy - (yiAux-y)*Fijx) * 12.0 *mu_perp / ((l + d)^2)
      
    end

    return Fasx, Fasy, Wij

end

# function attractiveForces_rods(
#     x,y,d,l,theta,
#     x2,y2,d2,l2,theta2,eta, eps)

#     Fijx = 0.
#     Fijy = 0.
#     Wij = 0.

#     #Function that finds the virtual spheres of contact between both rods
#     xiAux,yiAux,xjAux,yjAux = CBMMetrics.rodIntersection(x,y,l,theta,x2,y2,l2,theta2)
   

#     #Compute distance between virtual spheres
#     rij = sqrt((xiAux-xjAux)^2 +(yiAux-yjAux)^2)
#     if rij > (d+d2)/2. && rij < 5*(d+d2)/2. #If it is smaller than a diameter compute forces
       
#         #Compute direction
#         nijx = (xiAux-xjAux)/rij
#         nijy = (yiAux-yjAux)/rij
    
#         FnAux = -24*eps* ((d+d2)/2)^6/rij^7 / (eta * (l + d))
           
#         Fijx = FnAux * nijx
#         Fijy = FnAux * nijy
#         Wij = ((xiAux-x)*Fijy - (yiAux-y)*Fijx) * 12. / ((l + d)^2)
      
#     end

#     return Fijx, Fijy, Wij

# end

# function attractiveForces_rods_yukawa(
#     x,y,d,l,theta,
#     x2,y2,d2,l2,theta2,
#     eta, eps, gamma)

#     Fijx = 0.0
#     Fijy = 0.0
#     Wij  = 0.0

#     # Centros de contacto virtual
#     xiAux, yiAux, xjAux, yjAux = CBMMetrics.rodIntersection(x,y,l,theta, x2,y2,l2,theta2)

#     # Distancia entre "esferas" virtuales
#     rij = sqrt((xiAux - xjAux)^2 + (yiAux - yjAux)^2)
#     σ = 0.5*(d + d2)
#     rcut = 5σ

#     if rij > σ && rij < rcut
#         # Dirección ij
#         nijx = (xiAux - xjAux)/rij
#         nijy = (yiAux - yjAux)/rij

#         # Fuerza Yukawa atractiva (eps>0 → atracción)
#         # u(r) = eps * exp(-γ (r-σ)) / r   (desplazo para que r≈σ no sea enorme)
#         # F = -du/dr * n  ⇒ magnitud:
#         # du/dr = eps * exp(-γ (r-σ)) * (-γ*r - 1) / r^2
#         # ⇒ -du/dr = eps * exp(-γ (r-σ)) * (γ*r + 1) / r^2  (repulsivo si eps>0)
#         # Hacemos atractiva invirtiendo el signo:
#         Fmag = - eps * exp(-gamma*(rij - σ)) * (gamma*rij + 1.0) / (rij^2)

#         # (opcional) reescala como en tu versión:
#         Fmag /= (eta * (l + d))

#         Fijx = Fmag * nijx
#         Fijy = Fmag * nijy

#         # Mismo torque por palanca
#         Wij = ((xiAux - x)*Fijy - (yiAux - y)*Fijx) * 12.0 / ((l + d)^2)
#     end

#     return Fijx, Fijy, Wij
# end


function attractiveForces_rods_yukawa3d(
    x,y,z,d,l,theta,phi,
    x2,y2,z2,d2,l2,theta2,phi2,
    eta, eps, gamma)

    Fijx = 0.0; Fijy = 0.0; Fijz = 0.0
    Wijtheta = 0.0; Wijphi = 0.0

    # 1) Puntos "más cercanos" sobre cada segmento (ajusta el namespace si toca)
    xiAux, yiAux, ziAux,  xjAux, yjAux, zjAux =
        CBMMetrics.rodIntersection3d_(x, y, z, l, theta, phi,
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


# function alignmentTorque_rods3d(
#     x,y,z,d,l,theta,phi,
#     x2,y2,z2,d2,l2,theta2,phi2,
#     k_align, r_align, p; eps=0.7)

#     # salida (sin fuerzas traslacionales)
#     Wijtheta = 0.0; Wijphi = 0.0

#     # puntos más cercanos
#     xiAux, yiAux, ziAux,  xjAux, yjAux, zjAux =
#         CBMMetrics.rodIntersection3d_(x, y, z, l, theta, phi,
#                                      x2, y2, z2, l2, theta2, phi2)

#     # sanity
#     if !(isfinite(xiAux) & isfinite(yiAux) & isfinite(ziAux) &
#          isfinite(xjAux) & isfinite(yjAux) & isfinite(zjAux))
#         return 0.0, 0.0, 0.0, 0.0, 0.0
#     end

#     dx = xiAux - xjAux; dy = yiAux - yjAux; dz = ziAux - zjAux
#     rij2 = dx*dx + dy*dy + dz*dz
#     if !(rij2 ≥ 0.0)            # cubre NaN
#         return 0.0,0.0,0.0,0.0,0.0
#     end
#     rij = sqrt(rij2)

#     # vecindad
#     if !(isfinite(rij)) || rij > r_align
#         return 0.0,0.0,0.0,0.0,0.0
#     end

#     # direcciones (φ = elevación)
#     ux =  cos(theta)*cos(phi)
#     uy =  sin(theta)*cos(phi)
#     uz =  sin(phi)

#     ux2 =  cos(theta2)*cos(phi2)
#     uy2 =  sin(theta2)*cos(phi2)
#     uz2 =  sin(phi2)

#     # productos
#     c   = ux*ux2 + uy*uy2 + uz*uz2
#     c   = clamp(c, -1.0, 1.0)     # por si hay drift numérico

#     cx = uy*uz2 - uz*uy2
#     cy = uz*ux2 - ux*uz2
#     cz = ux*uy2 - uy*ux2
#     normc2 = cx*cx + cy*cy + cz*cz
#     if !(normc2 > eps)            # paralelas/antiparalelas → no torque
#         return 0.0,0.0,0.0,0.0,0.0
#     end
#     inv_normc = inv(sqrt(normc2))

#     # peso radial suave (evita que actúe lejos)
#     wr = exp(-(rij/r_align)^2)

#     # magnitud (nemática: abs(c)^p)
#     mag = k_align * wr * (abs(c)^p)

#     # par vectorial
#     τx = mag * cx * inv_normc
#     τy = mag * cy * inv_normc
#     τz = mag * cz * inv_normc

#     # proyección a (θ, φ)
#     duθx, duθy, duθz = du_dtheta(theta, phi)
#     duφx, duφy, duφz = du_dphi(theta,  phi)

#     sθx = uy*duθz - uz*duθy
#     sθy = uz*duθx - ux*duθz
#     sθz = ux*duθy - uy*duθx

#     sφx = uy*duφz - uz*duφy
#     sφy = uz*duφx - ux*duφz
#     sφz = ux*duφy - uy*duφx

#     scale = 12.0 / ((l + d)^2)

#     Wijtheta = scale * (τx*sθx + τy*sθy + τz*sθz)
#     Wijphi   = scale * (τx*sφx + τy*sφy + τz*sφz)

#     # limpia salidas
#     if !(isfinite(Wijtheta) & isfinite(Wijphi))
#         Wijtheta = 0.0; Wijphi = 0.0
#     end
#     return 0.0, 0.0, 0.0, Wijtheta, Wijphi
# end



# function repulsiveForces_rods3dline_asym(
#     x,y,z,d,l,theta,phi,                     # rod i
#     x2,y2,z2,d2,l2,theta2,phi2,              # rod j
#     eta, E, A, K)

#     # Acumuladores de fuerza y par (sobre i por j)
#     Fijx = 0.0; Fijy = 0.0; Fijz = 0.0
#     τx_sum = 0.0; τy_sum = 0.0; τz_sum = 0.0

#     # Dirección del rod i (φ = elevación)
#     ux =  cos(theta)*cos(phi)
#     uy =  sin(theta)*cos(phi)
#     uz =  sin(phi)

#     # Derivadas para proyección (asumiendo tu convención)
#     duθx, duθy, duθz = du_dtheta(theta, phi)
#     duφx, duφy, duφz = du_dphi(theta,  phi)

#     # sθ = u × du/dθ ; sφ = u × du/dφ
#     sθx = uy*duθz - uz*duθy
#     sθy = uz*duθx - ux*duθz
#     sθz = ux*duθy - uy*duθx

#     sφx = uy*duφz - uz*duφy
#     sφy = uz*duφx - ux*duφz
#     sφz = ux*duφy - uy*duφx

#     # Movilidad anisótropa (como en tu versión actual)
#     mu_par  = 1.0 / eta           # μ∥
#     mu_perp = 1.0 / (eta * A)     # μ⊥
#     delta_s = l / (K - 1)        # separación entre puntos muestreados

#     for s in LinRange(-l/2, l/2, K)
#         # Punto sobre el eje del rod i
#         xpi = x + s*ux
#         ypi = y + s*uy
#         zpi = z + s*uz

#         # Punto más cercano entre (punto i) y (rod j completo)
#         xiAux, yiAux, ziAux,  xjAux, yjAux, zjAux =
#             CBMMetrics.rodIntersection3d_(xpi, ypi, zpi, 0.0, theta, phi,
#                                          x2,  y2,  z2,  l2, theta2, phi2)

#         dx = xiAux - xjAux
#         dy = yiAux - yjAux
#         dz = ziAux - zjAux
#         rij = sqrt(dx*dx + dy*dy + dz*dz)

#         if rij > 0.0 && rij < (d + d2)/2
#             hAux = (d + d2)/2 - rij
#             nijx = dx / rij
#             nijy = dy / rij
#             nijz = dz / rij

#             # Fuerza cruda de contacto (SIN /eta)
#             Fraw_mag = delta_s*E * sqrt(d2 * hAux^3) / (l + d)
#             Fx_raw = Fraw_mag * nijx
#             Fy_raw = Fraw_mag * nijy
#             Fz_raw = Fraw_mag * nijz

#             # TRASLACIÓN: movilidad anisótropa (tu estilo actual)
#             u_dot_Fraw = ux*Fx_raw + uy*Fy_raw + uz*Fz_raw
#             Fasx = mu_perp*Fx_raw + (mu_par - mu_perp)*u_dot_Fraw*ux
#             Fasy = mu_perp*Fy_raw + (mu_par - mu_perp)*u_dot_Fraw*uy
#             Fasz = mu_perp*Fz_raw + (mu_par - mu_perp)*u_dot_Fraw*uz

#             Fijx += Fasx
#             Fijy += Fasy
#             Fijz += Fasz

#             # ROTACIÓN: usa la fuerza "isotrópica" de tu código (aquí igual que la cruda)
#             Fx_iso = Fx_raw
#             Fy_iso = Fy_raw
#             Fz_iso = Fz_raw

#             # Brazo desde CM_i al punto de aplicación (xpi,ypi,zpi)
#             rix = xpi - x
#             riy = ypi - y
#             riz = zpi - z

#             τx = riy*Fz_iso - riz*Fy_iso
#             τy = riz*Fx_iso - rix*Fz_iso
#             τz = rix*Fy_iso - riy*Fx_iso

#             τx_sum += τx
#             τy_sum += τy
#             τz_sum += τz
#         end
#     end

#     # Proyección del par total a (θ, φ)
#     # Escala rotacional: tu versión llevaba *mu_perp; la mantenemos
#     scale = 12.0 / ((l + d)^2) * mu_perp
#     Wijtheta = scale * (τx_sum*sθx + τy_sum*sθy + τz_sum*sθz)
#     Wijphi   = scale * (τx_sum*sφx + τy_sum*sφy + τz_sum*sφz)

#     return Fijx, Fijy, Fijz, Wijtheta, Wijphi
# end

# # Perfil de radio con suavizado C¹ en las puntas
# # s ∈ [-L/2, L/2], a = radio del cuerpo, a_end = radio en la punta, l_tip = longitud de la zona de punta
# # p = exponente de suavizado (p=4 o 6 suele ir muy bien)
# @inline function radius_profile_C1(s, L, a, a_end, l_tip, p=1)
#     sabs = abs(s)
#     s0   = L/2 - l_tip
#     if sabs <= s0
#         return a
#     else
#         # t en [0,1] a lo largo de la punta
#         t = (sabs - s0) / l_tip
#         # Suavizado polinómico C¹: valor y derivada nulas en el empalme
#         # r = a_end + (a - a_end)*(1 - t^p)
#         return a_end + (a - a_end) * (1.0 - t^p)
#     end
# end


# function repulsiveForces_rods3dline_asym_cutC1(
#     x,y,z,d,l,theta,phi,                     # rod i (centro, diámetro, largo, orientación)
#     x2,y2,z2,d2,l2,theta2,phi2,              # rod j
#     eta, Eeff, A, K;                         # Eeff = E* (módulo efectivo de Hertz)
#     a_end_ratio = 0.8,                       # r_tip = a_end_ratio * a
#     l_tip  = 0.2,                      # l_tip = l_tip_frac * L
#     p_smooth    = 3                          # exponente de suavizado (C¹)
# )

#     # Acumuladores de fuerza y par (sobre i por j)
#     Fijx = 0.0; Fijy = 0.0; Fijz = 0.0
#     τx_sum = 0.0; τy_sum = 0.0; τz_sum = 0.0

#     # Dirección del rod i (φ = elevación)
#     ux =  cos(theta)*cos(phi)
#     uy =  sin(theta)*cos(phi)
#     uz =  sin(phi)

#     # Dirección del rod j
#     ux2 =  cos(theta2)*cos(phi2)
#     uy2 =  sin(theta2)*cos(phi2)
#     uz2 =  sin(phi2)

#     # Derivadas para proyección (asumiendo tu convención)
#     duθx, duθy, duθz = du_dtheta(theta, phi)
#     duφx, duφy, duφz = du_dphi(theta,  phi)

#     # sθ = u × du/dθ ; sφ = u × du/dφ
#     sθx = uy*duθz - uz*duθy
#     sθy = uz*duθx - ux*duθz
#     sθz = ux*duθy - uy*duθx

#     sφx = uy*duφz - uz*duφy
#     sφy = uz*duφx - ux*duφz
#     sφz = ux*duφy - uy*duφx

#     # Movilidad anisótropa (como en tu versión actual)
#     mu_par  = 1.0 / eta           # μ∥
#     mu_perp = 1.0 / (eta * A)     # μ⊥
#     delta_s = l / (K - 1)         # separación entre puntos muestreados

#     # Radios base
#     a_i  = 0.5*d
#     a_j  = 0.5*d2
    
#     a_end_i = a_end_ratio * a_i
#     a_end_j = a_end_ratio * a_j

#     for s in LinRange(-l/2, l/2, K)
#         # Punto sobre el eje del rod i
#         xpi = x + s*ux
#         ypi = y + s*uy
#         zpi = z + s*uz

#         # Punto más cercano sobre el rod j (como hacías con los polos)
#         xiAux, yiAux, ziAux,  xjAux, yjAux, zjAux =
#             CBMMetrics.rodIntersection3d_(xpi, ypi, zpi, 0.0, theta, phi,
#                                          x2,  y2,  z2,  l2, theta2, phi2)

#         # Vector j->i y distancia
#         dx = xiAux - xjAux
#         dy = yiAux - yjAux
#         dz = ziAux - zjAux
#         rij = sqrt(dx*dx + dy*dy + dz*dz)

#         # --- NUEVO: radios locales con suavizado C¹ ---
#         # r_i en s (ya lo tenemos)
#         r_i = radius_profile_C1(s, l, a_i, a_end_i, l_tip, p_smooth)

#         # Necesitamos s_j (parámetro axial del punto de j)
#         # Proyección del vector desde el centro de j hasta (xjAux,yjAux,zjAux) sobre el eje u2
#         sx = xjAux - x2; sy = yjAux - y2; sz = zjAux - z2
#         s_j = sx*ux2 + sy*uy2 + sz*uz2
#         # Clampear a [-l2/2, l2/2] por robustez numérica
#         if s_j >  l2/2; s_j =  l2/2; elseif s_j < -l2/2; s_j = -l2/2; end

#         r_j = radius_profile_C1(s_j, l2, a_j, a_end_j, l_tip, p_smooth)

#         # --- solapamiento local (Hertz) ---
#         # contacto si rij < r_i + r_j
#         if rij > 0.0
#             hAux = (r_i + r_j) - rij
#             if hAux > 0.0
#                 nijx = dx / rij
#                 nijy = dy / rij
#                 nijz = dz / rij

#                 # Curvatura efectiva local (esfera-esfera con radios locales)
#                 # Robusto y estándar para ABM con cápsulas:
#                 Rstar = (r_i * r_j) / (r_i + r_j)

#                 # Hertz local: F = (4/3) Eeff sqrt(R*) delta^(3/2)
#                 # factor delta_s para integrar a lo largo del eje de i (como haces tú)
#                 # Fuerza cruda de contacto (SIN /eta)

               
#                 F_mag = delta_s * Eeff * sqrt(Rstar) * (hAux^(1.5))

#                 # Fuerza normal
#                 Fx_raw = F_mag * nijx
#                 Fy_raw = F_mag * nijy
#                 Fz_raw = F_mag * nijz

#                 # TRASLACIÓN: movilidad anisótropa (tu estilo actual)
#                 u_dot_Fraw = ux*Fx_raw + uy*Fy_raw + uz*Fz_raw
#                 Fasx = mu_perp*Fx_raw + (mu_par - mu_perp)*u_dot_Fraw*ux
#                 Fasy = mu_perp*Fy_raw + (mu_par - mu_perp)*u_dot_Fraw*uy
#                 Fasz = mu_perp*Fz_raw + (mu_par - mu_perp)*u_dot_Fraw*uz

#                 Fijx += Fasx
#                 Fijy += Fasy
#                 Fijz += Fasz

#                 # ROTACIÓN: usa la fuerza isotrópica para el brazo
#                 Fx_iso = Fx_raw; Fy_iso = Fy_raw; Fz_iso = Fz_raw

#                 # Brazo desde CM_i al punto de aplicación (xpi,ypi,zpi)
#                 rix = xpi - x
#                 riy = ypi - y
#                 riz = zpi - z


#                 τx = riy*Fz_iso - riz*Fy_iso
#                 τy = riz*Fx_iso - rix*Fz_iso
#                 τz = rix*Fy_iso - riy*Fx_iso

#                 τx_sum += τx
#                 τy_sum += τy
#                 τz_sum += τz
#             end
#         end
#     end

#     # Proyección del par total a (θ, φ)
#     # Escala rotacional: tu versión llevaba *mu_perp; la mantenemos
#     scale = 12.0 / ((l + d)^2) * mu_perp
#     Wijtheta = scale * (τx_sum*sθx + τy_sum*sθy + τz_sum*sθz)
#     Wijphi   = scale * (τx_sum*sφx + τy_sum*sφy + τz_sum*sφz)

#     return Fijx, Fijy, Fijz, Wijtheta, Wijphi
# end



function repulsiveForces_rods3d_asym(
    x,y,z,d,l,theta,phi,
    x2,y2,z2,d2,l2,theta2,phi2,
    eta, E, A)

    Fijx = 0.0; Fijy = 0.0; Fijz = 0.0
    Wijtheta = 0.0; Wijphi = 0.0

    # Puntos más cercanos
    xiAux, yiAux, ziAux,  xjAux, yjAux, zjAux =
        CBMMetrics.rodIntersection3d_(x, y, z, l, theta, phi,
                                     x2, y2, z2, l2, theta2, phi2)

    dx = xiAux - xjAux
    dy = yiAux - yjAux
    dz = ziAux - zjAux
    rij = sqrt(dx*dx + dy*dy + dz*dz)

    # Dirección del rod i
    ux =  cos(theta)*cos(phi)
    uy =  sin(theta)*cos(phi)
    uz =  sin(phi)
    # (por si acaso) — debería estar ya normalizado
    # un = inv(sqrt(ux^2+uy^2+uz^2)); ux*=un; uy*=un; uz*=un

    if rij > 0.0 && rij < (d + d2)/2
        hAux = (d + d2)/2 - rij
        nijx = dx / rij
        nijy = dy / rij
        nijz = dz / rij

        # Fuerza "cruda" de contacto SIN /eta
        Fraw_mag = E * sqrt(d2 * hAux^3) / (l + d)
        Fx_raw = Fraw_mag * nijx
        Fy_raw = Fraw_mag * nijy
        Fz_raw = Fraw_mag * nijz

        # === 1) TRASLACIÓN con movilidad anisótropa ===
        mu_par  = 1.0 / eta  # μ∥
        mu_perp = 1.0 / (eta*A)         # μ⊥
        u_dot_Fraw = ux*Fx_raw + uy*Fy_raw + uz*Fz_raw
        Fasx = mu_perp*Fx_raw + (mu_par - mu_perp)*u_dot_Fraw*ux
        Fasy = mu_perp*Fy_raw + (mu_par - mu_perp)*u_dot_Fraw*uy
        Fasz = mu_perp*Fz_raw + (mu_par - mu_perp)*u_dot_Fraw*uz

        Fijx = Fasx;  Fijy = Fasy;  Fijz = Fasz

        # === 2) ROTACIÓN con fuerza ISOTRÓPICA (no anisotropa) ===
        # fuerza isotrópica = Fraw / eta
        Fx_iso = Fx_raw
        Fy_iso = Fy_raw
        Fz_iso = Fz_raw 

        rx = (xiAux - x); ry = (yiAux - y); rz = (ziAux - z)
        τx = ry*Fz_iso - rz*Fy_iso
        τy = rz*Fx_iso - rx*Fz_iso
        τz = rx*Fy_iso - ry*Fx_iso

        # Proyección a (θ, φ)
        duθx, duθy, duθz = du_dtheta(theta, phi)
        duφx, duφy, duφz = du_dphi(theta,  phi)

        sθx = uy*duθz - uz*duθy
        sθy = uz*duθx - ux*duθz
        sθz = ux*duθy - uy*duθx

        sφx = uy*duφz - uz*duφy
        sφy = uz*duφx - ux*duφz
        sφz = ux*duφy - uy*duφx

        # Escala rotacional del caso 3D isotrópico
        scale = 12.0 / ((l + d)^2)*mu_perp

        Wijtheta = scale * (τx*sθx + τy*sθy + τz*sθz)
        Wijphi   = scale * (τx*sφx + τy*sφy + τz*sφz)
    end

    return Fijx, Fijy, Fijz, Wijtheta, Wijphi
end

# function attractiveForces_rodline_yukawa3d(
#     x,y,z,d,l,theta,phi,                  # rod i
#     x2,y2,z2,d2,l2,theta2,phi2,           # rod j
#     eta, eps, gamma,K; rcut_factor=3.0
# )
#     # Acumuladores fuerza/par (sobre i por j)
#     Fijx = 0.0; Fijy = 0.0; Fijz = 0.0
#     τx = 0.0; τy = 0.0; τz = 0.0

#     # Dirección del rod i
#     ux =  cos(theta)*cos(phi)
#     uy =  sin(theta)*cos(phi)
#     uz =  sin(phi)

#     # Geometría
#     σ    = 0.5*(d + d2)
#     rcut = rcut_factor * σ

#     # Muestreo uniforme a lo largo del rod i
#     for s in LinRange(-l/2, l/2, K)
#         # Punto sobre el eje del rod i
#         xpi = x + s*ux
#         ypi = y + s*uy
#         zpi = z + s*uz

#         # Punto más cercano sobre el rod j (como hacías con los polos)
#         xiAux, yiAux, ziAux,  xjAux, yjAux, zjAux =
#             CBMMetrics.rodIntersection3d_(xpi, ypi, zpi, 0.0, theta, phi,
#                                          x2,  y2,  z2,  l2, theta2, phi2)

#         # Vector j->i y distancia
#         dx = xiAux - xjAux
#         dy = yiAux - yjAux
#         dz = ziAux - zjAux
#         rij = sqrt(dx*dx + dy*dy + dz*dz)
#         h = rij - σ

#         if rij <= rcut && h > 0
#             # Yukawa atractivo desplazado por σ:
#             # u(r) = -eps * exp(-gamma*(r-σ)) / r
#             # F = -du/dr * n̂
#             Fmag = - eps * exp(-gamma*h) * (gamma*rij + 1.0) / (rij^2)

#             # Misma normalización que tú
#             Fmag /= (eta * (l + d))

#             inv_r = 1.0/rij
#             nijx = dx*inv_r;  nijy = dy*inv_r;  nijz = dz*inv_r

#             Fx = Fmag * nijx
#             Fy = Fmag * nijy
#             Fz = Fmag * nijz

#             Fijx += Fx;  Fijy += Fy;  Fijz += Fz

#             # Par respecto al CM de i
#             rix = xpi - x;  riy = ypi - y;  riz = zpi - z
#             τx += riy*Fz - riz*Fy
#             τy += riz*Fx - rix*Fz
#             τz += rix*Fy - riy*Fx
#         end
#     end

#     # Proyección a (θ, φ) (igual que en tu código)
#     duθ = du_dtheta(theta, phi)
#     duφ = du_dphi(theta,  phi)

#     sθx = uy*duθ[3] - uz*duθ[2]
#     sθy = uz*duθ[1] - ux*duθ[3]
#     sθz = ux*duθ[2] - uy*duθ[1]

#     sφx = uy*duφ[3] - uz*duφ[2]
#     sφy = uz*duφ[1] - ux*duφ[3]
#     sφz = ux*duφ[2] - uy*duφ[1]

#     scale = 12.0 / ((l + d)^2)
#     Wijtheta = scale * (τx*sθx + τy*sθy + τz*sθz)
#     Wijphi   = scale * (τx*sφx + τy*sφy + τz*sφz)

#     return Fijx, Fijy, Fijz, Wijtheta, Wijphi
# end

# function attractiveForces_rodline_spring(
#     x,y,z,d,l,theta,phi,                  # rod i
#     x2,y2,z2,d2,l2,theta2,phi2,           # rod j
#     eta, eps,K; rcut_factor=2
# )
#     # Acumuladores fuerza/par (sobre i por j)
#     Fijx = 0.0; Fijy = 0.0; Fijz = 0.0
#     τx = 0.0; τy = 0.0; τz = 0.0

#     # Dirección del rod i
#     ux =  cos(theta)*cos(phi)
#     uy =  sin(theta)*cos(phi)
#     uz =  sin(phi)

#     # Geometría
#     σ    = 0.5*(d + d2)
#     rcut = rcut_factor * σ
#     delta_s = l / (K - 1)
#     # Muestreo uniforme a lo largo del rod i
#     for s in LinRange(-l/2, l/2, K)
#         # Punto sobre el eje del rod i
#         xpi = x + s*ux
#         ypi = y + s*uy
#         zpi = z + s*uz

#         # Punto más cercano sobre el rod j (como hacías con los polos)
#         xiAux, yiAux, ziAux,  xjAux, yjAux, zjAux =
#             CBMMetrics.rodIntersection3d_(xpi, ypi, zpi, 0.0, theta, phi,
#                                          x2,  y2,  z2,  l2, theta2, phi2)

#         # Vector j->i y distancia
#         dx = xiAux - xjAux
#         dy = yiAux - yjAux
#         dz = ziAux - zjAux
#         rij = sqrt(dx*dx + dy*dy + dz*dz)
#         h = rij - σ


#         if rij <= rcut && h > 0
#             # Yukawa atractivo desplazado por σ:
#             # u(r) = -eps * exp(-gamma*(r-σ)) / r
#             # F = -du/dr * n̂
#             Fmag = - eps * h

#             # Misma normalización que tú
#             Fmag *= delta_s/(eta * (l + d))

#             inv_r = 1.0/rij
#             nijx = dx*inv_r;  nijy = dy*inv_r;  nijz = dz*inv_r

#             Fx = Fmag * nijx
#             Fy = Fmag * nijy
#             Fz = Fmag * nijz

#             Fijx += Fx;  Fijy += Fy;  Fijz += Fz

#            # Par respecto al CM de i
#             rix = xpi - x;  riy = ypi - y;  riz = zpi - z
#             τx += riy*Fz - riz*Fy
#             τy += riz*Fx - rix*Fz
#             τz += rix*Fy - riy*Fx
#         end
#     end

#     # Proyección a (θ, φ) (igual que en tu código)
#     duθ = du_dtheta(theta, phi)
#     duφ = du_dphi(theta,  phi)

#     sθx = uy*duθ[3] - uz*duθ[2]
#     sθy = uz*duθ[1] - ux*duθ[3]
#     sθz = ux*duθ[2] - uy*duθ[1]

#     sφx = uy*duφ[3] - uz*duφ[2]
#     sφy = uz*duφ[1] - ux*duφ[3]
#     sφz = ux*duφ[2] - uy*duφ[1]

#     scale = 12.0 / ((l + d)^2)
#     Wijtheta = scale * (τx*sθx + τy*sθy + τz*sθz)
#     Wijphi   = scale * (τx*sφx + τy*sφy + τz*sφz)

#     return Fijx, Fijy, Fijz, Wijtheta, Wijphi
# end

# function attractiveForces_rods_yukawa3d(
#     x,y,z,d,l,theta,phi,                  # rod i
#     x2,y2,z2,d2,l2,theta2,phi2,           # rod j
#     eta, eps, gamma; rcut_factor=1.8, hmin_factor=1e-3
# )
#     # Fuerza y par (sobre i por j)
#     Fijx = 0.0; Fijy = 0.0; Fijz = 0.0

    
#     # Parámetros geométricos
#     σ    = 0.5*(d + d2)
#     rcut = rcut_factor * σ

#     # 2) Vector y distancia j->i
#     dx = x - x2
#     dy = y - y2
#     dz = z - z2
#     rij = sqrt(dx*dx + dy*dy + dz*dz)

#     h = rij - σ

#     if rij <= rcut && h>0
#         # Yukawa atractivo desplazado por σ:
#         # u(r) = -eps * exp(-gamma*(r-σ)) / r
#         # F = -du/dr * n̂ ⇒ magnitud (negativa: hacia el otro polo)
#         Fmag = - eps * exp(-gamma*h) * (gamma*rij + 1.0) / (rij^2)

#         # Reescalado consistente
#         Fmag /= (eta * (l + d))

#         inv_r = 1.0/rij
#         nijx = dx*inv_r;  nijy = dy*inv_r;  nijz = dz*inv_r

#         Fx_p = Fmag * nijx
#         Fy_p = Fmag * nijy
#         Fz_p = Fmag * nijz

#         Fijx += Fx_p;  Fijy += Fy_p;  Fijz += Fz_p


#     end
      

#     return Fijx, Fijy, Fijz, 0, 0
# end

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
# function substrateAttraction_rods_yukawa3d(
#     x,y,z,d,l,theta,phi,
#     eta, eps_cs, gamma; rcut_h = 5.0
# )
#     a = 0 # plano en z=0
#     (x1,y1,z1), (x2,y2,z2) = rod_poles3d_ordered(x,y,z,l,theta,phi)

#     # Acumuladores
#     Fx = 0.0; Fy = 0.0; Fz = 0.0
#     τx = 0.0; τy = 0.0; τz = 0.0

#     # Recorremos polos (aproximación por extremos del esferocilindro)
#     for pole in 1:2
#         xp, yp, zp = pole == 1 ? (x1,y1,z1) : (x2,y2,z2)

#         # gap h = distancia "libre" polo–plano menos el radio a
#         #   s = zp  (altura del polo)
#         #   h = s - a
#         s  = zp
#         h  = s - a + d/2

#         # Solo aporta si está por encima del plano efectivo (h > h0)
#         # y dentro de un corte razonable
#         if h > 0
#             Δh = h 
#             if Δh < rcut_h
#                 # f(h) = exp(gamma*(h-h0)) / (h-h0)
#                 # dU/dh = eps_cs * f * (gamma - 1/(h-h0))
#                 # F = -dU/dh  (hacia -z si dU/dh > 0)
#                 f     = exp(-gamma * Δh) / (Δh)
#                 dUdh  = eps_cs * f * (-gamma - 1.0/Δh)
#                 Fmag  = - dUdh / (eta * (l + d))   # reescalado consistente

#                 # Fuerza solo en z (hacia el sustrato si Fmag < 0)
#                 Fx_p, Fy_p, Fz_p = 0.0, 0.0, Fmag
#                 Fx += Fx_p; Fy += Fy_p; Fz += Fz_p

#                 # Par respecto al CM: r × F
#                 rx = xp - x; ry = yp - y; rz = zp - z
#                 τx += ry*Fz_p - rz*Fy_p
#                 τy += rz*Fx_p - rx*Fz_p
#                 τz += rx*Fy_p - ry*Fx_p
#             end
#         end
#     end

#     # --- Proyección del par a (θ, φ) ---
#     ux =  cos(theta)*cos(phi)
#     uy =  sin(theta)*cos(phi)
#     uz =  sin(phi)

#     duθ = du_dtheta(theta, phi)   # (-cosφ sinθ,  cosφ cosθ, 0)
#     duφ = du_dphi(theta,  phi)    # (-sinφ cosθ, -sinφ sinθ, cosφ)

#     # sθ = u × du/dθ
#     sθx = uy*duθ[3] - uz*duθ[2]
#     sθy = uz*duθ[1] - ux*duθ[3]
#     sθz = ux*duθ[2] - uy*duθ[1]

#     # sφ = u × du/dφ
#     sφx = uy*duφ[3] - uz*duφ[2]
#     sφy = uz*duφ[1] - ux*duφ[3]
#     sφz = ux*duφ[2] - uy*duφ[1]

#     # Mismo factor de escala que usas en el resto del código
#     scale = 12.0 / ((l + d)^2)
#     Wθ = scale * (τx*sθx + τy*sθy + τz*sθz)
#     Wφ = scale * (τx*sφx + τy*sφy + τz*sφz)

#     return Fx, Fy, Fz, Wθ, Wφ
# end

function substrateAttraction_rods_yukawa3d(
    x,y,z,d,l,theta,phi,
    eta, eps_cs, gamma; rcut_h = 5.0
)
    a = 0.0  # plano en z = 0
    (x1,y1,z1), (x2,y2,z2) = rod_poles3d_ordered(x,y,z,l,theta,phi)

    Fx = 0.0; Fy = 0.0; Fz = 0.0
    τx = 0.0; τy = 0.0; τz = 0.0

    for pole in 1:2
        xp, yp, zp = pole == 1 ? (x1,y1,z1) : (x2,y2,z2)

        # Gap libre esfera-plano
        h = zp - a - d/2

        if 0.0 < h < rcut_h
            # Potencial atractivo:
            # U(h) = -eps_cs * exp(-gamma*h) / h
            #
            # dU/dh = eps_cs * exp(-gamma*h) * (gamma/h + 1/h^2)
            # Fz = -dU/dh
            expterm = exp(-gamma * h)
            dUdh = eps_cs * expterm * (gamma/h + 1.0/h^2)

            Fmag = -dUdh / (eta * (l + d))

            Fx_p = 0.0
            Fy_p = 0.0
            Fz_p = Fmag

            Fx += Fx_p
            Fy += Fy_p
            Fz += Fz_p

            rx = xp - x
            ry = yp - y
            rz = zp - z

            τx += ry*Fz_p - rz*Fy_p
            τy += rz*Fx_p - rx*Fz_p
            τz += rx*Fy_p - ry*Fx_p
        end
    end

    ux = cos(theta)*cos(phi)
    uy = sin(theta)*cos(phi)
    uz = sin(phi)

    duθ = du_dtheta(theta, phi)
    duφ = du_dphi(theta, phi)

    sθx = uy*duθ[3] - uz*duθ[2]
    sθy = uz*duθ[1] - ux*duθ[3]
    sθz = ux*duθ[2] - uy*duθ[1]

    sφx = uy*duφ[3] - uz*duφ[2]
    sφy = uz*duφ[1] - ux*duφ[3]
    sφz = ux*duφ[2] - uy*duφ[1]

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

function padRepulsion_rods_point_3d(
    x,y,z, xp, yp, zp, d,l,theta,phi,
    eta, Ebp, hPad
)
    a = hPad
    

    # Acumuladores netos
    Fx = 0.0; Fy = 0.0; Fz = 0.0
    τx = 0.0; τy = 0.0; τz = 0.0

    
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


function padRepulsion_rods_point_3d_perp(
    x, y, z,
    xs, ys, zs,
    nx, ny, nz,
    d, l, theta, phi,
    eta, Ebp, hPad
)
    Fx = 0.0; Fy = 0.0; Fz = 0.0
    τx = 0.0; τy = 0.0; τz = 0.0

    δ = zs - hPad

    if δ > 0.0
        Fmag = Ebp * sqrt(d * δ^3) / (eta * (l + d))

        Fx = -Fmag * nx
        Fy = -Fmag * ny
        Fz = -Fmag * nz

        rx = xs - x
        ry = ys - y
        rz = zs - z

        τx = ry*Fz - rz*Fy
        τy = rz*Fx - rx*Fz
        τz = rx*Fy - ry*Fx
    end

    ux = cos(theta)*cos(phi)
    uy = sin(theta)*cos(phi)
    uz = sin(phi)

    duθ = du_dtheta(theta, phi)
    duφ = du_dphi(theta, phi)

    sθx = uy*duθ[3] - uz*duθ[2]
    sθy = uz*duθ[1] - ux*duθ[3]
    sθz = ux*duθ[2] - uy*duθ[1]

    sφx = uy*duφ[3] - uz*duφ[2]
    sφy = uz*duφ[1] - ux*duφ[3]
    sφz = ux*duφ[2] - uy*duφ[1]

    scale = 12.0 / ((l + d)^2)
    Wθ = scale * (τx*sθx + τy*sθy + τz*sθz)
    Wφ = scale * (τx*sφx + τy*sφy + τz*sφz)

    return Fx, Fy, Fz, Wθ, Wφ
end



# function apply_pad_contact_along_rod(
#     # estado bacteria
#     x, y, z, d, l, theta, phi,
#     # parámetros físicos
#     eta, Ebp,
#     # pad
#     HPad, simBox, NMedium,
#     # paso sobre el rod
#     Δs::Float64 = 1.0; F_pad = nothing
# )
    
#     if F_pad === nothing
#         F_pad = zeros(size(HPad))
#     end
#     # orientación del rod (tu convención)
#     ux =  cos(theta)*cos(phi)
#     uy =  sin(theta)*cos(phi)
#     uz =  sin(phi)

#     # grid spacing
#     hx = (simBox[1,2] - simBox[1,1]) / NMedium[1]
#     hy = (simBox[2,2] - simBox[2,1]) / NMedium[2]

#     Δs_eff = min(Δs, 0.5*min(hx, hy))  # para no saltar celdas


#     # Acumuladores totales bacteria
#     Fx_tot = 0.0; Fy_tot = 0.0; Fz_tot = 0.0
#     Wθ_tot = 0.0; Wφ_tot = 0.0

#     # recorrer s en [-l/2, l/2]
#     K = max(2, ceil(Int, l / Δs_eff))  # aseguras al menos 2
#     for s in LinRange(-l/2, l/2, K)

#         # punto en el rod
#         xp = x + s*ux
#         yp = y + s*uy
#         zp = z + s*uz

#         # celda del pad
#         idx = Int(floor(Int, xp/(simBox[1,2]-simBox[1,1])*NMedium[1]) + NMedium[1]/2)
#         if mod(xp, (simBox[1,2]-simBox[1,1])/NMedium[1])>0
#             idx += 1
#         end
#         idy = Int(floor(Int, yp/(simBox[2,2]-simBox[2,1])*NMedium[2])+ NMedium[2]/2)
#         if mod(yp, (simBox[2,2]-simBox[2,1])/NMedium[2])>0
#             idy += 1
#         end
        

#         # h local del pad (tu altura / “a” local)
#         hPad = HPad[idx, idy, 1]

#         # fuerza en ese punto (tu función)
#         Fx, Fy, Fz, Wθ, Wφ = padRepulsion_rods_point_3d(
#             x, y, z,  # CM bacteria
#             xp, yp, zp,  # punto sobre el rod
#             d, l, theta, phi,
#             eta, Ebp, hPad
#         )

#         # acumula bacteria
#         Fx_tot += Fx; Fy_tot += Fy; Fz_tot += Fz
#         Wθ_tot += Wθ; Wφ_tot += Wφ
#         F_pad[idx, idy, 1] += -Fz*eta*(l+d)

#     end

#     return Fx_tot, Fy_tot, Fz_tot, Wθ_tot, Wφ_tot, F_pad
# end




# function padRepulsion_rods3d(
#     x,y,z,d,l,theta,phi,
#     eta, Ebp, hPad
# )
#     a = hPad
#     (x1,y1,z1), (x2,y2,z2) = rod_poles3d_ordered(x,y,z,l,theta,phi)

#     # Acumuladores netos
#     Fx = 0.0; Fy = 0.0; Fz = 0.0
#     τx = 0.0; τy = 0.0; τz = 0.0

#     for pole in 1:2
#         xp, yp, zp = pole == 1 ? (x1,y1,z1) : (x2,y2,z2)

#         s = zp
#         δ = s + d/2 -a

#         # --- Repulsión pad (esfera–plano) ---
#         if δ > 0.0
#             Fh = -Ebp * sqrt(d * δ^3) / (eta * (l + d)) # negativo: hacia -z
#             # Fuerza en el polo: (0,0,Fh)
#             Fx += 0.0; Fy += 0.0; Fz +=  Fh
#             # Par respecto al CM: r × F
#             rx = xp - x; ry = yp - y; rz = zp - z
#             τx += ry*Fh - rz*0.0
#             τy += rz*0.0 - rx*Fh
#             τz += rx*0.0 - ry*0.0
#         end

#     end

#     # --- Proyección correcta del par a (θ, φ) ---
#     ux =  cos(theta)*cos(phi)
#     uy =  sin(theta)*cos(phi)
#     uz =  sin(phi)

#     duθ = du_dtheta(theta, phi)   # (-cosφ sinθ,  cosφ cosθ, 0)
#     duφ = du_dphi(theta,  phi)    # (-sinφ cosθ, -sinφ sinθ, cosφ)

#     # sθ = u × du/dθ
#     sθx = uy*duθ[3] - uz*duθ[2]
#     sθy = uz*duθ[1] - ux*duθ[3]
#     sθz = ux*duθ[2] - uy*duθ[1]

#     # sφ = u × du/dφ
#     sφx = uy*duφ[3] - uz*duφ[2]
#     sφy = uz*duφ[1] - ux*duφ[3]
#     sφz = ux*duφ[2] - uy*duφ[1]

#     scale = 12.0 / ((l + d)^2)
#     Wθ = scale * (τx*sθx + τy*sθy + τz*sθz)
#     Wφ = scale * (τx*sφx + τy*sφy + τz*sφz)

#     return Fx, Fy, Fz, Wθ, Wφ
# end

# # --- Fuerzas de pared y pili en 3D para un rod ---

# function wallForces_rods_bond3d(
#     x,y,z,d,l,theta,phi,
#     eta, Ebv,
#     kb, Δ0, Δmax, T_eng,hax,
#     xpili1::Float64, ypili1::Float64, t_pili1::Int,
#     xpili2::Float64, ypili2::Float64, t_pili2::Int
# )
#     a = 0.5*d
#     (x1,y1,z1), (x2,y2,z2) = rod_poles3d_ordered(x,y,z,l,theta,phi)

#     # Acumuladores netos
#     Fx = 0.0; Fy = 0.0; Fz = 0.0
#     τx = 0.0; τy = 0.0; τz = 0.0

#     for pole in 1:2
#         xp, yp, zp = pole == 1 ? (x1,y1,z1) : (x2,y2,z2)
#         xpili = pole == 1 ? xpili1 : xpili2
#         ypili = pole == 1 ? ypili1 : ypili2
#         t_pili = pole == 1 ? t_pili1 : t_pili2

#         s = zp
#         δ = a - s

#         # --- Repulsión suelo (esfera–plano) ---
#         if δ > 0.0
#             Fh = Ebv * sqrt(d * δ^3) / (eta * (l + d))
#             # Fuerza en el polo: (0,0,Fh)
#             Fx += 0.0; Fy += 0.0; Fz += Fh
#             # Par respecto al CM: r × F
#             rx = xp - x; ry = yp - y; rz = zp - z
#             τx += ry*Fh - rz*0.0
#             τy += rz*0.0 - rx*Fh
#             τz += rx*0.0 - ry*0.0
#         end

#         # --- PILI ---
#         if t_pili < 0 # Enganchado
#             xa, ya, za = xpili, ypili, 0.0
#             dx = xp - xa; dy = yp - ya; dz = zp - za
#             Δ = sqrt(dx*dx + dy*dy + dz*dz)
#             if Δ > 1e-12
#                 Fspr = -kb * (Δ - Δ0)
#                 if (Δ > Δmax) 
#                     t_pili = 0; xpili = NaN; ypili = NaN
#                 else
#                     Fx_p = Fspr * (dx/Δ) / (eta * (l + d))
#                     Fy_p = Fspr * (dy/Δ) / (eta * (l + d))
#                     Fz_p = Fspr * (dz/Δ) / (eta * (l + d))
#                     Fx += Fx_p; Fy += Fy_p; Fz += Fz_p
#                     rx = xp - x; ry = yp - y; rz = zp - z
#                     τx += ry*Fz_p - rz*Fy_p
#                     τy += rz*Fx_p - rx*Fz_p
#                     τz += rx*Fy_p - ry*Fx_p
#                 end
#             end

#         else  # No enganchado, contando tiempo
#             if (s - a) < hax
#                 if t_pili >= T_eng
#                     t_pili = -1; xpili = xp; ypili = yp
#                 else
#                     t_pili += 1
#                 end
#             else
#                 t_pili = 0
#             end
#         end

#         if pole==1
#             xpili1=xpili; ypili1=ypili; t_pili1=t_pili
#         else
#             xpili2=xpili; ypili2=ypili; t_pili2=t_pili
#         end
#     end

#     # --- Proyección correcta del par a (θ, φ) ---
#     ux =  cos(theta)*cos(phi)
#     uy =  sin(theta)*cos(phi)
#     uz =  sin(phi)

#     duθ = du_dtheta(theta, phi)   # (-cosφ sinθ,  cosφ cosθ, 0)
#     duφ = du_dphi(theta,  phi)    # (-sinφ cosθ, -sinφ sinθ, cosφ)

#     # sθ = u × du/dθ
#     sθx = uy*duθ[3] - uz*duθ[2]
#     sθy = uz*duθ[1] - ux*duθ[3]
#     sθz = ux*duθ[2] - uy*duθ[1]

#     # sφ = u × du/dφ
#     sφx = uy*duφ[3] - uz*duφ[2]
#     sφy = uz*duφ[1] - ux*duφ[3]
#     sφz = ux*duφ[2] - uy*duφ[1]

#     scale = 12.0 / ((l + d)^2)
#     Wθ = scale * (τx*sθx + τy*sθy + τz*sθz)
#     Wφ = scale * (τx*sφx + τy*sφy + τz*sφz)

#     return Fx, Fy, Fz, Wθ, Wφ,
#            xpili1, ypili1, t_pili1,
#            xpili2, ypili2, t_pili2
# end
