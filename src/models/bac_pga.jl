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
