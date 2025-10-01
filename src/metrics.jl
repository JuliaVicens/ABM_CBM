module CBMMetrics

    using CUDA
    import CellBasedModels: AGENT

    """
        function cellInMesh(edge,x,xMin,xMax,nX) 

    Give the integer position in a regular discrete mesh with poits separamtions of `edge`, given a position in `x`. The simulation domain being (`xMin`, `xMax`) and maximum number of mesh points `nX`.

    e.g Grid with 6 points at [0,.1,.2,.3,.4,.5]
    ```
    >>> cellInMesh(.1,.29,0.,.5,6)
    4
    ```
    which if the closest point in the mesh.
    """
    function cellInMesh(edge,x,xMin,xMax,nX) 
        return if x > xMax nX elseif x < xMin 1 else Int((x-xMin)÷edge)+1 end
    end

    ##############################################################################################################################
    # Distance metrics
    ##############################################################################################################################
    """
        euclidean(x1,x2)
        euclidean(x1,x2,y1,y2)
        euclidean(x1,x2,y1,y2,z1,z2)

    Euclidean distance metric between two positions.

    d = (x₁-x₂)²
    """
    euclidean(x1,x2) = sqrt((x1-x2)^2)
    euclidean(x1,x2,y1,y2) = sqrt((x1-x2)^2+(y1-y2)^2)
    euclidean(x1,x2,y1,y2,z1,z2) = sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)

    """
        macro euclidean(it2)
        macro euclidean(it1,it2)

    Macro that given the iterator symbos it1 and it2, give the corresponding euclidean distanve in the correct dimentions. 
    If it1 is not provided it asumes the default iteration index of agents (i1_).
    """
    macro euclidean(it2)

        abm = AGENT
    
        args = [:(x[i1_]), :(x[$it2])]
        if abm.dims > 1
            args = [args;[:(y[i1_]), :(y[$it2])]]
        end
        if abm.dims > 2
            args = [args;[:(z[i1_]), :(z[$it2])]]
        end
    
        return esc(:(CBMMetrics.euclidean($(args...))))
    
    end
    
    macro euclidean(it1,it2)
    
        abm = AGENT
    
        args = [:(x[$it1]), :(x[$it2])]
        if abm.dims > 1
            args = [args;[:(y[$it1]), :(y[$it2])]]
        end
        if abm.dims > 2
            args = [args;[:(z[$it1]), :(z[$it2])]]
        end
    
        return esc(:(CBMMetrics.euclidean($(args...))))
    
    end

    """
        manhattan(x1,x2)
        manhattan(x1,x2,y1,y2)
        manhattan(x1,x2,y1,y2,z1,z2)

    Manhattan distance metric between two positions.

    d = |x₁-x₂|
    """
    manhattan(x1,x2) = abs(x1-x2)
    manhattan(x1,x2,y1,y2) = abs(x1-x2)+abs(y1-y2)
    manhattan(x1,x2,y1,y2,z1,z2) = abs(x1-x2)+abs(y1-y2)+abs(z1-z2)
    
    """
        macro manhattan(it2)
        macro manhattan(it1,it2)

    Macro that given the iterator symbos it1 and it2, give the corresponding manhattan distanve in the correct dimentions. 
    If it1 is not provided it asumes the default iteration index of agents (i1_).
    """
    macro manhattan(it2)
    
        abm = AGENT
    
        args = [:(x[i1_]), :(x[$it2])]
        if abm.dims > 1
            args = [args;[:(y[i1_]), :(y[$it2])]]
        end
        if abm.dims > 2
            args = [args;[:(z[i1_]), :(z[$it2])]]
        end
    
        return esc(:(CBMMetrics.manhattan($(args...))))
    
    end

    macro manhattan(it1,it2)
    
        abm = AGENT
    
        args = [:(x[$it1]), :(x[$it2])]
        if abm.dims > 1
            args = [args;[:(y[$it1]), :(y[$it2])]]
        end
        if abm.dims > 2
            args = [args;[:(z[$it1]), :(z[$it2])]]
        end
    
        return esc(:(CBMMetrics.manhattan($(args...))))
    
    end

    """
        function intersection2lines(x1,y1,theta1,x2,y2,theta2,inf_eff=100000)
    
    Finds the point of intersection of two lines. You have to provide a point in space and and angle for eachline: (x1,y1,theta1) and (x2,y2,theta2).
    
    If the lines are parallel, it returns a point effectively in the infinite. The effective distance is described by `inf_eff``.
    
    Returns the point of intersection.
    """
    function intersection2lines(x1,y1,theta1,x2,y2,theta2,inf_eff=100000)

        pxIntersect = 0.
        pyIntersect = 0.
    
        cxAux = (x1-x2)
        cyAux = (y1-y2)
        normAux = cos(theta1)*sin(theta2)-sin(theta1)*cos(theta2)
        if abs(normAux) > 0.0000001
            scaleAux = (-sin(theta2)*cxAux+cos(theta2)*cyAux)/normAux
            pxIntersect = scaleAux*cos(theta1)+x1
            pyIntersect = scaleAux*sin(theta1)+y1
        else #if parallel send  point to an infinite
            pxIntersect = (x1+x2)/2
            pyIntersect = (y1+y2)/2
        end
    
        return pxIntersect,pyIntersect
    end

    """
        function point2line(x1,y1,x2,y2,theta2)
    
    Given a point (x1,x2), finds the closest point projected over a line described by a point in the line and the angle: (x2,y2,theta2).
    
    Returns the coordinates of the closest point over the line.
    """
    function point2line(x1,y1,x2,y2,theta2)
        #Compute closest point over the line axis of the other rod
        cx = x2-x1
        cy = y2-y1
        xjAux = (sin(theta2)*cx-cos(theta2)*cy)*sin(theta2)+x1
        yjAux = -(sin(theta2)*cx-cos(theta2)*cy)*cos(theta2)+y1
    
        return xjAux,yjAux
    end

    """
        function pointInsideRod(x1,y1,l1,theta1,pxAux,pyAux,separation)
    
    Given a line segment described by it central point (x1,y1), its angle in the plate theta1 and its length l1;
    and given a point over the save line (pxAux,pyAux), returns the point if inside the segment or the closes extreme of the segment. 
    If provided a separation (0,1), it moves the point that separation.
            
    Returns the coordinates of the closest point over the segment.
    """
    function pointInsideRod(x1,y1,l1,theta1,pxAux,pyAux,separation)
    
        di = max(sqrt((x1-pxAux)^2+(y1-pyAux)^2),0.00000001)
    
        dxi = (pxAux-x1)/di
        dyi = (pyAux-y1)/di
        return separation*min(di,l1/2)*dxi+x1, separation*min(di,l1/2)*dyi+y1
    
    end

    """
        function rodIntersection(x1,y1,l1,theta1,x2,y2,l2,theta2;separation=0.99)
    
    Given a two line segment described by it central point (x1,y2), its angle in the plate theta1 and its length l1;
    finds the closest spheres of both segments.
    
    Returns the coordinates of the closest spheres (x1Aux,y1Aux), (x2Aux,y2Aux).
    """
    function rodIntersection(x1,y1,l1,theta1,x2,y2,l2,theta2;separation=0.99)
            
        #Compute distance between centers of mass
        x1Aux = x1; x2Aux = x2; y1Aux = y1; y2Aux = y2; #Declare them in the global scope
    
        #Compute intersecting point of the extended direction
        pxAux, pyAux = intersection2lines(x1,y1,theta1,x2,y2,theta2)
    
        #Compute distance from mass center of both rods
        di = sqrt((x1-pxAux)^2+(y1-pyAux)^2)
        dj = sqrt((x2-pxAux)^2+(y2-pyAux)^2)
        normAux = cos(theta1)*sin(theta2)-sin(theta1)*cos(theta2)
        if abs(normAux) < 0.000001
            x1Aux,y1Aux= point2line(pxAux,pyAux,x1,y1,theta1)
            x1Aux,y1Aux = pointInsideRod(x1,y1,l1,theta1,x1Aux,y1Aux,separation)
    
            x2Aux,y2Aux= point2line(pxAux,pyAux,x2,y2,theta2)
            x2Aux,y2Aux = pointInsideRod(x2,y2,l2,theta2,x2Aux,y2Aux,separation) 
        elseif di<=l1/2 && dj<=l2/2 #Case that the intersecting point lies inside both rods
            x1Aux,y1Aux = pointInsideRod(x1,y1,l1,theta1,pxAux,pyAux,separation)
            x2Aux,y2Aux = pointInsideRod(x2,y2,l2,theta2,pxAux,pyAux,separation) 
        elseif di<=l1/2 && dj>l2/2
            x2Aux,y2Aux = pointInsideRod(x2,y2,l2,theta2,pxAux,pyAux,separation)    
            x1Aux,y1Aux= point2line(x2Aux,y2Aux,x1,y1,theta1)
        elseif di>l1/2 && dj<=l2/2
            x1Aux,y1Aux = pointInsideRod(x1,y1,l1,theta1,pxAux,pyAux,separation)    
            x2Aux,y2Aux= point2line(x1Aux,y1Aux,x2,y2,theta2)  
        else
            x2Aux,y2Aux = pointInsideRod(x2,y2,l2,theta2,pxAux,pyAux,separation)    
            x1Aux,y1Aux= point2line(x2Aux,y2Aux,x1,y1,theta1)
            dj = sqrt((x1-x1Aux)^2+(y1-y1Aux)^2)
            if dj > l1/2
                x1Aux,y1Aux = pointInsideRod(x1,y1,l1,theta1,pxAux,pyAux,separation)
                x2Aux_,y2Aux_= point2line(x1Aux,y1Aux,x2,y2,theta2)
    
                dj = sqrt((x2-x2Aux_)^2+(y2-y2Aux_)^2)
                if dj < l2/2
                    x2Aux = x2Aux_
                    y2Aux = y2Aux_
                end
            end
        end
    
        return x1Aux,y1Aux,x2Aux,y2Aux
    end
    



    # ===== utilidades comunes (simples) =====
    @inline dot3(a,b) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
    @inline sub3(a,b) = (a[1]-b[1], a[2]-b[2], a[3]-b[3])
    @inline add3(a,b) = (a[1]+b[1], a[2]+b[2], a[3]+b[3])
    @inline mul3(s,a) = (s*a[1], s*a[2], s*a[3])
    @inline norm3(a) = sqrt(dot3(a,a))

    # acimut-elevación -> vector unitario
    @inline function dir_from_angles(theta, phi)
        cφ = cos(phi)
        (cφ*cos(theta), cφ*sin(theta), sin(phi))
    end

    # proyección de P sobre la recta (O + t u), u unitario
    @inline function project_point_on_line(P, O, u)
        t = dot3(sub3(P,O), u)
        add3(O, mul3(t,u)), t
    end

    # ===== 1) intersection2lines -> 3D (devuelve UN punto como el original) =====
    """
        intersection2lines3d(x1,y1,z1,theta1,phi1, x2,y2,z2,theta2,phi2; tol=1e-7, inf_eff=1e5)

    Intersección de dos rectas en 3D siguiendo la filosofía 2D original:
    - Si se cruzan, devuelve el punto de cruce.
    - Si son alabeadas (no se cortan), devuelve el **punto medio** de los puntos más cercanos.
    - Si son (casi) paralelas, devuelve el **punto medio de los centros** (igual que tu 2D).

    Retorna: (px, py, pz)
    """
    function intersection2lines3d(x1,y1,z1,theta1,phi1, x2,y2,z2,theta2,phi2; tol=1e-7, inf_eff=1e5)
        O1 = (x1,y1,z1); u1 = dir_from_angles(theta1, phi1)
        O2 = (x2,y2,z2); u2 = dir_from_angles(theta2, phi2)

        r  = sub3(O1,O2)
        a  = 1.0
        b  = dot3(u1,u2)
        c  = 1.0
        d  = dot3(u1,r)
        e  = dot3(u2,r)
        den = a*c - b*b

        if abs(den) > tol
            t = (b*e - c*d)/den
            s = (a*e - b*d)/den
            P1 = add3(O1, mul3(t,u1))
            P2 = add3(O2, mul3(s,u2))
            # si realmente se cruzan, P1≈P2; devolvemos el medio igualmente
            Pm = mul3(0.5, add3(P1,P2))
            return Pm[1], Pm[2], Pm[3]
        else
            # paralelas o casi: como en tu 2D, el "punto al infinito" lo sustituimos por el medio de centros
            px = (x1 + x2)/2
            py = (y1 + y2)/2
            pz = (z1 + z2)/2
            return px, py, pz
        end
    end

    # ===== 2) point2line -> 3D =====
    """
        point2line3d(x1,y1,z1, x2,y2,z2, theta2,phi2)

    Proyecta el punto (x1,y1,z1) sobre la recta que pasa por (x2,y2,z2)
    con dirección (theta2,phi2). Devuelve las coords del punto proyectado.
    """
    function point2line3d(x1,y1,z1, x2,y2,z2, theta2,phi2)
        P = (x1,y1,z1)
        O = (x2,y2,z2)
        u = dir_from_angles(theta2, phi2)
        Q,_ = project_point_on_line(P, O, u)
        return Q[1], Q[2], Q[3]
    end

    # ===== 3) pointInsideRod -> 3D (misma lógica que tu 2D) =====
    """
        pointInsideRod3d(x1,y1,z1, l1, theta1,phi1, px,py,pz, separation)

    Igual que tu versión 2D: asume que (px,py,pz) está sobre la misma recta (o cerca),
    coge el vector desde el centro al punto y lo recorta a longitud l1/2, aplicando `separation`.
    """
    function pointInsideRod3d(x1,y1,z1, l1, theta1,phi1, px,py,pz, separation)
        # nota: theta1,phi1 no se usan aquí, igual que en tu 2D
        di = max(norm3((x1-px, y1-py, z1-pz)), 1e-8)
        dxi = (px - x1)/di
        dyi = (py - y1)/di
        dzi = (pz - z1)/di
        scale = separation * min(di, l1/2)
        return x1 + scale*dxi, y1 + scale*dyi, z1 + scale*dzi
    end

    # ===== 4) rodIntersection -> 3D (estructura de casos como el original) =====
    """
        rodIntersection3d(x1,y1,z1,l1,theta1,phi1, x2,y2,z2,l2,theta2,phi2; separation=0.99, tol=1e-6)

    Mismo flujo que tu 2D:
    - Calcula un "punto de intersección" de direcciones extendidas (en 3D: el medio de los más cercanos).
    - Comprueba si cae dentro de cada segmento por distancia al centro (di,dj).
    - Maneja paralelas con proyecciones y recortes usando `pointInsideRod3d`.

    Devuelve: (x1Aux,y1Aux,z1Aux, x2Aux,y2Aux,z2Aux)
    """

    function rodIntersection3d(x1,y1,z1,l1,theta1,phi1,
                            x2,y2,z2,l2,theta2,phi2; separation=0.99, tol=1e-6)

        # por compatibilidad con tu patrón
        x1Aux = x1; y1Aux = y1; z1Aux = z1
        x2Aux = x2; y2Aux = y2; z2Aux = z2

        # 1) "intersección" (medio de puntos más cercanos o medio de centros si paralelas)
        pxAux, pyAux, pzAux = intersection2lines3d(x1,y1,z1,theta1,phi1, x2,y2,z2,theta2,phi2; tol=tol)

        # 2) distancias desde cada centro a ese punto
        di = sqrt((x1 - pxAux)^2 + (y1 - pyAux)^2 + (z1 - pzAux)^2)
        dj = sqrt((x2 - pxAux)^2 + (y2 - pyAux)^2 + (z2 - pzAux)^2)

        # "normAux" 3D análogo a determinante 2D: si direcciones casi paralelas -> 0
        u1 = dir_from_angles(theta1, phi1)
        u2 = dir_from_angles(theta2, phi2)
        parallelish = (1 - abs(dot3(u1,u2))) < 1e-7  # |u1·u2|≈1

        if parallelish
            # proyecta pxAux en cada eje y recorta dentro del segmento
            x1Aux,y1Aux,z1Aux = point2line3d(pxAux,pyAux,pzAux, x1,y1,z1, theta1,phi1)
            x1Aux,y1Aux,z1Aux = pointInsideRod3d(x1,y1,z1, l1, theta1,phi1, x1Aux,y1Aux,z1Aux, separation)

            x2Aux,y2Aux,z2Aux = point2line3d(pxAux,pyAux,pzAux, x2,y2,z2, theta2,phi2)
            x2Aux,y2Aux,z2Aux = pointInsideRod3d(x2,y2,z2, l2, theta2,phi2, x2Aux,y2Aux,z2Aux, separation)

        elseif di <= l1/2 && dj <= l2/2
            # el punto medio de las rectas cae "dentro" de ambos rods por distancia a centros
            x1Aux,y1Aux,z1Aux = pointInsideRod3d(x1,y1,z1, l1, theta1,phi1, pxAux,pyAux,pzAux, separation)
            x2Aux,y2Aux,z2Aux = pointInsideRod3d(x2,y2,z2, l2, theta2,phi2, pxAux,pyAux,pzAux, separation)

        elseif di <= l1/2 && dj > l2/2
            x2Aux,y2Aux,z2Aux = pointInsideRod3d(x2,y2,z2, l2, theta2,phi2, pxAux,pyAux,pzAux, separation)
            x1Aux,y1Aux,z1Aux = point2line3d(x2Aux,y2Aux,z2Aux, x1,y1,z1, theta1,phi1)

        elseif di > l1/2 && dj <= l2/2
            x1Aux,y1Aux,z1Aux = pointInsideRod3d(x1,y1,z1, l1, theta1,phi1, pxAux,pyAux,pzAux, separation)
            x2Aux,y2Aux,z2Aux = point2line3d(x1Aux,y1Aux,z1Aux, x2,y2,z2, theta2,phi2)

        else
            # ambos fuera: sigue tu lógica de "intenta con uno, comprueba el otro"
            x2Aux,y2Aux,z2Aux = pointInsideRod3d(x2,y2,z2, l2, theta2,phi2, pxAux,pyAux,pzAux, separation)
            x1Aux,y1Aux,z1Aux = point2line3d(x2Aux,y2Aux,z2Aux, x1,y1,z1, theta1,phi1)

            dj2 = sqrt((x1 - x1Aux)^2 + (y1 - y1Aux)^2 + (z1 - z1Aux)^2)
            if dj2 > l1/2
                x1Aux,y1Aux,z1Aux = pointInsideRod3d(x1,y1,z1, l1, theta1,phi1, pxAux,pyAux,pzAux, separation)
                x2Aux_,y2Aux_,z2Aux_ = point2line3d(x1Aux,y1Aux,z1Aux, x2,y2,z2, theta2,phi2)

                dj3 = sqrt((x2 - x2Aux_)^2 + (y2 - y2Aux_)^2 + (z2 - z2Aux_)^2)
                if dj3 < l2/2
                    x2Aux,y2Aux,z2Aux = x2Aux_,y2Aux_,z2Aux_
                end
            end
        end

        return x1Aux,y1Aux,z1Aux, x2Aux,y2Aux,z2Aux
    end

end
