"""
    function saveRAM!(community::Community)

Function that stores the present configuration of community in the field `Community.pastTimes`.
"""
function saveRAM!(community::Community)
    
    com = Community()

    # Transform to the correct platform the parameters
    for (sym,prop) in pairs(community.abm.parameters)
        p = community.parameters[sym]
        if prop.scope == :agent
            com.parameters[sym] = copy(Array(p[1:community.N]))
        elseif prop.scope in [:model,:medium]
            com.parameters[sym] = copy(Array(p))
        end
    end

    for i in [:dt,:t,:N]
        setfield!(com,i,copy(getfield(community,i)))
    end
    #id
    setfield!(com,:id,copy(Array{Int64}(getfield(community,:id))[1:community.N]))
    setfield!(com,:NMedium,copy(Array{Int64}(getfield(community,:NMedium))))
    setfield!(com,:simBox,copy(Array{Float64}(getfield(community,:simBox))))
    setfield!(com,:abm,community.abm)

    push!(community.pastTimes,com)

    return

end

"""
    saveJLD2(file::String, community::Community; overwrite=false)

Save in `file` the current instance of the `community`. If some other community was being saved here, an error will raise. If `overwrite=true` is specified, it will remove the previous community.   
"""
# function saveJLD2(file::String, community::Community; overwrite=false)

#     #Check if file is open and assotiated with our community
#     if file in keys(SAVING) && isfile(file)
#         if SAVING[file].uuid == community.uuid.value
#             if !haskey(JLD2.OPEN_FILES,SAVING[file].file.path)
#                 SAVING[file].file = jldopen(file, "a+")
#             end
#         elseif !overwrite
#             error("File $file contains other community information inside it. If you want to overwrite it, set the key argument overwrite to true.")
#         else
#             if haskey(JLD2.OPEN_FILES,SAVING[file].file.path)
#                 close(SAVING[file].file)
#             end
#             f = jldopen(file, "w")   
#             SAVING[file] = SavingFile(community.uuid.value,f)
#         end
#     elseif isfile(file)
#         f = jldopen(file, "w")   
#         if "uuid" in keys(f)
#             if f["uuid"] == community.uuid.value
#                 close(file, "w")   
#                 f = jldopen(file, "a+")   
#                 SAVING[file] = SavingFile(community.uuid.value,f)
#             elseif !overwrite
#                 error("File $file contains other community information inside it. If you want to overwrite it, set the key argument overwrite to true.")
#             else
#                 SAVING[file] = SavingFile(community.uuid.value,f)
#             end
#         else
#             SAVING[file] = SavingFile(community.uuid.value,f)
#         end
#     else
#         f = jldopen(file, "w") 
#         SAVING[file] = SavingFile(community.uuid.value,f)    
#     end

#     f = SAVING[file].file

#     if !haskey(f,"uuid")
#         f["uuid"] = community.uuid
#     end

#     if !( "abm" in keys(f) )
#         f["abm/dims"] = community.abm.dims
#         f["abm/parameters"] = community.abm.parameters
#         f["abm/declaredUpdates"] = community.abm.declaredUpdates
#         f["abm/removalOfAgents_"] = community.abm.removalOfAgents_
#         f["abm/agentAlg/alg"] = community.abm.agentAlg
#         f["abm/agentAlg/args"] = community.abm.agentSolveArgs
#         f["abm/modelAlg/alg"] = community.abm.modelAlg
#         f["abm/modelAlg/args"] = community.abm.modelSolveArgs
#         f["abm/mediumAlg/alg"] = community.abm.mediumAlg
#         f["abm/mediumAlg/args"] = community.abm.mediumSolveArgs
#         f["abm/platform"] = community.abm.platform
#     end

#     t = 1
#     if "times" in keys(f)
#         t = length(f["times"]) + 1
#     end
#     # Transform to the correct platform the parameters
#     for (sym,prop) in pairs(community.abm.parameters)
#         p = community.parameters[sym]
#         if prop.scope == :agent
#             f["times/$t/parameters/$sym"] = copy(Array(p[1:community.N]))
#         elseif prop.scope in [:model,:medium]
#             f["times/$t/parameters/$sym"] = copy(Array(p))
#         end
#     end

#     for sym in [:N,:t,:dt]
#         f["times/$t/$sym"] = copy(getfield(community,sym))
#     end
#     for sym in [:NMedium,:simBox]
#         f["times/$t/$sym"] = copy(Array(getfield(community,sym)))
#     end
#     f["times/$t/id"] = copy(Array{Int64}(getfield(community,:id))[1:community.N])

#     return

# end



using JLD2

function saveJLD2(file::String, community::Community; overwrite=false, include_past::Bool=true)
    # --- abrir sin truncar salvo overwrite ---
    mode = (isfile(file) && !overwrite) ? "r+" : "w"
    f = jldopen(file, mode)

    try
        # uuid
        if haskey(f, "uuid")
            stored_uuid = f["uuid"]
            if stored_uuid != community.uuid.value
                overwrite || error("UUID distinta. Usa overwrite=true para sobrescribir.")
                close(f); f = jldopen(file, "w"); f["uuid"] = community.uuid
            end
        else
            f["uuid"] = community.uuid
        end

        # metadatos ABM una sola vez
        if !haskey(f, "abm")
            f["abm/dims"]             = community.abm.dims
            f["abm/parameters"]       = community.abm.parameters
            f["abm/declaredUpdates"]  = community.abm.declaredUpdates
            f["abm/removalOfAgents_"] = community.abm.removalOfAgents_
            f["abm/agentAlg/alg"]     = community.abm.agentAlg
            f["abm/agentAlg/args"]    = community.abm.agentSolveArgs
            f["abm/modelAlg/alg"]     = community.abm.modelAlg
            f["abm/modelAlg/args"]    = community.abm.modelSolveArgs
            f["abm/mediumAlg/alg"]    = community.abm.mediumAlg
            f["abm/mediumAlg/args"]   = community.abm.mediumSolveArgs
            f["abm/platform"]         = community.abm.platform
        end

        # grupo times
        if !haskey(f, "times"); JLD2.Group(f, "times"); end

        # siguiente índice base = maxclave+1
        keys_times = collect(keys(f["times"]))
        num_keys   = [parse(Int,k) for k in keys_times if !isnothing(tryparse(Int,k))]
        base = isempty(num_keys) ? 1 : maximum(num_keys)+1

        # helper para bajar a CPU
        to_cpu(x) = x
        to_cpu(x::AbstractArray) = Array(x)

        # construye la lista que vamos a guardar
        hist = include_past ? vcat(community.pastTimes, [community]) : [community]

        # guarda TODO el historial en una sola llamada
        for (i, com) in enumerate(hist)
            key = "times/$(base + i - 1)"
            JLD2.Group(f, key)
            JLD2.Group(f, "$key/parameters")

            # parámetros por scope
            for (sym, prop) in pairs(com.abm.parameters)
                p = com.parameters[sym]
                if prop.scope == :agent
                    f["$key/parameters/$sym"] = to_cpu(p[1:com.N])
                elseif prop.scope in (:model, :medium)
                    f["$key/parameters/$sym"] = to_cpu(p)
                end
            end

            # escalares y arrays básicos de ese snapshot
            f["$key/N"]   = com.N
            f["$key/t"]   = com.t
            f["$key/dt"]  = com.dt
            f["$key/id"]  = Array{Int64}(getfield(com,:id))[1:com.N]
            f["$key/NMedium"] = to_cpu(getfield(com,:NMedium))
            f["$key/simBox"]  = to_cpu(getfield(com,:simBox))
        end

        return base:base + length(hist) - 1  # rango de claves escritas
    finally
        close(f)
    end
end



"""
    function loadJLD2(file::String)

Load the Community structure saved in file.
"""
# function loadJLD2(file::String)

#     if file in keys(SAVING) && isfile(file)
#         close(SAVING[file].file)
#         delete!(SAVING,file)
#     end

#     jldopen(file, "r") do f

#         #Agent
#         abm = ABM()
#         abm.dims = f["abm/dims"] 
#         abm.parameters = f["abm/parameters"] 
#         abm.declaredUpdates = f["abm/declaredUpdates"] 
#         abm.removalOfAgents_ = f["abm/removalOfAgents_"] 

#         #Assign abm
#         t = length(f["times"])
#         community = Community(abm,
#             id = f["times/$t/id"],
#             N = f["times/$t/N"],
#             NMedium = f["times/$t/NMedium"],
#             t = f["times/$t/t"],
#             dt = f["times/$t/dt"],
#             simBox = f["times/$t/simBox"],
#             platform = f["abm/platform"],#eval(Meta.parse("CellBasedModels.$(f["platform/platform"])()")),
#             agentAlg = f["abm/agentAlg/alg"],#eval(Meta.parse("CellBasedModels.$(f["agentAlg/alg"])()")),
#             agentSolveArgs = f["abm/agentAlg/args"],
#             modelAlg = f["abm/modelAlg/alg"],#eval(Meta.parse("CellBasedModels.$(f["modelAlg/alg"])()")),
#             modelSolveArgs = f["abm/modelAlg/args"],
#             mediumAlg = f["abm/mediumAlg/alg"],#eval(Meta.parse("CellBasedModels.$(f["mediumAlg/alg"])()")),
#             mediumSolveArgs = f["abm/mediumAlg/args"],
#         )
#         setfield!(community,:uuid,f["uuid"])
#         for (sym,prop) in pairs(community.abm.parameters)
#             community.parameters[sym] = f["times/$t/parameters/$sym"]
#         end

#         #Base parameters
#         times = sort([Meta.parse(i) for i in keys(f["times"])])[1:end-1]
#         for t in times
#             com = Community(abm,
#                 id = f["times/$t/id"],
#                 N = f["times/$t/N"],
#                 NMedium = f["times/$t/NMedium"],
#                 t = f["times/$t/t"],
#                 dt = f["times/$t/dt"],
#                 simBox = f["times/$t/simBox"],
#                 platform = f["platform/platform"],#eval(Meta.parse("CellBasedModels.$(f["platform/platform"])()")),
#                 agentAlg = f["agentAlg/alg"],#eval(Meta.parse("CellBasedModels.$(f["agentAlg/alg"])()")),
#                 agentSolveArgs = f["agentAlg/args"],
#                 modelAlg = f["modelAlg/alg"],#eval(Meta.parse("CellBasedModels.$(f["modelAlg/alg"])()")),
#                 modelSolveArgs = f["modelAlg/args"],
#                 mediumAlg = f["mediumAlg/alg"],#eval(Meta.parse("CellBasedModels.$(f["mediumAlg/alg"])()")),
#                 mediumSolveArgs = f["mediumAlg/args"],
#             )
#             for (sym,prop) in pairs(community.abm.parameters)
#                 com.parameters[sym] = f["times/$t/parameters/$sym"]
#             end
    
#             push!(community.pastTimes, com)
#         end

#         setfield!(community,:loaded,false)

#         return community
#     end

# end


# function loadJLD2(file::String)
#     if file in keys(SAVING) && isfile(file)
#         close(SAVING[file].file)
#         delete!(SAVING, file)
#     end

#     return jldopen(file, "r") do f
#         # 1) Reconstruir ABM
#         abm = ABM()
#         abm.dims              = f["abm/dims"]
#         abm.parameters        = f["abm/parameters"]
#         abm.declaredUpdates   = f["abm/declaredUpdates"]
#         abm.removalOfAgents_  = f["abm/removalOfAgents_"]

#         # --- Compatibilidad hacia atrás: si no existen en abm/, busca legacy ---
#         has_abm_agent  = haskey(f, "abm/agentAlg/alg")
#         has_abm_model  = haskey(f, "abm/modelAlg/alg")
#         has_abm_medium = haskey(f, "abm/mediumAlg/alg")
#         has_abm_plat   = haskey(f, "abm/platform")

#         abm.agentAlg        = has_abm_agent  ? f["abm/agentAlg/alg"]    : get(f, "agentAlg/alg", abm.agentAlg)
#         abm.agentSolveArgs  = has_abm_agent  ? f["abm/agentAlg/args"]   : get(f, "agentAlg/args", Dict{Symbol,Any}())
#         abm.modelAlg        = has_abm_model  ? f["abm/modelAlg/alg"]    : get(f, "modelAlg/alg", abm.modelAlg)
#         abm.modelSolveArgs  = has_abm_model  ? f["abm/modelAlg/args"]   : get(f, "modelAlg/args", Dict{Symbol,Any}())
#         abm.mediumAlg       = has_abm_medium ? f["abm/mediumAlg/alg"]   : get(f, "mediumAlg/alg", abm.mediumAlg)
#         abm.mediumSolveArgs = has_abm_medium ? f["abm/mediumAlg/args"]  : get(f, "mediumAlg/args", Dict{Symbol,Any}())
#         abm.platform        = has_abm_plat   ? f["abm/platform"]         : get(f, "platform", abm.platform)

#         # 2) Último tiempo
#         tlast = length(f["times"])
#         community = Community(
#             abm;
#             id      = f["times/$tlast/id"],
#             N       = f["times/$tlast/N"],
#             NMedium = f["times/$tlast/NMedium"],
#             t       = f["times/$tlast/t"],
#             dt      = f["times/$tlast/dt"],
#             simBox  = f["times/$tlast/simBox"],
#         )

#         # UUID (guardado como valor escalar)
#         setfield!(community, :uuid, f["uuid"])
#         for (sym, prop) in pairs(community.abm.parameters)
#             community.parameters[sym] = f["times/$tlast/parameters/$sym"]
#         end

#         # 3) Reconstruir tiempos pasados
#         times_sym = sort([Meta.parse(i) for i in keys(f["times"])])[1:end-1]
#         for tt in times_sym
#             com = Community(
#                 abm;
#                 id      = f["times/$tt/id"],
#                 N       = f["times/$tt/N"],
#                 NMedium = f["times/$tt/NMedium"],
#                 t       = f["times/$tt/t"],
#                 dt      = f["times/$tt/dt"],
#                 simBox  = f["times/$tt/simBox"],
#             )
#             for (sym, prop) in pairs(community.abm.parameters)
#                 com.parameters[sym] = f["times/$tt/parameters/$sym"]
#             end
#             push!(community.pastTimes, com)
#         end

#         setfield!(community, :loaded, false)
#         return community
#     end
# end

function loadJLD2(file::String)
    if file in keys(SAVING) && isfile(file)
        close(SAVING[file].file)
        delete!(SAVING, file)
    end

    return jldopen(file, "r") do f
        # 1) Reconstruir ABM
        abm = ABM()
        abm.dims              = f["abm/dims"]
        abm.parameters        = f["abm/parameters"]
        abm.declaredUpdates   = f["abm/declaredUpdates"]
        abm.removalOfAgents_  = f["abm/removalOfAgents_"]

        has_abm_agent  = haskey(f, "abm/agentAlg/alg")
        has_abm_model  = haskey(f, "abm/modelAlg/alg")
        has_abm_medium = haskey(f, "abm/mediumAlg/alg")
        has_abm_plat   = haskey(f, "abm/platform")

        abm.agentAlg        = has_abm_agent  ? f["abm/agentAlg/alg"]   : get(f, "agentAlg/alg", abm.agentAlg)
        abm.agentSolveArgs  = has_abm_agent  ? f["abm/agentAlg/args"]  : get(f, "agentAlg/args", Dict{Symbol,Any}())
        abm.modelAlg        = has_abm_model  ? f["abm/modelAlg/alg"]   : get(f, "modelAlg/alg", abm.modelAlg)
        abm.modelSolveArgs  = has_abm_model  ? f["abm/modelAlg/args"]  : get(f, "modelAlg/args", Dict{Symbol,Any}())
        abm.mediumAlg       = has_abm_medium ? f["abm/mediumAlg/alg"]  : get(f, "mediumAlg/alg", abm.mediumAlg)
        abm.mediumSolveArgs = has_abm_medium ? f["abm/mediumAlg/args"] : get(f, "mediumAlg/args", Dict{Symbol,Any}())
        abm.platform        = has_abm_plat   ? f["abm/platform"]       : get(f, "platform", abm.platform)

        # 2) Claves de tiempos (robusto)
        times_group = f["times"]
        raw_keys = collect(keys(times_group))              # => Vector{String}
        time_keys = sort([k for k in raw_keys if !isnothing(tryparse(Int, k))];
                         by = x -> parse(Int, x))
        @assert !isempty(time_keys) "No se encontraron tiempos en el archivo"

        last_key = last(time_keys)                         # p.ej. "120"
        # 3) Construir community en el último tiempo real
        community = Community(
            abm;
            id      = f["times/$last_key/id"],
            N       = f["times/$last_key/N"],
            NMedium = f["times/$last_key/NMedium"],
            t       = f["times/$last_key/t"],
            dt      = f["times/$last_key/dt"],
            simBox  = f["times/$last_key/simBox"],
        )

        # UUID (si se guardó como escalar)
        if haskey(f, "uuid")
            setfield!(community, :uuid, f["uuid"])
        end

        # Cargar parámetros del último tiempo
        for (sym, _) in pairs(community.abm.parameters)
            path = "times/$last_key/parameters/$(String(sym))"
            if haskey(f, path)
                community.parameters[sym] = f[path]
            end
        end

        # 4) Reconstruir TODOS los tiempos pasados (en orden ascendente)
        past_keys = time_keys[1:end-1]
        for key in past_keys
            com = Community(
                abm;
                id      = f["times/$key/id"],
                N       = f["times/$key/N"],
                NMedium = f["times/$key/NMedium"],
                t       = f["times/$key/t"],
                dt      = f["times/$key/dt"],
                simBox  = f["times/$key/simBox"],
            )
            for (sym, _) in pairs(community.abm.parameters)
                path = "times/$key/parameters/$(String(sym))"
                if haskey(f, path)
                    com.parameters[sym] = f[path]
                end
            end
            push!(community.pastTimes, com)
        end

        setfield!(community, :loaded, false)
        return community
    end
end
