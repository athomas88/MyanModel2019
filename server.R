# server.r
#adding text
shinyServer(function(input, output) {
  values <- reactiveValues(
    forwardProjectionOutputs  = list(),
    combinedAggregatedOutput = data.frame(),
    speciesAggregatedOutput = data.frame(),
    params_belize = params_belize_default,
    params_temp = params_belize_default,
    virginTidy = data.frame(),
    sizeSpectrum = data.frame()
  )

  output$customControlF <- renderUI({
    if (2 %in% input$interventionGroup) {
      sliderInput(inputId = "custF", label = "What multiplier of natural mortality should F be? (F = M at 1)", 0, 2, 1, step = 0.1)
    }
  })
  output$customControlFsize <- renderUI({
    if (4 %in% input$interventionGroup) {
      sliderInput(inputId = "custFsize", label = "(Size and F limits) What multiplier of natural mortality should F be? (F = M at 1)", 0, 2, 1, step = 0.1)
    }
  })
  output$customControlNTZ <- renderUI({
    if (5 %in% input$interventionGroup) {
      sliderInput(inputId = "ntzSize", label = "What percentage of the fishery should the no-take zone protect? /ငါးမဖမ်းရဇုံ ရာခိုင်နှုန်းအဖုံး", 0, 100, 20, step = 1)
    }
  })
  output$customControlNTseason <- renderUI({
    if (6 %in% input$interventionGroup) {
      sliderTextInput(inputId = "ntSeason", label = "What months should be in the closed season?", choices = month.name, selected = month.name[c(6,8)])
    }
  })
  output$customControlNTF <- renderUI({
    if (6 %in% input$interventionGroup) {
      sliderInput(inputId = "ntF", label = "What percentage of fishing will be allowed during the closed season?", 0, 100, 20, step = 5)
    }
  })
  output$speciesSelect <- renderUI({
    if (input$siteAnalysis == "Rakhine"){
      selectInput(inputId = "sppSelect",label="Select species to plot / စံပြုပုံစံအကွက်တွင် ငါးမျိုးစိတ်များအားရွေးချယ်ရန်",choices=sort(params_belize_default$sciname[ which (params_belize_default$region == "Rakhine")]))
    }
    else if (input$siteAnalysis == "Delta"){
      selectInput(inputId = "sppSelect",label="Select species to plot / စံပြုပုံစံအကွက်တွင် ငါးမျိုးစိတ်များအားရွေးချယ်ရန်",choices=sort(params_belize_default$sciname[ which (params_belize_default$region == "Delta")]))
    }
    else if (input$siteAnalysis == "Tanintharyi"){
      selectInput(inputId = "sppSelect",label="Select species to plot / စံပြုပုံစံအကွက်တွင် ငါးမျိုးစိတ်များအားရွေးချယ်ရန်",choices=sort(params_belize_default$sciname[ which (params_belize_default$region == "Tanintharyi")]))
    } 
  })
  
  #If user clicks "Run Simulation" button, do:
  observeEvent(input$runSim, {
    values$params_belize<-values$params_temp
    progress <- shiny::Progress$new()
    progress$set(message = "Running simulation:", value = 0.01)
    on.exit(progress$close())
    progressCounter <- 0
    progressTotal <- 2 + length(input$interventionGroup)
    if (input$siteAnalysis == "Rakhine"){
      params_belize_analysis <- values$params_belize %>%
        filter(region == "Rakhine")
    }
    if (input$siteAnalysis == "Delta"){
      params_belize_analysis <- values$params_belize %>%
        filter(region == "Delta")
    }
    if (input$siteAnalysis == "Tanintharyi"){
      params_belize_analysis <- values$params_belize %>%
        filter(region == "Tanintharyi")
    }
    nsS <- nrow(params_belize_analysis)
    inter_belize <- matrix(1,nrow=nrow(params_belize_analysis),ncol=nrow(params_belize_analysis))
    rownames(inter_belize) <- params_belize_analysis$species
    colnames(inter_belize) <- params_belize_analysis$species
    values$params_belize<-params_belize_analysis
    
    MizerParamsBelize <- MizerParams(params_belize_analysis,no_w=nwS,no_w_pp=npS,interaction = inter_belize,kappa = kappa)
    MizerParamsBelize@species_params$M <- params_belize_analysis$M
    tRes <- 0.1
    fisheryAge <- 50
    t0 <-tNow - fisheryAge
    tVirgin <- t0 - 50
    tManagement <- input$managementStart
    tEnd <- tManagement + input$tProjection
    ## Define virgin unfished period
    virginInput <- list(managementName = "Virgin Unfished",
                        interactions = inter_belize,
                        params_data = params_belize_analysis,
                        sel_func = "knife_edge",
                        knife_edge_size = 0,
                        time = seq(tVirgin,t0-1),
                        e0 = rep(0,length(MizerParamsBelize@species_params$M)),
                        eDelta = 0,
                        initN = get_initial_n(MizerParamsBelize),
                        initNPP = MizerParamsBelize@cc_pp,
                        nws = nwS,
                        nps = npS,
                        nss = nsS,
                        res = tRes)
    progress$set(detail = virginInput$managementName)
    virginOutput <- tidySim(virginInput)
    virginN <- virginOutput$sim@n[dim(virginOutput$sim@n)[1],,]
    virginNPP <- virginOutput$sim@n_pp[dim(virginOutput$sim@n_pp)[1],]
    progressCounter <- progressCounter + 1
    #progress$set(value = progressCounter, detail = paste(round(progressCounter*100/progressTotal,0),"% done.",sep=""))
    ## Define historic fishing period
    historicInput <- list(managementName = "Historic Fishing",
                          interactions = inter_belize,
                          params_data = params_belize_analysis,
                          sel_func = "knife_edge",
                          knife_edge_size = 0,
                          time = seq(t0-1,tNow-1),
                          e0 = params_belize_analysis$FvM/fisheryAge,
                          eDelta = (fisheryAge)^(1/fisheryAge)-1,
                          initN = virginN,
                          initNPP = virginNPP,
                          nws = nwS,
                          nps = npS,
                          nss = nsS,
                          res = tRes)
    progress$set(value = progressCounter/progressTotal,detail = historicInput$managementName)
    historicOutput <- tidySim(historicInput)
    baselineN <- historicOutput$sim@n[dim(historicOutput$sim@n)[1],,]
    baselineNPP <- historicOutput$sim@n_pp[dim(historicOutput$sim@n_pp)[1],]
    progressCounter <- progressCounter + 1
    #progressCounter <- progressCounter + 1
    #progress$set(value = progressCounter, detail = paste(round(progressCounter*100/progressTotal,0),"% done.",sep=""))
    if (6 %in% input$interventionGroup){
      m1 <- match(input$ntSeason[1],month.name)
      m2 <- match(input$ntSeason[2],month.name)
      openF <- 1 - sum(monthlyCatch$CatchProp[m1:m2])
      closeF <- sum(monthlyCatch$CatchProp[m1:m2]) * (input$ntF/100)
    } else {
      closeF <- vector()
      openF <- vector()
    }
    comp <- rep(input$compliance/100, nsS)
    eDelta <- input$effortCreep/100
    forwardProjectionInputs <- list(
      statusQuo = list(managementName = "No Change",
                       interactions = inter_belize,
                       params_data = params_belize_analysis,
                       sel_func = historicInput$sel_func,
                       knife_edge_size = historicInput$knife_edge_size,
                       time = seq(tNow-1,tEnd),
                       e0 = params_belize_analysis$FvM,
                       eDelta = eDelta,
                       initN = baselineN,
                       initNPP = baselineNPP,
                       nws = nwS,
                       nps = npS,
                       nss = nsS,
                       res = tRes),
      FInput = list(managementName = "Limit catch relative to natural mortality",
                    interactions = inter_belize,
                    params_data = params_belize_analysis,
                    sel_func = "knife_edge",
                    knife_edge_size = historicInput$knife_edge_size,
                    time = seq(tManagement-1,tEnd),
                    e0 = (params_belize_analysis$M * input$custF) * comp + (params_belize_analysis$FvM * (1-comp)),
                    eDelta = 0,
                    initN = baselineN,
                    initNPP = baselineNPP,
                    nws = nwS,
                    nps = npS,
                    nss = nsS,
                    res = tRes),
      sizeLimitInput = list(managementName = "Minimum size Limit (for each species)",
                            interactions = inter_belize,
                            params_data = params_belize_analysis,
                            sel_func = "knife_edge",
                            knife_edge_size = (params_belize_analysis$w_mat * comp) + (historicInput$knife_edge_size * (1-comp)),
                            time = seq(tManagement-1,tEnd),
                            e0 = params_belize_analysis$FvM,
                            eDelta = eDelta,
                            initN = baselineN,
                            initNPP = baselineNPP,
                            nws = nwS,
                            nps = npS,
                            nss = nsS,
                            res = tRes),
      sizeAndFInput = list(managementName = "Size limits and catch limits",
                           interactions = inter_belize,
                           params_data = params_belize_analysis,
                           sel_func = "knife_edge",
                           knife_edge_size = (params_belize_analysis$w_mat * comp) + (historicInput$knife_edge_size * (1-comp)),
                           time = seq(tManagement-1,tEnd),
                           e0 = (params_belize_analysis$M * input$custFsize * comp) + (params_belize_analysis$FvM * (1-comp)),
                           eDelta = 0,
                           initN = baselineN,
                           initNPP = baselineNPP,
                           nws = nwS,
                           nps = npS,
                           nss = nsS,
                           res = tRes),
      ntz = list(managementName = "No-take zone",
                 interactions = inter_belize,
                 params_data = params_belize_analysis,
                 sel_func = historicInput$sel_func,
                 knife_edge_size = historicInput$knife_edge_size,
                 time = seq(tManagement-1,tEnd),
                 e0 = (params_belize_analysis$FvM*(1 - input$ntzSize/100) * comp) + (params_belize_analysis$FvM * (1-comp)),
                 eDelta = eDelta,
                 initN = baselineN,
                 initNPP = baselineNPP,
                 nws = nwS,
                 nps = npS,
                 nss = nsS,
                 res = tRes),
      seasonalClosure = list(managementName = "Seasonal closure",
                             interactions = inter_belize,
                             params_data = params_belize_analysis,
                             sel_func = historicInput$sel_func,
                             knife_edge_size = historicInput$knife_edge_size,
                             time = seq(tManagement-1,tEnd),
                             e0 = (params_belize_analysis$FvM*openF + params_belize_analysis$FvM*closeF) * comp + (params_belize_analysis$FvM * (1-comp)),
                             eDelta = eDelta,
                             initN = baselineN,
                             initNPP = baselineNPP,
                             nws = nwS,
                             nps = npS,
                             nss = nsS,
                             res = tRes))
    
    forwardProjectionInputs <- forwardProjectionInputs[as.numeric(input$interventionGroup)]
    values$forwardProjectionOutputs <- map(forwardProjectionInputs,function(x){
      progress$set(value = progressCounter/progressTotal,detail = x$managementName)
      progressCounter <<- progressCounter + 1
      tidySimCC(x)
    })
    
    values$combinedSpeciesOutput <- values$forwardProjectionOutputs %>%
      map("simResultsSpecies") %>% 
      bind_rows() %>%
      rbind(historicOutput$simResultsSpecies,
            virginOutput$simResultsSpecies) %>%
      filter(!is.na(Year))
    
    values$combinedAggregatedOutput <- values$forwardProjectionOutputs %>%
      map("simResultsAggregated") %>% 
      bind_rows() %>%
      rbind(historicOutput$simResultsAggregated,
            virginOutput$simResultsAggregated) %>%
      filter(!is.na(Year))
    
    speciesVirgin <- rownames(virginN)
    values$virginTidy <- virginN %>% as_data_frame()
    values$virginTidy$Species <- speciesVirgin
    values$virginTidy <- values$virginTidy %>%
      gather(weightClass,virginAbundance,-Species)
    values$virginTidy$weightClass <- as.numeric(values$virginTidy$weightClass)
    
    simS<-values$forwardProjectionOutputs %>% map("sim")
    timeStep <- input$tProjection +1
    simNames <- names(simS)
    mapNames <- map(forwardProjectionInputs,function(x){x$managementName}) %>% as.vector()
    values$sizeSpectrum <- map(simNames,function(x){
      temp <- get(x,simS)
      mName <- get(x,mapNames)
      temp <- temp@n[timeStep,,]
      speciesName <- rownames(temp)
      temp <- temp %>% as_data_frame()
      temp$Species <- speciesName
      temp <- temp %>%
        gather(weightClass,n,-Species)
      temp$weightClass <- as.numeric(temp$weightClass)
      temp$Management <- mName
      temp <- temp %>%
        mutate(Biomass = weightClass * n) %>%
        filter(n>0) %>%
        left_join(values$virginTidy) %>%
        mutate(relativeAbundance = n/virginAbundance)
    }) %>%
      bind_rows()
  })

  output$PlotProjectionsAggregated <- renderPlot({
    if(nrow(req(values$combinedAggregatedOutput)) > 0) {
      bmin <- min(values$combinedAggregatedOutput$Biomass)
      bmax <- max(values$combinedAggregatedOutput$Biomass)
      ymin <- min(values$combinedAggregatedOutput$Yield)
      ymax <- max(values$combinedAggregatedOutput$Yield)
      t0 <-tNow - fisheryAge
      tVirgin <- t0 - 100
      tManagement <- input$managementStart
      tEnd <- tManagement + input$tProjection
      bPlot <- values$combinedAggregatedOutput %>%
        mutate(Biomass = (Biomass - bmin) / (bmax - bmin)) %>% ## Convert from grams to kg
        filter(Year >= tVirgin) %>%
        ggplot(aes(x=Year, y=Biomass,color=Management,group=Management)) +
        scale_colour_manual(values=c("Virgin Unfished"="black", 
                                    "Historic Fishing"="blue",
                                    "No Change"="red",
                                    "Limit catch relative to natural mortality" = "#009E73",
                                    "Minimum size Limit (for each species)"= "slateblue",
                                    "Size limits and catch limits"= "yellow",
                                    "No-take zone" = "orange",
                                    "Seasonal closure" = "pink")) +
        annotate("rect", xmin = -Inf, xmax = t0-1, ymin = -Inf, ymax = Inf, fill = "lightskyblue", alpha = 0.2) +
        annotate("rect", xmin = t0-1, xmax = tNow - 1, ymin = -Inf, ymax = Inf, fill = "deepskyblue", alpha = 0.2) +
        geom_line(size=1) +
        geom_vline(xintercept = tManagement-1) +
        xlab("Year") +
        ylab("Relative Biomass in Water") +
        ggtitle("Projected relative biomass")  #+
        #theme(axis.text.y=element_blank())
      
      yPlot <- values$combinedAggregatedOutput %>%
        filter(Year >= tVirgin) %>%
        mutate(Yield = (Yield - ymin) / (ymax - ymin)) %>% ## Convert from grams to kg
        ggplot(aes(x=Year, y=Yield,color=Management,group=Management)) +
        scale_colour_manual(values=c("Virgin Unfished"="black", 
                                     "Historic Fishing"="blue",
                                     "No Change"="red",
                                     "Limit catch relative to natural mortality" = "#009E73",
                                     "Minimum size Limit (for each species)"= "slateblue",
                                     "Size limits and catch limits"= "yellow",
                                     "No-take zone" = "orange",
                                     "Seasonal closure" = "pink")) +
        annotate("rect", xmin = -Inf, xmax = t0-1, ymin = -Inf, ymax = Inf, fill = "lightskyblue", alpha = 0.2) +
        annotate("rect", xmin = t0-1, xmax = tNow - 1, ymin = -Inf, ymax = Inf, fill = "deepskyblue", alpha = 0.2) +
        geom_line(size=1) +
        geom_vline(xintercept = tManagement-1) +
        xlab("Year") +
        ylab("Relative Yield") +
        ggtitle("Projected relative yield") #+
        #theme(axis.text.y=element_blank())
      
      tradeoffPlot <- values$combinedAggregatedOutput %>%
        filter(Management != "Historic Fishing" & Management != "Virgin Unfished") %>%
        filter(Year == tEnd) %>%
        ggplot(aes(x=((Biomass - min(values$combinedAggregatedOutput$Biomass)) / (max(values$combinedAggregatedOutput$Biomass) - min(values$combinedAggregatedOutput$Biomass))),
                   y=((Yield - min(values$combinedAggregatedOutput$Yield)) / (max(values$combinedAggregatedOutput$Yield) - min(values$combinedAggregatedOutput$Yield))),
                   color=Management)) +
        scale_colour_manual(values=c("Virgin Unfished"="black", 
                                     "Historic Fishing"="blue",
                                     "No Change"="red",
                                     "Limit catch relative to natural mortality" = "#009E73",
                                     "Minimum size Limit (for each species)"= "slateblue",
                                     "Size limits and catch limits"= "yellow",
                                     "No-take zone" = "orange",
                                     "Seasonal closure" = "pink")) +
        geom_point(size = 5) +
        xlab("Relative Biomass in Water\n[kg]") +
        ylab("Relative Yield\n[kg]") +
        ggtitle(paste("Tradeoff plot of projected biomass and yield in year",tEnd,sep=" ")) #+
        #theme(axis.text.y=element_blank(),
         #     axis.text.x=element_blank())
      
      grid.arrange(tradeoffPlot,
                   bPlot,
                   yPlot,
                   ncol=1)
    }
  })
  output$PlotProjectionsSpecies <- renderPlot({
    if(nrow(req(values$combinedSpeciesOutput)) > 0) {
      t0 <-tNow - fisheryAge
      tVirgin <- t0 - 100
      tManagement <- input$managementStart
      tEnd <- tManagement + input$tProjection
      bPlot<-values$combinedSpeciesOutput %>%
        #filter(Management != "Historic Fishing" & Management != "Virgin Unfished") %>%
        filter(Species == paste(input$sppSelect,"_",input$siteAnalysis,sep="")) %>%
        mutate(Biomass = (Biomass - min(values$combinedSpeciesOutput$Biomass)) / (max(values$combinedSpeciesOutput$Biomass) - min(values$combinedSpeciesOutput$Biomass))) %>% ## Convert to normalized values
        ggplot(aes(x=Year, y=Biomass,color=Management,group=Management)) +
        scale_colour_manual(values=c("Virgin Unfished"="black", 
                                     "Historic Fishing"="blue",
                                     "No Change"="red",
                                     "Limit catch relative to natural mortality" = "#009E73",
                                     "Minimum size Limit (for each species)"= "slateblue",
                                     "Size limits and catch limits"= "yellow",
                                     "No-take zone" = "orange",
                                     "Seasonal closure" = "pink")) +
        #facet_grid(Management~.) +
        geom_line(size=1) +
        geom_vline(xintercept = tManagement-1) +
        ylab("Relative biomass") +
        scale_y_log10() +
        annotate("rect", xmin = -Inf, xmax = t0-1, ymin = 10 ^ -Inf, ymax = 10 ^ Inf, fill = "lightskyblue", alpha = 0.2) +
        annotate("rect", xmin = t0-1, xmax = tNow - 1, ymin = 10 ^ -Inf, ymax = 10 ^ Inf, fill = "deepskyblue", alpha = 0.2) +
        ggtitle("Projected relative biomass") #+
        #theme(axis.text.y=element_blank())
      
      yPlot<-values$combinedSpeciesOutput %>%
        #filter(Management != "Historic Fishing" & Management != "Virgin Unfished") %>%
        filter(Species == paste(input$sppSelect,"_",input$siteAnalysis,sep="")) %>%
        mutate(Yield = (Yield - min(values$combinedSpeciesOutput$Yield)) / (max(values$combinedSpeciesOutput$Yield) - min(values$combinedSpeciesOutput$Yield))) %>% ## Convert from grams to kg
        ggplot(aes(x=Year, y=Yield,color=Management,group=Management)) +
        scale_colour_manual(values=c("Virgin Unfished"="black", 
                                     "Historic Fishing"="blue",
                                     "No Change"="red",
                                     "Limit catch relative to natural mortality" = "#009E73",
                                     "Minimum size Limit (for each species)"= "slateblue",
                                     "Size limits and catch limits"= "yellow",
                                     "No-take zone" = "orange",
                                     "Seasonal closure" = "pink")) +
        annotate("rect", xmin = -Inf, xmax = t0-1, ymin = 10 ^ -Inf, ymax = 10 ^ Inf, fill = "lightskyblue", alpha = 0.2) +
        annotate("rect", xmin = t0-1, xmax = tNow - 1, ymin = 10 ^ -Inf, ymax = 10 ^ Inf, fill = "deepskyblue", alpha = 0.2) +
        #facet_grid(Management~.) +
        geom_line(size=1) +
        geom_vline(xintercept = tManagement-1) +
        ylab("Relative yield") +
        scale_y_log10() +
        ggtitle("Projected relative yield") #+
        #theme(axis.text.y=element_blank())
      
      tradeoffPlot <- values$combinedSpeciesOutput %>%
        filter(Management != "Historic Fishing" & Management != "Virgin Unfished") %>%
        filter(Species == paste(input$sppSelect,"_",input$siteAnalysis,sep="")) %>%
        filter(Year == tEnd) %>%
        ggplot(aes(x=((Biomass - min(values$combinedSpeciesOutput$Biomass)) / (max(values$combinedSpeciesOutput$Biomass) - min(values$combinedSpeciesOutput$Biomass))),
                   y=((Yield - min(values$combinedSpeciesOutput$Yield)) / (max(values$combinedSpeciesOutput$Yield) - min(values$combinedSpeciesOutput$Yield))),
                   color=Management)) +
        scale_colour_manual(values=c("Virgin Unfished"="black", 
                                     "Historic Fishing"="blue",
                                     "No Change"="red",
                                     "Limit catch relative to natural mortality" = "#009E73",
                                     "Minimum size Limit (for each species)"= "slateblue",
                                     "Size limits and catch limits"= "yellow",
                                     "No-take zone" = "orange",
                                     "Seasonal closure" = "pink")) +
        geom_point(size = 4) +
        xlab("Relative biomass in Water") +
        ylab("Relative yield") +
        ggtitle(paste("Tradeoff plot of projected biomass and yield in year",tEnd,sep=" ")) #+
        # theme(axis.text.y=element_blank(),
              # axis.text.x=element_blank())
      
      grid.arrange(tradeoffPlot,
                   bPlot,
                   yPlot,
                   ncol=1)
    }
  })
  
  output$PlotSizeSpectrum <- renderPlot({
    if(nrow(req(values$sizeSpectrum)) > 0) {
      t0 <-tNow - fisheryAge
      tVirgin <- t0 - 100
      tManagement <- input$managementStart
      tEnd <- tManagement + input$tProjection
      values$sizeSpectrum %>%
        group_by(Management,weightClass) %>%
        summarize(virginAbundance = sum(virginAbundance),
                  fishedAbundance = sum(n),
                  relativeAbundance = fishedAbundance/virginAbundance) %>%
        ggplot(aes(x=weightClass, y=relativeAbundance,color=Management)) +
        scale_colour_manual(values=c("Virgin Unfished"="black", 
                                     "Historic Fishing"="blue",
                                     "No Change"="red",
                                     "Limit catch relative to natural mortality" = "#009E73",
                                     "Minimum size Limit (for each species)"= "slateblue",
                                     "Size limits and catch limits"= "yellow",
                                     "No-take zone" = "orange",
                                     "Seasonal closure" = "pink")) +
        geom_line(size=1) +
        ylab(paste("Relative abundance in ",tEnd," compared to unfished",sep="")) +
        xlab("Size class [g]") +
        geom_hline(yintercept=1,linetype=2) +
        scale_x_log10() +
        ggtitle("Plot of size spectrum")
    }
  })
  observeEvent(input$revertParams,{
    values$params_belize <- params_belize_analysis
  })
  observeEvent(input$saveParams,{
    values$params_belize <- isolate(hot_to_r(input$hot))
  })
  output$hot <- renderRHandsontable({
    rhandsontable(values$params_belize)
    
  })
  
})