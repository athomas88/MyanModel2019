# ui.R
shinyUI(fluidPage(
  navbarPage(("Myanmar Fishery Simulation / မြန်မာ့ငါးလုပ်ငန်းစံပြုပုံစံ "),
             tabPanel("Model Inputs and Outputs / စံပြုပုံစံ သွင်းအားစု နှင့် ရလဒ်",
                      sidebarLayout(
                        sidebarPanel(
                          h2("Model Inputs / စံပြုပုံစံ သွင်း"),
                          selectInput("siteAnalysis",label="Select site / နေရာရွေးချယ်ရန်",choices=c("Rakhine", "Delta", "Tanintharyi"),selected="Delta"),
                          sliderInput(inputId = "managementStart", label = "When should new management start? / မည်သည့်အချိန်တွင် စီမံခန့်ခွဲမှုအသစ်အညွှန်းကို စတင်မည် ", 2019, 2039, 2019, step = 1, sep=""),
                          sliderInput(inputId = "tProjection", label = "Years projected into future / အနာဂတ်မျှော်မှန်းခုနှစ်", 1, 100, 50, step = 1),
                          sliderInput(inputId = "effortCreep", label = "Increase in fishing over time [% per year] / ငါးဖမ်းစီးမှုမြင့်တက်မည့်နှုန်း(တနှစ်အတွင်းတိုးလာမည့်
ရာခိုင်နှုန်း)", 0, 5, 1, step = 0.1),
                          sliderInput(inputId = "compliance", label = "What percentage of the fishery will obey management rules?", 0, 100, 80, step = 5),
                          checkboxGroupInput("interventionGroup", label = "Select potential management options / အလားအလာရှိသေစီမံခန့်ခွဲမှုပုံစံအားရွေးချယ်ခြင်း", 
                                             choices = list("No Change" =  1, 
                                                            "Limit catch relative to natural mortality" = 2,
                                                            "Minimum size limit (for each species)" = 3,
                                                            "Size limits and catch limits" = 4,
                                                            "No-take zone" = 5,
                                                            "Seasonal closure" = 6),
                                             selected = c(1)),
                          h5("No Change = ပြောင်းလဲမှုမရှိ"),
                          h5("Limit catch relative to natural mortality = သဘာဝအလျောက်သေဆုံးမှုအားရွေးချယ"),
                          h5("Minimum size limit (for each species) = အငယ်ဆုံးအရွယ်အစားကန့်သတ်ချက်"),
                          h5("Size limits and catch limits = အငယ်ဆုံးအရွယ်အစားကန့်သတ်ချက် နှင့် ငါးဖမ်းစီးမှုပမာဏကန့်သတ်ချက်"),
                          h5("No-take zone = ငါးမဖမ်းရဇုံ"),
#                          h5("Seasonal closure = TRANSLATE"),
                          hr(),
                          uiOutput("customControlF"),
                          uiOutput("customControlFsize"),
                          hr(),
                          uiOutput("customControlNTZ"),
                          hr(),
                          uiOutput("customControlNTseason"),
                          uiOutput("customControlNTF"),
                          hr(),
                          actionButton(inputId = "runSim", label = "Run / စံပြုပုံစံအားလုပ်ဆောင်စေခြင်")
                        ),
                        mainPanel(
                          h2("Model Outputs / စံပြုပုံစံရလဒ်များ"),
                          selectInput("plotSelect", label = "Select plot types / စံပြုပုံစံအမျိုးအစားရွေးချယ်ခြင်", 
                                      choices = list("Total catch and population of all species" = 1, "Species-specific catch and population" = 2,"Projected size spectrum" = 3), 
                                      selected = 1),
                          h5("Total catch and population of all species = ငါးမျိုးစိတ်အားလုံး၏ စုစုပေါင်းဖမ်းစီးမှုနှင့်ငါးပမာဏ"),
                          h5(" Species specific catch and population = သတ်မှတ်ရွေးချယ်ထားသောငါးမျိုးစိတ်၏ဖမ်းစီးမှုနှင့် ငါးပမာဏ"),
                          h5("Projected size spectrum = အရွယ်အစားသတ်မှတ်ချက်"),
                          
                          conditionalPanel(
                            condition = "input.plotSelect == 1",
                            plotOutput("PlotProjectionsAggregated",height="1200px",width="1000px")
                          ),
                          conditionalPanel(
                            condition = "input.plotSelect == 2",
                            uiOutput("speciesSelect"),
                            plotOutput("PlotProjectionsSpecies",height="1200px",width="1000px")
                          ),
                          conditionalPanel(
                            condition = "input.plotSelect == 3",
                            plotOutput("PlotSizeSpectrum",height="400px",width="1000px")
                          ),
                          hr(),
                          print("contact: Alice Thomas-Smyth - athomas@edf.org")
                        )
                      )
             ),
             tabPanel("Species-specific inputs / သတ်မှတ်ရွေးချယ်ထားသောငါးမျိုးစိတ်အတွက်သွင်းအားစုမ",
                      h5("သင်အနေဖြင့်ငါးမျိုးစိတ်သွင်းအားစုများအားပြုပြင်နိုင်ပါသည်။
ကျွနုပ်တို့၏သုတေသနရလဒ်မှသွင်းအားစုများမှာအောက်ပါ အတိုင်းဖြစ်သည်။စုများ"),
                      h5("ဇယားကွက်တွင်ပြောင်းလဲမှုပြုလုပ်ပြီး “Save” အားနှိပ်ပါ"),
                      h5("သွင်းအားစုများအားသင့်အနေဖြင့် ”Revert” အားနှိပ်ပြီးပြောင်းလဲနိုင်သည်။"),
                      rHandsontableOutput("hot"),
                      actionButton(inputId = "saveParams", label = "Save / အညွှန်းကိန်းများအားသိမ်းဆည်းမည်"),
                      actionButton(inputId = "revertParams", label = "Revert / အညွှန်းကိန်းများအားပြောင်းလဲမည်")
             )
  )
)
)