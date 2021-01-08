source('functions.R')
options(shiny.reactlog = FALSE)
options(shiny.fullstacktrace = TRUE)

# Define UI for application
ui <- dashboardPage(skin = "black",
                    dashboardHeader(title = 'xtal-sttl'),
                    dashboardSidebar(sidebarMenu(
                        menuItem(
                            'Magma Density',
                            tabName = 'density',
                            icon = icon('cube')
                        ),
                        menuItem(
                            'Magma Viscosity',
                            tabName = 'viscosity',
                            icon = icon('vial')
                        ),
                        menuItem(
                            'Stokes Velocity',
                            tabName = 'stokes',
                            icon = icon('tachometer-alt')
                        ),
                        menuItem('Options',
                                 tabName = 'options',
                                 icon = icon('sliders-h')),
                        menuItem('Acknowledgements',
                                 tabName = 'ack',
                                 icon = icon('book-open'))
                    )),
                    body
)

server <- function(input, output, session) {
    # Initialize reactive values (reactive values store values that need to change)
    vals <- reactiveValues()
    # Locations of peak positions (output from isoplotR function 'peakfit')
    vals$dens <- NULL
    vals$data <- tibble(ID = rep(NA_character_, 15),
                        SiO2 = rep(NA_real_, 15),
                        TiO2 = rep(NA_real_, 15),
                        Al2O3 = rep(NA_real_, 15),
                        Fe2O3 = rep(NA_real_, 15),
                        FeO = rep(NA_real_, 15),
                        MgO = rep(NA_real_, 15),
                        CaO = rep(NA_real_, 15),
                        Na2O = rep(NA_real_, 15),
                        K2O = rep(NA_real_, 15),
                        H2O = rep(NA_real_, 15))
    vals$dens <- tibble(ID = rep(NA_character_, 15),
                        SiO2 = rep(NA_real_, 15),
                        H2O = rep(NA_real_, 15),
                        T = rep(NA_real_, 15),
                        P = rep(NA_real_, 15),
                        rho.melt = rep(NA_real_, 15),
                        sig.rho = rep(NA_real_, 15))
    vals$log.n <- tibble(ID = rep(NA_character_, 15),
                         SiO2 = rep(NA_real_, 15),
                         H2O = rep(NA_real_, 15),
                         T = rep(NA_real_, 15),
                         P = rep(NA_real_, 15),
                         rho.melt = rep(NA_real_, 15),
                         sig.rho = rep(NA_real_, 15),
                         n.melt = rep(NA_real_, 15),
                         n.mush = rep(NA_real_, 15))
    vals$v.stokes <- tibble(ID = rep(NA_character_, 15),
                            SiO2 = rep(NA_real_, 15),
                            H2O = rep(NA_real_, 15),
                            T = rep(NA_real_, 15),
                            P = rep(NA_real_, 15),
                            rho.xtal = rep(NA_real_, 15),
                            sig.rho = rep(NA_real_, 15),
                            n.melt = rep(NA_real_, 15),
                            n.mush = rep(NA_real_, 15),
                            rho.melt = rep(NA_real_, 15),
                            vel = rep(NA_real_, 15))
    vals$x1 <- NULL
    vals$y1 <- NULL
    vals$x2 <- NULL
    vals$y2 <- NULL
    # Density___________________________________________________________________
    # Data table
    output$table <- renderRHandsontable({
        d <- vals$data
        rhandsontable(d, stretchH = 'all') %>%
            hot_cols(columnSorting = TRUE)
    })
    # Converts data table into R object
    observe({
        if (!is.null(input$table)) {vals$data <- hot_to_r(input$table)} 
    })
    # Handle toggle input
    observeEvent(input$toggle, {
        if(input$toggle == FALSE) {
            vals$data$toggle <- FALSE
        } else {
            vals$data$toggle <- TRUE
        }
    })
    # Handle density calculation (code modified from Iacovino and Till (2019))
    observeEvent(input$calc1, {
        # Tidy data
        d <- vals$data %>% drop_na() %>% replace(is.na(.), 0) %>%
            mutate(across(where(is.logical), as.numeric)) %>% filter(toggle == TRUE)
        # Make grid for P-T ranges
        grid <- expand.grid(T = seq(input$Tlow1, input$Thigh1, length.out = 5),
                            P = seq(input$Plow1, input$Phigh1, length.out = 5))
        # Calculate magma density for each composition for each P-T point
        withProgress(message = 'Calculating Density', value = 0, {
            vals$dens <- purrr::map2_df(grid$P, grid$T, ~{
                # Increment the progress bar
                incProgress(1/nrow(grid))
                # Calculate density
                magma_density(SiO2 = d$SiO2,
                              TiO2 = d$TiO2,
                              Al2O3 = d$Al2O3,
                              Fe2O3 = d$Fe2O3,
                              FeO = d$FeO,
                              MgO = d$MgO,
                              CaO = d$MgO,
                              Na2O = d$Na2O,
                              K2O = d$K2O,
                              H2O = d$H2O,
                              P = .x,
                              T = .y,
                              ID = d$ID)
            }) %>% group_by(ID) %>% arrange(ID, T, P)
        })
        i <<- vals$data
        j <<- vals$dens
    })
    # Results table
    output$table2 <- renderRHandsontable({
        d <- vals$dens
        rhandsontable(d, stretchH = 'all') %>%
            hot_cols(columnSorting = TRUE)
    })
    # Reactive plot function updates the data plot
    p1 <- reactive({
        d <- vals$data %>% drop_na() %>% replace(is.na(.), 0) %>%
            mutate(across(where(is.logical), as.numeric)) %>% filter(toggle == TRUE)
        if (nrow(d) != 0) {
            if(input$type1 == 'Summary'){
                p <- d %>% select(where(is.numeric), -toggle) %>% pivot_longer(everything(), 'variable') %>%
                    ggplot() + geom_histogram(aes(x = value)) + labs(x = NULL, y = NULL) + facet_wrap(~variable, scales = 'free') + theme_classic(base_size = 14) + theme(strip.background = element_rect(color = NA))
            } else if(input$type1 == 'Summary Density'){
                req(!is.null(vals$dens))
                p1 <- vals$dens %>% filter(ID %in% vals$data$ID[vals$data$toggle == TRUE]) %>% ggplot() + geom_histogram(aes(x = rho.melt)) + labs(x = 'Liquid Density (g/cc)', y = NULL) + theme_classic(base_size = 14)
                p2 <- vals$dens %>% filter(ID %in% vals$data$ID[vals$data$toggle == TRUE]) %>% ggplot() + geom_boxplot(aes(x = SiO2, y = rho.melt, group = ID), width = 0.5) + labs(x = 'SiO2 (wt%)', y = 'Liquid Density (g/cc)') + theme_classic(base_size = 14)
                p <- p1 + p2
            } else if(input$type1 == 'Density P' && !is.null(vals$dens)){
                req(!is.null(vals$dens))
                p <- vals$dens %>% filter(ID %in% vals$data$ID[vals$data$toggle == TRUE]) %>% ggplot() + geom_point(aes(x = P, y = rho.melt, group = interaction(ID, T)), alpha = 0.2, shape = 20) + geom_line(aes(x = P, y = rho.melt, group = interaction(ID, T), color = T), size = 0.8, se = F) + labs(x = 'P (kbar)', y = 'Liquid Density (g/cc)', color = 'Celcius') + scale_color_viridis_c(option = 'D') + theme_classic(base_size = 14)
            } else if(input$type1 == 'Density T'){
                req(!is.null(vals$dens))
                p <- vals$dens %>% filter(ID %in% vals$data$ID[vals$data$toggle == TRUE]) %>% ggplot() + geom_point(aes(x = T, y = rho.melt, group = interaction(ID, P)), alpha = 0.2, shape = 20) + geom_line(aes(x = T, y = rho.melt, group = interaction(ID, P), color = P), size = 0.8, se = F) + labs(x = 'T (Celcius)', y = 'Liquid Density (g/cc)', color = 'kbar') + scale_color_viridis_c(option = 'D') + theme_classic(base_size = 14)
            }
        } else {
            p <- NULL
        }
        return(p)
    })
    # Render plot
    output$p1 <- renderPlot({
        p <- p1 %>% debounce(200)
        p()
    })
    # Viscosity_______________________________________________________________
    # Handle viscosity calculation (VFT model from Hess and Dingwell (1996))
    observeEvent(input$calc2, {
        # Require data
        req(!is.null(vals$dens))
        # Constants
        # VFT model
        a1 <- -3.54
        a2 <- 0.83
        b1 <- 9601
        b2 <- -2366
        c1 <- 196
        c2 <- 32
        # Crystal packing
        phi0 <- 0.6
        # Modify using xtal fraction -> n_magma = n_melt * (1 - (phi/phi0))^(-5/2)
        vals$log.n <- vals$dens %>%
            mutate(n.melt = (a1 + a2*log(ifelse(H2O == 0, 1e-15, H2O))) + ((b1 + b2*log(ifelse(H2O == 0, 1e-15, H2O))) / ((T + 273.15) - (c1 + c2*log(ifelse(H2O == 0, 1e-15, H2O))))),
                   n.mush = n.melt * (1 - (input$phi/phi0))^(-5/2))
        k <<- vals$log.n
    })
    # Density results table
    output$table3 <- renderRHandsontable({
        d <- vals$dens
        rhandsontable(d, stretchH = 'all') %>%
            hot_cols(columnSorting = TRUE)
    })
    # Results table
    output$table4 <- renderRHandsontable({
        d <- vals$log.n
        rhandsontable(d, stretchH = 'all') %>%
            hot_cols(columnSorting = TRUE)
    })
    # Viscosity plots
    p2 <- reactive({
        if(input$type2 == 'Summary'){
            p1 <- vals$log.n %>% filter(ID %in% vals$data$ID[vals$data$toggle == TRUE]) %>% ggplot() + geom_histogram(aes(x = n.melt)) + labs(x = 'Melt Viscosity (Pa s)', y = NULL) + theme_classic(base_size = 14)
            p2 <- vals$log.n %>% filter(ID %in% vals$data$ID[vals$data$toggle == TRUE]) %>% ggplot() + geom_histogram(aes(x = n.mush)) + labs(x = 'Mush Viscosity (Pa s)', y = NULL) + theme_classic(base_size = 14)
            p <- p1 + p2
        } else if(input$type2 == 'Viscosity T'){
            p <- vals$log.n %>% filter(ID %in% vals$data$ID[vals$data$toggle == TRUE]) %>% ggplot() + geom_point(aes(x = T, y = n.melt), alpha = 0.2) + geom_line(aes(x = T, y = n.melt, color = H2O, group = H2O), size = 0.8) + labs(x = 'T (Celcius)', y = 'Melt Viscosity (Pa s)', color = 'H2O wt%') + scale_color_viridis_c(option = 'D') + theme_classic(base_size = 14)
        } else if(input$type2 == 'Viscosity H2O'){
            p <- vals$log.n %>% filter(ID %in% vals$data$ID[vals$data$toggle == TRUE]) %>% ggplot() + geom_point(aes(x = H2O, y = n.melt), alpha = 0.2) + geom_line(aes(x = H2O, y = n.melt, color = T, group = T), size = 0.8) + labs(x = 'H2O (wt%)', y = 'Melt Viscosity (Pa s)', color = 'Celcius') + scale_color_viridis_c(option = 'D') + theme_classic(base_size = 14)
        }
        return(p)
    })
    # Render plot
    output$p2 <- renderPlot({
        p <- p2 %>% debounce(200)
        p()
    })
    # Stokes_______________________________________________________________
    # Handle stokes settling velocity calculation
    observeEvent(input$calc3, {
        # Calculate stokes velocity = (2 * g * r^2 * (rho_crystal - rho_melt)) / (9 * n_melt)
        # Constants
        g <- 9.81
        vals$v.stokes <- vals$log.n %>%
            mutate(rho.xtal = input$rho/1000,
                   vel = ((2 * g * input$rad^2 * (input$rho - rho.melt * 1000)) / (9 * n.melt * 10)) / 1e2 * 3.154e7)
        l <<- vals$v.stokes
    })
    # Results table
    output$table5 <- renderRHandsontable({
        d <- vals$v.stokes
        rhandsontable(d, stretchH = 'all') %>%
            hot_cols(columnSorting = TRUE)
    })
    # Stokes plots
    p3 <- reactive({
        p1 <- vals$v.stokes %>% filter(ID %in% vals$data$ID[vals$data$toggle == TRUE]) %>% ggplot() + geom_raster(aes(x = T, y = P, fill = vel)) + labs(x = 'T (Celcius)', y = 'Pressure (kbar)', fill = 'cm/yr') + scale_y_reverse() + scale_fill_viridis_c(option = 'D') + theme_classic(base_size = 14)
        p2 <- vals$v.stokes %>% filter(ID %in% vals$data$ID[vals$data$toggle == TRUE]) %>% ggplot() + geom_point(aes(x = T, y = vel), alpha = 0.2) + geom_line(aes(x = T, y = vel, group = interaction(ID, P), color = H2O)) + labs(x = 'T (Celcius)', y = 'Velocity (cm/yr)') + scale_color_viridis_c(option = 'D') + theme_classic(base_size = 14)
        p3 <- vals$v.stokes %>% filter(ID %in% vals$data$ID[vals$data$toggle == TRUE]) %>% ggplot() + geom_point(aes(x = P, y = vel), alpha = 0.2) + geom_line(aes(x = P, y = vel, group = interaction(ID, T), color = H2O)) + labs(x = 'Pressure (kbar)', y = 'Velocity (cm/yr)') + scale_color_viridis_c(option = 'D') + theme_classic(base_size = 14)
        p <- p1 + p2 + p3
        return(p)
    })
    # Render plot
    output$p3 <- renderPlot({
        p <- p3 %>% debounce(200)
        p()
    })
    # Acknowledgements___________________________________________
    # Render acknoledgments verbatim text output
    output$acknowledgements <- renderText({
        paste0('Author: Buchanan Kerswell\nOpen source MIT license\nhttps://github.com/buchanankerswell/xtal-sttl\nDensity calculation uses code translated directly from the python script of Iacovino and Till (2019)\nViscosity calculation is based on the VFT model of Hess and Dingwell (1996)')
    })
    # This observer handles the plot download
    observe({
        # Needs a filename input by the user
        req(input$filename)
        # Download button
        output$dlplot <-
            # Download the current plot as a pdf, with a filename, width, and height defined
            # by user inputs
            downloadHandler(filename = paste0(input$filename, '.pdf'),
                            content =
                                function(filename) {if(input$p.select == 'Magma Density') {
                                        p <- p1()
                                    } else if(input$p.select == 'Magma Viscosity'){
                                        p <- p2()
                                    } else if(input$p.select == 'Stokes Velocity'){
                                        p <- p3()
                                    }
                                    ggsave(
                                        filename,
                                        plot = p,
                                        device = pdf(
                                            width = input$pdf_width,
                                            height = input$pdf_height
                                            ))})
        })
}
# Run the application
shinyApp(ui = ui, server = server, options = 'quiet')
