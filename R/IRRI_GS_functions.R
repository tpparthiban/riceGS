################################################################################
#############                 GS related functions                  ############
#############         IRRI Irrigated breeding Program               ############
#############  Parthiban Thathapalli Prakash, Jerome Bartholome     ############
#############                      July 2021                        ############
################################################################################

# Training set selection --------------------------------------------------

# This is a function is derived from the STPGA package documentation
# It allows for multiple starting points for the genetic algorithm and more robust results

#' Title
#'
#' @param K
#' @param sTS
#' @param repi
#' @param repo
#'
#' @return
#' @export
#'
#' @examples

OptiTS <- function(K,
                   sTS = 100,
                   rep = 1) {
  LambdaTrait <- 1 / ncol(K)
  svdK <- svd(K, nu = 10, nv = 10)
  PCAs <- K %*% svdK$v
  rownames(PCAs) <- rownames(K)
  colnames(PCAs) <- paste("PCA.", seq(1:ncol(PCAs)), sep = "")

  l_TS <- vector(mode = "list", length = rep)
  for (i in 1:rep) {
    print(paste("Rep:", i, sep = ""))
    cat("\n")
    TS <-
      GenAlgForSubsetSelectionNoTest(
        P = PCAs,
        ntoselect = sTS,
        InitPop = NULL,
        npop = 100,
        nelite = 10,
        mutprob = .5,
        mutintensity = rpois(1,4),
        niterations = 200,
        minitbefstop = 50,
        plotiters = F
      )
    l_TS[[i]] <- TS[[1]]
  }

  allTS <- sort(table(unlist(l_TS)), decreasing = T)
  fTS <- names(allTS)[1:sTS]

  return(fTS)
}




# pruner_asreml -----------------------------------------------------------

#' Title
#'
#' @param model
#' @param dataset
#'
#' @return
#' @export
#'
#' @examples

pruner_asreml <-
  function(model, dataset) {
    stdresid <- model$aom$R[, "stdCondRes"]
    ndf <- model$nedf
    studRes <- stdresid / sqrt((ndf - stdresid ^ 2) / (ndf - 1))

    #outlier points according to alpha level of 0.001
    out <-
      c(which(pt(studRes, df = ndf) > 0.99), which(pt(studRes, df = ndf) < 0.01))

    #if there are outliers remove them and re-fit model
    NA_before <- sum(is.na(dataset$trait))

    # Remove the points and uodate the model
    if (length(out) > 0) {
      #Convert outliers to NA
      dataset[out, 'trait'] <- NA
      model <- update(model, data = dataset)
    }

    NA_after <- sum(is.na(dataset$trait))
    Outliers_removed <- NA_after - NA_before
    print(paste("Number of outliers removed =", Outliers_removed))
    # plot(model)
  }


# pruner_sommer -----------------------------------------------------------


#' Title
#'
#' @param model
#' @param dataset
#'
#' @return
#' @export
#'
#' @examples

pruner_sommer <-
  function(model, dataset) {
    stdresid <- residuals(object = model)$trait.residuals
    ndf <- length(stdresid)
    studRes <- stdresid / sqrt((ndf - stdresid ^ 2) / (ndf - 1))

    #outlier points according to alpha level of 0.001
    out <-
      c(which(pt(studRes, df = ndf) > 0.99), which(pt(studRes, df = ndf) < 0.01))

    #if there are outliers remove them and re-fit model
    NA_before <- sum(is.na(dataset$trait))

    # Remove the points and update the model
    if (length(out) > 0) {
      #Convert outliers to NA
      dataset$trait[out] <- NA
      if (unique(dataset$DESIGN) == "RCBD")
      {
        model <-
          mmer(
            fixed = trait ~ 1 + REP,
            random = ~ DESIGN_X + DESIGN_Y + GID,
            rcov = ~ units,
            data = dataset,
            verbose = FALSE
          )
      }
      else {
        model <-
          mmer(
            fixed = trait ~ 1,
            random = ~ DESIGN_X + DESIGN_Y + GID,
            rcov = ~ units,
            data = dataset,
            verbose = FALSE
          )
      }
    }

    NA_after <- sum(is.na(dataset$trait))
    Outliers_removed <- NA_after - NA_before
    print(paste("Number of outliers removed =", Outliers_removed))
    # plot(model)
    return(model)
  }




# single_trial_asremlr ------------------------------------------------------------


#' Title
#'
#' @param dataset
#'
#' @return
#' @export
#'
#' @examples

single_trial_asreml <-
  function(dataset) {
    # fit model with asreml
    if (unique(dataset$DESIGN) == "RCBD"){
      asreml.options(aom = T)
      model <-
        asreml(
          fixed = trait ~ REP ,
          random = ~ DESIGN_X + DESIGN_Y + GID,
          na.action = na.method(x = "include"),
          data = dataset
        )
      pruner_asreml(model, dataset)
    } else {
      asreml.options(aom = T)
      model <-
        asreml(
          fixed = trait ~ 1 ,
          random = ~ DESIGN_X + DESIGN_Y + GID,
          na.action = na.method(x = "include"),
          data = dataset
        )
      pruner_asreml(model, dataset)
    }

    # plot(model, title = paste(loc, trait, sep = " "))

    #Get blups and pevs
    q <- predict(model, classify = 'GID')
    blup <- q$pvals$predicted.value
    gids <- q$pvals$GID
    names(blup) <- gids
    pev <- c(q$pvals$std.error) ^ 2

    #var comps
    sm <- summary(model)$varcomp
    vg <- sm[match('GID', row.names(sm)), 1]
    ve <- sm[match('units!R', row.names(sm)), 1]
    h2 <- vpredict(model, h2 ~ V3 / (V3 + V4))$Estimate

    #reliability
    rel <- 1 - pev / vg
    ### Deregress BLUP
    dblup <- (blup - mean(blup)) / rel
    dblup <- dblup + mean(blup)

    rslts <-
      tibble::tibble(
        location = unique(na.omit(dataset$LOCATION)),
        trait = trait,
        gid = gids,
        blups = blup,
        dblup = dblup,
        AIC = summary(model)$aic,
        rel = rel,
        pev = pev,
        Vg = vg,
        Ve = ve,
        h2 = h2
      )
    return(list(model, rslts))
  }


# single_trial_sommer ------------------------------------------------------------


#' single_trial_sommer
#'
#' @param dataset
#'
#' @return
#' @export
#'
#' @examples

single_trial_sommer <-
  function(dataset) {
    # fit model with asreml
    if (unique(dataset$DESIGN) == "RCBD"){
      model <-
        mmer(
          fixed = trait ~ REP,
          random = ~ DESIGN_X + DESIGN_Y + GID,
          rcov = ~ units,
          data = dataset,
          verbose = FALSE
        )
      # model <- pruner_sommer(model, dataset)
    } else {
      model <-
        mmer(
          fixed = trait ~ 1,
          random = ~ DESIGN_X + DESIGN_Y + GID,
          rcov = ~ units,
          data = dataset,
          verbose = FALSE
        )
      # model <- pruner_sommer(model, dataset)
    }

    # plot(model, title = paste(loc, trait, sep = " "))

    #Get blups and pevs
    q <- randef(model)
    blup <- q$GID$trait
    gids <- names(blup)

    pev <- diag(model$PevU$GID$trait)

    #var comps
    sm <- summary(model)$varcomp
    vg <- sm[grepl('GID', row.names(sm)), 1]
    ve <- sm[grepl('units', row.names(sm)), 1]
    h2 <- vpredict(model, h2 ~ V3 / (V3 + V4))$Estimate

    #reliability
    rel <- 1 - pev / vg
    ### Deregress BLUP
    dblup <- (blup) / rel
    # dblup <- dblup + mean(blup)

    mu <- model$Beta$Estimate

    blup <- blup + mu
    dblup <- dblup + mu

    rslts <-
      tibble::tibble(
        location = unique(na.omit(dataset$LOCATION)),
        trait = trait,
        gid = gids,
        blups = blup,
        dblup = dblup,
        AIC = model$AIC,
        rel = rel,
        pev = pev,
        Vg = vg,
        Ve = ve,
        h2 = h2
      )
    return(list(model, rslts))
  }





# gblup_asreml ------------------------------------------------------------

#' gblup_asreml
#'
#' @param dataset
#' @param invG
#'
#' @return
#' @export
#'
#' @examples

gblup_asreml <-
  function(dataset, invG) {
    # fit model with asreml

    asreml.options(aom = T)
    asreml.options(aom = T)
    model <-
      asreml(
        fixed = trait ~ 1 + location,
        random = ~ vm(gid, invG),
        data = dataset
      )


    #Get blups and pevs
    p <- predict(model, classify = 'gid', only = c('gid'))

    BV <- p$pvals$predicted.value
    gid <- p$pvals$gid
    se <- p$pvals$std.error
    names(BV) <- gid
    avsed <- p$avsed
    pev <- c(p$pvals$std.error) ^ 2

    #var comps
    sm <- summary(model)$varcomp
    nvar <- nrow(sm) + length(as.vector(coefficients(model)$fixed))
    vg <- sm[match('vm(gid, invG)', row.names(sm)), 1]
    ve <- sm[match('units!R', row.names(sm)), 1]

    #reliability
    rel <- 1 - (pev / vg)
    h2 <- vpredict(model, h2 ~ V1 / (V1 + V2))$Estimate

    rslts <-
      tibble::tibble(
        breedingzone = unique(na.omit(dataset$breedingzone)),
        trait = gsub("_dblup", "", trait),
        gid = gid,
        BV = BV,
        se = se,
        pev = pev,
        Vg = vg,
        h2 = h2,
        rel = rel
      )
    return(list(model, rslts))
  }



# gblup_sommer ------------------------------------------------------------


#' gblup_sommer
#'
#' @param dataset
#' @param G
#'
#' @return
#' @export
#'
#' @examples
#'

gblup_sommer <-
  function(dataset, G) {
    # fit model with asreml

    model <- mmer(
      fixed = trait ~ 1 + location,
      random = ~ vs(gid, Gu = G),
      rcov = ~ vs(units),
      data = dataset,
      verbose = FALSE
    )

    #Get blups and pevs
    BV <- randef(model)
    BV <- BV$`u:gid`$trait
    gid <- names(BV)

    pev <- diag(model$PevU$`u:gid`$trait)

    #var comps
    sm <- summary(model)$varcomp
    vg <- sm[grepl('gid', row.names(sm)), 1]
    ve <- sm[grepl('units', row.names(sm)), 1]
    h2 <- vpredict(model, h2 ~ V1 / (V1 + V2))$Estimate

    #reliability
    rel <- 1 - pev / vg
    ### Deregress BLUP

    mu <- model$Beta$Estimate[1]

    BV <- BV + mu

    rslts <-
      tibble::tibble(
        breedingzone = unique(na.omit(dataset$breedingzone)),
        trait = gsub("_dblup", "", trait),
        gid = gid,
        BV = BV,
        rel = rel,
        # se = se,
        pev = pev,
        Vg = vg,
        Ve = ve,
        h2 = h2
      )
    return(list(model, rslts))
  }



