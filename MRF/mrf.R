#library.path <- .libPaths()
#print(library.path)
#library(scatterplot3d, lib.loc = library.path)
#library(scatterplot3d)
# Plot RMSD against another variable in the data of either a singleton or pairwise# potential.
# self: assumed to be a singleton or pairwise "instance"
# secondary.column: string to identify the column in the second axis
# tertiary.column: another string that, when provided, triggers a 3D plot
# exclude.zeroes: true if the zeroes in the secondary and tertiary are ignored
# rmsd.max: rows in data with RMSD < rmsd.max are ignored
rmsd.plot.potential <- function(self, secondary.column, tertiary.column=NA, exclude.zeroes=TRUE, rmsd.max=NA) {
    selected.rows <- self$data
    if (exclude.zeroes) {
        selected.rows <- selected.rows[selected.rows[secondary.column] != 0,]
        if (!is.na(tertiary.column)) {
            selected.rows <- selected.rows[selected.rows[tertiary.column] != 0,]
        }
    }
    if (!is.na(rmsd.max)) {
        selected.rows <- selected.rows[selected.rows["RMSD"] < rmsd.max,]
    }
    # Doing [string] indexing returns a dataframe with just one column (the one specified in square brackts).
    # The additional [,1] extracts that single column which should be numeric. Otherwise plot wouldn't work because of the data type
    if (is.na(tertiary.column)) {
        if ("singleton.potential" %in% class(self)) {
            unit.column <- paste("Unit", self$label, sep=".")
            same.unit.factor <- selected.rows[unit.column] == self$label
        } else {
            unit.left.column <- paste("Unit", self$label.left, sep=".")
            unit.right.column <- paste("Unit", self$label.right, sep=".")
            same.unit.factor <- selected.rows[unit.left.column] == self$label.left & selected.rows[unit.right.column] == self$label.right
        }
        split.by.unit.data <- split(selected.rows, same.unit.factor)
        same.unit.data <- split.by.unit.data[["TRUE"]]
        swapped.unit.data <- split.by.unit.data[["FALSE"]]
        low.ylim <- min(selected.rows["Belief"])
        high.ylim <- max(selected.rows["Belief"])
        low.xlim <- min(selected.rows["RMSD"])
        high.xlim <- max(selected.rows["RMSD"])
        plot(same.unit.data["RMSD"][,1], same.unit.data[secondary.column][,1],
        main=paste("Unit(s):", self$label, "Potential:", secondary.column),
        xlab="RMSD", ylab=secondary.column, pch=16, cex=0.7,
        xlim=c(low.xlim, high.xlim), ylim=c(low.ylim, high.ylim))
        par(new=T)
        plot(swapped.unit.data["RMSD"][,1], swapped.unit.data[secondary.column][,1],
        xlab="", ylab="", pch=3, cex=0.7, col="red", axes=F,
        xlim=c(low.xlim, high.xlim), ylim=c(low.ylim, high.ylim))
        par(new=F)
    } else {
        scatterplot3d(selected.rows["RMSD"][,1],
        selected.rows[secondary.column][,1],
        selected.rows[tertiary.column][,1],
        main=paste("Unit(s):", self$label, "Potential:", secondary.column, "/", tertiary.column),
        xlab="RMSD", ylab=secondary.column, zlab=tertiary.column,
        pch=16, cex.symbols=0.7, angle=60, highlight.3d=TRUE)
    }
}

# wraps a png export around the regular rmsd.plot function by redirecting the output dev
export.rmsd.plot.potential <- function(self, output.prefix, secondary.column, tertiary.column=NA, exclude.zeroes=TRUE, rmsd.max=NA) {
    filename <- paste(output.prefix, self$label, secondary.column, sep="_")
    if (!is.na(tertiary.column)) {
        filename <- paste(filename, tertiary.column, sep="_");
    }
    if (!is.na(rmsd.max)) {
        filename <- paste(filename, "_rmsd-", rmsd.max, sep="")
    }
    if (!exclude.zeroes) {
        filename <- paste(filename, "with-zeroes", sep="_")
    }
    filename <- paste(filename, "png", sep=".")
    png(filename=filename, width=800, height=800, pointsize=18);
    rmsd.plot.potential(self, secondary.column, tertiary.column, exclude.zeroes, rmsd.max)
    dev.off()
}

# Auxiliary function used by both potential constructors in order to scale the feature columns so that they have std dev of 1
# and they all have positive values. It returns a data frame with the adjusted subset of columns indicated in column.names.
# The called is in charge of overwriting them if desired.
standardize.columns <- function(original.data, column.names) {
    # Subtract the mean and divide by std dev
    initial.columns <- as.data.frame(original.data[column.names])
    
    scalable <- apply(initial.columns, MARGIN=2, sd)
    scalable[scalable == 0] <- 1

    avg <- apply(initial.columns, MARGIN=2, mean)
   
    standardized.columns <- as.data.frame(initial.columns)
    standardized.minima <- sapply(standardized.columns, min)

    standardized.range <- sapply(standardized.columns, max) - standardized.minima
    standardized.range[standardized.range == 0] <- 1
    
    # Since values will be input into a log function, we want to avoid -Inf, thus we add a small value to avoid zero
    positive.columns <- scale(standardized.columns, center=standardized.minima - 0.000001, scale=standardized.range)
    
    return(as.data.frame(positive.columns))
    #return(standardized.columns)
}

# Constants for feature column names
feature.columns.singleton.potential <- c("CC", "Overlap")
feature.columns.pairwise.potential <- c("PhysicsScore", "vdw", "vdw_a", "vdw_r", "elec", "elec_sr_a", "elec_lr_a", "elec_sr_r", "elec_lr_r", "hb", "solv", "sasa", "acp","no_clashes")

# Creates a wrapper for ?.mrf file data
singleton.potential <- function(filename, prefix=NA, use.standardization=TRUE) {
    # Assume files end with .mrf extension
    label <- sub("\\.mrf$","", filename)
    if (!is.na(prefix)) {
        filename <- paste(prefix, filename, sep="")
    }
    data <- read.table(filename, header=TRUE, colClasses=c("numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","character"))
    # kept for file writing purposes, in case they are renamed by other functions
    original.columns <- names(data)
   
    if (use.standardization) {
        feature.columns <- feature.columns.singleton.potential
        adjusted.columns <- standardize.columns(data, feature.columns)
        
        for (column.name in names(adjusted.columns)) {
            data[[column.name]] <- adjusted.columns[[column.name]]
        }
    }
    
    out <- list(label=label, data=data, original.columns=original.columns)
    class(out) <- c("singleton.potential", "potential")
    invisible(out)
}

# Creates a wrapper for ?-?.mrf file data
pairwise.potential <- function(filename, prefix=NA, use.standardization=TRUE) {
    # Assume the filename ends in .mrf and that there's a single - separating both vertex labels in the filename.
    # e.g. A-B.mrf would represent the potential between vertex A and vertex B/
    both.labels <- sub("\\.mrf$","", filename)
    if (!is.na(prefix)) {
        filename <- paste(prefix, filename, sep="")
    }
    separate.labels <- strsplit(both.labels, "-")[[1]]
    label.left <- separate.labels[1]
    label.right <- separate.labels[2]
    data <- read.table(filename, header=TRUE,colClasses=c("numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","character","character"))
    
    # kept for file writing purposes, in case they are renamed by other functions
    original.columns <- names(data)
    original.clashes <- data[["no_clashes"]]
    data <- cbind(data,clashes = data[["no_clashes"]])

    if (use.standardization) {
        feature.columns <- feature.columns.pairwise.potential
        # These are lower-is-better columns that need to be inverted to be higher-is-better instead
        inverted.columns <- c("PhysicsScore", "vdw", "vdw_r", "elec", "elec_sr_r", "elec_lr_r", "solv", "sasa","no_clashes")
        
        for (column.name in inverted.columns) {
            data[[column.name]] <- -data[[column.name]]
        }
        
        adjusted.columns <- standardize.columns(data, feature.columns)
        
        for (column.name in names(adjusted.columns)) {
            data[[column.name]] <- adjusted.columns[[column.name]]
        }
    }
    
    out <- list(label=paste(label.left, "-", label.right, sep=""), label.left=label.left, label.right=label.right, data=data, original.columns=original.columns,clashes=original.clashes)
    class(out) <- c("pairwise.potential", "potential")
    invisible(out)
}

# A potential.collection contains two lists, one with singleton.potentials and another one with pairwise.potentials,
# created by using the filename vectors provided as parameters to this function.
# Additionally, it contains the description of the MRF graph that can be used for estimation purposes.
# Sample creation:
# potential.collection(c("A.mrf", "B.mrf", "C.mrf"), c("A-B.mrf", "B-C.mrf"))
# The two internal lists would contain keys "A", "B" and "C", in the example, for the singleton.potential list and "A-B",
# "B-C" would be the keys for the pairwise.potential list
potential.collection <- function(singleton.filenames, pairwise.filenames, use.standardization=TRUE) {
    singleton.potentials <- list()
    pairwise.potentials <- list()
    
    singleton.names <- c()
    pairwise.names <- c()
    
    for (singleton.file in singleton.filenames) {
        print(paste("Loading ", singleton.file))
        current.singleton <- singleton.potential(singleton.file, use.standardization=use.standardization)
        singleton.potentials[[current.singleton$label]] <- current.singleton
        # append the singleton name to the list used to create the graph
        singleton.names <- c(singleton.names, current.singleton$label)
    }
    for (pairwise.file in pairwise.filenames) {
        print(paste("Loading ", pairwise.file))
        current.pairwise <- pairwise.potential(pairwise.file, use.standardization=use.standardization)
        pairwise.potentials[[current.pairwise$label]] <- current.pairwise
        pairwise.names <- c(pairwise.names, current.pairwise$label)
    }
    factor.cliques <- factor.clique.collection(singleton.names, pairwise.names)
    graph <- factor.graph(factor.cliques)
    out <- list(singleton.potentials=singleton.potentials, pairwise.potentials=pairwise.potentials, graph=graph)
    class(out) <- "potential.collection"
    out <- rename.transformation.columns.potential.collection(out)
    invisible(out)
}

write.potential <- function(any.potential, output.prefix) {
    output.filename <- paste(output.prefix, "-", any.potential$label, ".mrf", sep="")
    sink(output.filename)
    is.first <- TRUE
    for (column.name in any.potential$original.columns) {
        if (!is.first) {
            cat("\t")
        }
        cat(column.name)
        is.first <- FALSE
    }
    cat("\n")
    sink()
    write.table(any.potential$data, file=output.filename, sep="\t",
    eol="\t\n", col.names=FALSE, row.names=FALSE, append=TRUE)
}

write.potential.collection <- function(self, output.prefix) {
    for (singleton in self$singleton.potentials) {
        write.potential(singleton, output.prefix)
    }
    for (pairwise in self$pairwise.potentials) {
        write.potential(pairwise, output.prefix)
    }
}

# Create a reduced version of a potential.collection keeping "samples" number of rows.
# A number of elements below low.rmsd.percentile and above high.rmsd.percentile are kept and the rest
# of the samples are drawn randomly from the complete potential tables.
# The samples are drawn from the singleton potentials and the rows in pairwise potential that do not
# contain the values sampled are removed from the potential.collection returned.
# Note: This function assumes that transformation columns have suffixes added by rename.transformation.columns.potential.collection
sample.potential.collection <- function(self, samples, low.samples, high.samples, low.rmsd.percentile, high.rmsd.percentile) {
    print("sample.potential.collection ")
    for (singleton in self$singleton.potentials) {
        # divide into < low.rmsd.threshold (TRUE) and > high.rmsd.threshold (FALSE)
        
        print("RMSD")
        print(singleton$data$RMSD)
        print("low.rmsd.percentile ")
        print(low.rmsd.percentile)
        print("high.rmsd.percentile ")
        print(high.rmsd.percentile)
        
        percentile.thresholds <- quantile(singleton$data$RMSD, c(low.rmsd.percentile, high.rmsd.percentile))
        
        low.rmsds <- split(singleton$data, singleton$data["RMSD"] < percentile.thresholds[1])[["TRUE"]]
        high.rmsds <- split(singleton$data, singleton$data["RMSD"] > percentile.thresholds[2])[["TRUE"]]
        # Get the samples from low and high RMSD's first
        reduced <- low.rmsds[sample(nrow(low.rmsds), low.samples),]
        reduced <- rbind(reduced, high.rmsds[sample(nrow(high.rmsds), high.samples),])
        while (nrow(reduced) < samples) {
            new.sample <- singleton$data[sample(nrow(singleton$data), 1),]
            if (!(row.names(new.sample) %in% row.names(reduced))) {
                reduced <- rbind(reduced, new.sample)
            }
        }
        self$singleton.potentials[[singleton$label]]$data <- reduced
    }
    for (pairwise in self$pairwise.potentials) {
        left.singleton <- self$singleton.potentials[[pairwise$label.left]]$data
        right.singleton <- self$singleton.potentials[[pairwise$label.right]]$data
        
        # merge by the transformation columns, which is all columns - features - RMSD
        left.columns <- setdiff(colnames(left.singleton), c("RMSD", feature.columns.singleton.potential))
        right.columns <- setdiff(colnames(right.singleton), c("RMSD", feature.columns.singleton.potential))
        
        reduced <- merge(x=left.singleton[left.columns], y=pairwise$data, by=left.columns)
        reduced <- merge(x=right.singleton[right.columns], y=reduced, by=right.columns)
        self$pairwise.potentials[[pairwise$label]]$data <- reduced
    }
    return(self)
}

# Call repeatedly export.rmsd.plot.potential for all singleton and pairwise potentials, with a manually
# curated set of columns that can be of interest for further analysis.
# self: an instance of potential.collection. If beliefs have been set, then plots related to beliefs are also generated.
# output.prefix: prefix added to all files generated
# beliefs.only: True if we just generate the images for beliefs, otherwise all relevant potentials are plotted.
export.plot.potential.collection <- function(self, output.prefix, beliefs.only=TRUE) {
    for (singleton.potential in self$singleton.potentials) {
        if (!beliefs.only) {
            export.rmsd.plot.potential(singleton.potential, output.prefix, "CC")
            export.rmsd.plot.potential(singleton.potential, output.prefix, "Overlap")
        }
        if (!is.null(singleton.potential$data$Belief)) {
            export.rmsd.plot.potential(singleton.potential, output.prefix, "Belief")
        }
    }
    for (pairwise.potential in self$pairwise.potentials) {
        if (!beliefs.only) {
            export.rmsd.plot.potential(pairwise.potential, output.prefix, "PhysicsScore")
            export.rmsd.plot.potential(pairwise.potential, output.prefix, "vdw")
            export.rmsd.plot.potential(pairwise.potential, output.prefix, "elec")
        }
        if (!is.null(pairwise.potential$data$Belief)) {
            export.rmsd.plot.potential(pairwise.potential, output.prefix, "Belief")
        }
    }
}

# Constants used to identify transformation columns
transformation.columns.singleton.potential <- c("Alpha_z", "Beta_y", "Gamma_x", "tx", "ty", "tz", "Unit")
left.transformation.columns.pairwise.potential <- c("Alpha_z_l", "Beta_y_l", "Gamma_x_l", "tx_l", "ty_l", "tz_l", "Unit_l")
right.transformation.columns.pairwise.potential <- c("Alpha_z_r", "Beta_y_r", "Gamma_x_r", "tx_r", "ty_r", "tz_r", "Unit_r")

# Every potential in the collection has a data frame with generic transformation column names like Alpha_z, tx, ty, etc.
# To easily merge different tables, without worrying about name clashes, this function will rename all columns to have suffixes like Alpha_z.A, tx.B, ty.C, etc.
# Each suffix will be the name of the node it represents. In singleton potentials name, that's simply the node name while in pairwise potentials,
# like in "A-B", it will replace _l by the first component in the edge name (e.g.. _l by .A) and _r suffixes using the second component (e.g. _r by .B).
# A copy of the potential.collection provided is returned, after renaming.
rename.transformation.columns.potential.collection <- function(self) {
    for (singleton.potential in self$singleton.potentials) {
        colnames(self$singleton.potentials[[singleton.potential$label]]$data) <- replace.single.column.names(colnames(singleton.potential$data),
        transformation.columns.singleton.potential, singleton.potential$label, FALSE)
    }
    for (pairwise.potential in self$pairwise.potentials) {
        # Two replacements are needed: left and right suffixes
        new.names <- replace.single.column.names(colnames(pairwise.potential$data), left.transformation.columns.pairwise.potential, pairwise.potential$label.left, TRUE)
        new.names <- replace.single.column.names(new.names, right.transformation.columns.pairwise.potential, pairwise.potential$label.right, TRUE)
        colnames(self$pairwise.potentials[[pairwise.potential$label]]$data) <- new.names
    }
    return(self)
}

# Auxiliary function that generates the new names for a single set of columns from just one potential. is.pairwise determines whether or not the suffix is
# just added to the end or if it replaces a _l or _r suffix. Called multiple times by rename.transformation.columns.potential.collection
replace.single.column.names <- function(column.names, template.column.names, suffix, is.pairwise) {
    return(as.vector(sapply(column.names,
    function(x) {
        if (x  %in% template.column.names) {
            if (is.pairwise) {
                # They should have a _l or _r suffix that we want removed
                paste(substr(x, 1, nchar(x)-2), suffix, sep=".")
            } else {
                paste(x, suffix, sep=".")
            }
        } else {
            x
        }
    })))
}

# Calculates the product of all potentials across all assignments of "X" variables in the MRF. This is used to compute the Z normalization constant
# in the MRF. For each X we compute exp(sum_i{phi_i(X)}), in other words, all potential values added and then exp applied to them
# The overall sum for all X is returned as a scalar.
# self: a potential collection
# potential.weights: model weights for CC, vdw, elec, etc. applied to each phi
# max.merged.rows: for efficiency purposes the values for each X are not computed one by one. max.merged.rows values of X are
#                  generated at a time and then they are merged with the other potential tables internally. Conceptually the same result
#                  would be obtained by evaluating each X separately.
compute.normalization.constant.potential.collection <- function(self, potential.weights, max.merged.rows) {
    total.nodes <- length(self$singleton.potentials)
    
    rows.in.nodes <- sapply(self$singleton.potentials, function(x) {nrow(x$data)})
   
    next.index.to.update <- total.nodes
    current.row.indices <- rep.int(1, total.nodes)
    
    all.x <- NULL
    counter <- 0
    z <- 0
    raw.sum <- 0
    while (next.index.to.update > 0) {
        # Add the X value to the current set and only compute the potential values if it has reached the merged.rows.threshold
        if (is.null(all.x)) {
            all.x <- get.transformation.from.row.indices(self, current.row.indices)
        } else {
            all.x <- rbind(all.x, get.transformation.from.row.indices(self, current.row.indices))
        }
        if (nrow(all.x) >= max.merged.rows) {
            counter <- counter + nrow(all.x)
            sums <- sum.beliefs.complete.transformations(self, all.x, potential.weights)
            z <- z + sums$z
            raw.sum <- raw.sum + sums$raw.sum
            all.x <- NULL
        }
        # prepare indices for next value of X
        if (current.row.indices[next.index.to.update] >= rows.in.nodes[next.index.to.update]) {
            # reset when we're past, many reset levels could be necessary
            while (next.index.to.update > 0) {
                if (current.row.indices[next.index.to.update] < rows.in.nodes[next.index.to.update]) {
                    current.row.indices[next.index.to.update] <- current.row.indices[next.index.to.update] + 1
                    next.index.to.update <- total.nodes
                    break;
                } else {
                    current.row.indices[next.index.to.update] <- 1
                    next.index.to.update <- next.index.to.update - 1
                }
            }
        } else { # update didn't end iteration over current index
            current.row.indices[next.index.to.update] <- current.row.indices[next.index.to.update] + 1
            # reset the inner loops to the min value and start from the innermost
            while (next.index.to.update < total.nodes) {
                next.index.to.update <- next.index.to.update + 1
                current.row.indices[next.index.to.update] <- 1
            }
        }
    }
    # The last subset may be of size < max.merged.rows so the loop will exit before processing it.
    if (!is.null(all.x)) {
        counter <- counter + nrow(all.x)
        sums <- sum.beliefs.complete.transformations(self, all.x, potential.weights)
        z <- z + sums$z
        raw.sum <- raw.sum + sums$raw.sum
    }
    return(list(m=counter, z=z, raw.sum=raw.sum))
}

# Given the transformations data frame that describes several translations and rotations (for each node in the graph), this function will go through
# every potential selecting the values that match the transformations and accumulate the beliefs.
# Returns a scalar with the sum of all potential values.
sum.beliefs.complete.transformations <- function(potential.collection, transformations, potential.weights) {
    # Create a new column that starts at zero values to accumulate
    transformations <- cbind(transformations, data.frame(Accumulated=rep(0, nrow(transformations))))
    for (singleton in potential.collection$singleton.potentials) {
        with.belief <- add.belief.column(singleton$data, potential.weights$singleton.weights)
        updated.sum <- merge(transformations, with.belief)
        updated.sum$Accumulated <- updated.sum$Accumulated + updated.sum$Belief
        # Only keep the transformation columns and accumulated values
        columns.kept <- setdiff(colnames(updated.sum), c("RMSD", "Belief", feature.columns.singleton.potential))
        transformations <- updated.sum[columns.kept]
    }
    for (pairwise in potential.collection$pairwise.potentials) {
        with.belief <- add.belief.column(pairwise$data, potential.weights$pairwise.weights)
        updated.sum <- merge(transformations, with.belief)
        updated.sum$Accumulated <- updated.sum$Accumulated + updated.sum$Belief
        # Only keep the transformation columns and accumulated values
        columns.kept <- setdiff(colnames(updated.sum), c("RMSD", "Belief", feature.columns.pairwise.potential))
        transformations <- updated.sum[columns.kept]
    }
    # Z is the normalization constant portion computed by this function.
    # After adding all potentials for a particular X it is necessary to apply exp to it because the original form of the potential is
    # exp(sum(phi_i * weight_i)) where each phi is a potential type raw.sum is also returned because it is used to compute the likelihood
    return(list(z=sum(exp(transformations$Accumulated)), raw.sum=sum(transformations$Accumulated)))
}

# Get a single row with all x,y,z and angles from singleton potentials, using row.indices as a guide to determine what rows are used.
get.transformation.from.row.indices <- function(potential.collection, row.indices) {
    node.labels <- as.vector(sapply(potential.collection$singleton.potentials, function(x) {x$label}))
    transformation.columns <- sapply(potential.collection$singleton.potentials,
    function(x) {
        setdiff(colnames(x$data), c("RMSD", feature.columns.singleton.potential))})
    current.x <- NULL
    for (node.index in 1:length(potential.collection$singleton.potentials)) {
        node.row <- row.indices[node.index]
        node.label <- node.labels[node.index]
        node <- potential.collection$singleton.potentials[[node.label]]
        node.x <- node$data[node.row, transformation.columns[,node.label]]
        if (is.null(current.x)) {
            current.x <- node.x
        } else {
            current.x <- cbind(current.x, node.x)
        }
    }
    return(current.x)
}

# Default values to be used as weighting factors in the MRF. Any of them can be replaced by a new value to produce a new model.
# These weights are fed into the inference functions.
default.potential.collection.weights <- function() {
    return(potential.collection.weights(CC=1.0, Overlap=1.0, PhysicsScore=1.0, vdw=0.1, elec=0.1))
}

# All weights are set to zero except the ones provided as parameters. All zero weights passed as single parameters will be removed
# from the attributes list of this object.
# override.weights.names and override.weights.values are an alternate way of initializing the weights used by the gradient optimization procedure.
# Since R's optim function only expects the values, separating the names and numeric values is more convenient. Zero Weights set through
# this mechanism are not deleted from the attributes list.
potential.collection.weights <- function(CC=0.0, Overlap=0.0, PhysicsScore=0.0, vdw=0.0, vdw_a=0.0, vdw_r=0.0, elec=0.0, elec_sr_a=0.0, elec_lr_a=0.0,
elec_sr_r=0.0, elec_lr_r=0.0, hb=0.0, solv=0.0, sasa=0.0, acp=0.0, no_clashes=0.0,
override.weights.names=NULL, override.weights.values=NULL) {
    singleton.weights <- list(CC=CC, Overlap=Overlap)
    pairwise.weights <- list(PhysicsScore=PhysicsScore, vdw=vdw, vdw_a=vdw_a, vdw_r=vdw_r, elec=elec, elec_sr_a=elec_sr_a, elec_lr_a=elec_lr_a,
    elec_sr_r=elec_sr_r, elec_lr_r=elec_lr_r, hb=hb, solv=solv, sasa=sasa, acp=acp, no_clashes=no_clashes)
    
    # Remove zero-value weights automatically to avoid computing those
    # Since their names will not be present the corresponding columns will be ignored by other functions
    for (name in names(singleton.weights)) {
        if (singleton.weights[[name]] == 0) {
            singleton.weights[[name]] <- NULL
        }
    }
    for (name in names(pairwise.weights)) {
        if (pairwise.weights[[name]] == 0) {
            pairwise.weights[[name]] <- NULL
        }
    }
    
    if (!is.null(override.weights.names)) {
        for (weight.index in 1:length(override.weights.names)) {
            name <- override.weights.names[weight.index]
            if (name %in% c("CC", "Overlap")) {
                singleton.weights[[name]] <- override.weights.values[weight.index]
            } else if (name %in% c("PhysicsScore", "vdw", "vdw_a", "vdw_r", "elec", "elec_sr_a", "elec_lr_a", "elec_sr_r", "elec_lr_r", "hb", "solv", "sasa", "acp","no_clashes")) {
                pairwise.weights[[name]] <- override.weights.values[weight.index]
            }
        }
    }
    
    out <-(list(singleton.weights=singleton.weights, pairwise.weights=pairwise.weights))
    class(out) <- "potential.collection.weights"
    invisible(out)
}

# self: a potential collection
# weights.used: array of string that contains the weights to be used
# initial.weights: 1-to-1 relationship with weights.used which provides the numeric values
# max.merged.rows: size of a buffer used for optimization purposes. Represents the number of rows kept in memory before calling merge operations
optimize.weights.potential.collection <- function(self, weights.used, initial.weights, max.merged.rows=250) {
    print("optimize.weights.potential.collection ")
    
    likelihood.closure <- function(weights) {
        calculate.likelihood(weights, weights.used, self, max.merged.rows)
    }
    optim.results <- optim(initial.weights, likelihood.closure, method="L-BFGS-B", control=list(fnscale=-1), lower=-5, upper=5)
    return(optim.results)
}

calculate.likelihood <- function(weights, weights.used, a.potential.collection, max.merged.rows) {
    # Assume that the weights are provided in the same order as the function defined.
    potential.weights <- potential.collection.weights(override.weights.names=weights.used, override.weights.values=weights)
    print(Sys.time())
    print(potential.weights)
    result <- compute.normalization.constant.potential.collection( a.potential.collection, potential.weights, max.merged.rows)
    likelihood.value <- result$raw.sum - result$m * log(result$z)
    print(list(Likelihood.Value=likelihood.value,Likelihood.Components=result))
    return(likelihood.value)
}

# Returns a copy of the potential.collection provided as the first parameter with an additional column added to each potential dataframe, called "Belief".
# Initially they are set to the local value given by applying the weights to each potential (using log values).
# Then, message passing is used to distribute the maximum values across the MRF which, in the end, should have enough information to generate a map assignment
# from the Belief values added to the potentials.
#
# self: an instance of potential.collection
# potential.weights: an instance of potential.collection.weights. If null beliefs are not recomputed and are assumed to be in the potentials already
map.beliefs.potential.collection <- function(self, potential.weights=NULL) {
    if (!is.null(potential.weights)) {
        # Add belief columns to each potential, based on the weights
        for (singleton.name in names(self$singleton.potentials)) {
            self$singleton.potentials[[singleton.name]]$data <- add.belief.column(self$singleton.potentials[[singleton.name]]$data, potential.weights$singleton.weights)
        }
        for (pairwise.name in names(self$pairwise.potentials)) {
            self$pairwise.potentials[[pairwise.name]]$data <- add.belief.column(self$pairwise.potentials[[pairwise.name]]$data, potential.weights$pairwise.weights)
        }
    }
    
    # Store all forward messages because, on the backward pass, they need to be subtracted to avoid double counting.
    all.messages <- list()
    
    for (index in 1:length(self$graph$edges)) {
        edge <- self$graph$edges[[index]]
        message.id <- paste(edge$source, edge$destination, sep="_")
    }
        
    # Pass forward messages
    for (forward.index in 1:length(self$graph$edges)) {
        forward.edge <- self$graph$edges[[forward.index]]
        message.id <- paste(forward.edge$source, forward.edge$destination, sep="_")
        print(paste("Sending forward message", message.id))
        message.result <- get.updated.beliefs.from.message(self, forward.edge)
        self$pairwise.potentials[[forward.edge$destination]] <- message.result$updated.factor
        all.messages[[message.id]] <- message.result$message.passed
    }
    
    # Pass backward ones. We can be updating singletons on the way back (check that we update the correct potential list)
    for (backward.index in length(self$graph$edges):1) {
        backward.edge <- self$graph$edges[[backward.index]]
        print(paste("Sending backward message ", backward.edge$destination, "_", backward.edge$source, sep=""))
        if (backward.edge$is.singleton) {
            # Note that on the backward pass we only care about the $updated.factor part of the result, since we're not storing the backward messages
            self$singleton.potentials[[backward.edge$source]] <- get.updated.beliefs.from.message(self, backward.edge, TRUE, all.messages)$updated.factor
        } else {
            self$pairwise.potentials[[backward.edge$source]] <- get.updated.beliefs.from.message(self, backward.edge, TRUE, all.messages)$updated.factor
        }
    }
    
    return(self)
}

add.belief.column <- function(potential.data.frame, weights.list) {
    print("add.belief.column ")
    
    if (length(weights.list) == 1) {
        potential.data.frame[["Belief"]] <- potential.data.frame[, names(weights.list)] * unlist(weights.list)
    } else {
        weighted.columns <- sweep(potential.data.frame[, names(weights.list)], 2, unlist(weights.list), "*")
        # Add the new column as the sum of the weighted columns
        potential.data.frame[["Belief"]] <- Reduce("+", weighted.columns[names(weights.list)])
    }
    return(potential.data.frame)
}

# Returns a subset of the potential data frame contained in self, keeping only the maximum Belief values for individual combinations of
# alpha, beta, gamma, tx, ty, tz. In singleton.potential, it means that we just remove a subset of the columns, but in pairwise.potential, we will select
# either the left or the right transformation and return the max for each of the transformation values.
# For example, if a potential represents A-B and we want to maximize by="B", we will select the maximum Belief for each combination of Alpha_z_r,
# Beta_y_r, Gamma_z_r, tx_r, ty_r, tz_r, Unit_r. This effectively removes all _l columns that represent the "A" part in this potential.
# The data frame returned will have a Belief column and the 7 columns mentioned, but the _l or _r suffix will be stripped from the name.
#
# self: a singleton.potential or pairwise.potential
# by: the name of the node that will be kept (e.g. "A" or "B" in an "A-B" factor) ignored in the case of singleton potentials
# returns: a data.frame with a Belief column and the 6 transformation columns
get.factor.maximization <- function(self, by) {
    transformation.columns <- add.transformation.column.suffix(by)
    maximized.factor <- aggregate(self$data[c("Belief")], by=as.list(self$data[transformation.columns]), FUN=max)
    return(maximized.factor)
}

# Returns a list of size 2 with:
# 1) an updated version of the factor that the message is passed to, where the Belief column values are updated by adding the beliefs from the message.
# 2) a copy of the Belief column passed from source to destination, which represents the "message" passed
#
# a.potential.collection: An instance of potential.collection where all potential data frames have their beliefs set
# a.factor.graph.edge: One of the edges already contained in the potential.collection passed to the function as a convenience
# is.direction.inverted: All edges are stored in the "forward direction" but they also need to be passed back eventually. When this is true,
#                       the potential returned is the modified source based on the message passed back
# previous.messages: list with all the messages forward pass, expected only when is.direction.inverted=TRUE
#
# returns a list with 2 elements:
# list$updated.factor, should replace the factor on the function caller side list$message.passed, belief values passed from source to destination
get.updated.beliefs.from.message <- function(a.potential.collection, a.factor.graph.edge, is.direction.inverted=FALSE, previous.messages=NULL) {
    source.name <- a.factor.graph.edge$source
    destination.name <- a.factor.graph.edge$destination
    message.variable <- a.factor.graph.edge$shared.variable
    # In singleton edges the source is always the singleton and the destination is pairwise.
    if (a.factor.graph.edge$is.singleton) {
        source.factor <- a.potential.collection$singleton.potentials[[source.name]]
    } else {
        source.factor <- a.potential.collection$pairwise.potentials[[source.name]]
    }
    destination.factor <- a.potential.collection$pairwise.potentials[[destination.name]]
    
    if (is.direction.inverted) {
        tmp.factor <- source.factor
        source.factor <- destination.factor
        destination.factor <- tmp.factor
    }
    # Figure out the names used to join the message to the destination
    message.columns <- add.transformation.column.suffix(message.variable)
    
    if (!is.null(previous.messages) && is.direction.inverted) {
        # Backward message, we need to subtract the message in the original direction to avoid double counting
        original.message.id <- paste(source.name, destination.name, sep="_")
        original.message <- previous.messages[[original.message.id]]
        adjusted.data <- merge(y=original.message, x=source.factor$data, by=message.columns)
        adjusted.data[["Belief.x"]] <- adjusted.data[["Belief.x"]] - adjusted.data[["Belief.y"]]
        names(adjusted.data)[names(adjusted.data) == "Belief.x"] <- "Belief"
        adjusted.data[["Belief.y"]] <- NULL
        source.factor$data <- adjusted.data
    }
    
    maximized.message <- get.factor.maximization(source.factor, message.variable)
    joint <- merge(y=maximized.message, x=destination.factor$data, by=message.columns)
    
    # Store the message, which is just the transformation columns plus the beliefs
    #message.passed <- joint[c(destination.columns, "Belief.y")]
    #names(message.passed)[names(message.passed) == "Belief.y"] <- "Belief"
    # Reformat joint to have both Beliefs summed names as Belief, and remove both temporary Belief columns
    joint[["Belief.y"]] <- joint[["Belief.x"]] + joint [["Belief.y"]]
    names(joint)[names(joint) == "Belief.y"] <- "Belief"
    joint[["Belief.x"]] <- NULL
    destination.factor$data <- joint
    
    return(list(updated.factor=destination.factor, message.passed=maximized.message))
}

# Add a node suffix to the standards transformation columns.
# e.g. Alpha_z will be returned as Alpha_z.A if variable.name is "A"
add.transformation.column.suffix <- function(variable.name) {
    column.suffix <- paste(".", variable.name, sep="")
    column.names <- c("Alpha_z", "Beta_y", "Gamma_x", "tx", "ty", "tz", "Unit")
    return(as.vector(sapply(column.names, function(x) paste(x, column.suffix, sep=""))))
}

# Prints top n rows sorted by the column specified. By default it returns the highest belief values from each potential table.
# Both the column to sort by and whether it's decreasing or not can be altered via parameters
#
# self: a potential.collection instance with Belief values set
# top.n: positive integer stating how many of the top elements are returned
# column: Column to sort by. It's intended to be either Belief or RMSD
# decreasing.order: It should be TRUE for Belief or FALSE for RMSD
# z : Optional parameter used when Belief is the target column and we want to show the probability
print.top.by.column.potential.collection <- function(self, top.n, output.file=NULL, column="Belief", decreasing.order=TRUE, z=NULL) {
    print("Write output file.")
    #options(max.print = 1000000)
    if(!is.null(output.file)) {
        sink(file=output.file);
    }
    #print("+++Singleton potentials+++")
    for (singleton in self$singleton.potentials) {
        cat(singleton$label)
        if (!is.null(z)) {
            singleton$data$Belief <- exp(singleton$data$Belief) / z
        }
        positions <- order(singleton$data[, column], decreasing=decreasing.order)#[1:top.n]
        options("width"=200)
        
        print(singleton$data[positions, c("RMSD", add.transformation.column.suffix(singleton$label), "Belief")],row.names = FALSE)
    }
    for (pairwise in self$pairwise.potentials) {
        cat(pairwise$label)
        if (!is.null(z)) {
            pairwise$data$Belief <- exp(pairwise$data$Belief) / z
        }
        positions <- order(pairwise$data[, column], decreasing=decreasing.order)#[1:(top.n*top.n)]
        print(pairwise$data[positions, c("RMSD", add.transformation.column.suffix(pairwise$label.left), add.transformation.column.suffix(pairwise$label.right)
        ,"clashes","no_clashes","Belief")],,row.names = FALSE)
    }
    if(!is.null(output.file)) {
        sink()
    }
}

factor.clique.collection <- function(singleton.names, pairwise.names) {
    neighbor.counts <- list()
    # Initialize counts to zero
    for (singleton.name in singleton.names) {
        neighbor.counts[[singleton.name]] = 0
    }
    for (pairwise.name in pairwise.names) {
        vertices <- strsplit(pairwise.name, "-")[[1]]
        neighbor.counts[[vertices[1]]] = neighbor.counts[[vertices[1]]] + 1
        neighbor.counts[[vertices[2]]] = neighbor.counts[[vertices[2]]] + 1
    }
    sorted.singletons <- neighbor.counts[order(as.numeric(neighbor.counts))]
    singleton.names <- names(sorted.singletons)
    cliques <- list()
    
    for (var.index in 1:length(sorted.singletons)) {
        variable.eliminated <- singleton.names[var.index]
        # Look for the pairwise terms that contain the variable
        clique.variables <- c()
        clique.pairwise.names <- c()
        for (pairwise.name.index in 1:length(pairwise.names)) {
            pairwise.name <- pairwise.names[pairwise.name.index]
            if (pairwise.name != "") { # i.e. it hasn't been used already
                vertices <- strsplit(pairwise.name, "-")[[1]]
                if (variable.eliminated %in% vertices) {
                    clique.pairwise.names = c(clique.pairwise.names, pairwise.name)
                    for (vertex in vertices) {
                        if (vertex != variable.eliminated) {
                            clique.variables = c(clique.variables, vertex)
                        }
                    }
                    # Mark as used
                    pairwise.names[pairwise.name.index] = ""
                }
            }
        }
        # Look for already found cliques that use variable.eliminated
        other.factor.indices = c()
        if (length(cliques) > 0) {
            for (factor.index in 1:length(cliques)) {
                if (variable.eliminated %in% cliques[[factor.index]]$clique.variables) {
                    other.factor.indices = c(other.factor.indices, factor.index)
                }
            }
        }
        # Create a new clique object with the information gathered
        current.clique <- factor.clique(variable.eliminated, clique.variables, c(variable.eliminated), clique.pairwise.names, other.factor.indices)
        cliques[[var.index]] <- current.clique
    }
    out <- list(cliques=cliques)
    class(out) <- "factor.clique.collection"
    invisible(out)
}

# Every tau clique in factor.clique.collection$cliques will be connected to other cliques upon creation of the factor graph. Connections in the factor
# graph will only be made between singleton-pairwise or pairwise-pairwise but NOT between clique-clique, singleton-clique or pairwise-clique.
# Thus, for the tau cases where cliques need to be connected, a representative pairwise factor will be selected from the clique, in order to make connections.
# This can either be the first pairwise term in a given clique that contains variable.name or a factor from other tau cliques contained as
# subfunctions of this factor clique.
#
# clique.collection: an instance of factor.clique.collection
# position: the position in clique.collection$cliques where the potentially recursive search starts from
# variable.name: the pairwise factor returned must contain this variable
get.first.pairwise.factor <- function (clique.collection, position, variable.name) {
    current.clique <- clique.collection$cliques[[position]]
    for (pairwise.name in current.clique$pairwise.names) {
        if (variable.name %in% strsplit(pairwise.name, "-")[[1]]) {
            return(pairwise.name)
        }
    }
    # Look for a factor in the other tau cliques associated
    for (other.clique in current.clique$other.factor.indices) {
        returned.factor <- get.first.pairwise.factor(clique.collection, other.clique, variable.name)
        if (!is.null(returned.factor)) {
            return(returned.factor)
        }
    }
    # return NULL otherwise to indicate that nothing was found
    return(NULL)
}

# Factor cliques represent cluster graph nodes that are created through the Variable Elimination algorithm, based on a series of
# input potentials given to potential.collection.
# However, these will only be symbolic cliques, i.e. they will only store the names of their singleton or pairwise components, and not the actual data.
# In addition, since the Variable Elimination algorithm creates intermediate factor cliques, each instance of factor.clique will also store indices for
# any other factor.clique that is required to compute the current factor.clique.
# The main purpose of factor.clique instances are to help in the construction of factor graphs, by storing the intermediate results for the VE algorithm.
# variable.eliminated: A string that identifies what variable was eliminated this clique was created
# clique.variables: Other than variable.eliminated, list all the variables that appear in this clique's potentials (including the other tau factors linked
#                   through other.factor.indices)
# singleton.names: List of strings that identify singleton.potentials (e.g. "A", "B", "C")
# pairwise.names: List of strings that identify pairwise.potentials (e.g. "A-B", "B-C")
# other.factor.indices: List of indices to an outer list that holds all factor.clique instances. These numbers represent the other factors (tau) that should
# be connected to this factor.clique
factor.clique <- function(variable.eliminated, clique.variables, singleton.names, pairwise.names, other.factor.indices) {
    out <- list(variable.eliminated=variable.eliminated, clique.variables=clique.variables, singleton.names=singleton.names,
    pairwise.names=pairwise.names, other.factor.indices=other.factor.indices)
    class(out) <- "factor.clique"
    invisible(out)
}

# Represents edges in a factor graph.
#
# source/destination: Messages first flow from source to destination. These are just the string labels for each potential
# shared.variable: The intersection variable between source and destination. The message passed marginalizes over this variable
# is.singleton: Marks that one of the two vertices is a singleton potential. Edges coming from singletons will be the first ones passed.
factor.graph.edge <- function(source, destination, shared.variable, is.singleton) {
    out <- list(source=source, destination=destination, shared.variable=shared.variable, is.singleton=is.singleton)
    class(out) <- "factor.graph.edge"
    invisible(out)
}

# Iterates through the cliques to infer what the graph should look like, expressed as a collection of edges. The edges are stored in the order
# that should be followed for message passing purposes
factor.graph <- function(clique.collection) {
    edges <- list()
    # temporary list of pairwise-pairwise edges that will be added to edges later, based on the lowest number of neighbors per pairwise factor
    unsorted.edges <- list()
    neighbor.counts <- list()
    for (clique in clique.collection$cliques) {
        # Since variable elimination was used to create the cliques, we know that each clique will have one (and only one) singleton
        singleton <- clique$singleton.names[1]
        all.pairwise <- clique$pairwise.names
        for (factor.position in clique$other.factor.indices) {
            all.pairwise <- c(all.pairwise, get.first.pairwise.factor(clique.collection, factor.position, singleton))
        }
        # The first edge starts with the singleton
        source <- singleton
        is.singleton <- TRUE
        
        for (pairwise.name in all.pairwise) {
            new.edge <- factor.graph.edge(source, pairwise.name, singleton, is.singleton)
            if (is.singleton) {
                edges[[length(edges) + 1]] <- new.edge
                # is.singleton <- FALSE
            }
            #else {
            #              unsorted.edges[[length(unsorted.edges) + 1]] <- new.edge
               source.count <- if(is.null(neighbor.counts[[source]]))
               0 else neighbor.counts[[source]]
               destination.count <- if(is.null(neighbor.counts[[pairwise.name]]))
               0 else neighbor.counts[[pairwise.name]]
               neighbor.counts[[source]] <- source.count + 1
               neighbor.counts[[pairwise.name]] <- destination.count + 1
                # }
                #  source <- pairwise.name
        }
    }
   
    out <- list(edges=edges, neighbor.counts=neighbor.counts)
    class(out) <- "factor.graph"
    invisible(out)
}

main.map <- function(singleton.files, pairwise.files, model.weights, prefix, top.output.per.potential) {
    mrf.collection <- potential.collection(singleton.files, pairwise.files)
    mrf.beliefs <- map.beliefs.potential.collection(mrf.collection, model.weights)

    export.plot.potential.collection(mrf.beliefs, prefix)
    sink(file=paste(prefix, "_weights.txt", sep=""))
    cat("Weights used for MAP estimation\n")
    print(model.weights)
    sink()
    print.top.by.column.potential.collection(mrf.beliefs, top.output.per.potential, output.file=paste(prefix, "_top100.txt", sep=""))
}

main.learn <- function(singleton.files, pairwise.files, weight.names, initial.model.weights, number.of.samples, low.samples,
high.samples, low.rmsd.percentile, high.rmsd.percentile, prefix, top.output.per.potential) {
    print(low.rmsd.percentile)
    print(high.rmsd.percentile)
    # standarize col, compute factor graph
    mrf.collection <- potential.collection(singleton.files, pairwise.files)
    sampled.collection <- sample.potential.collection(mrf.collection, number.of.samples, low.samples, high.samples, low.rmsd.percentile, high.rmsd.percentile)
    optimization.results <- optimize.weights.potential.collection(sampled.collection, weight.names, initial.model.weights)
    optimized.weights <- potential.collection.weights(override.weights.names=weight.names, override.weights.values=optimization.results$par)
    sampled.beliefs <- map.beliefs.potential.collection(sampled.collection, optimized.weights)
    all.beliefs <- map.beliefs.potential.collection(mrf.collection, optimized.weights)
    # Output the optimization results
    sink(file=paste(prefix, "_optimization.txt", sep=""))
    print(optimization.results)
    print(sampled.beliefs$singleton.potentials)
    print(sampled.beliefs$pairwise.potentials)
    sink()
    # Output beliefs for the sampled set
    export.plot.potential.collection(sampled.beliefs, paste(prefix, "_sampled", sep =""))
    print.top.by.column.potential.collection(sampled.beliefs, top.output.per.potential, output.file=paste(prefix, "_sampled_top.txt", sep=""))
    # Output beliefs for the complete set
    export.plot.potential.collection(all.beliefs, paste(prefix, "_complete", sep =""))
    print.top.by.column.potential.collection(all.beliefs, top.output.per.potential, output.file=paste(prefix, "_complete_top.txt", sep=""))
}

main.plot.singleton.rmsds <- function (singleton.files, image.prefix) {
    print("main.plot.singleton.rmsds ")
    
    percentiles.list <- c(0.05,0.1,0.15,0.2,0.25,0.75,0.8,0.85,0.9,0.95)
    if (is.na(image.prefix)) {
        sink(file="rmsd_percentiles.txt")
    } else {
        sink(file=paste(image.prefix, "_rmsd_percentiles.txt", sep=""))
    }
    cat("Unit", percentiles.list, "\n")
    for (singleton.file in singleton.files) {
        potential <- singleton.potential(singleton.file)
        if (is.na(image.prefix)) {
            output.filename <- paste(potential$label, "-hist.png", sep="")
        } else {
            output.filename <- paste(image.prefix, "-", potential$label, "-hist.png", sep="")
        }
        png(filename=output.filename, width=800, height=800, pointsize=18)
        hist(potential$data$RMSD, xlab="RMSD", main=potential$label)
        dev.off()
        # Concatenate some of the RMSDs at various percentile thresholds
        cat(potential$label, quantile(potential$data$RMSD, percentiles.list), "\n")
    }
    sink()
}

# Main code that triggers either the MAP assignment or learning
# operations
params <- commandArgs(trailingOnly=TRUE);
map.usage <-
paste("Usage: R --slave --vanilla --file=mrf.R --args ",
"\"operation='map'\" \"singletons=c('A.mrf', 'B.mrf', 'C.mrf')\"",
"\"pairwise=c('A-C.mrf', 'B-C.mrf')\"",
"\"weights=potential.collection.weights(CC=1,Overlap=1,PhysicsScore=1)\"",
"\"output.prefix='output-prefix'\"\n");
learning.usage <-
paste("Usage: R --slave --vanilla --file=mrf.R --args ",
"\"operation='learn'\" \"singletons=c('A.mrf', 'B.mrf', 'C.mrf')\"",
"\"pairwise=c('A-C.mrf', 'B-C.mrf')\"",
"\"features=c('CC', 'Overlap', 'PhysicsScore', 'vdw', 'elec')\"",
"\"initial.weights=c(1,1,1,0.1,0.1)\"",
"\"samples=10\" \"low.samples=3\" \"high.samples=3\"",
"\"low.rmsd.percentile=0.1\" \"high.rmsd.percentile=0.9\"",
"\"output.prefix='output-prefix'\"\n")
rmsd.histograms.usage <-
paste("Usage: R --slave --vanilla --file=mrf.R --args ",
"\"operation='rmsd-histograms'\" \"singletons=c('A.mrf', 'B.mrf', 'C.mrf')\"",
"\"output.prefix='output-prefix'\"\n")
if (length(params) == 2 || length(params) == 3 || length(params) == 5 || length(params) == 11) {
    for(param in params)
    {
        eval(parse(text=param));
    }
    if (operation == "map") {
        main.map(singletons, pairwise, weights, output.prefix, 100)
    } else if(operation == "learn") {
        main.learn(singletons, pairwise, features, initial.weights, samples, low.samples, high.samples, low.rmsd.percentile, high.rmsd.percentile, output.prefix, 10)
    } else if(operation == "rmsd-histograms") {
        if (exists("output.prefix")) {
            main.plot.singleton.rmsds(singletons, output.prefix)
        } else {
            main.plot.singleton.rmsds(singletons, NA)
        }
    } else {
        cat(paste(map.usage, learning.usage, rmsd.histograms.usage, sep="\n"))
    }
} else {
    cat(paste(map.usage, learning.usage, rmsd.histograms.usage, sep="\n"))
}
