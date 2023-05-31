
carstm_saveRDS = function (object, file = "", ascii = FALSE, version = NULL, compress = TRUE, compression_level = 3, refhook = NULL)  {
    # copied from base and modified to allow other compression_level
    if (is.character(file)) {
        if (file == "") 
            stop("'file' must be non-empty string")
        object <- object
        mode <- ifelse(ascii %in% FALSE, "wb", "w")
        if (is.logical(compress)) { 
            if (compress) {
                con <- gzfile(file, mode, compression = compression_level)
            } else {
                con <- file(file, mode)
            }
        } else {
                con <- switch(compress, 
                    bzip2 = bzfile(file, mode, compression = compression_level), 
                    xz = xzfile(file, mode, compression = compression_level), 
                    gzip = gzfile(file, mode, compression = compression_level), 
                    stop("invalid 'compress' argument: ", compress)
                )
        }
        on.exit(close(con))
    } else if (inherits(file, "connection")) {
        if (!missing(compress))  {
            warning("'compress' is ignored unless 'file' is a file name")
        }
        con <- file
    } else {
        stop("bad 'file' argument")
    }
    .Internal(serializeToConn(object, con, ascii, version, refhook))
}
