#!/usr/bin/env Rscript
# Run from the repository root. Generated files are reproducible from this script,
# the Rust adapter, and the pinned rextendr release; never edit FFI glue manually.
if (!file.exists("Cargo.toml") || !dir.exists("withinr")) {
  stop("Run this script from the within repository root.")
}
if (!requireNamespace("rextendr", quietly = TRUE) ||
    packageVersion("rextendr") != "0.5.0") {
  stop("Install rextendr 0.5.0 to reproduce the committed build scaffolding.")
}
if (!requireNamespace("roxygen2", quietly = TRUE)) stop("Install roxygen2 first.")
root <- normalizePath(".", winslash = "/")
pkg <- file.path(root, "withinr")
args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% c("--vendor", "--no-document"))) stop("Unknown argument.")

# Keep local path dependencies inside the source package. The small nested
# workspace retains inherited manifest fields without rewriting dependencies or
# substituting the (potentially older) crates.io release of within.
snapshot <- file.path(pkg, "src/rust/within-core")
unlink(snapshot, recursive = TRUE)
dir.create(file.path(snapshot, "crates"), recursive = TRUE)
manifest <- readLines(file.path(root, "Cargo.toml"))
manifest <- manifest[!grepl('"crates/within-py"', manifest, fixed = TRUE)]
writeLines(manifest, file.path(snapshot, "Cargo.toml"))
stopifnot(file.copy(file.path(root, "LICENSE"), snapshot, overwrite = TRUE))
for (crate in c("within", "schwarz-precond")) {
  from <- file.path(root, "crates", crate)
  to <- file.path(snapshot, "crates", crate)
  dir.create(to)
  # Runtime source only. Strip test/benchmark manifests from the packaged copy;
  # the repository remains the sole source of truth for these crates.
  manifest <- readLines(file.path(from, "Cargo.toml"))
  manifest <- sub("^\\[package\\]$", '[package]\nworkspace = "../.."', manifest)
  section <- cumsum(grepl("^\\[", manifest))
  headers <- manifest[grepl("^\\[", manifest)]
  drop <- which(grepl("^\\[dev-dependencies|^\\[\\[bench", headers))
  writeLines(manifest[!section %in% drop], file.path(to, "Cargo.toml"))
  stopifnot(file.copy(file.path(from, "src"), to, recursive = TRUE))
  stopifnot(file.copy(file.path(from, "README.md"), to))
}

# Generate in a disposable package so use_extendr(overwrite = TRUE) cannot
# overwrite our handwritten API adapter or its dependency manifest.
scaffold <- tempfile("withinr-scaffold-")
dir.create(scaffold)
stopifnot(file.copy(file.path(pkg, "DESCRIPTION"), scaffold))
rextendr::use_extendr(scaffold, crate_name = "withinr", lib_name = "withinr",
                      overwrite = TRUE, quiet = TRUE)
generated <- c("src/entrypoint.c", "src/Makevars.in", "src/Makevars.win.in",
               "src/withinr-win.def", "src/rust/document.rs", "tools/msrv.R",
               "tools/config.R", "configure", "configure.win", "cleanup", "cleanup.win")
for (name in generated) {
  dest <- file.path(pkg, name)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(file.path(scaffold, name), dest, overwrite = TRUE))
}
Sys.chmod(file.path(pkg, c("configure", "configure.win", "cleanup", "cleanup.win")), "0755")
if (!file.exists(file.path(pkg, "R/extendr-wrappers.R"))) {
  stopifnot(file.copy(file.path(scaffold, "R/extendr-wrappers.R"), file.path(pkg, "R")))
}
unlink(scaffold, recursive = TRUE)

# Policy additions to rextendr's build templates, applied reproducibly here.
# Keep dependency resolution locked and compiler concurrency bounded even in
# development builds; the FFI and registration code remain helper-generated.
for (name in c("src/Makevars.in", "src/Makevars.win.in")) {
  dest <- file.path(pkg, name)
  lines <- readLines(dest)
  lines <- gsub("cargo build @CRAN_FLAGS@", "CARGO_BUILD_JOBS=2 cargo build --locked @CRAN_FLAGS@", lines, fixed = TRUE)
  lines <- gsub("cargo run @CRAN_FLAGS@", "CARGO_BUILD_JOBS=2 cargo run --locked @CRAN_FLAGS@", lines, fixed = TRUE)
  writeLines(lines, dest)
}

if ("--vendor" %in% args) {
  rextendr::vendor_crates(pkg, clean = FALSE)
  # Some WASI-only transitive dependencies ship precompiled WebAssembly/object
  # fixtures. CRAN source packages cannot contain executable binaries. Native R
  # builds do not use these files; retain the crate sources and licenses and
  # update Cargo's per-file checksums after removing the binary artifacts.
  vendor <- file.path(pkg, "src/rust/vendor")
  for (crate in list.dirs(vendor, recursive = FALSE)) {
    files <- list.files(crate, recursive = TRUE, all.files = TRUE)
    binaries <- files[grepl("\\.(wasm|o|a|so|dll|dylib|exe)$", files)]
    if (length(binaries)) {
      checksum_path <- file.path(crate, ".cargo-checksum.json")
      checksum <- jsonlite::read_json(checksum_path)
      for (name in binaries) {
        unlink(file.path(crate, name))
        checksum$files[[name]] <- NULL
      }
      jsonlite::write_json(checksum, checksum_path, auto_unbox = TRUE)
    }
  }
  processx::run("tar", c("-cJ", "--no-xattrs", "-f", "vendor.tar.xz", "vendor"),
                wd = file.path(pkg, "src/rust"))
  unlink(vendor, recursive = TRUE)
  # Include attribution for every bundled dependency, including target-specific
  # crates. Original copyright notices and license files remain in the archive.
  metadata <- processx::run("cargo", c("metadata", "--locked", "--offline",
    "--format-version=1", "--manifest-path", file.path(pkg, "src/rust/Cargo.toml")))
  crates <- jsonlite::fromJSON(metadata$stdout, simplifyVector = FALSE)$packages
  crates <- crates[order(vapply(crates, function(x) x$name, character(1)))]
  lines <- c("Bundled Rust source attribution (generated; do not edit).",
    "Original license and copyright notices are retained in src/rust/vendor.tar.xz",
    "and src/rust/within-core. The source manifest of each crate is authoritative.",
    "Precompiled WASI objects/fixtures are omitted. This package targets native R.", "")
  for (crate in crates) {
    lines <- c(lines, paste(crate$name, crate$version),
      paste("Authors:", paste(unlist(crate$authors), collapse = "; ")),
      paste("License:", crate$license), paste("Repository:", crate$repository), "")
  }
  dir.create(file.path(pkg, "inst"), showWarnings = FALSE)
  writeLines(lines, file.path(pkg, "inst/COPYRIGHTS"), useBytes = TRUE)
}
if (!"--no-document" %in% args) {
  # Bootstrap NAMESPACE without loading a DLL when every generated file has
  # been deleted. The placeholder wrappers above come from use_extendr itself.
  roxygen2::roxygenise(pkg, load_code = roxygen2::load_source)
  # pkgload compiles the library and runs the rextendr-generated document binary;
  # roxygen then creates NAMESPACE and Rd pages from R source and those wrappers.
  roxygen2::roxygenise(pkg)
}
message("R package sources prepared in ", pkg)
