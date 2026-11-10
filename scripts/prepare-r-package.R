#!/usr/bin/env Rscript
# Generate every non-handwritten file of the withinr R package.
#
# Run from the repository root:
#
#   Rscript scripts/prepare-r-package.R            # glue + docs (online cargo build)
#   Rscript scripts/prepare-r-package.R --vendor   # ... plus offline vendored deps (CRAN)
#   Rscript scripts/prepare-r-package.R --no-document   # glue only, no compile
#   Rscript scripts/prepare-r-package.R --clean    # delete everything generated
#   Rscript scripts/prepare-r-package.R --check    # fail if the core snapshot is stale
#
# `--no-document` leaves rextendr's placeholder wrappers and no NAMESPACE: enough
# for `cargo check`/`clippy` on the adapter, not for loading the package.
#
# Handwritten sources are `withinr/R/*.R` (except extendr-wrappers.R),
# `withinr/src/rust/{Cargo.toml,Cargo.lock,src/}`, DESCRIPTION, tests, and docs.
# Everything listed in `generated` below is produced here and ignored by Git.
if (!file.exists("Cargo.toml") || !dir.exists("withinr")) {
  stop("Run this script from the within repository root.")
}
args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% c("--vendor", "--no-document", "--clean", "--check"))) {
  stop("Unknown argument; use --vendor, --no-document, --clean, or --check.")
}
if (any(c("--clean", "--check") %in% args) && length(args) > 1L) {
  stop("--clean and --check cannot be combined with other arguments.")
}
root <- normalizePath(".", winslash = "/")
pkg <- file.path(root, "withinr")

# extendr/rextendr scaffolding: FFI entrypoint, build templates, configure
# scripts, the Rust `document` binary, and the R `.Call` wrappers it writes.
scaffolding <- c(
  "src/entrypoint.c", "src/Makevars.in", "src/Makevars.win.in", "src/withinr-win.def",
  "src/rust/document.rs", "tools/msrv.R", "tools/config.R",
  "configure", "configure.win", "cleanup", "cleanup.win"
)
vendored <- c("src/rust/vendor.tar.xz", "src/rust/vendor-config.toml", "inst/COPYRIGHTS")
documentation <- c("NAMESPACE", "man")
# Build leftovers; a stale extracted `src/vendor` would shadow the tarball.
leftovers <- c(
  "src/vendor", "src/rust/vendor", "src/.cargo", "src/Makevars", "src/Makevars.win",
  file.path("src", list.files(file.path(pkg, "src"), "\\.(o|so|dll|lib)$"))
)
generated <- c(scaffolding, "R/extendr-wrappers.R", "src/rust/within-core", leftovers)

# Copy the workspace's Rust core (runtime sources only) into `dest`, so the
# source package builds without this repository or a crates.io release of
# within. The small nested workspace keeps inherited manifest fields.
write_core_snapshot <- function(dest) {
  dir.create(file.path(dest, "crates"), recursive = TRUE)
  manifest <- readLines(file.path(root, "Cargo.toml"))
  python <- grepl("crates/within-py", manifest, fixed = TRUE)
  if (sum(python) != 1L || !grepl('^\\s*"crates/within-py",\\s*$', manifest[python])) {
    stop("Expected one `\"crates/within-py\",` line in Cargo.toml's workspace members; ",
         "update write_core_snapshot() in this script.")
  }
  writeLines(manifest[!python], file.path(dest, "Cargo.toml"))
  stopifnot(file.copy(file.path(root, "LICENSE"), dest, overwrite = TRUE))
  for (crate in c("within", "schwarz-precond")) {
    from <- file.path(root, "crates", crate)
    to <- file.path(dest, "crates", crate)
    dir.create(to)
    # Drop test/benchmark manifest sections.
    manifest <- readLines(file.path(from, "Cargo.toml"))
    manifest <- sub("^\\[package\\]$", '[package]\nworkspace = "../.."', manifest)
    section <- cumsum(grepl("^\\[", manifest))
    headers <- manifest[grepl("^\\[", manifest)]
    drop <- which(grepl("^\\[dev-dependencies|^\\[\\[bench", headers))
    writeLines(manifest[!section %in% drop], file.path(to, "Cargo.toml"))
    stopifnot(file.copy(file.path(from, "src"), to, recursive = TRUE))
    stopifnot(file.copy(file.path(from, "README.md"), to))
  }
}

snapshot <- file.path(pkg, "src/rust/within-core")

if ("--check" %in% args) {
  # The snapshot is a copy: edits to crates/ reach the package only after this
  # script reruns. Compare it with a fresh snapshot, file by file.
  fresh <- tempfile("within-core-")
  write_core_snapshot(fresh)
  digest <- function(dir) {
    files <- sort(list.files(dir, recursive = TRUE, all.files = TRUE))
    stats::setNames(unname(tools::md5sum(file.path(dir, files))), files)
  }
  current <- if (dir.exists(snapshot)) digest(snapshot) else character()
  expected <- digest(fresh)
  unlink(fresh, recursive = TRUE)
  files <- union(names(current), names(expected))
  changed <- files[is.na(current[files]) | is.na(expected[files]) |
                     current[files] != expected[files]]
  if (length(changed)) {
    stop("withinr/src/rust/within-core is stale (", length(changed), " file(s) differ, e.g. ",
         paste(utils::head(changed, 3), collapse = ", "),
         "); rerun Rscript scripts/prepare-r-package.R")
  }
  message("withinr's core snapshot matches the workspace.")
  quit(save = "no")
}

if ("--clean" %in% args) {
  unlink(file.path(pkg, c(generated, vendored, documentation, "src/rust/target")),
         recursive = TRUE)
  for (dir in c("tools", "inst")) {
    path <- file.path(pkg, dir)
    if (dir.exists(path) && !length(dir(path, all.files = TRUE, no.. = TRUE))) {
      unlink(path, recursive = TRUE)
    }
  }
  message("Removed generated withinr files.")
  quit(save = "no")
}

for (tool in c("rextendr", "roxygen2", "pkgload")) {
  if (!requireNamespace(tool, quietly = TRUE)) stop("Install the R package '", tool, "' first.")
}
if (packageVersion("rextendr") < "0.5.0") stop("rextendr >= 0.5.0 is required.")
if (!nzchar(Sys.which("cargo"))) stop("cargo is not on PATH; install Rust (>= 1.85).")

# Start from a clean slate so no stale artifact can leak into the build. The
# Cargo target directory is left alone: pkgload's debug builds (`load_all()`)
# reuse `target/debug` incrementally, while the generated Makevars deletes the
# target directory after every release build (`R CMD INSTALL`) unless DEBUG is
# set. Without --vendor a stale vendor archive is removed too, and cargo
# resolves dependencies online.
unlink(file.path(pkg, c(generated, if (!"--no-document" %in% args) documentation,
                        if (!"--vendor" %in% args) vendored)), recursive = TRUE)

# 1. Snapshot the workspace's Rust core into the package.
write_core_snapshot(snapshot)

# 2. Generate the extendr scaffolding in a disposable package, so
#    use_extendr(overwrite = TRUE) cannot touch the handwritten bridge.
scaffold <- tempfile("withinr-scaffold-")
dir.create(scaffold)
stopifnot(file.copy(file.path(pkg, "DESCRIPTION"), scaffold))
rextendr::use_extendr(scaffold, crate_name = "withinr", lib_name = "withinr",
                      overwrite = TRUE, quiet = TRUE)
for (name in c(scaffolding, "R/extendr-wrappers.R")) {
  dest <- file.path(pkg, name)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(file.path(scaffold, name), dest, overwrite = TRUE))
}
Sys.chmod(file.path(pkg, c("configure", "configure.win", "cleanup", "cleanup.win")), "0755")
unlink(scaffold, recursive = TRUE)

# Policy additions to rextendr's build templates: locked dependency resolution
# and at most two compiler jobs (CRAN), even in development builds.
for (name in c("src/Makevars.in", "src/Makevars.win.in")) {
  dest <- file.path(pkg, name)
  lines <- readLines(dest)
  patched <- gsub("cargo (build|run) @CRAN_FLAGS@",
                  "CARGO_BUILD_JOBS=2 cargo \\1 --locked @CRAN_FLAGS@", lines)
  if (sum(patched != lines) != 2L) {
    stop("The rextendr template in ", name, " changed; update the patch in this script.")
  }
  writeLines(patched, dest)
}

# 3. Optionally vendor the locked crates.io dependencies for offline builds.
if ("--vendor" %in% args) {
  rextendr::vendor_crates(pkg, clean = FALSE)
  # Some WASI-only transitive dependencies ship precompiled WebAssembly/object
  # fixtures. CRAN source packages cannot contain executable binaries, and
  # native R builds never use them: drop them and update Cargo's checksums.
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
  # Attribution for every bundled crate; original notices stay in the archive.
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

# 4. Compile, write the real `.Call` wrappers, and document.
if (!"--no-document" %in% args) {
  # Bootstrap NAMESPACE (useDynLib) from source without loading a DLL.
  roxygen2::roxygenise(pkg, load_code = roxygen2::load_source)
  # pkgload compiles the crate and runs its `document` binary, which writes
  # R/extendr-wrappers.R from extendr's metadata; roxygen then writes the
  # final NAMESPACE and man/ pages.
  roxygen2::roxygenise(pkg)
}
message("withinr sources prepared in ", pkg)
