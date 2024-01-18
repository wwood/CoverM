#!/bin/bash -e

set -o pipefail

echo "Building normally .."
cargo build --release

export VERSION=`cargo run -- --version |awk '{print $2}'`

# For minimap header fix binary
export PATH=target/release:$PATH

echo "Found version $VERSION .."



echo "Testing release version .."
cargo test --release

echo "Building musl static binary .."
# Use cross not cargo here so htslib is OK - see https://github.com/rust-bio/rust-htslib
cross build --target x86_64-unknown-linux-musl --release

echo "Making static dist .."
mkdir -p dist/coverm-x86_64-unknown-linux-musl-$VERSION
cp \
 target/x86_64-unknown-linux-musl/release/coverm \
 INSTALL.md \
 dist/coverm-x86_64-unknown-linux-musl-$VERSION/
cd dist
tar czf coverm-x86_64-unknown-linux-musl-$VERSION.tar.gz coverm-x86_64-unknown-linux-musl-$VERSION
cd ..

echo "Building HTML versions of man pages .."
for SUBCOMMAND in genome cluster contig filter make
do
    echo "Documenting $SUBCOMMAND .."
    cargo run -- $SUBCOMMAND --full-help-roff |pandoc - -t markdown -f man |sed 's/\\\[/[/g; s/\\\]/]/g' |cat <(sed s/SUBCOMMAND/$SUBCOMMAND/ prelude) - >docs/coverm-$SUBCOMMAND.Rmd
    echo "library(prettydoc); setwd('docs'); rmarkdown::render('coverm-$SUBCOMMAND.Rmd','prettydoc::html_pretty','coverm-$SUBCOMMAND.html')" |R --no-save
    rm docs/coverm-$SUBCOMMAND.Rmd
    echo "Finished documenting $SUBCOMMAND"
done


echo "Now make sure git is up to date, the documentation HTML has been properly generated and run cargo publish"
