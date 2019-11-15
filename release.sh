#!/bin/bash -e

set -o pipefail

export LIBCLANG_PATH=~/.guix-profile/lib

export VERSION=`cargo run -- --version |awk '{print $2}'`

echo "Found version $VERSION .."

echo "Building normally .."
cargo build --release

echo "Testing release version .."
cargo test --release

echo "Building musl static binary .."
cargo build --target x86_64-unknown-linux-musl --release

echo "Making static dist .."
mkdir dist/coverm-x86_64-unknown-linux-musl-$VERSION
cp target/x86_64-unknown-linux-musl/release/coverm dist/coverm-x86_64-unknown-linux-musl-$VERSION/
cd dist
tar czf coverm-x86_64-unknown-linux-musl-$VERSION.tar.gz coverm-x86_64-unknown-linux-musl-$VERSION
cd ..

echo "Now make sure git is up to date, and run LIBCLANG_PATH=~/.guix-profile/lib cargo publish"