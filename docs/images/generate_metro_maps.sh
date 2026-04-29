#!/bin/sh

## render static maps
nf-metro render genomeassembly.mmd \
    -o sanger-tol-genomeassembly_metro_map_dark.svg \
    --logo sanger-tol-genomeassembly_logo_dark.png
nf-metro render genomeassembly.mmd \
    -o sanger-tol-genomeassembly_metro_map_light.svg \
    --theme light \
    --logo sanger-tol-genomeassembly_logo_light.png

## render animated maps
nf-metro render genomeassembly.mmd \
    --animate \
    -o sanger-tol-genomeassembly_metro_map_dark_animated.svg \
    --logo sanger-tol-genomeassembly_logo_dark.png
nf-metro render genomeassembly.mmd \
    --animate \
    -o sanger-tol-genomeassembly_metro_map_light_animated.svg \
    --theme light \
    --logo sanger-tol-genomeassembly_logo_light.png
