#!/bin/bash

if [[ $# -ne 1 ]]
then
  echo "Usage: $0 /path/to/diagram.mmd"
  exit 1
fi

set -e

MMD="$1"
NAME="${MMD//.mmd}"
render () {
  echo "Rendering $1"
  nf-metro render "${MMD}" -o "${NAME}_$1.svg" --theme $2 --format svg
  nf-metro render "${MMD}" -o "${NAME}_$1_animated.svg" --theme $2 --format svg --animate
  nf-metro render "${MMD}" -o "${NAME}_$1.png" --theme $2 --format png
}

render dark nfcore-dark
render light nfcore-light
