#!/bin/bash

echo "Building Apptainer images..."

# --- Function to build Apptainer images ---
build_apptainer_image() {
  container_sif="$1"
  container_def="$2"
  algo_name="$3"
  force_rebuild="$4"

  # Check if the .sif image already exists
  if [ -f "$container_sif" ]; then
    if [ "$force_rebuild" = true ]; then
      echo "Forcing rebuild of $algo_name image..."
      apptainer build --force "$container_sif" "$container_def"
    else
      # Prompt for force rebuild
      read -p "A .sif image for $algo_name already exists. Force rebuild? (y/N) " -n 1 -r
      echo    # Move to a new line after input
      if [[ $REPLY =~ ^[Yy]$ ]]; then
        # Force rebuild
        apptainer build --force "$container_sif" "$container_def"
      else
        echo "Skipping rebuild for $algo_name image."
      fi
    fi
  else
    # Build the image (no existing image found)
    echo "Building $algo_name image..."
    apptainer build "$container_sif" "$container_def"
  fi
}

# Get the optional algorithm argument and --force flag
requested_algo=""
force_rebuild=false

for arg in "$@"; do
  if [ "$arg" == "--force" ]; then
    force_rebuild=true
  elif [ -z "$requested_algo" ]; then
    requested_algo="$arg"
  fi
done

# --- Build algorithm images ---
# Loop through each directory (algorithm) in the "algorithms" folder
for algo_name in algorithms/*; do
  # Check if it's a directory and not the "base" folder
  if [ -d "$algo_name" ] && [ "${algo_name##*/}" != "base" ]; then 
    current_algo="${algo_name##*/}"

    # If an algorithm is provided, only build that one
    if [ -z "$requested_algo" ] || [ "$requested_algo" == "$current_algo" ]; then
      # Construct the full paths for the command
      container_sif="algorithms/$current_algo/container.sif"
      container_def="algorithms/$current_algo/container.def"

      # Build the algorithm image
      build_apptainer_image "$container_sif" "$container_def" "$current_algo" "$force_rebuild"
    fi
  fi
done

# --- Build evaluation image if no algorithm was provided ---
if [ -z "$requested_algo" ] || [ "$requested_algo" == "evaluation" ]; then
  # Construct paths for the evaluation image
  evaluation_sif="evaluation.sif"
  evaluation_def="evaluation.def"

  # Build the evaluation image
  build_apptainer_image "$evaluation_sif" "$evaluation_def" "evaluation" "$force_rebuild"
fi

echo "Building Apptainer images finished."