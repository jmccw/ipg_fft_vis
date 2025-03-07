#!/bin/bash

# Bif diagram plotter for dynamics models
# Jordan Walsh Integrated Photonics

# Set the range and step size for the parameters
param1_start=-2.5
param1_end=2.5
param1_step=0.01

param2_start=0.0
param2_end=2.0
param2_steps=100

# Set constants
gamma_=0.65 # qw laser on the order of ~0.01
alpha=2.0
P=0.5

# Set upper bounds on variable initial conditions (these are randomised between 0 and X_upper to capture more information ~ at a glance)
N0_upper=1.0
R0_upper=1.0
phi0_upper=1.0

# ------------------------------- #

silence_outputs=1
output_file="output.dat"

# Calculate the total number of iterations
param1_steps=$(awk "BEGIN { print (($param1_end - $param1_start) / $param1_step) + 1 }")
total_iterations=$(awk "BEGIN { print $param1_steps }")

# Function to display the progress bar
show_progress() {
  local current=$1
  local total=$2
  local width=50  # Width of the progress bar
  local progress=$((current * width / total))
  local remaining=$((width - progress))

  # Build the progress bar string
  printf "\r["
  printf "%0.s#" $(seq 1 $progress)
  printf "%0.s " $(seq 1 $remaining)
  printf "] %d%% (%d/%d)" $((current * 100 / total)) $current $total
}

iteration=0  # Initialize iteration counter

# clear previous data
> $output_file

# Run classB over the range of parameter values
for param1 in $(seq $param1_start $param1_step $param1_end); do
  # Run classB with the current parameters
  ./main_threaded $param2_start $param2_end $param2_steps $param1 $gamma_ $alpha $R0_upper $phi0_upper $N0_upper $P $silence_outputs

  # Update progress
  iteration=$((iteration + 1))
  show_progress $iteration $total_iterations
done

gnuplot -p plot2.gp

printf '\nFinished!\n'
