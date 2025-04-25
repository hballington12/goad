#!/bin/bash

# ANSI color codes
CYAN='\033[0;36m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
RED='\033[0;31m'
RESET='\033[0m'

CONFIG_FILE="./config/local.toml"

clear
echo -e "${CYAN}┌────────────────────────────────────────────────────────┐${RESET}"
echo -e "${CYAN}│                                                        │${RESET}"
echo -e "${CYAN}│  ${BOLD}GOAD Configuration Editor${RESET}${CYAN}                      │${RESET}"
echo -e "${CYAN}│                                                        │${RESET}"
echo -e "${CYAN}└────────────────────────────────────────────────────────┘${RESET}"
echo

# Check if the config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo -e "${RED}Configuration file not found at $CONFIG_FILE${RESET}"
    echo -e "${YELLOW}Please run the setup script first.${RESET}"
    exit 1
fi

echo -e "${BLUE}This tool will help you customize important GOAD parameters.${RESET}"
echo -e "${YELLOW}Configuration file: $CONFIG_FILE${RESET}"
echo

# Function to update a simple scalar TOML parameter
update_param() {
    local param=$1
    local description=$2
    local current_value=$3
    local default_value=$4
    
    echo -e "${CYAN}$description${RESET}"
    echo -e "${BLUE}Current value: ${GREEN}$current_value${RESET}"
    echo -e "${BLUE}Default value: ${YELLOW}$default_value${RESET}"
    
    read -p "Enter new value (or press Enter to keep current): " new_value
    
    if [ -n "$new_value" ]; then
        sed -i "s/^$param = .*/$param = $new_value/" "$CONFIG_FILE"
        echo -e "${GREEN}Updated $param to $new_value${RESET}"
    else
        echo -e "${BLUE}Keeping current value: $current_value${RESET}"
    fi
    echo
}

# Function to update a string parameter (with quotes)
update_string_param() {
    local param=$1
    local description=$2
    local current_value=$3
    local default_value=$4
    
    echo -e "${CYAN}$description${RESET}"
    echo -e "${BLUE}Current value: ${GREEN}$current_value${RESET}"
    echo -e "${BLUE}Default value: ${YELLOW}$default_value${RESET}"
    
    read -p "Enter new value (or press Enter to keep current): " new_value
    
    if [ -n "$new_value" ]; then
        # Use a different pattern for strings that need quotes
        sed -i "s|^$param = \".*\"|$param = \"$new_value\"|" "$CONFIG_FILE"
        echo -e "${GREEN}Updated $param to \"$new_value\"${RESET}"
    else
        echo -e "${BLUE}Keeping current value: $current_value${RESET}"
    fi
    echo
}

# Function to update geometry file path with validation
update_geometry_path() {
    local param="geom_name"
    local current_value=$1
    local default_value=$2
    
    echo -e "${CYAN}Geometry file path (OBJ format)${RESET}"
    echo -e "${BLUE}Current value: ${GREEN}$current_value${RESET}"
    echo -e "${BLUE}Default value: ${YELLOW}$default_value${RESET}"
    
    while true; do
        read -p "Enter new path (or press Enter to keep current): " new_value
        
        # Keep current if empty
        if [ -z "$new_value" ]; then
            echo -e "${BLUE}Keeping current value: $current_value${RESET}"
            break
        fi
        
        # Check if file exists
        if [ -f "$new_value" ]; then
            sed -i "s|^$param = \".*\"$|$param = \"$new_value\"|" "$CONFIG_FILE"
            echo -e "${GREEN}Updated $param to \"$new_value\"${RESET}"
            break
        else
            echo -e "${RED}Warning: File '$new_value' does not exist!${RESET}"
            read -p "Try again? (y/n, n will use this path anyway): " retry
            if [[ $retry != [Yy]* ]]; then
                sed -i "s|^$param = \".*\"$|$param = \"$new_value\"|" "$CONFIG_FILE"
                echo -e "${YELLOW}Using non-existent path: \"$new_value\"${RESET}"
                echo -e "${YELLOW}Make sure to create this file before running GOAD${RESET}"
                break
            fi
        fi
    done
    echo
}

# Function to update the particle refractive index field
update_particle_refr_index() {
    local param="particle_refr_index"
    echo -e "${CYAN}Refractive index of particle (real and imaginary parts)${RESET}"
    
    # Extract current values from the file using grep/sed
    local current_real=$(sed -n '/particle_refr_index/,/]/p' "$CONFIG_FILE" | grep -Eo '[0-9.]+,' | head -1 | tr -d ',')
    local current_imag=$(sed -n '/particle_refr_index/,/]/p' "$CONFIG_FILE" | grep -Eo '[0-9.]+,' | tail -1 | tr -d ',')

    echo -e "${BLUE}Current real part: ${GREEN}$current_real${RESET}"
    echo -e "${BLUE}Current imaginary part: ${GREEN}$current_imag${RESET}"
    
    read -p "Enter new real part (or press Enter to keep current): " new_real
    read -p "Enter new imaginary part (or press Enter to keep current): " new_imag

    new_real=${new_real:-$current_real}
    new_imag=${new_imag:-$current_imag}

    # Replace the full block (assumes 1 entry in nested list)
    sed -i "/^particle_refr_index = \[/,/]/c\particle_refr_index = [\n    [\n        $new_real,\n        $new_imag,\n    ]," "$CONFIG_FILE"
    echo -e "${GREEN}Updated particle_refr_index to [[$new_real, $new_imag]]${RESET}"
    echo
}

# Function to update the medium refractive index field
update_medium_refr_index() {
    local param="medium_refr_index"
    echo -e "${CYAN}Refractive index of the medium (real and imaginary parts)${RESET}"
    
    # Extract current values
    local current_real=$(sed -n '/^medium_refr_index = \[/,/]/p' "$CONFIG_FILE" | grep -Eo '[0-9.]+,' | head -1 | tr -d ',')
    local current_imag=$(sed -n '/^medium_refr_index = \[/,/]/p' "$CONFIG_FILE" | grep -Eo '[0-9.]+' | tail -1)

    echo -e "${BLUE}Current real part: ${GREEN}$current_real${RESET}"
    echo -e "${BLUE}Current imaginary part: ${GREEN}$current_imag${RESET}"
    
    read -p "Enter new real part (or press Enter to keep current): " new_real
    read -p "Enter new imaginary part (or press Enter to keep current): " new_imag

    new_real=${new_real:-$current_real}
    new_imag=${new_imag:-$current_imag}

    # Replace the entire line
    # Replace the full block (handle multi-line array)
    sed -i "/^medium_refr_index = \[/,/]/c\medium_refr_index = [\n    $new_real,\n    $new_imag,\n]" "$CONFIG_FILE"

    echo -e "${GREEN}Updated medium_refr_index to [$new_real, $new_imag]${RESET}"
    echo
}

# Function to extract value for a parameter
extract_value() {
    local param=$1
    grep -E "^$param = " "$CONFIG_FILE" | sed -E "s/^$param = (.*)/\1/" | sed -E 's/".+"/"&"/; s/#.*$//'
}

# Extract current values
wavelength=$(extract_value "wavelength")
geom_name=$(extract_value "geom_name" | tr -d '"')
beam_power_threshold=$(extract_value "beam_power_threshold")
beam_area_threshold_fac=$(extract_value "beam_area_threshold_fac")
cutoff=$(extract_value "cutoff")
max_rec=$(extract_value "max_rec")
max_tir=$(extract_value "max_tir")
distortion=$(extract_value "distortion")
directory=$(extract_value "directory" | tr -d '"')

# Main editing flow
echo -e "${BOLD}Basic Settings:${RESET}"
update_param "wavelength" "Wavelength of incident light" "$wavelength" "0.532"
update_medium_refr_index
update_particle_refr_index
update_geometry_path "$geom_name" "./examples/data/hex.obj"

echo -e "${BOLD}Simulation Settings:${RESET}"
update_param "beam_power_threshold" "Minimum beam power threshold" "$beam_power_threshold" "0.005"
update_param "beam_area_threshold_fac" "Beam area threshold factor" "$beam_area_threshold_fac" "0.1"
update_param "cutoff" "Total power cutoff fraction" "$cutoff" "0.99"
update_param "max_rec" "Maximum number of recursions" "$max_rec" "10"
update_param "max_tir" "Maximum number of total internal reflections" "$max_tir" "10"

echo -e "${BOLD}Advanced Settings:${RESET}"
update_param "distortion" "Geometry distortion factor" "$distortion" "0.0"
update_string_param "directory" "Output directory" "$directory" "goad_run"

echo -e "${GREEN}Configuration updated successfully!${RESET}"
echo -e "${BLUE}You can manually edit $CONFIG_FILE for more detailed configuration options.${RESET}"
