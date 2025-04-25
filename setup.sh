#!/bin/bash

# ANSI color codes
CYAN='\033[0;36m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
RED='\033[0;31m'
RESET='\033[0m'

clear
echo -e "${CYAN}┌────────────────────────────────────────────────────────┐${RESET}"
echo -e "${CYAN}│                                                        │${RESET}"
echo -e "${CYAN}│  ${BOLD}GOAD - Geometric Optics with Aperture Diffraction${RESET}${CYAN}     │${RESET}"
echo -e "${CYAN}│                                                        │${RESET}"
echo -e "${CYAN}└────────────────────────────────────────────────────────┘${RESET}"
echo

echo -e "${GREEN}Welcome to the GOAD interactive setup!${RESET}"
echo -e "${BLUE}This tool helps you simulate light scattering using geometric optics${RESET}"
echo -e "${BLUE}with aperture diffraction for accurate Mueller matrix computation.${RESET}"
echo

echo -e "${YELLOW}What GOAD can do:${RESET}"
echo -e " • Calculate full Mueller matrix elements"
echo -e " • Handle complex geometries including concavities and layered media"
echo -e " • Compute scattering patterns in fixed or random orientations"
echo

# Check if Cargo is installed and accessible
echo -e "${CYAN}Checking for Rust and Cargo installation...${RESET}"

if command -v cargo >/dev/null 2>&1; then
    echo -e "${GREEN}✓ Cargo is installed and available in your PATH${RESET}"
    cargo_version=$(cargo --version)
    echo -e "  ${BLUE}$cargo_version${RESET}"
    
    # Check if GOAD is already compiled
    echo
    echo -e "${CYAN}Checking for GOAD executable...${RESET}"
    
    if [ -f "./target/release/goad" ]; then
        echo -e "${GREEN}✓ GOAD is already compiled${RESET}"
        
        echo -e "${YELLOW}Would you like to recompile GOAD to ensure you have the latest version?${RESET}"
        read -p "Recompile? (y/n) " -n 1 -r
        echo
        
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            echo -e "${CYAN}Recompiling GOAD...${RESET}"
            cargo build --release
            echo -e "${GREEN}✓ GOAD has been recompiled successfully${RESET}"
        fi
    else
        echo -e "${YELLOW}⚠ GOAD executable not found${RESET}"
        echo -e "${YELLOW}Would you like to compile GOAD now?${RESET}"
        read -p "Compile? (y/n) " -n 1 -r
        echo
        
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            echo -e "${CYAN}Compiling GOAD...${RESET}"
            echo -e "${BLUE}This may take a few minutes for the first build.${RESET}"
            cargo build --release
            
            # Check if compilation was successful
            if [ -f "./target/release/goad" ]; then
                echo -e "${GREEN}✓ GOAD has been compiled successfully${RESET}"
            else
                echo -e "${RED}✗ Compilation failed. Please check for errors above.${RESET}"
                exit 1
            fi
        else
            echo -e "${YELLOW}Skipping compilation. You'll need to compile GOAD before using it.${RESET}"
        fi
    fi
else
    echo -e "${YELLOW}⚠ Cargo is not found in your PATH${RESET}"
    
    # Check common locations where Cargo might be installed
    COMMON_CARGO_LOCATIONS=(
        "$HOME/.cargo/bin/cargo"
        "/usr/local/bin/cargo"
        "/usr/bin/cargo"
        "/opt/cargo/bin/cargo"
    )
    
    cargo_found=false
    cargo_location=""
    for location in "${COMMON_CARGO_LOCATIONS[@]}"; do
        if [ -x "$location" ]; then
            cargo_found=true
            cargo_location=$(dirname "$location")
            echo -e "${YELLOW}Cargo was found at:${RESET} $location"
            echo -e "${YELLOW}But it's not in your PATH.${RESET}"
            
            # Detect user's shell
            current_shell=$(basename "$SHELL")
            echo -e "${BLUE}Detected shell:${RESET} $current_shell"
            
            case "$current_shell" in
                bash)
                    echo -e "${YELLOW}To add Cargo to your PATH in bash, run:${RESET}"
                    echo -e "  echo 'export PATH=\"$cargo_location:\$PATH\"' >> ~/.bashrc"
                    echo -e "  source ~/.bashrc"
                    ;;
                zsh)
                    echo -e "${YELLOW}To add Cargo to your PATH in zsh, run:${RESET}"
                    echo -e "  echo 'export PATH=\"$cargo_location:\$PATH\"' >> ~/.zshrc"
                    echo -e "  source ~/.zshrc"
                    ;;
                tcsh|csh)
                    echo -e "${YELLOW}To add Cargo to your PATH in tcsh/csh, run:${RESET}"
                    echo -e "  echo 'set path = ($cargo_location \$path)' >> ~/.tcshrc"
                    echo -e "  source ~/.tcshrc"
                    ;;
                fish)
                    echo -e "${YELLOW}To add Cargo to your PATH in fish, run:${RESET}"
                    echo -e "  fish_add_path $cargo_location"
                    echo -e "  # Or for older fish versions:"
                    echo -e "  # set -U fish_user_paths $cargo_location \$fish_user_paths"
                    ;;
                *)
                    echo -e "${YELLOW}To add Cargo to your PATH, add this line to your shell's config file:${RESET}"
                    echo -e "  export PATH=\"$cargo_location:\$PATH\""
                    echo -e "${YELLOW}Then restart your terminal or source your config file.${RESET}"
                    ;;
            esac
            
            # Ask if user wants to add Cargo to PATH automatically
            echo
            read -p "Would you like to try adding Cargo to PATH automatically? (y/n) " -n 1 -r
            echo
            if [[ $REPLY =~ ^[Yy]$ ]]; then
                case "$current_shell" in
                    bash)
                        echo 'export PATH="'"$cargo_location"':$PATH"' >> ~/.bashrc
                        echo -e "${GREEN}Added Cargo to PATH in ~/.bashrc${RESET}"
                        echo -e "${YELLOW}To apply changes in current session, run: source ~/.bashrc${RESET}"
                        ;;
                    zsh)
                        echo 'export PATH="'"$cargo_location"':$PATH"' >> ~/.zshrc
                        echo -e "${GREEN}Added Cargo to PATH in ~/.zshrc${RESET}"
                        echo -e "${YELLOW}To apply changes in current session, run: source ~/.zshrc${RESET}"
                        ;;
                    tcsh|csh)
                        echo "set path = ($cargo_location \$path)" >> ~/.tcshrc
                        echo -e "${GREEN}Added Cargo to PATH in ~/.tcshrc${RESET}"
                        echo -e "${YELLOW}To apply changes in current session, run: source ~/.tcshrc${RESET}"
                        ;;
                    fish)
                        if type -q fish_add_path 2>/dev/null; then
                            fish -c "fish_add_path $cargo_location"
                            echo -e "${GREEN}Added Cargo to PATH using fish_add_path${RESET}"
                        else
                            fish -c "set -U fish_user_paths $cargo_location \$fish_user_paths"
                            echo -e "${GREEN}Added Cargo to PATH using fish_user_paths${RESET}"
                        fi
                        ;;
                    *)
                        echo -e "${YELLOW}Unable to automatically add to PATH for your shell.${RESET}"
                        echo -e "${YELLOW}Please follow the manual instructions above.${RESET}"
                        ;;
                esac
                
                echo -e "${YELLOW}Please restart this script after restarting your terminal or sourcing your config.${RESET}"
                exit 0
            fi
            break
        fi
    done
    
    if [ "$cargo_found" = false ]; then
        echo -e "${RED}✗ Cargo is not installed on this system${RESET}"
        echo -e "${YELLOW}To install Rust and Cargo, run the following command:${RESET}"
        echo -e "  ${GREEN}curl https://sh.rustup.rs -sSf | sh${RESET}"
        echo -e "${YELLOW}After installation, restart this setup script.${RESET}"
        echo
        read -p "Do you want to install Rust and Cargo now? (y/n) " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            echo -e "${CYAN}Installing Rust and Cargo...${RESET}"
            curl https://sh.rustup.rs -sSf | sh
            echo -e "${GREEN}Please restart your terminal and run this setup script again.${RESET}"
            exit 0
        fi
    fi
fi

echo

echo -e "${BOLD}We'll help you get started with running your first simulation.${RESET}"
echo -e "${BLUE}Let's set up your configuration for GOAD.${RESET}"

# Check if a local config file exists
echo
echo -e "${CYAN}Checking for local configuration file...${RESET}"

if [ -f "./config/local.toml" ]; then
    echo -e "${GREEN}✓ Local configuration file already exists at ./config/local.toml${RESET}"
    echo -e "${YELLOW}Would you like to reset it by creating a new one from the default?${RESET}"
    read -p "Create new local.toml? (y/n) " -n 1 -r
    echo
    
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        # Create backup of existing local.toml
        backup_name="./config/local.toml.backup.$(date +%Y%m%d%H%M%S)"
        echo -e "${BLUE}Creating backup of existing local.toml as ${backup_name}${RESET}"
        cp "./config/local.toml" "${backup_name}"
        
        # Copy default config to local
        echo -e "${CYAN}Creating new local configuration file...${RESET}"
        cp "./config/default.toml" "./config/local.toml"
        echo -e "${GREEN}✓ New local configuration created at ./config/local.toml${RESET}"
        echo -e "${BLUE}Your previous configuration was backed up to ${backup_name}${RESET}"
    fi
else
    echo -e "${YELLOW}⚠ No local configuration file found${RESET}"
    echo -e "${YELLOW}Would you like to create one from the default template? (Recommended)${RESET}"
    read -p "Create local.toml? (y/n) " -n 1 -r
    echo
    
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        # Check if the config directory exists
        if [ ! -d "./config" ]; then
            echo -e "${YELLOW}⚠ Config directory not found. Creating it...${RESET}"
            mkdir -p ./config
        fi
        
        # Check if the default config exists
        if [ -f "./config/default.toml" ]; then
            echo -e "${CYAN}Creating local configuration file...${RESET}"
            cp "./config/default.toml" "./config/local.toml"
            echo -e "${GREEN}✓ Local configuration created at ./config/local.toml${RESET}"
            echo -e "${BLUE}You can customize this file to change simulation parameters.${RESET}"
        else
            echo -e "${RED}✗ Default configuration file not found at ./config/default.toml${RESET}"
            echo -e "${YELLOW}Creating an empty local.toml file. You'll need to populate it manually.${RESET}"
            touch "./config/local.toml"
        fi
    else
        echo -e "${YELLOW}Skipping local configuration creation. GOAD will use default settings.${RESET}"
    fi
fi

# Ask if user wants to customize the configuration
if [ -f "./config/local.toml" ]; then
    echo
    echo -e "${YELLOW}Would you like to customize the configuration file now? (Recommended)${RESET}"
    read -p "Customize configuration? (y/n) " -n 1 -r
    echo
    
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        # Check if the config editor script exists
        if [ -f "./config_editor.sh" ]; then
            echo -e "${CYAN}Launching configuration editor...${RESET}"
            chmod +x ./config_editor.sh
            ./config_editor.sh
        else
            # Create a basic config editor script if it doesn't exist
            echo -e "${YELLOW}⚠ Configuration editor not found. Creating a simple editor script...${RESET}"
            
            cat > ./config_editor.sh << 'EOF'
#!/bin/bash
# Simple configuration editor for GOAD
echo "This is a placeholder for the configuration editor."
echo "Edit this script to add interactive configuration options."
echo "For now, please manually edit ./config/local.toml with your text editor."
EOF
            
            chmod +x ./config_editor.sh
            echo -e "${CYAN}Launching simple configuration editor...${RESET}"
            ./config_editor.sh
            echo -e "${BLUE}For more options, please edit ./config/local.toml directly with your preferred text editor.${RESET}"
        fi
    else
        echo -e "${BLUE}Skipping configuration customization. You can manually edit ./config/local.toml anytime.${RESET}"
        echo -e "${BLUE}You can also place a copy of local.toml in another working directory, which will be prioritised over the copy in ./config.${RESET}"
    fi
fi

# More detailed setup steps will be added here in future updates
