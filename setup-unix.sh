#!/usr/bin/env bash
set -e

echo "=== TPMS setup (Linux / macOS) ==="

# --- Detect OS and architecture ---
OS="$(uname -s)"
ARCH="$(uname -m)"

# --- Check if conda already exists ---
if command -v conda >/dev/null 2>&1; then
    echo "Found existing conda installation."
    CONDA_AVAILABLE=true
else
    echo "Conda not found. Will install Miniconda locally."
    CONDA_AVAILABLE=false
fi

# --- Install Miniconda only if needed ---
if [[ "$CONDA_AVAILABLE" == false ]]; then
    if [[ "$OS" == "Linux" ]]; then
        if [[ "$ARCH" == "x86_64" ]]; then
            MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
        elif [[ "$ARCH" == "aarch64" ]]; then
            MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh"
        else
            echo "Unsupported Linux architecture: $ARCH"
            exit 1
        fi
    elif [[ "$OS" == "Darwin" ]]; then
        if [[ "$ARCH" == "arm64" ]]; then
            MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh"
        else
            echo "Unsupported macOS architecture: $ARCH"
            exit 1
        fi
    else
        echo "Unsupported OS: $OS"
        exit 1
    fi

    echo "Downloading Miniconda..."
    curl -L -o miniconda.sh "$MINICONDA_URL"

    echo "Installing Miniconda to \$HOME/miniconda3"
    bash miniconda.sh -b -p "$HOME/miniconda3"

    rm miniconda.sh

    source "$HOME/miniconda3/etc/profile.d/conda.sh"
else
    # Make sure conda is usable in non-interactive shell
    CONDA_BASE="$(conda info --base)"
    source "$CONDA_BASE/etc/profile.d/conda.sh"
fi

# --- Create conda environment if missing ---
if ! conda env list | grep -q modeling_studio; then
    echo "Creating conda environment..."
    conda env create -f requirements.yml
else
    echo "Conda environment already exists."
fi

# --- Final instructions ---
echo
echo "=== Setup complete ==="
echo "Next steps:"
echo
echo "1. Activate the environment:"
echo "   conda activate modeling_studio"
echo
echo "2. Change to the GUI directory:"
echo "   cd gui"
echo
echo "3. Run the application:"
echo "   python top6meta.py"
echo