#!/bin/bash

# Script to check if Python can detect CUDA and GPU, and reinstall conda packages if needed
# Usage: bash check_cuda_and_reinstall.sh [environment_name] [environment_yml_file]

set -e

ENV_NAME=${1:-"scTranslator_env"}
ENV_YML=${2:-"$HOME/workspace/analyses/protein_prediction_publication/repo/envs/CrossModalNet_env.yml"}

echo "Checking CUDA availability in Python..."

# Function to check CUDA availability
check_cuda() {
    python -c "
import torch
import sys
print(f'Python executable: {sys.executable}')
print(f'PyTorch version: {torch.__version__}')
print(f'CUDA available: {torch.cuda.is_available()}')
if torch.cuda.is_available():
    print(f'CUDA version: {torch.version.cuda}')
    print(f'GPU count: {torch.cuda.device_count()}')
    for i in range(torch.cuda.device_count()):
        print(f'GPU {i}: {torch.cuda.get_device_name(i)}')
    exit(0)
else:
    print('CUDA is NOT available!')
    exit(1)
"
}

# Function to reinstall PyTorch with CUDA
reinstall_pytorch_cuda() {
    echo "CUDA not detected. Reinstalling PyTorch with CUDA support..."
    
    # Load CUDA module first
    echo "Loading CUDA module..."
    module load linux-rhel8-x86_64/Core/cuda/12.3.0-cmkgc57
    
    # Check if we're in a conda environment
    if [[ -n "$CONDA_DEFAULT_ENV" ]]; then
        echo "Currently in conda environment: $CONDA_DEFAULT_ENV"
        
        # Force reinstall PyTorch with CUDA
        echo "Uninstalling existing PyTorch..."
        conda uninstall -y pytorch torchvision torchaudio pytorch-cuda || true
        
        echo "Installing PyTorch with CUDA 11.8..."
        conda install -y pytorch torchvision torchaudio pytorch-cuda=11.7 -c pytorch -c nvidia
        
    else
        echo "Not in a conda environment. Activating $ENV_NAME..."
        source activate $ENV_NAME
        
        # Force reinstall PyTorch with CUDA
        echo "Uninstalling existing PyTorch..."
        conda uninstall -y pytorch torchvision torchaudio pytorch-cuda || true
        
        echo "Installing PyTorch with CUDA 11.8..."
        conda install -y pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia
    fi
}

# Function to recreate environment from scratch
recreate_environment() {
    echo "Recreating conda environment from scratch with CUDA loaded..."
    
    # Load CUDA module first
    echo "Loading CUDA module..."
    module load linux-rhel8-x86_64/Core/cuda/12.3.0-cmkgc57
    
    # Remove existing environment
    echo "Removing existing environment $ENV_NAME..."
    conda env remove -n $ENV_NAME -y || true
    
    # Create new environment with CUDA available
    if [[ -f "$ENV_YML" ]]; then
        echo "Creating environment from $ENV_YML..."
        conda env create -f $ENV_YML -n $ENV_NAME
    else
        echo "Creating basic environment with PyTorch CUDA..."
        conda create -n $ENV_NAME python=3.9 -y
        source activate $ENV_NAME
        conda install -y pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia
    fi
    
    # Activate the new environment
    source activate $ENV_NAME
}

# Main logic
echo "=== CUDA Detection Check ==="
echo "Environment: $ENV_NAME"
echo "YAML file: $ENV_YML"
echo "Current conda env: $CONDA_DEFAULT_ENV"

# First attempt: Check if CUDA is already working
if check_cuda; then
    echo "✅ CUDA is working correctly!"
    exit 0
else
    echo "❌ CUDA not detected. Attempting fixes..."
    
    # Try reinstalling PyTorch first (faster)
    echo "=== Attempting PyTorch reinstall ==="
    reinstall_pytorch_cuda
    
    # Check again
    if check_cuda; then
        echo "✅ CUDA working after PyTorch reinstall!"
        exit 0
    else
        echo "❌ PyTorch reinstall didn't work. Recreating environment..."
        
        # Nuclear option: recreate environment
        echo "=== Recreating environment ==="
        recreate_environment
        
        # Final check
        if check_cuda; then
            echo "✅ CUDA working after environment recreation!"
            exit 0
        else
            echo "❌ CUDA still not working after environment recreation!"
            echo "Manual intervention required."
            exit 1
        fi
    fi
fi